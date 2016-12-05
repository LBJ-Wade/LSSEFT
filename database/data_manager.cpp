//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include <iostream>
#include <sstream>
#include <assert.h>

#include <set>
#include <unordered_set>

#include "data_manager.h"
#include "data_manager_impl.h"

#include "sqlite3_detail/utilities.h"
#include "sqlite3_detail/operations.h"

#include "utilities/formatter.h"

#include "defaults.h"

#include "boost/timer/timer.hpp"


data_manager::data_manager(const boost::filesystem::path& c)
  : container(c),
    handle(nullptr),   // try to catch handle-not-initialized errors
    policy(),
    FRW_model_tol(LSSEFT_DEFAULT_FRW_MODEL_PARAMETER_TOLERANCE),
    z_tol(LSSEFT_DEFAULT_REDSHIFT_CONFIGURATION_TOLERANCE),
    k_tol(LSSEFT_DEFAULT_WAVENUMBER_CONFIGURATION_TOLERANCE)
  {
    // check whether container already exists
    if(boost::filesystem::exists(container))
      {
        if(boost::filesystem::is_regular_file(container))
          {
            if(sqlite3_open_v2(container.string().c_str(), &handle, SQLITE_OPEN_READWRITE, nullptr) != SQLITE_OK)
              {
                std::ostringstream msg;
                msg << ERROR_DATABASE_SQLITE_OPEN_FAILED << " " << container;
                throw runtime_exception(exception_type::database_error, msg.str());
              }

            return;
          }
        else
          {
            std::ostringstream msg;
            msg << ERROR_DATABASE_IS_NOT_FILE_A << " " << container << " " << ERROR_DATABASE_IS_NOT_FILE_B;
            throw runtime_exception(exception_type::database_error, msg.str());
          }
      }

    // if we get to here, the container does not already exist so we should create it
    if(sqlite3_open_v2(container.string().c_str(), &handle, SQLITE_OPEN_CREATE | SQLITE_OPEN_READWRITE, nullptr) != SQLITE_OK)
      {
        std::ostringstream msg;
        msg << ERROR_DATABASE_SQLITE_CREATE_FAILED << " " << container;
        throw runtime_exception(exception_type::database_error, msg.str());
      }

    // set up tables
    sqlite3_operations::create_tables(handle, policy);
  }


data_manager::~data_manager()
  {
    assert(this->handle != nullptr);

    // perform routine maintenance on the container
    sqlite3_operations::exec(this->handle, "VACUUM;");

    sqlite3_close(this->handle);
  }


// TOKENIZATION

// TODO: consider simplifying this code via traits and templates


std::unique_ptr<FRW_model_token> data_manager::tokenize(const FRW_model& obj)
  {
    // open a new transaction on the database
    std::shared_ptr<transaction_manager> transaction = this->open_transaction();

    // lookup id for this model, or generate one if it does not already exist
    std::unique_ptr<FRW_model_token> id = this->tokenize(*transaction, obj);

    // commit the transaction
    transaction->commit();

    return std::move(id);
  }


std::unique_ptr<FRW_model_token> data_manager::tokenize(transaction_manager& mgr, const FRW_model& obj)
  {
    // lookup id for this model, or generate one if it does not already exist
    unsigned int id = this->lookup_or_insert(mgr, obj);
    return std::make_unique<FRW_model_token>(id);
  }


std::unique_ptr<z_token> data_manager::tokenize(double z)
  {
    // open a new transaction on the database
    std::shared_ptr<transaction_manager> transaction = this->open_transaction();

    // lookup id for this redshift, or generate one if it does not already exist
    std::unique_ptr<z_token> id = this->tokenize(*transaction, z);

    // commit the transaction
    transaction->commit();

    return std::move(id);
  }


std::unique_ptr<z_token> data_manager::tokenize(transaction_manager& mgr, double z)
  {
    // lookup id for this redshift, or generate one if it does not already exist
    unsigned int id = this->lookup_or_insert(mgr, z);
    return std::make_unique<z_token>(id);
  }


template <typename Token>
std::unique_ptr<Token> data_manager::tokenize(const Mpc_units::energy& k)
  {
    // open a new transaction on the database
    std::shared_ptr<transaction_manager> transaction = this->open_transaction();

    // lookup id for this wavenumber, or generate one if it does not already exist
    std::unique_ptr<Token> id = this->tokenize<Token>(*transaction, k);

    // commit the transaction
    transaction->commit();

    return std::move(id);
  }


template <typename Token>
std::unique_ptr<Token> data_manager::tokenize(transaction_manager& mgr, const Mpc_units::energy& k)
  {
    // lookup id for this wavenumber, or generate one if it does not already exist
    unsigned int id = this->lookup_or_insert<Token>(mgr, k);
    return std::make_unique<Token>(id);
  }


std::unique_ptr<linear_Pk_token>
data_manager::tokenize(const FRW_model_token& model, const linear_Pk& Pk_lin)
  {
    // open a new transaction on the database
    std::shared_ptr<transaction_manager> transaction = this->open_transaction();
    
    // lookup id for this power spectrum, or generate one if it doesn't already exist
    std::unique_ptr<linear_Pk_token> id = this->tokenize(*transaction, model, Pk_lin);
    
    // commit the transaction
    transaction->commit();
    
    return std::move(id);
  }


std::unique_ptr<linear_Pk_token>
data_manager::tokenize(transaction_manager& mgr, const FRW_model_token& model, const linear_Pk& Pk_lin)
  {
    // lookup id for this power spectrum, or generate one if it doesn't already exist
    unsigned int id = this->lookup_or_insert(mgr, model, Pk_lin);
    return std::make_unique<linear_Pk_token>(id);
  }


// TRANSACTIONS


std::shared_ptr<transaction_manager> data_manager::open_transaction()
  {
    // check whether a transaction is already in progress; if so, raise an exception
    std::shared_ptr<transaction_manager> check = this->current_transaction.lock();
    if(check) throw runtime_exception(exception_type::transaction_error, ERROR_TRANSACTION_IN_PROGRESS);
    check.reset();

    // create a new transaction manager
    transaction_manager::open_handler     do_open     = std::bind(&data_manager::begin_transaction, this);
    transaction_manager::commit_handler   do_commit   = std::bind(&data_manager::commit_transaction, this);
    transaction_manager::rollback_handler do_rollback = std::bind(&data_manager::rollback_transaction, this);
    transaction_manager::release_handler  do_release  = std::bind(&data_manager::release_transaction, this);

    std::shared_ptr<transaction_manager> transaction = std::make_shared<transaction_manager>(do_open, do_commit, do_rollback, do_release);

    // record this transaction
    this->current_transaction = transaction;

    return(transaction);
  }


void data_manager::begin_transaction()
  {
    assert(this->handle != nullptr);
    sqlite3_operations::exec(this->handle, "BEGIN TRANSACTION");
  }


void data_manager::commit_transaction()
  {
    assert(this->handle != nullptr);
    sqlite3_operations::exec(this->handle, "COMMIT");
  }


void data_manager::rollback_transaction()
  {
    assert(this->handle != nullptr);
    sqlite3_operations::exec(this->handle, "ROLLBACK");
  }


void data_manager::release_transaction()
  {
    assert(this->handle != nullptr);

    // check whether a transaction is already in progress
    std::shared_ptr<transaction_manager> check = this->current_transaction.lock();
    if(!check) throw runtime_exception(exception_type::transaction_error, ERROR_NO_TRANSACTION_IN_PROGRESS);
    check.reset();

    this->current_transaction.reset();
  }


// LOOKUP AND INSERT


unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, const FRW_model& obj)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_FRW_model(this->handle, mgr, obj, this->policy, this->FRW_model_tol);
    if(id) return(*id);

    return sqlite3_operations::insert_FRW_model(this->handle, mgr, obj, this->policy);
  }


unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, double z)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_redshift(this->handle, mgr, z, this->policy, this->z_tol);
    if(id) return(*id);

    return sqlite3_operations::insert_redshift(this->handle, mgr, z, this->policy);
  }


template <typename Token>
unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, const Mpc_units::energy& k)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_wavenumber<Token>(this->handle, mgr, k, this->policy, this->k_tol);
    if(id) return(*id);

    return sqlite3_operations::insert_wavenumber<Token>(this->handle, mgr, k, this->policy);
  }


unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, const FRW_model_token& model, const linear_Pk& Pk_lin)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_Pk_linear(this->handle, mgr, model, Pk_lin, this->policy);
    if(id) return(*id);
    
    return sqlite3_operations::insert_Pk_linear(this->handle, mgr, model, Pk_lin, this->policy);
  }


// GENERATE REDSHIFT AND WAVENUMBER DATABASES


std::unique_ptr<z_database> data_manager::build_redshift_db(range<double>& sample)
  {
    // construct an empty redshift database
    std::unique_ptr<z_database> z_db(new z_database);

    // grab the grid of redshift samples
    const std::vector<double>& z_samples = sample.grid();

    for(std::vector<double>::const_iterator t = z_samples.begin(); t != z_samples.end(); ++t)
      {
        std::unique_ptr<z_token> tok = this->tokenize(*t);
        z_db->add_record(*t, *tok);
      }

    return(z_db);
  }


template <typename Token>
std::unique_ptr< wavenumber_database<Token> > data_manager::build_wavenumber_db(range<Mpc_units::energy>& sample)
  {
    // construct an empty wavenumber database
    std::unique_ptr< wavenumber_database<Token> > k_db = std::make_unique< wavenumber_database<Token> >();

    // grab the grid of wavenumber samples
    const std::vector<Mpc_units::energy>& k_samples = sample.grid();

    for(std::vector<Mpc_units::energy>::const_iterator t = k_samples.begin(); t != k_samples.end(); ++t)
      {
        std::unique_ptr<Token> tok = this->tokenize<Token>(*t);
        k_db->add_record(*t, *tok);
      }

    return(k_db);
  }


std::unique_ptr<k_database> data_manager::build_k_db(range<Mpc_units::energy>& sample)
  {
    return this->build_wavenumber_db<k_token>(sample);
  }


std::unique_ptr<IR_cutoff_database> data_manager::build_IR_cutoff_db(range<Mpc_units::energy>& sample)
  {
    return this->build_wavenumber_db<IR_cutoff_token>(sample);
  }


std::unique_ptr<UV_cutoff_database> data_manager::build_UV_cutoff_db(range<Mpc_units::energy>& sample)
  {
    return this->build_wavenumber_db<UV_cutoff_token>(sample);
  }


std::unique_ptr<IR_resum_database> data_manager::build_IR_resum_db(range<Mpc_units::energy>& sample)
  {
    return this->build_wavenumber_db<IR_resum_token>(sample);
  }


std::unique_ptr<k_database> data_manager::build_k_db(transaction_manager& mgr, linear_Pk& Pk_lin)
  {
    // construct an empty wavenumber database
    std::unique_ptr<k_database> k_db = std::make_unique<k_database>();
    
    // get power spectrum database underlying this container
    const tree_Pk::database_type& Pk_db = Pk_lin.get_db();
    
    for(tree_Pk::database_type::const_record_iterator t = Pk_db.record_cbegin(); t != Pk_db.record_cend(); ++t)
      {
        // ask linear_Pk container whether this P(k) value is acceptable
        const Mpc_units::energy& k = t->get_wavenumber();
        if(Pk_lin.is_valid(k))
          {
            std::unique_ptr<k_token> tok = this->tokenize<k_token>(mgr, k);
            k_db->add_record(k, *tok);
          }
      }
    
    return k_db;
  }


std::unique_ptr<transfer_work_list>
data_manager::build_transfer_work_list(FRW_model_token& model, k_database& k_db, z_database& z_db)
  {
    // start timer
    boost::timer::cpu_timer timer;

    // construct an empty work list
    std::unique_ptr<transfer_work_list> work_list = std::make_unique<transfer_work_list>();

    // open a transaction on the database
    std::shared_ptr<transaction_manager> mgr = this->open_transaction();

    // set up temporary table of desired z identifiers
    std::string z_table = sqlite3_operations::z_table(this->handle, *mgr, this->policy, z_db);

    // for each wavenumber in k_db, find which z-values are missing
    for(k_database::record_iterator t = k_db.record_begin(); t != k_db.record_end(); ++t)
      {
//        std::cout << "lsseft: checking missing redshift values for k = " << (*(*t) * Mpc_units::Mpc) << " h/Mpc = " << (*(*t)) / Mpc_units::eV << " eV" << '\n';

        // get a database of missing redshifts for this k-value.
        // sqlite3_operations::missing_redshifts() returns a std::unique_ptr which transfers ownership,
        // but we want to convert that to a std::shared_ptr which is what transfer_work_item expects,
        // because it shares ownership with objects representing MPI messages
        // (see comments in transfer_work_item constructor)
        std::shared_ptr<z_database> missing(
          std::move(sqlite3_operations::missing_transfer_redshifts(this->handle, *mgr, this->policy, model,
                                                                   t->get_token(), z_db, z_table)));

        // if any redshifts were missing, set up a record in the work list
        if(missing)
          {
//            std::cout << "  -- " << missing_values->size() << " redshifts" << '\n';

            work_list->emplace_back(*(*t), t->get_token(), missing);
          }
      }

    // drop unneeded temporary tables
    sqlite3_operations::drop_temp(this->handle, *mgr, z_table);

    // commit the transaction before allowing it to go out of scope
    mgr->commit();

    timer.stop();
    std::cout << "lsseft: constructed transfer function work list (" << work_list->size() << " items) in time " << format_time(timer.elapsed().wall) << '\n';
    
    // release list if it contains no work
    if(work_list->empty()) work_list.release();

    return(work_list);
  }


std::unique_ptr<z_database> data_manager::build_oneloop_work_list(FRW_model_token& model, z_database& z_db)
  {
    // start timer
    boost::timer::cpu_timer timer;

    // open a transaction on the database
    std::shared_ptr<transaction_manager> mgr = this->open_transaction();

    // set up temporary table of desired z identifiers
    std::string z_table = sqlite3_operations::z_table(this->handle, *mgr, this->policy, z_db);
    
    std::unique_ptr<z_database> work_list =
      sqlite3_operations::missing_oneloop_growth_redshifts(this->handle, *mgr, this->policy, model,
                                                           z_db, z_table);

    // drop unneeded temporary tables
    sqlite3_operations::drop_temp(this->handle, *mgr, z_table);

    // close transaction
    mgr->commit();

    timer.stop();
    std::cout << "lsseft: constructed one-loop growth factor work list ("
              << (work_list ? work_list->size() : 0) << " items) in time "
              << format_time(timer.elapsed().wall) << '\n';

    return(work_list);
  }


loop_configs data_manager::tensor_product(k_database& k_db, IR_cutoff_database& IR_db, UV_cutoff_database& UV_db)
  {
    loop_configs tensor_prod;
    
    for(k_database::const_record_iterator t = k_db.record_cbegin(); t != k_db.record_cend(); ++t)
      {
        for(UV_cutoff_database::const_record_iterator u = UV_db.record_cbegin(); u != UV_db.record_cend(); ++u)
          {
            for(IR_cutoff_database::const_record_iterator v = IR_db.record_cbegin(); v != IR_db.record_cend(); ++v)
              {
                tensor_prod.emplace(t, u, v);
              }
          }
      }
    
    return tensor_prod;
  }


resum_Pk_configs data_manager::tensor_product(k_database& k_db, IR_cutoff_database& IR_cutoff_db,
                                              UV_cutoff_database& UV_cutoff_db, IR_resum_database& IR_resum_db)
  {
    resum_Pk_configs tensor_prod;
    
    for(k_database::const_record_iterator t = k_db.record_cbegin(); t != k_db.record_cend(); ++t)
      {
        for(UV_cutoff_database::const_record_iterator u = UV_cutoff_db.record_cbegin(); u != UV_cutoff_db.record_cend(); ++u)
          {
            for(IR_cutoff_database::const_record_iterator v = IR_cutoff_db.record_cbegin(); v != IR_cutoff_db.record_cend(); ++v)
              {
                for(IR_resum_database::const_record_iterator w = IR_resum_db.record_cbegin(); w != IR_resum_db.record_cend(); ++w)
                  {
                    tensor_prod.emplace(t, u, v, w);
                  }
              }
          }
      }
    
    return tensor_prod;
  }


std::unique_ptr<loop_momentum_work_list>
data_manager::build_loop_momentum_work_list(FRW_model_token& model, k_database& k_db,
                                            IR_cutoff_database& IR_db, UV_cutoff_database& UV_db,
                                            std::shared_ptr<tree_Pk>& Pk)
  {
    // start timer
    boost::timer::cpu_timer timer;

    // construct an empty work list
    std::unique_ptr<loop_momentum_work_list> work_list = std::make_unique<loop_momentum_work_list>();

    // open a transaction on the database
    std::shared_ptr<transaction_manager> mgr = this->open_transaction();

    // tensor together the desired k-values with the IR and UV cutoffs to obtain a set of
    // desired combinations
    loop_configs required_configs = this->tensor_product(k_db, IR_db, UV_db);
    
    // obtain set of configurations that actually need to be computed, ie. are not already present
    // in the database
    loop_configs missing = sqlite3_operations::missing_loop_integral_configurations(this->handle, *mgr, this->policy,
                                                                                    model, required_configs);
    
    // add these missing configurations to the work list
    for(const loop_configs::value_type& record : missing)
      {
        work_list->emplace_back(*(*record.k), record.k->get_token(), *(*record.UV_cutoff),
                                record.UV_cutoff->get_token(), *(*record.IR_cutoff), record.IR_cutoff->get_token(), Pk);
      }

    // close transaction
    mgr->commit();

    timer.stop();
    std::cout << "lsseft: constructed loop momentum work list (" << work_list->size() << " items) in time " << format_time(timer.elapsed().wall) << '\n';

    // release list if it contains no work
    if(work_list->empty()) work_list.release();

    return(work_list);
  }


template <>
oneloop_growth data_manager::find<oneloop_growth>(transaction_manager& mgr, const FRW_model_token& model, z_database& z_db)
  {
    // construct payload and ask SQLite backend to populate it
    oneloop_growth payload(std::move(sqlite3_operations::find(this->handle, mgr, this->policy, model, z_db)));
    
    return std::move(payload);
  }


template <>
loop_integral
data_manager::find<loop_integral>(transaction_manager& mgr, const FRW_model_token& model, const k_token& k,
                                  const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff)
  {
    // construct payload and ask SQLite backend to populate it
    loop_integral payload(std::move(
      sqlite3_operations::find(this->handle, mgr, this->policy, model, k, IR_cutoff, UV_cutoff))
    );
    
    return std::move(payload);
  }


template <>
oneloop_Pk
data_manager::find<oneloop_Pk>(transaction_manager& mgr, const FRW_model_token& model, const k_token& k, const z_token& z,
                               const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff)
  {
    // construct payload and ask SQLite backend to populate it
    oneloop_Pk payload(std::move(
      sqlite3_operations::find(this->handle, mgr, this->policy, model, k, z, IR_cutoff, UV_cutoff))
    );
    
    return std::move(payload);
  }


template <>
Matsubara_A
data_manager::find(transaction_manager& mgr, const FRW_model_token& model, const IR_resum_token& IR_resum)
  {
    // construct payload nad ask SQLite backend to populate it
    Matsubara_A payload(std::move(
      sqlite3_operations::find(this->handle, mgr, this->policy, model, IR_resum))
    );
    
    return std::move(payload);
  }


std::unique_ptr<one_loop_Pk_work_list>
data_manager::build_one_loop_Pk_work_list(FRW_model_token& model, z_database& z_db, k_database& k_db,
                                          IR_cutoff_database& IR_db, UV_cutoff_database& UV_db,
                                          std::shared_ptr<tree_Pk>& Pk)
  {
    // start timer
    boost::timer::cpu_timer timer;
    
    // construct an empty work list
    std::unique_ptr<one_loop_Pk_work_list> work_list = std::make_unique<one_loop_Pk_work_list>();
    
    // open a transaction on the database
    std::shared_ptr<transaction_manager> mgr = this->open_transaction();
    
    // set up temporary table of desired z identifiers
    std::string z_table = sqlite3_operations::z_table(this->handle, *mgr, this->policy, z_db);
    
    // tensor together the desired k-values with the UV and IR cutoffs
    loop_configs required_configs = this->tensor_product(k_db, IR_db, UV_db);
    
    for(const loop_configs::value_type& record : required_configs)
      {
        // find redshifts that are missing for this configuration, if any
        std::unique_ptr<z_database> missing_zs =
          sqlite3_operations::missing_one_loop_Pk_redshifts(this->handle, *mgr, this->policy, model,
                                                            z_table, z_db, record);
        
        // schedule a task to compute any missing redshifts
        if(missing_zs)
          {
            std::shared_ptr<oneloop_growth> g = std::make_shared<oneloop_growth>(
              this->find<oneloop_growth>(*mgr, model, z_db));
            std::shared_ptr<loop_integral> l = std::make_shared<loop_integral>(
              this->find<loop_integral>(*mgr, model, record.k->get_token(), record.IR_cutoff->get_token(), record.UV_cutoff->get_token()));

            work_list->emplace_back(*(*record.k), g, l, Pk);
          }
      }
    
    // drop unneeded temporary tables
    sqlite3_operations::drop_temp(this->handle, *mgr, z_table);
    
    timer.stop();
    std::cout << "lsseft: constructed one-loop P(k) work list (" << work_list->size() << " items) in time " << format_time(timer.elapsed().wall) << '\n';
    
    // close transaction
    mgr->commit();
    
    // release list if it contains no work
    if(work_list->empty()) work_list.release();
    
    return work_list;
  }


std::unique_ptr<multipole_Pk_work_list>
data_manager::build_multipole_Pk_work_list(FRW_model_token& model, z_database& z_db, k_database& k_db,
                                           IR_cutoff_database& IR_cutoff_db, UV_cutoff_database& UV_cutoff_db,
                                           IR_resum_database& IR_resum_db, std::shared_ptr<tree_Pk>& Pk)
  {
    // start timer
    boost::timer::cpu_timer timer;
    
    // construct an empty work list
    std::unique_ptr<multipole_Pk_work_list> work_list = std::make_unique<multipole_Pk_work_list>();
    
    // open a transaction on the database
    std::shared_ptr<transaction_manager> mgr = this->open_transaction();
    
    // set up temporary table of desired z identifiers
    std::string z_table = sqlite3_operations::z_table(this->handle, *mgr, this->policy, z_db);
    
    // tensor together the desired k-values with the UV and IR cutoffs
    resum_Pk_configs required_configs = this->tensor_product(k_db, IR_cutoff_db, UV_cutoff_db, IR_resum_db);
    
    for(const resum_Pk_configs::value_type& record : required_configs)
      {
        // find redshifts that are missing for this configuration, if any
        std::unique_ptr<z_database> missing_zs =
          sqlite3_operations::missing_multipole_Pk_redshifts(this->handle, *mgr, this->policy, model,
                                                             z_table, z_db, record);
        
        // schedule a task to compute any missing redshifts
        if(missing_zs)
          {
            oneloop_growth gf_data = this->find<oneloop_growth>(*mgr, model, *missing_zs);
            
            for(oneloop_growth::const_iterator t = gf_data.cbegin(); t != gf_data.cend(); ++t)
              {
                // lookup one-loop data for this redshift and loop configuration
                std::shared_ptr<oneloop_Pk> loop_data = std::make_shared<oneloop_Pk>(
                  this->find<oneloop_Pk>(*mgr, model, record.k->get_token(), (*t).first,
                                         record.IR_cutoff->get_token(), record.UV_cutoff->get_token())
                );
                
                // lookup Matsubara-A coefficient for this IR resummation scale
                Matsubara_A A_coeff = this->find<Matsubara_A>(*mgr, model, record.IR_resum->get_token());
    
                work_list->emplace_back(*(*record.k), A_coeff, loop_data, (*t).second, Pk);
              }
          }
      }
    
    // drop unneeded temporary tables
    sqlite3_operations::drop_temp(this->handle, *mgr, z_table);
    
    // close transaction
    mgr->commit();
    
    timer.stop();
    std::cout << "lsseft: constructed one-loop multipole P(k) work list (" << work_list->size() << " items) in time " << format_time(timer.elapsed().wall) << '\n';
    
    // release list if it contains no work
    if(work_list->empty()) work_list.release();
    
    return work_list;
  }


std::unique_ptr<Matsubara_A_work_list>
data_manager::build_Matsubara_A_work_list(FRW_model_token& model, IR_resum_database& IR_resum_db,
                                          std::shared_ptr<tree_Pk>& Pk)
  {
    // start timer
    boost::timer::cpu_timer timer;
    
    // construct an empty work list
    std::unique_ptr<Matsubara_A_work_list> work_list = std::make_unique<Matsubara_A_work_list>();
    
    // open a transaction on the database
    std::shared_ptr<transaction_manager> mgr = this->open_transaction();
    
    // obatain list of missing configurations
    Matsubara_configs missing = sqlite3_operations::missing_Matsubara_A_configurations(this->handle, *mgr, this->policy,
                                                                                       model, IR_resum_db);
    
    // add these configurations to the work list
    for(const Matsubara_configs::value_type& record : missing)
      {
        work_list->emplace_back(*(*record.IR_resum), record.IR_resum->get_token(), Pk);
      }
    
    // close transaction
    mgr->commit();
    
    timer.stop();
    std::cout << "lsseft: constructed Matsubara-A work list (" << work_list->size() << " items) in time " << format_time(timer.elapsed().wall) << '\n';
    
    // release list if it contains no work
    if(work_list->empty()) work_list.release();
    
    return work_list;
  }


std::unique_ptr<filter_Pk_work_list>
data_manager::build_filter_Pk_work_list(linear_Pk_token& token, std::shared_ptr<linear_Pk>& Pk_lin)
  {
    // start timer
    boost::timer::cpu_timer timer;
    
    // construct an empty work list
    std::unique_ptr<filter_Pk_work_list> work_list = std::make_unique<filter_Pk_work_list>();
    
    // open a transaction on the database
    std::shared_ptr<transaction_manager> mgr = this->open_transaction();
    
    // set up temporary table of desired wavenumber identifiers
    std::unique_ptr<k_database> k_db = this->build_k_db(*mgr, *Pk_lin);
    std::string k_table = sqlite3_operations::k_table(this->handle, *mgr, this->policy, *k_db);
    
    // obtain list of missing configurations
    std::unique_ptr<k_database> missing =
      sqlite3_operations::missing_filter_Pk_wavenumbers(this->handle, *mgr, this->policy, token, *k_db, k_table);
    
    if(missing)
      {
        // add these configurations to the work list
        for(k_database::const_record_iterator t = missing->record_cbegin(); t != missing->record_cend(); ++t)
          {
            work_list->emplace_back(*(*t), t->get_token(), Pk_lin, token);
          }
      }
    
    // drop temporary table
    sqlite3_operations::drop_temp(this->handle, *mgr, k_table);
    
    // close transaction
    mgr->commit();
    
    timer.stop();
    std::cout << "lsseft: constructed no-wiggle filter work list (" << work_list->size() << " items) in time " << format_time(timer.elapsed().wall) << '\n';
    
    // release list if it contains to work
    if(work_list->empty()) work_list.release();
    
    return work_list;
  }
