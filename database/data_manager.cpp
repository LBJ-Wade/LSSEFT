//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include <iostream>
#include <sstream>
#include <assert.h>

#include "data_manager.h"

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
    unsigned int id = this->lookup_or_insert(*transaction, obj);

    // commit the transaction
    transaction->commit();

    return std::unique_ptr<FRW_model_token>(new FRW_model_token(id));
  }


std::unique_ptr<redshift_token> data_manager::tokenize(double z)
  {
    // open a new transaction on the database
    std::shared_ptr<transaction_manager> transaction = this->open_transaction();

    // lookup id for this redshift, or generate one if it does not already exist
    unsigned int id = this->lookup_or_insert(*transaction, z);

    // commit the transaction
    transaction->commit();

    return std::unique_ptr<redshift_token>(new redshift_token(id));
  }


std::unique_ptr<wavenumber_token> data_manager::tokenize(const eV_units::energy& k)
  {
    // open a new transaction on the database
    std::shared_ptr<transaction_manager> transaction = this->open_transaction();

    // lookup id for this wavenumber, or generate one if it does not already exist
    unsigned int id = this->lookup_or_insert(*transaction, k);

    // commit the transaction
    transaction->commit();

    return std::unique_ptr<wavenumber_token>(new wavenumber_token(id));
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


unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, const eV_units::energy& k)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_wavenumber(this->handle, mgr, k, this->policy, this->z_tol);
    if(id) return(*id);

    return sqlite3_operations::insert_wavenumber(this->handle, mgr, k, this->policy);
  }


// GENERATE REDSHIFT AND WAVENUMBER DATABASES


std::unique_ptr<redshift_database> data_manager::build_db(range<double>& sample)
  {
    // construct an empty redshift database
    std::unique_ptr<redshift_database> z_db(new redshift_database);

    // grab the grid of redshift samples
    const std::vector<double>& z_samples = sample.grid();

    for(std::vector<double>::const_iterator t = z_samples.begin(); t != z_samples.end(); ++t)
      {
        std::unique_ptr<redshift_token> tok = this->tokenize(*t);
        z_db->add_record(*t, *tok);
      }

    return(z_db);
  }


std::unique_ptr<wavenumber_database> data_manager::build_db(range<eV_units::energy>& sample)
  {
    // construct an empty wavenumber database
    std::unique_ptr<wavenumber_database> k_db(new wavenumber_database);

    // grab the grid of wavenumber samples
    const std::vector<eV_units::energy>& k_samples = sample.grid();

    for(std::vector<eV_units::energy>::const_iterator t = k_samples.begin(); t != k_samples.end(); ++t)
      {
        std::unique_ptr<wavenumber_token> tok = this->tokenize(*t);
        k_db->add_record(*t, *tok);
      }

    return(k_db);
  }


std::unique_ptr<transfer_work_list> data_manager::build_transfer_work_list(FRW_model_token& model, wavenumber_database& k_db,
                                                                           redshift_database& z_db)
  {
    // start timer
    boost::timer::auto_cpu_timer timer;

    // construct an empty work list
    std::unique_ptr<transfer_work_list> work_list(new transfer_work_list);

    // open a transaction on the database
    std::shared_ptr<transaction_manager> transaction = this->open_transaction();

    // set up temporary table of desired z identifiers
    std::string z_table = sqlite3_operations::z_table(this->handle, *transaction, this->policy, z_db);

    // for each wavenumber in k_db, find which z-values are missing
    for(wavenumber_database::record_iterator t = k_db.record_begin(); t != k_db.record_end(); ++t)
      {
//        std::cout << "lsseft: checking missing redshift values for k = " << (*(*t) * eV_units::Mpc) << " h/Mpc = " << static_cast<double>(*(*t)) << " eV" << '\n';

        // get a database of missing redshifts for this k-value
        std::shared_ptr<redshift_database> missing_values = sqlite3_operations::missing_redshifts(this->handle, *transaction, this->policy, model, t->get_token(), z_db, z_table);

        // if any redshifts were missing, set up a record in the work list
        if(missing_values)
          {
//            std::cout << "  -- " << missing_values->size() << " redshifts" << '\n';

            work_list->emplace_back(*(*t), t->get_token(), missing_values);
          }
      }

    // drop unneeded temporary tables
    sqlite3_operations::drop_temp(this->handle, *transaction, z_table);

    // commit the transaction before allowing it to go out of scope
    transaction->commit();

    timer.stop();
    std::cout << "lsseft: constructed transfer function work list in time " << format_time(timer.elapsed().wall) << '\n';

    return(work_list);
  }


void data_manager::store(const FRW_model_token& model_token, const transfer_function& sample)
  {
    // open a transaction on the database
    std::shared_ptr<transaction_manager> transaction = this->open_transaction();

    sqlite3_operations::store(this->handle, *transaction, this->policy, model_token, sample);

    // commit the transaction before allowing it to go out of scope
    transaction->commit();
  }
