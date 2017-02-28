//
// Created by David Seery on 11/08/2015.
// --@@ // Copyright (c) 2017 University of Sussex. All rights reserved.
//
// This file is part of the Sussex Effective Field Theory for
// Large-Scale Structure platform (LSSEFT).
//
// LSSEFT is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// LSSEFT is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with LSSEFT.  If not, see <http://www.gnu.org/licenses/>.
//
// @license: GPL-2
// @contributor: David Seery <D.Seery@sussex.ac.uk>
// --@@
//

#ifndef LSSEFT_DATA_MANAGER_H
#define LSSEFT_DATA_MANAGER_H


#include <memory>

#include "tokens.h"
#include "transaction_manager.h"
#include "z_database.h"
#include "k_database.h"
#include "IR_cutoff_database.h"
#include "UV_cutoff_database.h"
#include "IR_resum_database.h"

#include "cosmology/types.h"
#include "cosmology/FRW_model.h"
#include "cosmology/concepts/transfer_function.h"
#include "cosmology/concepts/oneloop_growth.h"
#include "cosmology/concepts/range.h"
#include "cosmology/concepts/power_spectrum.h"
#include "cosmology/concepts/loop_integral.h"

#include "sqlite3_detail/sqlite3_policy.h"
#include "sqlite3_detail/operations.h"

#include "utilities/formatter.h"

#include "boost/filesystem/operations.hpp"

#include "sqlite3.h"


constexpr double FILTER_PK_DEFAULT_BOTTOM_CLEARANCE = 1.25;
constexpr double FILTER_PK_DEFAULT_TOP_CLEARANCE    = 0.75;


class data_manager
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor opens SQLite3 handle, and creates tables if this is a new database
    data_manager(const boost::filesystem::path& c);

    //! destructor closes SQLite3 handle
    ~data_manager();


    // GENERATE WAVENUMBER AND REDSHIFT DATABASES

  public:

    //! generate redshift database from a set of samples
    std::unique_ptr<z_database> build_redshift_db(range<double>& sample);

    //! generate wavenumber database from a set of samples
    std::unique_ptr<k_database> build_k_db(range<Mpc_units::energy>& sample);
    
    //! generate wavenumber database from a linear power spectrum container
    //! only contains wavenumbers that can actually be evaluated by the container
    template <typename PkContainer>
    std::unique_ptr<k_database> build_k_db(transaction_manager& mgr, const PkContainer& Pk_lin,
                                           double bottom_clearance=SPLINE_PK_DEFAULT_BOTTOM_CLEARANCE,
                                           double top_clearance=SPLINE_PK_DEFAULT_TOP_CLEARANCE);

    //! generate IR cutoff database from a set of samples
    std::unique_ptr<IR_cutoff_database> build_IR_cutoff_db(range<Mpc_units::energy>& sample);

    //! generate UV cutoff database from a set of samples
    std::unique_ptr<UV_cutoff_database> build_UV_cutoff_db(range<Mpc_units::energy>& sample);
    
    //! generate IR resummation scale database from a set of samples
    std::unique_ptr<IR_resum_database> build_IR_resum_db(range<Mpc_units::energy>& sample);

  protected:

    template <typename Token>
    std::unique_ptr< wavenumber_database<Token> > build_wavenumber_db(range<Mpc_units::energy>& sample);


    // SERVICES

  public:

    //! build a work list representing z-values that are missing from the SQLite backing store
    //! for each k-value in a wavenumber database representing transfer functions.
    //! generates a new transaction on the database; will fail if a transaction is in progress
    std::unique_ptr<transfer_work_list>
    build_transfer_work_list(FRW_model_token& model, k_database& k_db, z_database& z_db);

    //! build a work list representing z-values that are missing from the SQLite backing store
    //! for each z-value needed for the one-loop growth factors.
    //! generates a new transaction on the database; will fail if a transaction is in progress
    std::unique_ptr<z_database> build_oneloop_work_list(FRW_model_token& model, z_database& z_db);

    //! build a work list representing k-values that are missing from the SQLite backing store
    //! for each k-value in a wavenumber database representing the momentum loop integral.
    //! generates a new transaction on the database; will fail if a transaction is in progress
    std::unique_ptr<loop_momentum_work_list>
    build_loop_momentum_work_list(FRW_model_token& model, k_database& k_db, IR_cutoff_database& IR_db,
                                  UV_cutoff_database& UV_db, std::shared_ptr<initial_filtered_Pk>& Pk);
    
    //! build a work list representing (k, z, IR, UV) combinations of the one-loop power spectra
    //! that are missing from the SQLite backing store.
    //! generates a new transaction on the database; will fail if a transaction is in progress
    std::unique_ptr<one_loop_Pk_work_list>
    build_one_loop_Pk_work_list(FRW_model_token& model, z_database& z_db, k_database& k_db,
                                IR_cutoff_database& IR_db, UV_cutoff_database& UV_db,
                                std::shared_ptr<initial_filtered_Pk>& Pk_init,
                                std::shared_ptr<final_filtered_Pk>& Pk_final);
    
    //! build a work list represenitng (k, z, IR_cutoff, UV_cutoff, IR_resum) combinations of the one-loop
    //! power spectrum that are missing from the SQLite backing store
    //! generates a new transaction on the database; will fail if a transaction is in progress
    std::unique_ptr<one_loop_resum_Pk_work_list>
    build_one_loop_resum_Pk_work_list(FRW_model_token& model, z_database& z_db, k_database& k_db,
                                          IR_cutoff_database& IR_cutoff_db, UV_cutoff_database& UV_cutoff_db,
                                          IR_resum_database& IR_resum_db,
                                          std::shared_ptr<initial_filtered_Pk>& Pk_init,
                                          std::shared_ptr<final_filtered_Pk>& Pk_final);
    
    //! build a work list representing (k, z, IR_cutoff, UV_cutoff, IR_resum) combinations of the one-loop
    //! multipole power spectra that are missing from the SQLite backing store.
    //! generates a new transaction on the database; will fail if a transaction is in progress
    std::unique_ptr<multipole_Pk_work_list>
    build_multipole_Pk_work_list(FRW_model_token& model, z_database& z_db, k_database& k_db,
                                     IR_cutoff_database& IR_cutoff_db, UV_cutoff_database& UV_cutoff_db,
                                     IR_resum_database& IR_resum_db,
                                     std::shared_ptr<initial_filtered_Pk>& Pk_init,
                                     std::shared_ptr<final_filtered_Pk>& Pk_final);
    
    //! build a work list representing data for calculation of the Matsubara- X & Y coefficients
    std::unique_ptr<Matsubara_XY_work_list>
    build_Matsubara_XY_work_list(FRW_model_token& model, IR_resum_database& IR_resum_db, std::shared_ptr<initial_filtered_Pk>& Pk);
    
    //! build a work list representing k-modes for which we need to produce a filtered wiggle/no-wiggle power spectrum
    std::unique_ptr<filter_Pk_work_list>
    build_filter_Pk_work_list(linear_Pk_token& token, std::shared_ptr<filterable_Pk>& Pk_lin);

    //! exchange a linear power spectrum container for a wiggle-Pk container
    template <typename PkContainer>
    std::unique_ptr<typename PkContainer::filtered_Pk_type> build_wiggle_Pk(const linear_Pk_token& token, const PkContainer& Pk_lin);
    
    //! compute how to rescale a final power spectrum to the same amplitude as an initial power spectrum
    template <typename PkContainer>
    PkContainer& rescale_final_Pk(const FRW_model_token& model, PkContainer& Pk, const z_database& z_db);
    
  protected:
    
    //! tensor together (k, IR cutoff, UV cutoff) combinations for loop integrals
    loop_configs tensor_product(k_database& k_db, IR_cutoff_database& IR_db, UV_cutoff_database& UV_db);
    
    //! tensor together (k, IR cutoff, UV cutoff, IR resummation scale) combinations for loop integrals
    resum_Pk_configs tensor_product(k_database& k_db, IR_cutoff_database& IR_cutoff_db, UV_cutoff_database& UV_cutoff_db,
                                    IR_resum_database& IR_resum_db);


    // TOKENS
    // tokens are the basic unit of currency used when interacting with the database

  public:
    
    //! tokenize an FRW model
    std::unique_ptr<FRW_model_token> tokenize(const FRW_model& obj);

    //! tokenize a redshift
    std::unique_ptr<z_token> tokenize(double z);

    //! tokenize a wavenumber of the type specified in the template
    template <typename Token>
    std::unique_ptr<Token> tokenize(const Mpc_units::energy& k);
    
    //! tokenize a linear power spectrum
    template <typename PkContainer>
    std::unique_ptr<linear_Pk_token> tokenize(const FRW_model_token& model, const PkContainer& Pk_lin);
  
  protected:
    
    //! tokenize an FRW model
    //! generates a new transaction on the database
    std::unique_ptr<FRW_model_token> tokenize(transaction_manager& mgr, const FRW_model& obj);
    
    //! tokenize a redshift
    //! generates a new transaction on the database; will fail if a transaction is in progress
    std::unique_ptr<z_token> tokenize(transaction_manager& mgr, double z);
    
    //! tokenize a wavenumber of the type specified in the template
    //! generates a new transaction on the database; will fail if a transaction is in progress
    template <typename Token>
    std::unique_ptr<Token> tokenize(transaction_manager& mgr, const Mpc_units::energy& k);
    
    //! tokenize a linear power spectrum
    //! generates a new transaction on the database; will fail if a transaction is in progress
    template <typename PkContainer>
    std::unique_ptr<linear_Pk_token> tokenize(transaction_manager& mgr, const FRW_model_token& model, const PkContainer& Pk_lin);

    // DATA STORAGE

  public:

    //! store a sample of some kind (the exact behaviour is determined by template specialization)
    //! generates a new transaction on the database; will fail if a transaction is in progress
    template <typename SampleType>
    void store(const FRW_model_token& model, const SampleType& sample);
    
    
    // DATA EXTRACTION
    
  protected:
    
    //! extract a sample of a z-dependent but not k-dependent quantity, of the type specified by
    //! the payload
    template <typename PayloadType>
    std::unique_ptr<PayloadType> find(transaction_manager& mgr, const FRW_model_token& model, const z_database& z_db);
    
    //! extract a sample of a power spectrum-like quantity that depends on k but not z
    template <typename PayloadType>
    std::unique_ptr<PayloadType> find(transaction_manager& mgr, const linear_Pk_token& token, const k_database& k_db);
    
    //! extract a sample of a loop integral-like quantity that is k-dependent, UV and IR cutoff-dependent
    //! but not z-dependent
    template <typename PayloadType>
    std::unique_ptr<PayloadType>
    find(transaction_manager& mgr, const FRW_model_token& model, const k_token& k, const linear_Pk_token& Pk,
         const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff);
    
    //! extract a sample of a P(k)-like quantity that is k-dependent, z-dependent,
    //! and IR/UV-cutoff dependent
    template <typename PayloadType>
    std::unique_ptr<PayloadType>
    find(transaction_manager& mgr, const FRW_model_token& model, const k_token& k, const z_token& z,
         const linear_Pk_token& Pk_init, const boost::optional<linear_Pk_token>& Pk_final,
         const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff);
    
    //! extract a quantity of a IR-resummation-scale dependent quantity
    template <typename PayloadType>
    std::unique_ptr<PayloadType>
    find(transaction_manager& mgr, const FRW_model_token& model, const linear_Pk_token& Pk,
         const IR_resum_token& IR_resum);

    
    // TRANSACTIONS

  protected:

    //! open a transaction; throws an exception if a transaction is already held open
    //! note we have to use a std::shared_ptr<> here, rather than a std::unique_ptr<>,
    //! because we hold a std::weak_ptr<> internally to keep track of whether a transaction is open
    std::shared_ptr<transaction_manager> open_transaction();

    //! begin a new transaction on the database
    void begin_transaction();

    //! commit a transaction on the database
    void commit_transaction();

    //! rollback a transaction on the database
    void rollback_transaction();

    //! release a transaction
    void release_transaction();


    // LOOKUP OR INSERT RECORDS

  protected:

    //! lookup or insert a new FRW model
    unsigned int lookup_or_insert(transaction_manager& mgr, const FRW_model &obj);

    //! lookup or insert a redshift
    unsigned int lookup_or_insert(transaction_manager& mgr, double z);

    //! lookup or insert a wavenumber
    template <typename Token>
    unsigned int lookup_or_insert(transaction_manager& mgr, const Mpc_units::energy &k);
    
    //! lookup or insert a linear power spectrum identifier
    template <typename PkContainer>
    unsigned int lookup_or_insert(transaction_manager& mgr, const FRW_model_token& model, const PkContainer& Pk_lin);
    
    
    // INTERNAL DATA

  private:

    //! path to data container
    boost::filesystem::path container;

    //! SQLite3 handler for container
    sqlite3* handle;


    // TRANSACTIONS

    //! current transaction manager, if one exists
    std::weak_ptr<transaction_manager> current_transaction;


    // SQLite3 policies

    //! sqlite3_policy object
    sqlite3_policy policy;


    // SEARCH TOLERANCES

    //! tolerance to use when searching for FRW model parameters
    double FRW_model_tol;

    //! tolerance to use when searching for redshift configurations
    double z_tol;

    //! tolerance to use when searching for wavenumber configurations
    double k_tol;

  };


template <typename SampleType>
void data_manager::store(const FRW_model_token& model, const SampleType& sample)
  {
    // open a transaction on the database
    std::shared_ptr<transaction_manager> transaction = this->open_transaction();

    sqlite3_operations::store(this->handle, *transaction, this->policy, model, sample);

    // commit the transaction
    transaction->commit();
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


template <typename PkContainer>
std::unique_ptr<k_database> data_manager::build_k_db(transaction_manager& mgr, const PkContainer& Pk_lin,
                                                     double bottom_clearance, double top_clearance)
  {
    // construct an empty wavenumber database
    std::unique_ptr<k_database> k_db = std::make_unique<k_database>();
    
    // get power spectrum database underlying this container
    const tree_Pk::database_type& Pk_db = Pk_lin.get_db();
    
    for(tree_Pk::database_type::const_record_iterator t = Pk_db.record_cbegin(); t != Pk_db.record_cend(); ++t)
      {
        // ask initial_Pk container whether this P(k) value is acceptable
        const Mpc_units::energy& k = t->get_wavenumber();
        if(Pk_lin.is_valid(k, bottom_clearance, top_clearance))
          {
            std::unique_ptr<k_token> tok = this->tokenize<k_token>(mgr, k);
            k_db->add_record(k, *tok);
          }
      }
    
    return k_db;
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


template <typename PkContainer>
std::unique_ptr<linear_Pk_token>
data_manager::tokenize(const FRW_model_token& model, const PkContainer& Pk_lin)
  {
    // open a new transaction on the database
    std::shared_ptr<transaction_manager> transaction = this->open_transaction();
    
    // lookup id for this power spectrum, or generate one if it doesn't already exist
    std::unique_ptr<linear_Pk_token> id = this->tokenize(*transaction, model, Pk_lin);
    
    // commit the transaction
    transaction->commit();
    
    return std::move(id);
  }


template <typename PkContainer>
std::unique_ptr<linear_Pk_token>
data_manager::tokenize(transaction_manager& mgr, const FRW_model_token& model, const PkContainer& Pk_lin)
  {
    // lookup id for this power spectrum, or generate one if it doesn't already exist
    unsigned int id = this->lookup_or_insert(mgr, model, Pk_lin);
    return std::make_unique<linear_Pk_token>(id);
  }


template <typename Token>
unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, const Mpc_units::energy& k)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_wavenumber<Token>(this->handle, mgr, k, this->policy, this->k_tol);
    if(id) return(*id);
    
    return sqlite3_operations::insert_wavenumber<Token>(this->handle, mgr, k, this->policy);
  }


template <typename PkContainer>
unsigned int data_manager::lookup_or_insert(transaction_manager& mgr, const FRW_model_token& model, const PkContainer& Pk_lin)
  {
    boost::optional<unsigned int> id = sqlite3_operations::lookup_Pk_linear(this->handle, mgr, model, Pk_lin, this->policy);
    if(id) return(*id);
    
    return sqlite3_operations::insert_Pk_linear(this->handle, mgr, model, Pk_lin, this->policy);
  }


template <typename PkContainer>
std::unique_ptr<typename PkContainer::filtered_Pk_type> data_manager::build_wiggle_Pk(const linear_Pk_token& token, const PkContainer& Pk_lin)
  {
    // open a transaction on the database
    std::shared_ptr<transaction_manager> mgr = this->open_transaction();
    
    // extract database of wavenumber configurations from linear power spectrum container
    std::unique_ptr<k_database> k_db = this->build_k_db(*mgr, Pk_lin, FILTER_PK_DEFAULT_BOTTOM_CLEARANCE, FILTER_PK_DEFAULT_TOP_CLEARANCE);
    
    // extract initial_filtered_Pk container
    std::unique_ptr<typename PkContainer::filtered_Pk_type> payload = this->find<typename PkContainer::filtered_Pk_type>(*mgr, token, *k_db);
    
    // close transaction
    mgr->commit();
    
    return std::move(payload);
  }


template <typename PkContainer>
PkContainer& data_manager::rescale_final_Pk(const FRW_model_token& model, PkContainer& Pk, const z_database& z_db)
  {
    // open a transaction on the database
    std::shared_ptr<transaction_manager> mgr = this->open_transaction();
    
    // extract growth functions for the redshift database
    std::unique_ptr<oneloop_growth> data = this->find<oneloop_growth>(*mgr, model, z_db);
    
    oneloop_value z_init = *data->begin();
    oneloop_value z_final = *(--data->end());
    
    double rescale = z_init.second.g / z_final.second.g;
    
    // rescaling for power spectrum goes like the square of the growth factor
    Pk.set_rescaling(rescale*rescale);
    
    // close transaction
    mgr->commit();
    
    return Pk;
  }


#endif //LSSEFT_DATA_MANAGER_H
