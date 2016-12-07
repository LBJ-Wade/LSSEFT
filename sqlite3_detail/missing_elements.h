//
// Created by David Seery on 13/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_MISSING_ELEMENTS_H
#define LSSEFT_SQLITE3_MISSING_ELEMENTS_H


#include <set>
#include <unordered_set>

#include "sqlite3_policy.h"

#include "database/transaction_manager.h"
#include "database/tokens.h"
#include "database/z_database.h"
#include "database/k_database.h"
#include "database/data_manager_impl.h"

#include "sqlite3.h"


typedef std::unordered_set<data_manager_impl::loop_momentum_configuration> loop_configs;
typedef std::unordered_set<data_manager_impl::Matsubara_XY_configuration> Matsubara_configs;
typedef std::unordered_set<data_manager_impl::resummed_Pk_configuration> resum_Pk_configs;

namespace sqlite3_operations
  {

    //! construct a database of redshifts which need to be computed for the transfer function at a given wavenumber
    //! ownership of the resulting database is transferred via std::unique_ptr<>
    std::unique_ptr<z_database>
    missing_transfer_redshifts(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                               const FRW_model_token& model, const k_token& k, const z_database& z_db,
                               const std::string& z_table);


    //! construct a database of redshifts which need to be computed for the one-loop growth functions
    //! ownership of the resulting database is transferred via std::unique_ptr<>
    std::unique_ptr<z_database>
    missing_oneloop_growth_redshifts(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                     const FRW_model_token& model, const z_database& z_db, const std::string& z_table);


    //! process a list of configurations for loop momentum integrals;
    //! we detect which ones are already present in the database and avoid computing them
    loop_configs
    missing_loop_integral_configurations(sqlite3* db, transaction_manager& mgr,
                                         const sqlite3_policy& policy, const FRW_model_token& model,
                                         const loop_configs& required_configs);
    
    
    //! process a list of configurations for one-loop P(k) calculations;
    //! we detect which ones are already present in the database and avoid computing them
    std::unique_ptr<z_database>
    missing_one_loop_Pk_redshifts(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                  const FRW_model_token& model, const std::string& z_table, const z_database& z_db,
                                  const loop_configs::value_type& record);
    
    //! processs a list of configurations for the Matsubara-A coefficient;
    //! we detect which ones are already present in the database and avoid computing them
    Matsubara_configs
    missing_Matsubara_A_configurations(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                           const FRW_model_token& model, const linear_Pk_token& Pk,
                                           const IR_resum_database& IR_db);
    
    //! process a list of configurations for one-loop multipole P(k) calculations;
    //! we detect which ones are already present in the database and avoid computing them
    std::unique_ptr<z_database>
    missing_multipole_Pk_redshifts(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                   const FRW_model_token& model, const std::string& z_table, const z_database& z_db,
                                   const resum_Pk_configs::value_type& record);
    
    //! process a list of configurations for filtering the wiggle & no-wiggle part of Pk
    std::unique_ptr<k_database>
    missing_filter_Pk_wavenumbers(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                  const linear_Pk_token& token, const k_database& k_db, const std::string& k_table);

  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_MISSING_ELEMENTS_H
