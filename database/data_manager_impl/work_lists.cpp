//
// Created by David Seery on 09/12/2016.
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

#include <iostream>
#include <sstream>
#include <assert.h>

#include <set>
#include <unordered_set>

#include "database/data_manager.h"
#include "database/data_manager_impl/types.h"

#include "sqlite3_detail/utilities.h"
#include "sqlite3_detail/operations.h"

#include "utilities/formatter.h"

#include "defaults.h"

#include "boost/timer/timer.hpp"


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


std::unique_ptr<loop_momentum_work_list>
data_manager::build_loop_momentum_work_list(FRW_model_token& model, k_database& k_db,
                                            IR_cutoff_database& IR_db, UV_cutoff_database& UV_db,
                                            std::shared_ptr<wiggle_Pk>& Pk)
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
    loop_configs missing =
      sqlite3_operations::missing_loop_integral_configurations(this->handle, *mgr, this->policy, model, Pk->get_token(),
                                                               required_configs);
    
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


std::unique_ptr<one_loop_Pk_work_list>
data_manager::build_one_loop_Pk_work_list(FRW_model_token& model, z_database& z_db, k_database& k_db,
                                          IR_cutoff_database& IR_db, UV_cutoff_database& UV_db,
                                          std::shared_ptr<wiggle_Pk>& Pk)
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
          sqlite3_operations::missing_one_loop_Pk_redshifts(this->handle, *mgr, this->policy, model, Pk->get_token(),
                                                            z_table, z_db, record);
        
        // schedule a task to compute any missing redshifts
        if(missing_zs)
          {
            std::shared_ptr<oneloop_growth> g = this->find<oneloop_growth>(*mgr, model, z_db);
            
            std::shared_ptr<loop_integral> l =
              this->find<loop_integral>(*mgr, model, record.k->get_token(), Pk->get_token(),
                                        record.IR_cutoff->get_token(), record.UV_cutoff->get_token());
            
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


std::unique_ptr<one_loop_resum_Pk_work_list>
data_manager::build_one_loop_resum_Pk_work_list(FRW_model_token& model, z_database& z_db, k_database& k_db,
                                                IR_cutoff_database& IR_cutoff_db, UV_cutoff_database& UV_cutoff_db,
                                                IR_resum_database& IR_resum_db, std::shared_ptr<wiggle_Pk>& Pk)
  {
    // start timer
    boost::timer::cpu_timer timer;
    
    // construct an empty work list
    std::unique_ptr<one_loop_resum_Pk_work_list> work_list = std::make_unique<one_loop_resum_Pk_work_list>();
    
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
          sqlite3_operations::missing_one_loop_resum_Pk_redshifts(this->handle, *mgr, this->policy, model,
                                                                  Pk->get_token(), z_table, z_db, record);
        
        // schedule a task to compute any missing redshifts
        if(missing_zs)
          {
            std::unique_ptr<oneloop_growth> gf_data = this->find<oneloop_growth>(*mgr, model, *missing_zs);
            
            for(oneloop_growth::const_iterator t = gf_data->cbegin(); t != gf_data->cend(); ++t)
              {
                // lookup one-loop data for this redshift and loop configuration
                std::shared_ptr<oneloop_Pk> loop_data =
                  this->find<oneloop_Pk>(*mgr, model, record.k->get_token(), (*t).first, Pk->get_token(),
                                         record.IR_cutoff->get_token(), record.UV_cutoff->get_token());
                
                // lookup Matsubara X & Y coefficients for this IR resummation scale
                std::unique_ptr<Matsubara_XY> XY_coeffs =
                  this->find<Matsubara_XY>(*mgr, model, Pk->get_token(), record.IR_resum->get_token());
                
                work_list->emplace_back(*(*record.k), *XY_coeffs, loop_data, (*t).second, Pk);
              }
          }
      }
    
    // drop unneeded temporary tables
    sqlite3_operations::drop_temp(this->handle, *mgr, z_table);
    
    // close transaction
    mgr->commit();
    
    timer.stop();
    std::cout << "lsseft: constructed one-loop resummed P(k) work list (" << work_list->size() << " items) in time " << format_time(timer.elapsed().wall) << '\n';
    
    // release list if it contains no work
    if(work_list->empty()) work_list.release();
    
    return work_list;
  }


std::unique_ptr<multipole_Pk_work_list>
data_manager::build_multipole_Pk_work_list(FRW_model_token& model, z_database& z_db, k_database& k_db,
                                           IR_cutoff_database& IR_cutoff_db, UV_cutoff_database& UV_cutoff_db,
                                           IR_resum_database& IR_resum_db, std::shared_ptr<wiggle_Pk>& Pk)
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
          sqlite3_operations::missing_multipole_Pk_redshifts(this->handle, *mgr, this->policy, model, Pk->get_token(),
                                                             z_table, z_db, record);
        
        // schedule a task to compute any missing redshifts
        if(missing_zs)
          {
            std::unique_ptr<oneloop_growth> gf_data = this->find<oneloop_growth>(*mgr, model, *missing_zs);
            
            for(oneloop_growth::const_iterator t = gf_data->cbegin(); t != gf_data->cend(); ++t)
              {
                // lookup one-loop data for this redshift and loop configuration
                std::shared_ptr<oneloop_Pk> loop_data =
                  this->find<oneloop_Pk>(*mgr, model, record.k->get_token(), (*t).first, Pk->get_token(),
                                         record.IR_cutoff->get_token(), record.UV_cutoff->get_token());
                
                // lookup Matsubara X & Y coefficients for this IR resummation scale
                std::unique_ptr<Matsubara_XY> XY_coeffs =
                  this->find<Matsubara_XY>(*mgr, model, Pk->get_token(), record.IR_resum->get_token());
                
                work_list->emplace_back(*(*record.k), *XY_coeffs, loop_data, (*t).second, Pk);
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


std::unique_ptr<Matsubara_XY_work_list>
data_manager::build_Matsubara_XY_work_list(FRW_model_token& model, IR_resum_database& IR_resum_db,
                                           std::shared_ptr<wiggle_Pk>& Pk)
  {
    // start timer
    boost::timer::cpu_timer timer;
    
    // construct an empty work list
    std::unique_ptr<Matsubara_XY_work_list> work_list = std::make_unique<Matsubara_XY_work_list>();
    
    // open a transaction on the database
    std::shared_ptr<transaction_manager> mgr = this->open_transaction();
    
    // obatain list of missing configurations
    Matsubara_configs missing =
      sqlite3_operations::missing_Matsubara_XY_configurations(this->handle, *mgr, this->policy, model, Pk->get_token(),
                                                              IR_resum_db);
    
    // add these configurations to the work list
    for(const Matsubara_configs::value_type& record : missing)
      {
        work_list->emplace_back(*(*record.IR_resum), record.IR_resum->get_token(), Pk);
      }
    
    // close transaction
    mgr->commit();
    
    timer.stop();
    std::cout << "lsseft: constructed Matsubara XY work list (" << work_list->size() << " items) in time " << format_time(timer.elapsed().wall) << '\n';
    
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
    std::unique_ptr<k_database> k_db = this->build_k_db(*mgr, *Pk_lin, FILTER_PK_DEFAULT_BOTTOM_CLEARANCE, FILTER_PK_DEFAULT_TOP_CLEARANCE);
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
    std::cout << "lsseft: constructed wiggle/no-wiggle filter work list (" << work_list->size() << " items) in time " << format_time(timer.elapsed().wall) << '\n';
    
    // release list if it contains to work
    if(work_list->empty()) work_list.release();
    
    return work_list;
  }

