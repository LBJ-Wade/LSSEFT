//
// Created by David Seery on 13/08/2015.
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
#include <list>
#include <algorithm>
#include <assert.h>
#include <set>
#include <unordered_set>
#include <cosmology/types.h>

#include "utilities.h"
#include "missing_elements.h"


namespace sqlite3_operations
  {

    //! find missing redshifts for a k-dependent table, returned as a std::set<> of ids
    std::set<unsigned int> missing_redshifts_for_table(sqlite3* db, const FRW_model_token& model, const k_token& k,
                                                       const std::string table, const std::string z_table)
      {
        assert(db != nullptr);

        // set up SQL statement to create table of missing z-numbers
        // note order by statement, which is important; the results from this function may be used
        // in std::list::merge(), which assumes the lists to be sorted
        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id FROM " << z_table << " "
          << "WHERE id NOT IN "
          << "(SELECT zid FROM " << table << " WHERE " << table << ".mid=@mid AND " << table << ".kid=@kid) "
          << "ORDER BY id;";

        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));

        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), k.get_id()));

        std::set<unsigned int> results;

        int status = 0;
        while((status = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(status == SQLITE_ROW)
              {
                results.insert(static_cast<unsigned int>(sqlite3_column_int(stmt, 0)));
              }
          }

        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));

        return(results);
      }


    //! find missing redshifts for a named k-independent table, returned as a std::set<> of ids
    std::set<unsigned int> missing_redshifts_for_table(sqlite3* db, const FRW_model_token& model, const growth_params_token& params,
                                                           const std::string& table, const std::string& z_table)
      {
        assert(db != nullptr);

        // set up SQL statement to create table of missing z-numbers
        // note order by statement, which is important; the results from this function may be used
        // in std::list::merge(), which assumes the lists to be sorted
        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id FROM " << z_table << " "
          << "WHERE id NOT IN "
          << "(SELECT zid FROM " << table << " WHERE mid=@mid AND params_id=@params_id) "
          << "ORDER BY id;";

        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));

        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@params_id"), params.get_id()));

        std::set<unsigned int> results;

        int status = 0;
        while((status = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(status == SQLITE_ROW)
              {
                results.insert(static_cast<unsigned int>(sqlite3_column_int(stmt, 0)));
              }
          }

        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));

        return(results);
      }
    
    
    //! find missing redshifts for a named loop-k-dependent table, returned as a std::set<> of ints
    std::set<unsigned int>
    missing_redshifts_for_table(sqlite3* db, const FRW_model_token& model, const growth_params_token& growth_params,
                                const loop_integral_params_token& loop_params, const linear_Pk_token& init_Pk,
                                const boost::optional<linear_Pk_token>& final_Pk, const std::string& table,
                                const std::string& z_table, const loop_configs::value_type& record)
      {
        assert(db != nullptr);
        
        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id FROM " << z_table << " "
          << "WHERE id NOT IN "
          << "(SELECT zid FROM " << table << " WHERE mid=@mid AND growth_params=@growth_params AND loop_params=@loop_params "
          << "AND kid=@kid AND init_Pk_id=@init_Pk_id "
          << "AND ((@final_Pk_id IS NULL AND final_Pk_id IS NULL) OR final_Pk_id=@final_Pk_id) AND IR_id=@IR_id AND UV_id=@UV_id) "
          << "ORDER BY id;";
        
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));
    
        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@growth_params"), growth_params.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@loop_params"), loop_params.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), record.k->get_token().get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@init_Pk_id"), init_Pk.get_id()));
        if(final_Pk)
          {
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@final_Pk_id"), final_Pk->get_id()));
          }
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), record.IR_cutoff->get_token().get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), record.UV_cutoff->get_token().get_id()));
        
        std::set<unsigned int> results;
    
        int status = 0;
        while((status = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(status == SQLITE_ROW)
              {
                results.insert(static_cast<unsigned int>(sqlite3_column_int(stmt, 0)));
              }
          }
    
        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));
    
        return(results);
      }
    
    
    //! find missing redshifts for a named loop-k-dependent table, returned as a std::set<> of ints
    std::set<unsigned int>
    missing_redshifts_for_table(sqlite3* db, const FRW_model_token& model, const growth_params_token& growth_params,
                                const loop_integral_params_token& loop_params,
                                const MatsubaraXY_params_token& XY_params,
                                const linear_Pk_token& init_Pk, const boost::optional<linear_Pk_token>& final_Pk,
                                const std::string& table, const std::string& z_table,
                                const resum_Pk_configs::value_type& record)
      {
        assert(db != nullptr);
        
        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id FROM " << z_table << " "
          << "WHERE id NOT IN "
          << "(SELECT zid FROM " << table << " WHERE mid=@mid AND growth_params=@growth_params AND loop_params=@loop_params AND XY_params=@XY_params "
          << "AND kid=@kid AND init_Pk_id=@init_Pk_id AND ((@final_Pk_id IS NULL AND final_Pk_id IS NULL) OR final_Pk_id=@final_Pk_id) "
          << "AND IR_cutoff_id=@IR_cutoff_id AND UV_cutoff_id=@UV_cutoff_id AND IR_resum_id=@IR_resum_id) "
          << "ORDER BY id;";
        
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));
        
        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@growth_params"), growth_params.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@loop_params"), loop_params.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@XY_params"), XY_params.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), record.k->get_token().get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@init_Pk_id"), init_Pk.get_id()));
        if(final_Pk)
          {
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@final_Pk_id"), final_Pk->get_id()));
          }
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_cutoff_id"), record.IR_cutoff->get_token().get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_cutoff_id"), record.UV_cutoff->get_token().get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_resum_id"), record.IR_resum->get_token().get_id()));
        
        std::set<unsigned int> results;
        
        int status = 0;
        while((status = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(status == SQLITE_ROW)
              {
                results.insert(static_cast<unsigned int>(sqlite3_column_int(stmt, 0)));
              }
          }
        
        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));
        
        return(results);
      }


    //! find missing redshifts for a named loop-k-dependent table, returned as a std::set<> of ints
    std::set<unsigned int>
    missing_redshifts_for_table(sqlite3* db, const FRW_model_token& model, const growth_params_token& growth_params,
                                const MatsubaraXY_params_token& XY_params,
                                const linear_Pk_token& init_Pk, const boost::optional<linear_Pk_token>& final_Pk,
                                const std::string& table, const std::string& z_table,
                                const resum_Pk_configs::value_type& record)
      {
        assert(db != nullptr);

        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id FROM " << z_table << " "
          << "WHERE id NOT IN "
          << "(SELECT zid FROM " << table << " WHERE mid=@mid AND growth_params=@growth_params AND XY_params=@XY_params "
          << "AND kid=@kid AND init_Pk_id=@init_Pk_id AND ((@final_Pk_id IS NULL AND final_Pk_id IS NULL) OR final_Pk_id=@final_Pk_id) "
          << "AND IR_cutoff_id=@IR_cutoff_id AND UV_cutoff_id=@UV_cutoff_id AND IR_resum_id=@IR_resum_id) "
          << "ORDER BY id;";

        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));

        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@growth_params"), growth_params.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@XY_params"), XY_params.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), record.k->get_token().get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@init_Pk_id"), init_Pk.get_id()));
        if(final_Pk)
          {
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@final_Pk_id"), final_Pk->get_id()));
          }
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_cutoff_id"), record.IR_cutoff->get_token().get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_cutoff_id"), record.UV_cutoff->get_token().get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_resum_id"), record.IR_resum->get_token().get_id()));

        std::set<unsigned int> results;

        int status = 0;
        while((status = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(status == SQLITE_ROW)
              {
                results.insert(static_cast<unsigned int>(sqlite3_column_int(stmt, 0)));
              }
          }

        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));

        return(results);
      }


    //! find missing wavenumbers for a filtered power spectrum
    std::set<unsigned int> missing_filtered_wavenumbers(sqlite3* db, const linear_Pk_token& Pk_token,
                                                        const filter_params_token& params_token,
                                                        const std::string& table, const std::string& k_table)
      {
        assert(db != nullptr);

        // set up SQL statement to create table of missing k-numbers
        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id FROM " << k_table << " "
          << "WHERE id NOT IN "
          << "(SELECT kid FROM " << table << " WHERE " << table << ".Pk_id=@Pk_id AND " << table << ".params_id=@params_id) "
          << "ORDER BY id;";

        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));

        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), Pk_token.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@params_id"), params_token.get_id()));

        std::set<unsigned int> results;

        int status = 0;
        while((status = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(status == SQLITE_ROW)
              {
                results.insert(static_cast<unsigned int>(sqlite3_column_int(stmt, 0)));
              }
          }

        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));

        return(results);
      }


    //! find missing redshifts for a transfer function sample
    //! the required redshifts are passed in as the z_database z_db
    //! the return value is a database of redshifts for which values need to be computed
    std::unique_ptr<z_database> missing_transfer_redshifts(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                                           const FRW_model_token& model, const k_token& k,
                                                           const z_database& z_db, const std::string& z_table)
      {
        assert(db != nullptr);

        // set up null pointer; will be attached to an empty database later if needed
        std::unique_ptr<z_database> missing_db;

        // get list of missing z-values for this k-mode of the transfer function
        std::set<unsigned int> missing = missing_redshifts_for_table(db, model, k, policy.transfer_table(), z_table);

        // if any elements are missing, push them into a database
        if(missing.size() > 0)
          {
            missing_db = std::make_unique<z_database>();

            for(unsigned int t : missing)
              {
                // lookup record for this identifier
                z_database::const_record_iterator rec = z_db.lookup(z_token(t));

                // add a corresponding record to the missing database
                missing_db->add_record(*(*rec), rec->get_token());
              }
          }

        return(missing_db);
      }
    
    
    std::set<unsigned int>
    update_missing_oneloop_growth_redshifts(sqlite3* db, const FRW_model_token& model, const growth_params_token& params,
                                            const std::string& table, const std::string& z_table,
                                            std::set<unsigned int>& total_missing)
      {
        std::set<unsigned int> missing = missing_redshifts_for_table(db, model, params, table, z_table);
        
        if(!missing.empty()) total_missing.insert(missing.begin(), missing.end());
        
        return missing;
      }
    
    
    std::set<unsigned int>
    update_missing_one_loop_Pk(sqlite3* db, const FRW_model_token& model, const growth_params_token& growth_params,
                               const loop_integral_params_token& loop_params, const linear_Pk_token& init_Pk,
                               const boost::optional<linear_Pk_token>& final_Pk, const std::string& table,
                               const std::string& z_table, const loop_configs::value_type& record,
                               std::set<unsigned int>& total_missing)
      {
        std::set<unsigned int> missing = missing_redshifts_for_table(db, model, growth_params, loop_params, init_Pk, final_Pk, table, z_table, record);
        
        if(!missing.empty()) total_missing.insert(missing.begin(), missing.end());
        
        return missing;
      }
    
    
    std::set<unsigned int>
    update_missing_multipole_Pk(sqlite3* db, const FRW_model_token& model, const growth_params_token& growth_params,
                                const loop_integral_params_token& loop_params,
                                const MatsubaraXY_params_token& XY_params,
                                const linear_Pk_token& init_Pk, const boost::optional<linear_Pk_token>& final_Pk,
                                const std::string& table, const std::string& z_table,
                                const resum_Pk_configs::value_type& record, std::set<unsigned int>& total_missing)
      {
        std::set<unsigned int> missing = missing_redshifts_for_table(db, model, growth_params, loop_params, XY_params,
                                                                     init_Pk, final_Pk, table, z_table, record);
    
        total_missing.insert(missing.begin(), missing.end());
    
        return missing;
      }


    std::set<unsigned int>
    update_missing_counterterms(sqlite3* db, const FRW_model_token& model, const growth_params_token& growth_params,
                                const MatsubaraXY_params_token& XY_params,
                                const linear_Pk_token& init_Pk, const boost::optional<linear_Pk_token>& final_Pk,
                                const std::string& table, const std::string& z_table,
                                const resum_Pk_configs::value_type& record, std::set<unsigned int>& total_missing)
      {
        std::set<unsigned int> missing = missing_redshifts_for_table(db, model, growth_params, XY_params, init_Pk,
                                                                     final_Pk, table, z_table, record);

        total_missing.insert(missing.begin(), missing.end());

        return missing;
      }


    void drop_inconsistent_redshifts(sqlite3* db, const FRW_model_token& model, const growth_params_token& params,
                                     const std::string& table, const std::set<unsigned int>& missing,
                                     const std::set<unsigned int>& total_missing)
      {
        std::set<unsigned int> inconsistent_set;
        
        std::set_difference(total_missing.begin(), total_missing.end(),
                            missing.begin(), missing.end(), std::inserter(inconsistent_set, inconsistent_set.begin()));
        
        if(inconsistent_set.size() > 0)
          {
            std::ostringstream drop_stmt;
            drop_stmt
              << "DELETE FROM " << table << " WHERE mid=@mid AND params_id=@params_id AND zid=@zid;";
    
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, drop_stmt.str().c_str(), drop_stmt.str().length()+1, &stmt, nullptr));
            
            for(unsigned int t : inconsistent_set)
              {
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@params_id"), params.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), t));
                
                check_stmt(db, sqlite3_step(stmt), SQLITE_DONE);
    
                // release bindings and reset statement
                check_stmt(db, sqlite3_clear_bindings(stmt));
                check_stmt(db, sqlite3_reset(stmt));
              }
    
            // finalize statement to release resources
            check_stmt(db, sqlite3_finalize(stmt));
          }
      }


    void
    drop_inconsistent_redshifts(sqlite3* db, const FRW_model_token& model, const linear_Pk_token& init_Pk,
                                const growth_params_token& growth_params,
                                const loop_integral_params_token& loop_params,
                                const boost::optional<linear_Pk_token>& final_Pk, const std::string& table,
                                const loop_configs::value_type& record, const std::set<unsigned int>& missing,
                                const std::set<unsigned int>& total_missing)
      {
        std::set<unsigned int> inconsistent_set;
    
        std::set_difference(total_missing.begin(), total_missing.end(),
                            missing.begin(), missing.end(), std::inserter(inconsistent_set, inconsistent_set.begin()));
    
        if(inconsistent_set.size() > 0)
          {
            std::ostringstream drop_stmt;
            drop_stmt
              << "DELETE FROM " << table << " WHERE mid=@mid AND growth_params=@growth_params AND loop_params=@loop_params "
              << "AND zid=@zid AND kid=@kid AND init_Pk_id=@init_Pk_id "
              << "AND ((@final_Pk_id IS NULL AND final_Pk_id IS NULL) OR final_Pk_id=@final_Pk_id) "
              << "AND IR_id=@IR_id AND UV_id=@UV_id;";
        
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, drop_stmt.str().c_str(), drop_stmt.str().length()+1, &stmt, nullptr));
        
            for(unsigned int t : inconsistent_set)
              {
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@growth_params"), growth_params.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@loop_params"), loop_params.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), t));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), record.k->get_token().get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@init_Pk_id"), init_Pk.get_id()));
                if(final_Pk)
                  {
                    check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@final_Pk_id"), final_Pk->get_id()));
                  }
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), record.UV_cutoff->get_token().get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), record.IR_cutoff->get_token().get_id()));
            
                check_stmt(db, sqlite3_step(stmt), SQLITE_DONE);
            
                // release bindings and reset statement
                check_stmt(db, sqlite3_clear_bindings(stmt));
                check_stmt(db, sqlite3_reset(stmt));
              }
        
            // finalize statement to release resources
            check_stmt(db, sqlite3_finalize(stmt));
          }
      }
    
    
    void
    drop_inconsistent_redshifts(sqlite3* db, const FRW_model_token& model, const growth_params_token& growth_params,
                                const loop_integral_params_token& loop_params,
                                const MatsubaraXY_params_token& XY_params, const linear_Pk_token& init_Pk,
                                const boost::optional<linear_Pk_token>& final_Pk, const std::string& table,
                                const resum_Pk_configs::value_type& record, const std::set<unsigned int>& missing,
                                const std::set<unsigned int>& total_missing)
      {
        std::set<unsigned int> inconsistent_set;
        
        std::set_difference(total_missing.begin(), total_missing.end(),
                            missing.begin(), missing.end(), std::inserter(inconsistent_set, inconsistent_set.begin()));
        
        if(inconsistent_set.size() > 0)
          {
            std::ostringstream drop_stmt;
            drop_stmt
              << "DELETE FROM " << table << " WHERE mid=@mid AND growth_params=@growth_params AND loop_params=@loop_params AND XY_params=@XY_params "
              << "AND zid=@zid AND kid=@kid AND init_Pk_id=@init_Pk_id AND ((@final_Pk_id IS NULL AND final_Pk_id IS NULL) OR final_Pk_id=@final_Pk_id) "
              << "AND IR_cutoff_id=@IR_cutoff_id AND UV_cutoff_id=@UV_cutoff_id AND IR_resum_id=@IR_resum_id;";
            
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, drop_stmt.str().c_str(), drop_stmt.str().length()+1, &stmt, nullptr));
            
            for(unsigned int t : inconsistent_set)
              {
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@growth_params"), growth_params.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@loop_params"), loop_params.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@XY_params"), XY_params.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), t));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), record.k->get_token().get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@init_Pk_id"), init_Pk.get_id()));
                if(final_Pk)
                  {
                    check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@final_Pk_id"), final_Pk->get_id()));
                  }
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_cutoff_id"), record.IR_cutoff->get_token().get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_cutoff_id"), record.UV_cutoff->get_token().get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_resum_id"), record.IR_resum->get_token().get_id()));
                
                check_stmt(db, sqlite3_step(stmt), SQLITE_DONE);
                
                // release bindings and reset statement
                check_stmt(db, sqlite3_clear_bindings(stmt));
                check_stmt(db, sqlite3_reset(stmt));
              }
            
            // finalize statement to release resources
            check_stmt(db, sqlite3_finalize(stmt));
          }
      }


    void
    drop_inconsistent_redshifts(sqlite3* db, const FRW_model_token& model, const growth_params_token& growth_params,
                                const MatsubaraXY_params_token& XY_params, const linear_Pk_token& init_Pk,
                                const boost::optional<linear_Pk_token>& final_Pk, const std::string& table,
                                const resum_Pk_configs::value_type& record, const std::set<unsigned int>& missing,
                                const std::set<unsigned int>& total_missing)
      {
        std::set<unsigned int> inconsistent_set;

        std::set_difference(total_missing.begin(), total_missing.end(),
                            missing.begin(), missing.end(), std::inserter(inconsistent_set, inconsistent_set.begin()));

        if(inconsistent_set.size() > 0)
          {
            std::ostringstream drop_stmt;
            drop_stmt
              << "DELETE FROM " << table << " WHERE mid=@mid AND growth_params=@growth_params AND XY_params=@XY_params "
              << "AND zid=@zid AND kid=@kid AND init_Pk_id=@init_Pk_id AND ((@final_Pk_id IS NULL AND final_Pk_id IS NULL) OR final_Pk_id=@final_Pk_id) "
              << "AND IR_cutoff_id=@IR_cutoff_id AND UV_cutoff_id=@UV_cutoff_id AND IR_resum_id=@IR_resum_id;";

            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, drop_stmt.str().c_str(), drop_stmt.str().length()+1, &stmt, nullptr));

            for(unsigned int t : inconsistent_set)
              {
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@growth_params"), growth_params.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@XY_params"), XY_params.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), t));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), record.k->get_token().get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@init_Pk_id"), init_Pk.get_id()));
                if(final_Pk)
                  {
                    check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@final_Pk_id"), final_Pk->get_id()));
                  }
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_cutoff_id"), record.IR_cutoff->get_token().get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_cutoff_id"), record.UV_cutoff->get_token().get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_resum_id"), record.IR_resum->get_token().get_id()));

                check_stmt(db, sqlite3_step(stmt), SQLITE_DONE);

                // release bindings and reset statement
                check_stmt(db, sqlite3_clear_bindings(stmt));
                check_stmt(db, sqlite3_reset(stmt));
              }

            // finalize statement to release resources
            check_stmt(db, sqlite3_finalize(stmt));
          }
      }

    
    //! find missing redshifts for a one-loop growth function sample
    //! the required redshifts are passed in as the z_database z_db
    //! the return value is a database of redshifts for which values need to be computed
    std::unique_ptr<z_database>
    missing_oneloop_growth_redshifts(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                         const FRW_model_token& model, const growth_params_token& params,
                                         const z_database& z_db, const std::string& z_table)
      {
        assert(db != nullptr);

        // set up null pointer; will be attached to an empty database later if needed
        std::unique_ptr<z_database> missing_db;

        // get list of missing z-values from each relevant table
        std::set<unsigned int> missing;
        std::set<unsigned int> missing_g = update_missing_oneloop_growth_redshifts(db, model, params,
                                                                                   policy.D_factor_table(), z_table, missing);
        std::set<unsigned int> missing_f = update_missing_oneloop_growth_redshifts(db, model, params, policy.f_factor_table(), z_table, missing);
        
        // drop any inconsistent redshifts -- ie. those present in some tables but not others
        drop_inconsistent_redshifts(db, model, params, policy.D_factor_table(), missing_g, missing);
        drop_inconsistent_redshifts(db, model, params, policy.f_factor_table(), missing_f, missing);

        // if any elements are missing, push them into a database
        if(!missing.empty())
          {
            missing_db = std::make_unique<z_database>();

            for(unsigned int t : missing)
              {
                // lookup record for this identifier
                z_database::const_record_iterator rec = z_db.lookup(z_token(t));

                // add a corresponding record to the missing database
                missing_db->add_record(*(*rec), rec->get_token());
              }
          }

        return(missing_db);
      }

    
    //! search for missing wavenumber/UV limit/IR limit configurations for a named table
    loop_configs
    missing_loop_integral_configurations_for_table(sqlite3* db, const FRW_model_token& model,
                                                   const loop_integral_params_token& params,
                                                   const linear_Pk_token& Pk_lin, const std::string& table,
                                                   const std::unordered_set<data_manager_impl::loop_momentum_configuration>& required_configs)
      {
        assert(db != nullptr);
    
        loop_configs missing_configs;
    
        std::ostringstream count_stmt;
        count_stmt
          << "SELECT COUNT(*) FROM (SELECT * FROM " << table << " "
          << "WHERE mid=@mid AND params_id=@params_id AND kid=@kid AND Pk_id=@Pk_id AND UV_id=@UV_id AND IR_id=@IR_id);";
    
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, count_stmt.str().c_str(), count_stmt.str().length()+1, &stmt, nullptr));
        
        for(const loop_configs::value_type& t : required_configs)
          {
            // bind parameter values to search for an entry with this combination of model id, k id, UV limit id, IR limit id
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@params_id"), params.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), t.k->get_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), Pk_lin.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), t.UV_cutoff->get_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), t.IR_cutoff->get_token().get_id()));

            int status = 0;
            while((status = sqlite3_step(stmt)) != SQLITE_DONE)
              {
                if(status == SQLITE_ROW)
                  {
                    if(sqlite3_column_int(stmt, 0) == 0)
                      {
                        missing_configs.insert(t);
                      }
                  }
              }
    
            // release bindings and reset statement
            check_stmt(db, sqlite3_clear_bindings(stmt));
            check_stmt(db, sqlite3_reset(stmt));
          }
    
        // finalize statement to release resources
        check_stmt(db, sqlite3_finalize(stmt));
        
        return missing_configs;
      }
    
    
    loop_configs
    update_missing_loop_integral_configurations(sqlite3* db, const FRW_model_token& model, const loop_integral_params_token& params,
                                                const linear_Pk_token& Pk_lin, const std::string& table,
                                                const loop_configs& required_configs, loop_configs& total_missing)
      {
        loop_configs missing = missing_loop_integral_configurations_for_table(db, model, params, Pk_lin, table, required_configs);
        
        total_missing.insert(missing.begin(), missing.end());
        
        return missing;
      }


    void drop_inconsistent_configurations(sqlite3* db, const FRW_model_token& model, const loop_integral_params_token& params,
                                          const linear_Pk_token& Pk_lin, const std::string& table, const loop_configs& missing,
                                          const loop_configs& total_missing)
      {
        loop_configs inconsistent_set;
    
        // note std::set_difference does not work for unordered_sets, at least up to C++17
        // so, we have to implement an equivalent function ourselves
        std::copy_if(total_missing.begin(), total_missing.end(),
                     std::inserter(inconsistent_set, inconsistent_set.begin()),
                     [&missing] (const loop_configs::value_type& v) -> bool { return missing.find(v) == missing.end(); });
        
        if(inconsistent_set.size() > 0)
          {
            std::ostringstream drop_stmt;
            drop_stmt
              << "DELETE FROM " << table << " WHERE mid=@mid AND params_id=@params_id AND kid=@kid AND Pk_id=@Pk_id AND UV_id=@UV_id AND IR_id=@IR_id;";
            
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, drop_stmt.str().c_str(), drop_stmt.str().length()+1, &stmt, nullptr));
            
            for(const loop_configs::value_type& t : inconsistent_set)
              {
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@params_id"), params.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), t.k->get_token().get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), Pk_lin.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), t.UV_cutoff->get_token().get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), t.IR_cutoff->get_token().get_id()));
                
                check_stmt(db, sqlite3_step(stmt), SQLITE_DONE);
                
                // release bindings and reset statement
                check_stmt(db, sqlite3_clear_bindings(stmt));
                check_stmt(db, sqlite3_reset(stmt));
              }
            
            // finalize statement to release resources
            check_stmt(db, sqlite3_finalize(stmt));
          }
      }

    
    //! find missing wavenumber/UV limit/IR limit configurations for a loop integral sample
    loop_configs missing_loop_integral_configurations(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                                      const FRW_model_token& model, const loop_integral_params_token& params,
                                                      const linear_Pk_token& Pk_lin, const loop_configs& required_configs)
      {
        loop_configs total_missing;

#include "autogenerated/missing_kernel_stmts.cpp"

        return total_missing;
      }
    
    
    std::unique_ptr<z_database>
    missing_one_loop_Pk_redshifts(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                  const FRW_model_token& model, const growth_params_token& growth_params,
                                  const loop_integral_params_token& loop_params, const linear_Pk_token& init_Pk,
                                  const boost::optional<linear_Pk_token>& final_Pk, const std::string& z_table,
                                  const z_database& z_db, const loop_configs::value_type& record)
      {
        assert(db != nullptr);
        
        // set up a null pointer for the returned database; will be attached to an empty instance later if needed
        std::unique_ptr<z_database> missing_db;
        
        // get list of missing z-values from each relevant table
        std::set<unsigned int> missing;

#include "autogenerated/missing_Pk_stmts.cpp"

        // push any missing elements into the result database
        if(!missing.empty())
          {
            missing_db = std::make_unique<z_database>();
            
            for(unsigned int t : missing)
              {
                // lookup record for this identifier
                z_database::const_record_iterator rec = z_db.lookup(z_token(t));
                
                // add a corresponding record to the missing database
                missing_db->add_record(*(*rec), rec->get_token());
              }
          }
        
        return missing_db;
      }


    std::unique_ptr<z_database>
    missing_multipole_Pk_redshifts(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                       const FRW_model_token& model, const growth_params_token& growth_params,
                                       const loop_integral_params_token& loop_params, const MatsubaraXY_params_token& XY_params,
                                       const linear_Pk_token& init_Pk, const boost::optional<linear_Pk_token>& final_Pk,
                                       const std::string& z_table, const z_database& z_db,
                                       const resum_Pk_configs::value_type& record)
      {
        assert(db != nullptr);
    
        // set up a null pointer for the returned database; will be attached to an empty instance later if needed
        std::unique_ptr<z_database> missing_db;
    
        std::set<unsigned int> missing;

#include "autogenerated/missing_multipole_stmts.cpp"

        // push any missing elements into the result database
        if(!missing.empty())
          {
            missing_db = std::make_unique<z_database>();
        
            for(unsigned int t : missing)
              {
                // lookup record for this identifier
                z_database::const_record_iterator rec = z_db.lookup(z_token(t));
            
                // add a corresponding record to the missing database
                missing_db->add_record(*(*rec), rec->get_token());
              }
          }
    
        return missing_db;
      }
    
    
    Matsubara_configs
    missing_Matsubara_XY_configurations(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                            const FRW_model_token& model, const linear_Pk_token& Pk,
                                            const IR_resum_database& IR_db, const MatsubaraXY_params_token& params_tok)
      {
        assert(db != nullptr);
        
        Matsubara_configs missing_configs;
        
        std::ostringstream count_stmt;
        count_stmt
          << "SELECT COUNT(*) FROM (SELECT * FROM " << policy.Matsubara_XY_table() << " "
          << "WHERE mid=@mid AND params_id=@params_id AND Pk_id=@Pk_id AND IR_resum_id=@IR_resum_id);";
    
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, count_stmt.str().c_str(), count_stmt.str().length()+1, &stmt, nullptr));
        
        for(IR_resum_database::const_record_iterator t = IR_db.record_cbegin(); t != IR_db.record_cend(); ++t)
          {
            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@params_id"), params_tok.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), Pk.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_resum_id"), t->get_token().get_id()));
    
            int status = 0;
            while((status = sqlite3_step(stmt)) != SQLITE_DONE)
              {
                if(status == SQLITE_ROW)
                  {
                    if(sqlite3_column_int(stmt, 0) == 0)
                      {
                        missing_configs.emplace(t);
                      }
                  }
              }
    
            // release bindings and reset statement
            check_stmt(db, sqlite3_clear_bindings(stmt));
            check_stmt(db, sqlite3_reset(stmt));
          }
    
        // finalize statement to release resources
        check_stmt(db, sqlite3_finalize(stmt));
    
        return missing_configs;
      }
    
    
    std::unique_ptr<k_database>
    missing_filter_Pk_wavenumbers(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                  const linear_Pk_token& Pk_token, const filter_params_token& params_token,
                                  const k_database& k_db, const std::string& k_table)
      {
        assert(db != nullptr);
        
        // set up null pointer; will be attached to an empty database later if needed
        std::unique_ptr<k_database> missing_db;
        
        // get list of missing k-values for this linear power spectrum
        std::set<unsigned int> missing = missing_filtered_wavenumbers(db, Pk_token, params_token, policy.Pk_linear_table(), k_table);
        
        if(missing.size() > 0)
          {
            missing_db = std::make_unique<k_database>();
            
            for(unsigned int t : missing)
              {
                // lookup record for this identifier
                k_database::const_record_iterator rec = k_db.lookup(k_token(t));
                
                // add a corresponding element to the missing database
                missing_db->add_record(*(*rec), rec->get_token());
              }
          }
        
        return missing_db;
      }


    std::unique_ptr<z_database>
    missing_counterterm_redshifts(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                  const FRW_model_token& model, const growth_params_token& growth_params,
                                  const MatsubaraXY_params_token& XY_params, const linear_Pk_token& init_Pk,
                                  const boost::optional<linear_Pk_token>& final_Pk, const std::string& z_table,
                                  const z_database& z_db, const resum_Pk_configs::value_type& record)
      {
        assert(db != nullptr);

        // set up a null pointer for the returned database; will be attached to an empty instance later if needed
        std::unique_ptr<z_database> missing_db;

        std::set<unsigned int> missing;

        // find missing configurations for each of P0, P2, P4 and merge
        std::set<unsigned int> missing_P0 =
          update_missing_counterterms(db, model, growth_params, XY_params, init_Pk, final_Pk,
                                      policy.counterterms_c0_table(), z_table, record, missing);
        std::set<unsigned int> missing_P2 =
          update_missing_counterterms(db, model, growth_params, XY_params, init_Pk, final_Pk,
                                      policy.counterterms_c0_table(), z_table, record, missing);
        std::set<unsigned int> missing_P4 =
          update_missing_counterterms(db, model, growth_params, XY_params, init_Pk, final_Pk,
                                      policy.counterterms_c0_table(), z_table, record, missing);

        // bring database to consistent state by dropping any conflicting results
        drop_inconsistent_redshifts(db, model, growth_params, XY_params, init_Pk, final_Pk,
                                    policy.counterterms_c0_table(), record, missing_P0, missing);
        drop_inconsistent_redshifts(db, model, growth_params, XY_params, init_Pk, final_Pk,
                                    policy.counterterms_c2_table(), record, missing_P2, missing);
        drop_inconsistent_redshifts(db, model, growth_params, XY_params, init_Pk, final_Pk,
                                    policy.counterterms_c4_table(), record, missing_P4, missing);

        // push any missing elements into the result database
        if(!missing.empty())
          {
            missing_db = std::make_unique<z_database>();

            for(unsigned int t : missing)
              {
                // lookup record for this identifier
                z_database::const_record_iterator rec = z_db.lookup(z_token(t));

                // add a corresponding record to the missing database
                missing_db->add_record(*(*rec), rec->get_token());
              }
          }

        return missing_db;
      }


  }   // namespace sqlite3_operations
