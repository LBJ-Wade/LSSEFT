//
// Created by David Seery on 13/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <list>
#include <algorithm>
#include <assert.h>
#include <set>
#include <unordered_set>

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
    std::set<unsigned int> missing_redshifts_for_table(sqlite3* db, const FRW_model_token& model,
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
          << "(SELECT zid FROM " << table << " WHERE " << table << ".mid=@mid) "
          << "ORDER BY id;";

        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));

        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));

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


    //! find missing wavenumbers in a named table, returned as a std::set<> of ids
    std::set<unsigned int> missing_wavenumbers_for_table(sqlite3* db, const FRW_model_token& model,
                                                         const std::string table, const std::string k_table)
      {
        assert(db != nullptr);

        // set up SQL statement to create table of missing k-numbers
        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id FROM " << k_table << " "
          << "WHERE id NOT IN "
          << "(SELECT kid FROM " << table << " WHERE " << table << ".mid=@mid) "
          << "ORDER BY id;";

        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));

        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));

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
        std::set<unsigned int> missing_list = missing_redshifts_for_table(db, model, k, policy.transfer_table(), z_table);

        // if any elements are missing, push them into a database
        if(missing_list.size() > 0)
          {
            missing_db = std::make_unique<z_database>();

            for(unsigned int t : missing_list)
              {
                // lookup record for this identifier
                z_database::const_record_iterator rec = z_db.lookup(z_token(t));

                // add a corresponding record to the missing database
                missing_db->add_record(*(*rec), rec->get_token());
              }
          }

        return(missing_db);
      }


    //! find missing redshifts for a one-loop growth function sample
    //! the required redshifts are passed in as the z_database z_db
    //! the return value is a database of redshifts for which values need to be computed
    std::unique_ptr<z_database> missing_oneloop_growth_redshifts(sqlite3* db, transaction_manager& mgr,
                                                                 const sqlite3_policy& policy,
                                                                 const FRW_model_token& model,
                                                                 const z_database& z_db, const std::string& z_table)
      {
        assert(db != nullptr);

        // set up null pointer; will be attached to an empty database later if needed
        std::unique_ptr<z_database> missing_db;

        // get list of missing z-values
        std::set<unsigned int> missing_list = missing_redshifts_for_table(db, model, policy.oneloop_table(), z_table);

        // if any elements are missing, push them into a database
        if(missing_list.size() > 0)
          {
            missing_db = std::make_unique<z_database>();

            for(unsigned int t : missing_list)
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
    std::unordered_set<data_manager_impl::loop_momentum_configuration>
    missing_loop_integral_configurations_for_table(sqlite3* db, const FRW_model_token& model, const std::string& table,
                                                   const std::unordered_set<data_manager_impl::loop_momentum_configuration>& required_configs)
      {
        assert(db != nullptr);
    
        std::unordered_set<data_manager_impl::loop_momentum_configuration> missing_configs;
    
        std::ostringstream count_stmt;
        count_stmt
          << "SELECT COUNT(*) FROM (SELECT * FROM " << table << " "
          << "WHERE mid=@mid AND kid=@kid AND UV_id=@UV_id AND IR_id=@IR_id);";
    
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, count_stmt.str().c_str(), count_stmt.str().length()+1, &stmt, nullptr));
        
        for(const loop_configs::value_type& t : required_configs)
          {
            // bind parameter values to search for an entry with this combination of model id, k id, UV limit id, IR limit id
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), t.k->get_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), t.UV->get_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), t.IR->get_token().get_id()));

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
    
    
    loop_configs update_missing_loop_integral_configurations(sqlite3* db, const FRW_model_token& model,
                                                             const std::string& table, const loop_configs& required_configs,
                                                             loop_configs& total_missing)
      {
        loop_configs missing = missing_loop_integral_configurations_for_table(db, model, table, required_configs);
        
        total_missing.insert(missing.begin(), missing.end());
        
        return missing;
      }
    
    
    void drop_inconsistent_configurations(sqlite3* db, const FRW_model_token& model, const std::string& table,
                                          const loop_configs& missing, const loop_configs& total_missing)
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
              << "DELETE FROM " << table << " WHERE mid=@mid AND kid=@kid AND UV_id=@UV_id AND IR_id=@IR_id;";
            
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, drop_stmt.str().c_str(), drop_stmt.str().length()+1, &stmt, nullptr));
            
            for(const loop_configs::value_type& t : inconsistent_set)
              {
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), t.k->get_token().get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), t.UV->get_token().get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), t.IR->get_token().get_id()));
                
                check_stmt(db, sqlite3_step(stmt));
                
                // release bindings and reset statement
                check_stmt(db, sqlite3_clear_bindings(stmt));
                check_stmt(db, sqlite3_reset(stmt));
              }
            
            // finalize statement to release resources
            check_stmt(db, sqlite3_finalize(stmt));
          }
      }

    
    //! find missing wavenumber/UV limit/IR limit configurations for a loop integral sample
    loop_configs missing_loop_integral_configurations(sqlite3* db, transaction_manager& mgr,
                                                      const sqlite3_policy& policy, const FRW_model_token& model,
                                                      const loop_configs& required_configs)
      {
        loop_configs total_missing;

        // find configurations missing from each table, and merge into global set of missing configurations
        loop_configs missing_AA = update_missing_loop_integral_configurations(db, model, policy.AA_table(), required_configs, total_missing);
        loop_configs missing_AB = update_missing_loop_integral_configurations(db, model, policy.AB_table(), required_configs, total_missing);
        loop_configs missing_BB = update_missing_loop_integral_configurations(db, model, policy.BB_table(), required_configs, total_missing);
        loop_configs missing_D = update_missing_loop_integral_configurations(db, model, policy.D_table(), required_configs, total_missing);
        loop_configs missing_E = update_missing_loop_integral_configurations(db, model, policy.E_table(), required_configs, total_missing);
        loop_configs missing_F = update_missing_loop_integral_configurations(db, model, policy.F_table(), required_configs, total_missing);
        loop_configs missing_G = update_missing_loop_integral_configurations(db, model, policy.G_table(), required_configs, total_missing);
        loop_configs missing_J1 = update_missing_loop_integral_configurations(db, model, policy.J1_table(), required_configs, total_missing);
        loop_configs missing_J2 = update_missing_loop_integral_configurations(db, model, policy.J2_table(), required_configs, total_missing);
        
        // drop any inconsistent configurations -- ie. those present in some tables, but not others
        drop_inconsistent_configurations(db, model, policy.AA_table(), missing_AA, total_missing);
        drop_inconsistent_configurations(db, model, policy.AB_table(), missing_AB, total_missing);
        drop_inconsistent_configurations(db, model, policy.BB_table(), missing_BB, total_missing);
        drop_inconsistent_configurations(db, model, policy.D_table(), missing_D, total_missing);
        drop_inconsistent_configurations(db, model, policy.E_table(), missing_E, total_missing);
        drop_inconsistent_configurations(db, model, policy.F_table(), missing_F, total_missing);
        drop_inconsistent_configurations(db, model, policy.G_table(), missing_G, total_missing);
        drop_inconsistent_configurations(db, model, policy.J1_table(), missing_J1, total_missing);
        drop_inconsistent_configurations(db, model, policy.J2_table(), missing_J2, total_missing);

        return total_missing;
      }


  }   // namespace sqlite3_operations
