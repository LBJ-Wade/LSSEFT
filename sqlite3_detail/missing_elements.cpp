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
    std::set<unsigned int> missing_redshifts_for_table(sqlite3* db, const FRW_model_token& model,
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
          << "(SELECT zid FROM " << table << " WHERE mid=@mid) "
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
    
    
    //! find missing redshifts for a named loop-k-dependent table, returned as a std::set<> of ints
    std::set<unsigned int>
    missing_redshifts_for_table(sqlite3* db, const FRW_model_token& model, const linear_Pk_token& Pk_lin,
                                const std::string& table, const std::string& z_table,
                                const loop_configs::value_type& record)
      {
        assert(db != nullptr);
        
        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id FROM " << z_table << " "
          << "WHERE id NOT IN "
          << "(SELECT zid FROM " << table << " WHERE mid=@mid AND kid=@kid AND Pk_id=@Pk_id AND IR_id=@IR_id AND UV_id=@UV_id) "
          << "ORDER BY id;";
        
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));
    
        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), record.k->get_token().get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), Pk_lin.get_id()));
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
    missing_redshifts_for_table(sqlite3* db, const FRW_model_token& model, const linear_Pk_token& Pk_lin,
                                const std::string& table, const std::string& z_table,
                                const resum_Pk_configs::value_type& record)
      {
        assert(db != nullptr);
        
        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id FROM " << z_table << " "
          << "WHERE id NOT IN "
          << "(SELECT zid FROM " << table << " WHERE mid=@mid AND kid=@kid AND Pk_id=@Pk_id "
          << "AND IR_cutoff_id=@IR_cutoff_id AND UV_cutoff_id=@UV_cutoff_id AND IR_resum_id=@IR_resum_id) "
          << "ORDER BY id;";
        
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));
        
        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), record.k->get_token().get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), Pk_lin.get_id()));
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


    //! find missing wavenumbers in a named table, returned as a std::set<> of ids
    std::set<unsigned int> missing_wavenumbers_for_table(sqlite3* db, const linear_Pk_token& token,
                                                         const std::string& table, const std::string& k_table)
      {
        assert(db != nullptr);

        // set up SQL statement to create table of missing k-numbers
        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id FROM " << k_table << " "
          << "WHERE id NOT IN "
          << "(SELECT kid FROM " << table << " WHERE " << table << ".Pk_id=@Pk_id) "
          << "ORDER BY id;";

        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));

        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), token.get_id()));

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
    
    
    std::set<unsigned int> update_missing_oneloop_growth_redshifts(sqlite3* db, const FRW_model_token& model,
                                                                   const std::string& table, const std::string& z_table,
                                                                   std::set<unsigned int>& total_missing)
      {
        std::set<unsigned int> missing = missing_redshifts_for_table(db, model, table, z_table);
        
        if(!missing.empty()) total_missing.insert(missing.begin(), missing.end());
        
        return missing;
      }
    
    
    std::set<unsigned int>
    update_missing_one_loop_Pk(sqlite3* db, const FRW_model_token& model, const linear_Pk_token& Pk_lin,
                               const std::string& table, const std::string& z_table,
                               const loop_configs::value_type& record, std::set<unsigned int>& total_missing)
      {
        std::set<unsigned int> missing = missing_redshifts_for_table(db, model, Pk_lin, table, z_table, record);
        
        if(!missing.empty()) total_missing.insert(missing.begin(), missing.end());
        
        return missing;
      }
    
    
    std::set<unsigned int>
    update_missing_multipole_Pk(sqlite3* db, const FRW_model_token& model, const linear_Pk_token& Pk_lin,
                                const std::string& table, const std::string& z_table,
                                const resum_Pk_configs::value_type& record, std::set<unsigned int>& total_missing)
      {
        std::set<unsigned int> missing = missing_redshifts_for_table(db, model, Pk_lin, table, z_table, record);
    
        total_missing.insert(missing.begin(), missing.end());
    
        return missing;
      }
    
    
    void drop_inconsistent_redshifts(sqlite3* db, const FRW_model_token& model, const std::string& table,
                                     const std::set<unsigned int>& missing, const std::set<unsigned int>& total_missing)
      {
        std::set<unsigned int> inconsistent_set;
        
        std::set_difference(total_missing.begin(), total_missing.end(),
                            missing.begin(), missing.end(), std::inserter(inconsistent_set, inconsistent_set.begin()));
        
        if(inconsistent_set.size() > 0)
          {
            std::ostringstream drop_stmt;
            drop_stmt
              << "DELETE FROM " << table << " WHERE mid=@mid AND zid=@zid;";
    
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, drop_stmt.str().c_str(), drop_stmt.str().length()+1, &stmt, nullptr));
            
            for(unsigned int t : inconsistent_set)
              {
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
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
    
    
    void drop_inconsistent_redshifts(sqlite3* db, const FRW_model_token& model, const linear_Pk_token& Pk_lin,
                                     const std::string& table, const loop_configs::value_type& record,
                                     const std::set<unsigned int>& missing, const std::set<unsigned int>& total_missing)
      {
        std::set<unsigned int> inconsistent_set;
    
        std::set_difference(total_missing.begin(), total_missing.end(),
                            missing.begin(), missing.end(), std::inserter(inconsistent_set, inconsistent_set.begin()));
    
        if(inconsistent_set.size() > 0)
          {
            std::ostringstream drop_stmt;
            drop_stmt
              << "DELETE FROM " << table << " WHERE mid=@mid AND zid=@zid AND kid=@kid AND Pk_id=@Pk_id "
              << "AND IR_id=@IR_id AND UV_id=@UV_id;";
        
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, drop_stmt.str().c_str(), drop_stmt.str().length()+1, &stmt, nullptr));
        
            for(unsigned int t : inconsistent_set)
              {
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), t));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), record.k->get_token().get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), Pk_lin.get_id()));
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
    
    
    void drop_inconsistent_redshifts(sqlite3* db, const FRW_model_token& model, const linear_Pk_token& Pk_lin,
                                     const std::string& table, const resum_Pk_configs::value_type& record,
                                     const std::set<unsigned int>& missing, const std::set<unsigned int>& total_missing)
      {
        std::set<unsigned int> inconsistent_set;
        
        std::set_difference(total_missing.begin(), total_missing.end(),
                            missing.begin(), missing.end(), std::inserter(inconsistent_set, inconsistent_set.begin()));
        
        if(inconsistent_set.size() > 0)
          {
            std::ostringstream drop_stmt;
            drop_stmt
              << "DELETE FROM " << table << " WHERE mid=@mid AND zid=@zid AND kid=@kid AND Pk_id=@Pk_id "
              << "AND IR_cutoff_id=@IR_cutoff_id AND UV_cutoff_id=@UV_cutoff_id AND IR_resum_id=@IR_resum_id;";
            
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, drop_stmt.str().c_str(), drop_stmt.str().length()+1, &stmt, nullptr));
            
            for(unsigned int t : inconsistent_set)
              {
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), t));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), record.k->get_token().get_id()));
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), Pk_lin.get_id()));
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
                                     const FRW_model_token& model, const z_database& z_db, const std::string& z_table)
      {
        assert(db != nullptr);

        // set up null pointer; will be attached to an empty database later if needed
        std::unique_ptr<z_database> missing_db;

        // get list of missing z-values from each relevant table
        std::set<unsigned int> missing;
        std::set<unsigned int> missing_g = update_missing_oneloop_growth_redshifts(db, model, policy.g_factor_table(), z_table, missing);
        std::set<unsigned int> missing_f = update_missing_oneloop_growth_redshifts(db, model, policy.f_factor_table(), z_table, missing);
        
        // drop any inconsistent redshifts -- ie. those present in some tables but not others
        drop_inconsistent_redshifts(db, model, policy.g_factor_table(), missing_g, missing);
        drop_inconsistent_redshifts(db, model, policy.f_factor_table(), missing_f, missing);

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
                                                   const linear_Pk_token& Pk_lin, const std::string& table,
                                                   const std::unordered_set<data_manager_impl::loop_momentum_configuration>& required_configs)
      {
        assert(db != nullptr);
    
        loop_configs missing_configs;
    
        std::ostringstream count_stmt;
        count_stmt
          << "SELECT COUNT(*) FROM (SELECT * FROM " << table << " "
          << "WHERE mid=@mid AND kid=@kid AND Pk_id=@Pk_id AND UV_id=@UV_id AND IR_id=@IR_id);";
    
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, count_stmt.str().c_str(), count_stmt.str().length()+1, &stmt, nullptr));
        
        for(const loop_configs::value_type& t : required_configs)
          {
            // bind parameter values to search for an entry with this combination of model id, k id, UV limit id, IR limit id
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
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
    update_missing_loop_integral_configurations(sqlite3* db, const FRW_model_token& model,
                                                const linear_Pk_token& Pk_lin, const std::string& table,
                                                const loop_configs& required_configs, loop_configs& total_missing)
      {
        loop_configs missing = missing_loop_integral_configurations_for_table(db, model, Pk_lin, table, required_configs);
        
        total_missing.insert(missing.begin(), missing.end());
        
        return missing;
      }
    
    
    void drop_inconsistent_configurations(sqlite3* db, const FRW_model_token& model, const linear_Pk_token& Pk_lin,
                                          const std::string& table, const loop_configs& missing,
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
              << "DELETE FROM " << table << " WHERE mid=@mid AND kid=@kid AND Pk_id=@Pk_id AND UV_id=@UV_id AND IR_id=@IR_id;";
            
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, drop_stmt.str().c_str(), drop_stmt.str().length()+1, &stmt, nullptr));
            
            for(const loop_configs::value_type& t : inconsistent_set)
              {
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
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
                                                      const FRW_model_token& model, const linear_Pk_token& Pk_lin,
                                                      const loop_configs& required_configs)
      {
        loop_configs total_missing;

        // find configurations missing from each table, and merge into global set of missing configurations
        loop_configs missing_AA = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.AA_table(), required_configs, total_missing);
        loop_configs missing_AB = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.AB_table(), required_configs, total_missing);
        loop_configs missing_BB = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.BB_table(), required_configs, total_missing);

        loop_configs missing_D = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.D_table(), required_configs, total_missing);
        loop_configs missing_E = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.E_table(), required_configs, total_missing);
        loop_configs missing_F = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.F_table(), required_configs, total_missing);
        loop_configs missing_G = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.G_table(), required_configs, total_missing);
        loop_configs missing_J1 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.J1_table(), required_configs, total_missing);
        loop_configs missing_J2 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.J2_table(), required_configs, total_missing);

        loop_configs missing_RSD13_a = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD13_a_table(), required_configs, total_missing);
        loop_configs missing_RSD13_b = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD13_b_table(), required_configs, total_missing);
        loop_configs missing_RSD13_c = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD13_c_table(), required_configs, total_missing);
        loop_configs missing_RSD13_d = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD13_d_table(), required_configs, total_missing);
        loop_configs missing_RSD13_e = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD13_e_table(), required_configs, total_missing);
        loop_configs missing_RSD13_f = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD13_f_table(), required_configs, total_missing);
        loop_configs missing_RSD13_g = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD13_g_table(), required_configs, total_missing);
    
        loop_configs missing_RSD22_A1 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD22_A1_table(), required_configs, total_missing);
        loop_configs missing_RSD22_A2 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD22_A2_table(), required_configs, total_missing);
        loop_configs missing_RSD22_A3 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD22_A3_table(), required_configs, total_missing);
        loop_configs missing_RSD22_A4 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD22_A4_table(), required_configs, total_missing);
        loop_configs missing_RSD22_A5 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD22_A5_table(), required_configs, total_missing);
        loop_configs missing_RSD22_B2 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD22_B2_table(), required_configs, total_missing);
        loop_configs missing_RSD22_B3 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD22_B3_table(), required_configs, total_missing);
        loop_configs missing_RSD22_B6 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD22_B6_table(), required_configs, total_missing);
        loop_configs missing_RSD22_B8 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD22_B8_table(), required_configs, total_missing);
        loop_configs missing_RSD22_B9 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD22_B9_table(), required_configs, total_missing);
        loop_configs missing_RSD22_C1 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD22_C1_table(), required_configs, total_missing);
        loop_configs missing_RSD22_C2 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD22_C2_table(), required_configs, total_missing);
        loop_configs missing_RSD22_C4 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD22_C4_table(), required_configs, total_missing);
        loop_configs missing_RSD22_D1 = update_missing_loop_integral_configurations(db, model, Pk_lin, policy.RSD22_D1_table(), required_configs, total_missing);
    
        // drop any inconsistent configurations -- ie. those present in some tables, but not others
        drop_inconsistent_configurations(db, model, Pk_lin, policy.AA_table(), missing_AA, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.AB_table(), missing_AB, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.BB_table(), missing_BB, total_missing);

        drop_inconsistent_configurations(db, model, Pk_lin, policy.D_table(), missing_D, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.E_table(), missing_E, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.F_table(), missing_F, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.G_table(), missing_G, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.J1_table(), missing_J1, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.J2_table(), missing_J2, total_missing);
    
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD13_a_table(), missing_RSD13_a, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD13_b_table(), missing_RSD13_b, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD13_c_table(), missing_RSD13_c, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD13_d_table(), missing_RSD13_d, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD13_e_table(), missing_RSD13_e, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD13_f_table(), missing_RSD13_f, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD13_g_table(), missing_RSD13_g, total_missing);
    
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD22_A1_table(), missing_RSD22_A1, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD22_A2_table(), missing_RSD22_A2, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD22_A3_table(), missing_RSD22_A3, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD22_A4_table(), missing_RSD22_A4, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD22_A5_table(), missing_RSD22_A5, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD22_B2_table(), missing_RSD22_B2, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD22_B3_table(), missing_RSD22_B3, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD22_B6_table(), missing_RSD22_B6, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD22_B8_table(), missing_RSD22_B8, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD22_B9_table(), missing_RSD22_B9, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD22_C1_table(), missing_RSD22_C1, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD22_C2_table(), missing_RSD22_C2, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD22_C4_table(), missing_RSD22_C4, total_missing);
        drop_inconsistent_configurations(db, model, Pk_lin, policy.RSD22_D1_table(), missing_RSD22_D1, total_missing);
    
        return total_missing;
      }
    
    
    std::unique_ptr<z_database>
    missing_one_loop_Pk_redshifts(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                  const FRW_model_token& model, const linear_Pk_token& Pk_lin,
                                  const std::string& z_table, const z_database& z_db,
                                  const loop_configs::value_type& record)
      {
        assert(db != nullptr);
        
        // set up a null pointer for the returned database; will be attached to an empty instance later if needed
        std::unique_ptr<z_database> missing_db;
        
        // get list of missing z-values from each relevant table
        std::set<unsigned int> missing;
        std::set<unsigned int> missing_delta = update_missing_one_loop_Pk(db, model, Pk_lin, policy.dd_Pk_table(), z_table, record, missing);
        std::set<unsigned int> missing_rsd0 = update_missing_one_loop_Pk(db, model, Pk_lin, policy.dd_rsd_mu0_Pk_table(), z_table, record, missing);
        std::set<unsigned int> missing_rsd2 = update_missing_one_loop_Pk(db, model, Pk_lin, policy.dd_rsd_mu2_Pk_table(), z_table, record, missing);
        std::set<unsigned int> missing_rsd4 = update_missing_one_loop_Pk(db, model, Pk_lin, policy.dd_rsd_mu4_Pk_table(), z_table, record, missing);
        std::set<unsigned int> missing_rsd6 = update_missing_one_loop_Pk(db, model, Pk_lin, policy.dd_rsd_mu6_Pk_table(), z_table, record, missing);
        std::set<unsigned int> missing_rsd8 = update_missing_one_loop_Pk(db, model, Pk_lin, policy.dd_rsd_mu8_Pk_table(), z_table, record, missing);
        
        // drop any inconsistent redshifts that are present in one table but not the others
        drop_inconsistent_redshifts(db, model, Pk_lin, policy.dd_Pk_table(), record, missing_delta, missing);
        drop_inconsistent_redshifts(db, model, Pk_lin, policy.dd_rsd_mu0_Pk_table(), record, missing_rsd0, missing);
        drop_inconsistent_redshifts(db, model, Pk_lin, policy.dd_rsd_mu2_Pk_table(), record, missing_rsd2, missing);
        drop_inconsistent_redshifts(db, model, Pk_lin, policy.dd_rsd_mu4_Pk_table(), record, missing_rsd4, missing);
        drop_inconsistent_redshifts(db, model, Pk_lin, policy.dd_rsd_mu6_Pk_table(), record, missing_rsd6, missing);
        drop_inconsistent_redshifts(db, model, Pk_lin, policy.dd_rsd_mu8_Pk_table(), record, missing_rsd8, missing);
        
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
    missing_one_loop_resum_Pk_redshifts(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                        const FRW_model_token& model, const linear_Pk_token& Pk_lin,
                                        const std::string& z_table, const z_database& z_db,
                                        const resum_Pk_configs::value_type& record)
      {
        assert(db != nullptr);
        
        // set up a null pointer for the returned database; will be attached to an empty instance later if needed
        std::unique_ptr<z_database> missing_db;

        // get missing redshifts for this configuration
        std::set<unsigned int> missing =
          missing_redshifts_for_table(db, model, Pk_lin, policy.dd_Pk_resum_table(), z_table, record);
    
        if(!missing.empty())
          {
            missing_db = std::make_unique<z_database>();
            
            for(unsigned int t : missing)
              {
                // lookup record for this identifier
                z_database::const_record_iterator rec = z_db.lookup(z_token(t));
                
                // add a corresponding record to this missing database
                missing_db->add_record(*(*rec), rec->get_token());
              }
          }
    
        return missing_db;
      }
    
    
    std::unique_ptr<z_database>
    missing_multipole_Pk_redshifts(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                   const FRW_model_token& model, const linear_Pk_token& Pk_lin,
                                   const std::string& z_table, const z_database& z_db,
                                   const resum_Pk_configs::value_type& record)
      {
        assert(db != nullptr);
    
        // set up a null pointer for the returned database; will be attached to an empty instance later if needed
        std::unique_ptr<z_database> missing_db;
    
        // get list of missing z-values from each relevant table
        std::set<unsigned int> missing;
        std::set<unsigned int> missing_P0 = update_missing_multipole_Pk(db, model, Pk_lin, policy.P0_table(), z_table, record, missing);
        std::set<unsigned int> missing_P2 = update_missing_multipole_Pk(db, model, Pk_lin, policy.P2_table(), z_table, record, missing);
        std::set<unsigned int> missing_P4 = update_missing_multipole_Pk(db, model, Pk_lin, policy.P4_table(), z_table, record, missing);
    
        // drop any inconsistent redshifts that are present in one table but not the others
        drop_inconsistent_redshifts(db, model, Pk_lin, policy.P0_table(), record, missing_P0, missing);
        drop_inconsistent_redshifts(db, model, Pk_lin, policy.P2_table(), record, missing_P2, missing);
        drop_inconsistent_redshifts(db, model, Pk_lin, policy.P4_table(), record, missing_P4, missing);
        
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
                                        const IR_resum_database& IR_db)
      {
        assert(db != nullptr);
        
        Matsubara_configs missing_configs;
        
        std::ostringstream count_stmt;
        count_stmt
          << "SELECT COUNT(*) FROM (SELECT * FROM " << policy.Matsubara_XY_table() << " "
          << "WHERE mid=@mid AND Pk_id=@Pk_id AND IR_resum_id=@IR_resum_id);";
    
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, count_stmt.str().c_str(), count_stmt.str().length()+1, &stmt, nullptr));
        
        for(IR_resum_database::const_record_iterator t = IR_db.record_cbegin(); t != IR_db.record_cend(); ++t)
          {
            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
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
                                  const linear_Pk_token& token, const k_database& k_db, const std::string& k_table)
      {
        assert(db != nullptr);
        
        // set up null pointer; will be attached to an empty database later if needed
        std::unique_ptr<k_database> missing_db;
        
        // get list of missing k-values for this linear power spectrum
        std::set<unsigned int> missing = missing_wavenumbers_for_table(db, token, policy.Pk_linear_table(), k_table);
        
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
    
    
  }   // namespace sqlite3_operations
