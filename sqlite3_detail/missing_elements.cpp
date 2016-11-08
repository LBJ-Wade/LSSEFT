//
// Created by David Seery on 13/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <list>
#include <algorithm>
#include <assert.h>

#include "utilities.h"
#include "missing_elements.h"


namespace sqlite3_operations
  {

    //! find missing redshifts for a k-dependent table
    std::list<unsigned int> missing_redshifts_for_table(sqlite3* db, const FRW_model_token& model,
                                                        const k_token& k,
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

        std::list<unsigned int> results;

        int status = 0;
        while((status = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(status == SQLITE_ROW)
              {
                results.push_back(static_cast<unsigned int>(sqlite3_column_int(stmt, 0)));
              }
          }

        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));

        return(results);
      }


    //! find missing redshifts for a k-independent table
    std::list<unsigned int> missing_redshifts_for_table(sqlite3* db, const FRW_model_token& model,
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

        std::list<unsigned int> results;

        int status = 0;
        while((status = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(status == SQLITE_ROW)
              {
                results.push_back(static_cast<unsigned int>(sqlite3_column_int(stmt, 0)));
              }
          }

        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));

        return(results);
      }


    //! find missing wavenumbers
    std::list<unsigned int> missing_wavenumbers_for_table(sqlite3* db, const FRW_model_token& model,
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

        std::list<unsigned int> results;

        int status = 0;
        while((status = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(status == SQLITE_ROW)
              {
                results.push_back(static_cast<unsigned int>(sqlite3_column_int(stmt, 0)));
              }
          }

        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));

        return(results);
      }


    //! find missing redshifts for a transfer function sample
    std::unique_ptr<z_database> missing_transfer_redshifts(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                                           const FRW_model_token& model, const k_token& k,
                                                           const z_database& z_db, const std::string& z_table)
      {
        assert(db != nullptr);

        // set up null pointer; will be attached to an empty database later if needed
        std::unique_ptr<z_database> missing_db;

        // get list of missing z-values for this k-mode of the transfer function
        std::list<unsigned int> missing_list = missing_redshifts_for_table(db, model, k, policy.transfer_table(), z_table);

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
    std::unique_ptr<z_database> missing_oneloop_redshifts(sqlite3* db, transaction_manager& mgr,
                                                          const sqlite3_policy& policy, const FRW_model_token& model,
                                                          const z_database& z_db, const std::string& z_table)
      {
        assert(db != nullptr);

        // set up null pointer; will be attached to an empty database later if needed
        std::unique_ptr<z_database> missing_db;

        // get list of missing z-values
        std::list<unsigned int> missing_list = missing_redshifts_for_table(db, model, policy.oneloop_table(), z_table);

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


    //! find missing wavenumbers for a loop momentum sample
    void missing_loop_momentum(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                               const FRW_model_token& model, const k_database& k_db,
                               std::list<data_manager_impl::loop_momentum_configuration>& combinations)
      {
        assert(db != nullptr);

        std::ostringstream count_stmt;
        count_stmt
        << "SELECT COUNT(*) FROM (SELECT * FROM " << policy.loop_momentum_table() << " "
        << "WHERE mid=@mid AND kid=@kid AND UV_id=@UV_id AND IR_id=@IR_id);";

        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, count_stmt.str().c_str(), count_stmt.str().length()+1, &stmt, nullptr));

        for(data_manager_impl::loop_momentum_configuration& t : combinations)
          {
            // bind paramaeter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), t.k->get_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), t.UV->get_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), t.IR->get_token().get_id()));

            int status = 0;
            while((status = sqlite3_step(stmt)) != SQLITE_DONE)
              {
                if(status == SQLITE_ROW)
                  {
                    if(sqlite3_column_int(stmt, 0) > 0) t.include = false;
                  }
              }

            // release bindings and reset statement
            check_stmt(db, sqlite3_clear_bindings(stmt));
            check_stmt(db, sqlite3_reset(stmt));
          }

        // finalize stataement to release resources
        check_stmt(db, sqlite3_finalize(stmt));
      }


  }   // namespace sqlite3_operations
