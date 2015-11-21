//
// Created by David Seery on 17/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include <iostream>
#include <sstream>

#include "store.h"
#include "utilities.h"

#include "localizations/messages.h"


namespace sqlite3_operations
  {


    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const transfer_function& sample)
      {
        assert(db != nullptr);

        // construct SQL insert statement
        std::ostringstream insert_stmt;
        insert_stmt
          << "INSERT INTO " << policy.transfer_table() << " VALUES (@mid, @kid, @zid, @delta_m, @delta_r, @theta_m, @theta_r, @Phi);";

        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));

        // get wavenumber token
        const k_token& k_token = sample.get_k_token();

        // loop through sample, writing its values into the database
        for(transfer_function::const_token_iterator t = sample.token_begin(); t != sample.token_end(); ++t)
          {
            transfer_value val = *t;

            // bind values to the statement
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), k_token.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), val.first.get_id()));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@delta_m"), val.second.delta_m));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@delta_r"), val.second.delta_r));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@theta_m"), val.second.theta_m));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@theta_r"), val.second.theta_r));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Phi"), val.second.Phi));

            // perform insertion
            check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_TRANSFER_FAIL, SQLITE_DONE);

            // clear bindings and reset statement
            check_stmt(db, sqlite3_clear_bindings(stmt));
            check_stmt(db, sqlite3_reset(stmt));
          }

        // finalize statement and release resources
        check_stmt(db, sqlite3_finalize(stmt));
      }


    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const oneloop_growth& sample)
      {
        assert(db != nullptr);

        // construct SQL insert statement
        std::ostringstream insert_stmt;
        insert_stmt
          << "INSERT INTO " << policy.oneloop_table() << " VALUES (@mid, @zid, @g_linear, @A, @B, @D, @E, @F, @G);";

        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));

        // loop through sample, writing its values into the database
        for(oneloop_growth::const_token_iterator t = sample.token_begin(); t != sample.token_end(); ++t)
          {
            oneloop_value val = *t;

            // bind values to the statement
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), val.first.get_id()));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@g_linear"), val.second.g));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@A"), val.second.A));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@B"), val.second.B));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@D"), val.second.D));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@E"), val.second.E));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@F"), val.second.F));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@G"), val.second.G));

            // perform insertion
            check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_ONELOOP_FAIL, SQLITE_DONE);

            // clear bindings and reset statement
            check_stmt(db, sqlite3_clear_bindings(stmt));
            check_stmt(db, sqlite3_reset(stmt));
          }

        // finalize statement and release resources
        check_stmt(db, sqlite3_finalize(stmt));
      }


    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const loop_integral& sample)
      {
        assert(db != nullptr);

        // construct SQL insert statement
        std::ostringstream insert_stmt;
        insert_stmt
          << "INSERT INTO " << policy.loop_momentum_table() << " VALUES (@mid, @kid, @IR_id, @UV_id, @A, @B, @D, @E, @F, @G);";

        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));

        // bind values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), sample.get_k_token().get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), sample.get_IR_token().get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), sample.get_UV_token().get_id()));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@A"), sample.get_A()));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@B"), sample.get_B()));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@D"), sample.get_D()));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@E"), sample.get_E()));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@F"), sample.get_F()));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@G"), sample.get_G()));

        // perform insertion
        check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_LOOP_MOMENTUM_FAIL, SQLITE_DONE);

        // clear bindings and release
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));
      }


  }   // namespace sqlite3_operations
