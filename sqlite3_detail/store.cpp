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
    
    
    template <typename DimensionfulType>
    DimensionfulType dimensionful_unit();

    
    template <>
    double dimensionful_unit<double>()
      {
        return 1.0;
      }
    
    template <>
    Mpc_units::inverse_energy3 dimensionful_unit<Mpc_units::inverse_energy3>()
      {
        return Mpc_units::Mpc3;
      }
    
    
    namespace store_impl
      {
    
        template <typename KernelType>
        void store_loop_integral(sqlite3* db, const std::string table_name, const KernelType& kernel,
                                 const FRW_model_token& model, const loop_integral& sample)
          {
            std::ostringstream insert_stmt;
            insert_stmt << "INSERT INTO " << table_name << " VALUES (@mid, @kid, @IR_id, @UV_id, @value, @regions, @evals, @err, @time);";
            
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));
            
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), sample.get_k_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), sample.get_IR_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), sample.get_UV_token().get_id()));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@value"), kernel.value / dimensionful_unit<typename KernelType::value_type>()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@regions"), kernel.regions));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@evals"), kernel.evaluations));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@err"), kernel.error));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@time"), kernel.time));
    
            // perform insertion
            check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_LOOP_MOMENTUM_FAIL, SQLITE_DONE);
    
            // clear bindings and release
            check_stmt(db, sqlite3_clear_bindings(stmt));
            check_stmt(db, sqlite3_finalize(stmt));
          }
        
      }   // namespace store_impl


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
        for(const transfer_value& val : sample)
          {
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

        // construct SQL insert statements
        std::ostringstream insert_g_stmt;
        insert_g_stmt
          << "INSERT INTO " << policy.g_factor_table() << " VALUES (@mid, @zid, @g_linear, @A, @B, @D, @E, @F, @G, @J);";
        
        std::ostringstream insert_f_stmt;
        insert_f_stmt
          << "INSERT INTO " << policy.f_factor_table() << " VALUES (@mid, @zid, @f_linear, @fA, @fB, @fD, @fE, @fF, @fG, @fJ);";
        
        std::ostringstream insert_meta_stmt;
        insert_meta_stmt
          << "INSERT INTO " << policy.gf_metadata_table() << " VALUES (@mid, @time, @steps);";

        // prepare statements
        sqlite3_stmt* g_stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_g_stmt.str().c_str(), insert_g_stmt.str().length()+1, &g_stmt, nullptr));
    
        sqlite3_stmt* f_stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_f_stmt.str().c_str(), insert_f_stmt.str().length()+1, &f_stmt, nullptr));
        
        sqlite3_stmt* meta_stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_meta_stmt.str().c_str(), insert_meta_stmt.str().length()+1, &meta_stmt, nullptr));
        
        check_stmt(db, sqlite3_bind_int(meta_stmt, sqlite3_bind_parameter_index(meta_stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int64(meta_stmt, sqlite3_bind_parameter_index(meta_stmt, "@time"), sample.get_integration_time()));
        check_stmt(db, sqlite3_bind_int(meta_stmt, sqlite3_bind_parameter_index(meta_stmt, "@steps"), static_cast<int>(sample.get_integration_steps())));
        
        check_stmt(db, sqlite3_step(meta_stmt), ERROR_SQLITE3_INSERT_GROWTH_META_FAIL, SQLITE_DONE);
    
        // loop through sample, writing its values into the database
        for(const oneloop_value& val : sample)
          {
            // bind values to the g statement
            check_stmt(db, sqlite3_bind_int(g_stmt, sqlite3_bind_parameter_index(g_stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(g_stmt, sqlite3_bind_parameter_index(g_stmt, "@zid"), val.first.get_id()));
            check_stmt(db, sqlite3_bind_double(g_stmt, sqlite3_bind_parameter_index(g_stmt, "@g_linear"), val.second.g));
            check_stmt(db, sqlite3_bind_double(g_stmt, sqlite3_bind_parameter_index(g_stmt, "@A"), val.second.A));
            check_stmt(db, sqlite3_bind_double(g_stmt, sqlite3_bind_parameter_index(g_stmt, "@B"), val.second.B));
            check_stmt(db, sqlite3_bind_double(g_stmt, sqlite3_bind_parameter_index(g_stmt, "@D"), val.second.D));
            check_stmt(db, sqlite3_bind_double(g_stmt, sqlite3_bind_parameter_index(g_stmt, "@E"), val.second.E));
            check_stmt(db, sqlite3_bind_double(g_stmt, sqlite3_bind_parameter_index(g_stmt, "@F"), val.second.F));
            check_stmt(db, sqlite3_bind_double(g_stmt, sqlite3_bind_parameter_index(g_stmt, "@G"), val.second.G));
            check_stmt(db, sqlite3_bind_double(g_stmt, sqlite3_bind_parameter_index(g_stmt, "@J"), val.second.J));

            // perform g insertion
            check_stmt(db, sqlite3_step(g_stmt), ERROR_SQLITE3_INSERT_GROWTH_G_FAIL, SQLITE_DONE);
    
            // bind values to the f statement
            check_stmt(db, sqlite3_bind_int(f_stmt, sqlite3_bind_parameter_index(f_stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(f_stmt, sqlite3_bind_parameter_index(f_stmt, "@zid"), val.first.get_id()));
            check_stmt(db, sqlite3_bind_double(f_stmt, sqlite3_bind_parameter_index(f_stmt, "@f_linear"), val.second.f));
            check_stmt(db, sqlite3_bind_double(f_stmt, sqlite3_bind_parameter_index(f_stmt, "@fA"), val.second.fA));
            check_stmt(db, sqlite3_bind_double(f_stmt, sqlite3_bind_parameter_index(f_stmt, "@fB"), val.second.fB));
            check_stmt(db, sqlite3_bind_double(f_stmt, sqlite3_bind_parameter_index(f_stmt, "@fD"), val.second.fD));
            check_stmt(db, sqlite3_bind_double(f_stmt, sqlite3_bind_parameter_index(f_stmt, "@fE"), val.second.fE));
            check_stmt(db, sqlite3_bind_double(f_stmt, sqlite3_bind_parameter_index(f_stmt, "@fF"), val.second.fF));
            check_stmt(db, sqlite3_bind_double(f_stmt, sqlite3_bind_parameter_index(f_stmt, "@fG"), val.second.fG));
            check_stmt(db, sqlite3_bind_double(f_stmt, sqlite3_bind_parameter_index(f_stmt, "@fJ"), val.second.fJ));
    
            // perform f insertion
            check_stmt(db, sqlite3_step(f_stmt), ERROR_SQLITE3_INSERT_GROWTH_F_FAIL, SQLITE_DONE);

            // clear bindings and reset statement
            check_stmt(db, sqlite3_clear_bindings(g_stmt));
            check_stmt(db, sqlite3_reset(g_stmt));
            check_stmt(db, sqlite3_clear_bindings(f_stmt));
            check_stmt(db, sqlite3_reset(f_stmt));
          }

        // finalize statement and release resources
        check_stmt(db, sqlite3_finalize(g_stmt));
        check_stmt(db, sqlite3_finalize(f_stmt));
        check_stmt(db, sqlite3_finalize(meta_stmt));
      }


    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model, const loop_integral& sample)
      {
        assert(db != nullptr);

        const delta_22_integrals& delta22 = sample.get_delta22();
        const delta_13_integrals& delta13 = sample.get_delta13();
        const rsd_22_integrals& rsd22     = sample.get_rsd22();
        const rsd_13_integrals& rsd13     = sample.get_rsd13();

        if(!delta22.get_fail() && !delta13.get_fail() && !rsd22.get_fail() && !rsd13.get_fail())
          {
            store_impl::store_loop_integral(db, policy.AA_table(), delta22.get_AA(), model, sample);
            store_impl::store_loop_integral(db, policy.AB_table(), delta22.get_AB(), model, sample);
            store_impl::store_loop_integral(db, policy.BB_table(), delta22.get_BB(), model, sample);

            store_impl::store_loop_integral(db, policy.D_table(), delta13.get_D(), model, sample);
            store_impl::store_loop_integral(db, policy.E_table(), delta13.get_E(), model, sample);
            store_impl::store_loop_integral(db, policy.F_table(), delta13.get_F(), model, sample);
            store_impl::store_loop_integral(db, policy.G_table(), delta13.get_G(), model, sample);
            store_impl::store_loop_integral(db, policy.J1_table(), delta13.get_J1(), model, sample);
            store_impl::store_loop_integral(db, policy.J2_table(), delta13.get_J2(), model, sample);
          }
      }


  }   // namespace sqlite3_operations
