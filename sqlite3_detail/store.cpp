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
    
    namespace store_impl
      {
        
        template <typename ValueType>
        double make_dimensionless(const ValueType& value)
          {
            return value / dimensionful_unit<ValueType>();
          }
        
        
        template <typename ValueContainer>
        auto extract_value(const ValueContainer& container) -> decltype(container.value)
          {
            return container.value;
          }
    
    
        template <typename ValueContainer>
        auto extract_value(const ValueContainer& container) -> decltype(container.get_value())
          {
            return container.get_value();
          }
    
    
        template <typename ValueContainer>
        auto extract_error(const ValueContainer& container) -> decltype(container.error)
          {
            return container.error;
          }
    
    
        template <typename ValueContainer>
        auto extract_error(const ValueContainer& container) -> decltype(container.get_error())
          {
            return container.get_error();
          }

        
        template <typename ValueContainer>
        double dimensionless_value(const ValueContainer& container)
          {
            return make_dimensionless(extract_value(container));
          }
    
    
        template <typename ValueContainer>
        double dimensionless_error(const ValueContainer& container)
          {
            return make_dimensionless(extract_error(container));
          }
    
    
        template <typename KernelType>
        void store_loop_kernel(sqlite3* db, const std::string& table_name, const KernelType& kernel,
                               const FRW_model_token& model, const loop_integral& sample)
          {
            std::ostringstream insert_stmt;
            insert_stmt
              << "INSERT INTO " << table_name << " VALUES (@mid, @kid, @Pk_id, @IR_id, @UV_id, "
              << "@raw_value, @raw_regions, @raw_evals, @raw_err, @raw_time, "
              << "@wiggle_value, @wiggle_regions, @wiggle_evals, @wiggle_err, @wiggle_time);";
            
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));
            
            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), sample.get_k_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), sample.get_Pk_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), sample.get_IR_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), sample.get_UV_token().get_id()));

            auto raw = kernel.get_raw();
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@raw_value"), dimensionless_value(raw)));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@raw_regions"),raw.regions));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@raw_evals"), raw.evaluations));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@raw_err"), dimensionless_error(raw)));
            check_stmt(db, sqlite3_bind_int64(stmt, sqlite3_bind_parameter_index(stmt, "@raw_time"), raw.time));
    
            auto wiggle = kernel.get_wiggle();
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@wiggle_value"), dimensionless_value(wiggle)));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@wiggle_regions"),wiggle.regions));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@wiggle_evals"), wiggle.evaluations));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@wiggle_err"), dimensionless_error(wiggle)));
            check_stmt(db, sqlite3_bind_int64(stmt, sqlite3_bind_parameter_index(stmt, "@wiggle_time"), wiggle.time));
    
            // perform insertion
            check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_LOOP_MOMENTUM_FAIL, SQLITE_DONE);
    
            // clear bindings and release
            check_stmt(db, sqlite3_clear_bindings(stmt));
            check_stmt(db, sqlite3_finalize(stmt));
          }
        
        
        template <typename ValueType>
        void store_Pk_value(sqlite3* db, sqlite3_stmt* stmt, const std::string& value_raw, const std::string& error_raw,
                            const std::string& value_wiggle, const std::string& error_wiggle, const ValueType& value)
          {
            const auto& raw = value.get_raw();
            const auto& wiggle = value.get_wiggle();
            
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, value_raw.c_str()), dimensionless_value(raw)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, error_raw.c_str()), dimensionless_error(raw)));
    
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, value_wiggle.c_str()), dimensionless_value(wiggle)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, error_wiggle.c_str()), dimensionless_error(wiggle)));
          }
    
    
        template <typename ValueType>
        void store_Pk_value(sqlite3* db, sqlite3_stmt* stmt, const std::string& value_raw, const std::string& value_wiggle, const ValueType& value)
          {
            const auto& raw = value.get_raw();
            const auto& wiggle = value.get_wiggle();

            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, value_raw.c_str()), dimensionless_value(raw)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, value_wiggle.c_str()), dimensionless_value(wiggle)));
          }
    
    
        template <typename PkType>
        void store_one_loop_Pk(sqlite3* db, const std::string& table_name, const PkType& value,
                               const FRW_model_token& model, const oneloop_Pk& sample)
          {
            std::ostringstream insert_stmt;
            insert_stmt
              << "INSERT INTO " << table_name << " VALUES (@mid, @zid, @kid, @Pk_id, @IR_id, @UV_id, "
              << "@Ptree_raw, @err_tree_raw, @P13_raw, @err_13_raw, @P22_raw, @err_22_raw, @P1loopSPT_raw, @err_1loopSPT_raw, @Z2_delta_raw, "
              << "@Ptree_wiggle, @err_tree_wiggle, @P13_wiggle, @err_13_wiggle, @P22_wiggle, @err_22_wiggle, @P1loopSPT_wiggle, @err_1loopSPT_wiggle, @Z2_delta_wiggle"
              << ");";
    
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));
    
            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), sample.get_z_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), sample.get_k_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), sample.get_Pk_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), sample.get_IR_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), sample.get_UV_token().get_id()));
    
            store_Pk_value(db, stmt, "@Ptree_raw", "@err_tree_raw", "@Ptree_wiggle", "@err_tree_wiggle", value.get_tree());
            store_Pk_value(db, stmt, "@P13_raw", "@err_13_raw", "@P13_wiggle", "@err_13_wiggle", value.get_13());
            store_Pk_value(db, stmt, "@P22_raw", "@err_22_raw", "@P22_wiggle", "@err_22_wiggle", value.get_22());
            store_Pk_value(db, stmt, "@P1loopSPT_raw", "@err_1loopSPT_raw", "@P1loopSPT_wiggle", "@err_1loopSPT_wiggle", value.get_1loop_SPT());
            store_Pk_value(db, stmt, "@Z2_delta_raw", "@Z2_delta_wiggle", value.get_Z2_delta());
            
            // perform insertion
            check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_ONELOOP_PK_FAIL, SQLITE_DONE);
    
            // clear bindings and release
            check_stmt(db, sqlite3_clear_bindings(stmt));
            check_stmt(db, sqlite3_finalize(stmt));
          }
    
        
        template <typename PkType>
        void store_one_loop_rsd_Pk(sqlite3* db, const std::string& table_name, const PkType& value,
                                   const FRW_model_token& model, const oneloop_Pk& sample)
          {
            std::ostringstream insert_stmt;
            insert_stmt
              << "INSERT INTO " << table_name << " VALUES (@mid, @zid, @kid, @Pk_id, @IR_id, @UV_id, "
              << "@Ptree_raw, @err_tree_raw, @P13_raw, @err_13_raw, @P22_raw, @err_22_raw, @P1loopSPT_raw, @err_1loopSPT_raw, @Z2_delta_raw, @Z0_v_raw, @Z2_v_raw, @Z0_vdelta_raw, @Z2_vdelta_raw, @Z2_vv_raw, @Z2_vvdelta_raw, @Z2_vvv_raw, "
              << "@Ptree_wiggle, @err_tree_wiggle, @P13_wiggle, @err_13_wiggle, @P22_wiggle, @err_22_wiggle, @P1loopSPT_wiggle, @err_1loopSPT_wiggle, @Z2_delta_wiggle, @Z0_v_wiggle, @Z2_v_wiggle, @Z0_vdelta_wiggle, @Z2_vdelta_wiggle, @Z2_vv_wiggle, @Z2_vvdelta_wiggle, @Z2_vvv_wiggle"
              << ");";
        
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));
        
            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), sample.get_z_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), sample.get_k_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), sample.get_Pk_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), sample.get_IR_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), sample.get_UV_token().get_id()));

            store_Pk_value(db, stmt, "@Ptree_raw", "@err_tree_raw", "@Ptree_wiggle", "@err_tree_wiggle", value.get_tree());
            store_Pk_value(db, stmt, "@P13_raw", "@err_13_raw", "@P13_wiggle", "@err_13_wiggle", value.get_13());
            store_Pk_value(db, stmt, "@P22_raw", "@err_22_raw", "@P22_wiggle", "@err_22_wiggle", value.get_22());
            store_Pk_value(db, stmt, "@P1loopSPT_raw", "@err_1loopSPT_raw", "@P1loopSPT_wiggle", "@err_1loopSPT_wiggle", value.get_1loop_SPT());
            store_Pk_value(db, stmt, "@Z2_delta_raw", "@Z2_delta_wiggle", value.get_Z2_delta());
            store_Pk_value(db, stmt, "@Z0_v_raw", "@Z0_v_wiggle", value.get_Z0_v());
            store_Pk_value(db, stmt, "@Z2_v_raw", "@Z2_v_wiggle", value.get_Z2_v());
            store_Pk_value(db, stmt, "@Z0_vdelta_raw", "@Z0_vdelta_wiggle", value.get_Z0_vdelta());
            store_Pk_value(db, stmt, "@Z2_vdelta_raw", "@Z2_vdelta_wiggle", value.get_Z2_vdelta());
            store_Pk_value(db, stmt, "@Z2_vv_raw", "@Z2_vv_wiggle", value.get_Z2_vv());
            store_Pk_value(db, stmt, "@Z2_vvdelta_raw", "@Z2_vvdelta_wiggle", value.get_Z2_vvdelta());
            store_Pk_value(db, stmt, "@Z2_vvv_raw", "@Z2_vvv_wiggle", value.get_Z2_vvv());
            
            // perform insertion
            check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_ONELOOP_RSD_PK_FAIL, SQLITE_DONE);
        
            // clear bindings and release
            check_stmt(db, sqlite3_clear_bindings(stmt));
            check_stmt(db, sqlite3_finalize(stmt));
          }
    
    
        void store_multipole_Pk(sqlite3* db, const std::string& table_name, const Pk_ell& value,
                                const FRW_model_token& model, const multipole_Pk& sample)
          {
            std::ostringstream insert_stmt;
            insert_stmt
              << "INSERT INTO " << table_name << " VALUES (@mid, @zid, @kid, @Pk_id, @IR_cutoff_id, @UV_cutoff_id, @IR_resum_id, "
                                              << "@Ptree, @Ptree_resum, @P13, @P13_resum, @P22, @P22_resum, "
                                              << "@P1loopSPT, @P1loopSPT_resum, @Z2_delta, @Z0_v, @Z2_v, @Z0_vdelta, @Z2_vdelta, "
                                              << "@Z2_vv, @Z2_vvdelta, @Z2_vvv);";
    
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));
    
            // bind parameter values
            const auto& tree = value.get_tree();
            const auto& tree_resum = value.get_tree_resum();
            const auto& P13 = value.get_13();
            const auto& P13_resum = value.get_13_resum();
            const auto& P22 = value.get_22();
            const auto& P22_resum = value.get_22_resum();
            const auto& P1loopSPT = value.get_1loop_SPT();
            const auto& P1loopSPT_resum = value.get_1loop_SPT_resum();
            const auto& Z2_delta = value.get_Z2_delta();
            const auto& Z0_v = value.get_Z0_v();
            const auto& Z2_v = value.get_Z2_v();
            const auto& Z0_vdelta = value.get_Z0_vdelta();
            const auto& Z2_vdelta = value.get_Z2_vdelta();
            const auto& Z2_vv = value.get_Z2_vv();
            const auto& Z2_vvdelta = value.get_Z2_vvdelta();
            const auto& Z2_vvv = value.get_Z2_vvv();
            
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), sample.get_z_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), sample.get_k_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), sample.get_Pk_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_cutoff_id"), sample.get_IR_cutoff_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_cutoff_id"), sample.get_UV_cutoff_token().get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_resum_id"), sample.get_IR_resum_token().get_id()));

            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Ptree"), make_dimensionless(tree)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Ptree_resum"), make_dimensionless(tree_resum)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@P13"), make_dimensionless(P13)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@P13_resum"), make_dimensionless(P13_resum)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@P22"), make_dimensionless(P22)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@P22_resum"), make_dimensionless(P22_resum)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@P1loopSPT"), make_dimensionless(P1loopSPT)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@P1loopSPT_resum"), make_dimensionless(P1loopSPT_resum)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Z2_delta"), make_dimensionless(Z2_delta)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Z0_v"), make_dimensionless(Z0_v)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Z2_v"), make_dimensionless(Z2_v)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Z0_vdelta"), make_dimensionless(Z0_vdelta)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Z2_vdelta"), make_dimensionless(Z2_vdelta)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Z2_vv"), make_dimensionless(Z2_vv)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Z2_vvdelta"), make_dimensionless(Z2_vvdelta)));
            check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Z2_vvv"), make_dimensionless(Z2_vvv)));
    
            // perform insertion
            check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_MULTIPOLE_PK_FAIL, SQLITE_DONE);
    
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

        // prepare statements
        sqlite3_stmt* g_stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_g_stmt.str().c_str(), insert_g_stmt.str().length()+1, &g_stmt, nullptr));
    
        sqlite3_stmt* f_stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_f_stmt.str().c_str(), insert_f_stmt.str().length()+1, &f_stmt, nullptr));
    
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
            store_impl::store_loop_kernel(db, policy.AA_table(), delta22.get_AA(), model, sample);
            store_impl::store_loop_kernel(db, policy.AB_table(), delta22.get_AB(), model, sample);
            store_impl::store_loop_kernel(db, policy.BB_table(), delta22.get_BB(), model, sample);
    
            store_impl::store_loop_kernel(db, policy.D_table(), delta13.get_D(), model, sample);
            store_impl::store_loop_kernel(db, policy.E_table(), delta13.get_E(), model, sample);
            store_impl::store_loop_kernel(db, policy.F_table(), delta13.get_F(), model, sample);
            store_impl::store_loop_kernel(db, policy.G_table(), delta13.get_G(), model, sample);
            store_impl::store_loop_kernel(db, policy.J1_table(), delta13.get_J1(), model, sample);
            store_impl::store_loop_kernel(db, policy.J2_table(), delta13.get_J2(), model, sample);
            
            store_impl::store_loop_kernel(db, policy.RSD13_a_table(), rsd13.get_a(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD13_b_table(), rsd13.get_b(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD13_c_table(), rsd13.get_c(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD13_d_table(), rsd13.get_d(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD13_e_table(), rsd13.get_e(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD13_f_table(), rsd13.get_f(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD13_g_table(), rsd13.get_g(), model, sample);
    
            store_impl::store_loop_kernel(db, policy.RSD22_A1_table(), rsd22.get_A1(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD22_A2_table(), rsd22.get_A2(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD22_A3_table(), rsd22.get_A3(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD22_A4_table(), rsd22.get_A4(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD22_A5_table(), rsd22.get_A5(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD22_B2_table(), rsd22.get_B2(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD22_B3_table(), rsd22.get_B3(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD22_B6_table(), rsd22.get_B6(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD22_B8_table(), rsd22.get_B8(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD22_B9_table(), rsd22.get_B9(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD22_C1_table(), rsd22.get_C1(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD22_C2_table(), rsd22.get_C2(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD22_C4_table(), rsd22.get_C4(), model, sample);
            store_impl::store_loop_kernel(db, policy.RSD22_D1_table(), rsd22.get_D1(), model, sample);
          }
      }
    
    
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
               const oneloop_Pk& sample)
      {
        assert(db != nullptr);
        
        store_impl::store_one_loop_Pk(db, policy.dd_Pk_table(), sample.get_dd(), model, sample);

        store_impl::store_one_loop_rsd_Pk(db, policy.dd_rsd_mu0_Pk_table(), sample.get_dd_rsd_mu0(), model, sample);
        store_impl::store_one_loop_rsd_Pk(db, policy.dd_rsd_mu2_Pk_table(), sample.get_dd_rsd_mu2(), model, sample);
        store_impl::store_one_loop_rsd_Pk(db, policy.dd_rsd_mu4_Pk_table(), sample.get_dd_rsd_mu4(), model, sample);
        store_impl::store_one_loop_rsd_Pk(db, policy.dd_rsd_mu6_Pk_table(), sample.get_dd_rsd_mu6(), model, sample);
        store_impl::store_one_loop_rsd_Pk(db, policy.dd_rsd_mu8_Pk_table(), sample.get_dd_rsd_mu8(), model, sample);
      }
    
    
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
               const multipole_Pk& sample)
      {
        assert(db != nullptr);
        
        store_impl::store_multipole_Pk(db, policy.P0_table(), sample.get_P0(), model, sample);
        store_impl::store_multipole_Pk(db, policy.P2_table(), sample.get_P2(), model, sample);
        store_impl::store_multipole_Pk(db, policy.P4_table(), sample.get_P4(), model, sample);
      }
    
    
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
               const Matsubara_XY& sample)
      {
        assert(db != nullptr);
        
        std::ostringstream insert_stmt;
        insert_stmt
          << "INSERT INTO " << policy.Matsubara_XY_table() << " VALUES (@mid, @Pk_id, @IR_resum_id, @X, @Y);";
    
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));
    
        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), sample.get_Pk_token().get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_resum_id"), sample.get_IR_resum_token().get_id()));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@X"), store_impl::make_dimensionless(sample.get_X())));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Y"), store_impl::make_dimensionless(sample.get_Y())));
    
        // perform insertion
        check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_MATSUBARA_XY_FAIL, SQLITE_DONE);
    
        // clear bindings and release
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));
      }
    
    
    void store(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token&,
               const filtered_Pk& sample)
      {
        assert(db != nullptr);
        
        std::ostringstream insert_stmt;
        insert_stmt
          << "INSERT INTO " << policy.Pk_linear_table() << " VALUES (@Pk_id, @kid, @Pk_raw, @Pk_w, @Pk_ref);";
    
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));
    
        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), sample.get_Pk_token().get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), sample.get_k_token().get_id()));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_raw"), store_impl::make_dimensionless(sample.get_Pk_raw())));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_w"), store_impl::make_dimensionless(sample.get_Pk_w())));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_ref"), store_impl::make_dimensionless(sample.get_Pk_ref())));
    
        // perform insertion
        check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_MATSUBARA_XY_FAIL, SQLITE_DONE);
    
        // clear bindings and release
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));
      }
    
    
  }   // namespace sqlite3_operations
