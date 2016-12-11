//
// Created by David Seery on 10/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#include "find.h"
#include "utilities.h"

#include "exceptions.h"
#include "localizations/messages.h"

#include "temporary_tables.h"

namespace sqlite3_operations
  {
    
    namespace find_impl
      {
    
        template <typename KernelType>
        void read_loop_kernel(sqlite3* db, const std::string& table, const FRW_model_token& model, const k_token& k,
                              const linear_Pk_token& Pk, const UV_cutoff_token& UV_cutoff, KernelType& kernel,
                              const IR_cutoff_token& IR_cutoff)
          {
            std::ostringstream read_stmt;
            read_stmt
              << "SELECT raw_value, raw_regions, raw_evals, raw_err, raw_time, nw_value, nw_regions, nw_evals, nw_err, nw_time FROM "
              << table << " WHERE mid=@mid AND kid=@kid AND Pk_id=@Pk_id AND UV_id=@UV_id AND IR_id=@IR_id;";
            
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, read_stmt.str().c_str(), read_stmt.str().length()+1, &stmt, nullptr));
            
            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), k.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), Pk.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), UV_cutoff.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), IR_cutoff.get_id()));
            
            // perform read
            int result = 0;
            unsigned int count = 0;
            while((result = sqlite3_step(stmt)) != SQLITE_DONE)
              {
                if(result == SQLITE_ROW)
                  {
                    auto& raw = kernel.get_raw();
                    auto& nw = kernel.get_nowiggle();
                    
                    raw.value = sqlite3_column_double(stmt, 0) * dimensionful_unit<typename KernelType::value_type>();
                    raw.regions = sqlite3_column_int(stmt, 1);
                    raw.evaluations = sqlite3_column_int(stmt, 2);
                    raw.error = sqlite3_column_double(stmt, 3) * dimensionful_unit<typename KernelType::value_type>();
                    raw.time = sqlite3_column_int64(stmt, 4);
    
                    nw.value = sqlite3_column_double(stmt, 5) * dimensionful_unit<typename KernelType::value_type>();
                    nw.regions = sqlite3_column_int(stmt, 6);
                    nw.evaluations = sqlite3_column_int(stmt, 7);
                    nw.error = sqlite3_column_double(stmt, 8) * dimensionful_unit<typename KernelType::value_type>();
                    nw.time = sqlite3_column_int64(stmt, 9);
    
                    ++count;
                  }
                else
                  {
                    check_stmt(db, sqlite3_clear_bindings(stmt));
                    check_stmt(db, sqlite3_finalize(stmt));

                    throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_READ_LOOP_MOMENTUM_FAIL);
                  }
              }
            
            // clear bindings and release
            check_stmt(db, sqlite3_clear_bindings(stmt));
            check_stmt(db, sqlite3_finalize(stmt));
            
            if(count != 1) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_LOOP_MOMENTUM_MISREAD);
          }
        
        
        template <typename DataType>
        void read_Pk_value(sqlite3_stmt* stmt, unsigned int value_raw, unsigned int err_raw, unsigned int value_nw,
                           unsigned int err_wiggle, DataType& data)
          {
            using value_type = typename DataType::value_type;

            value_type v_raw = sqlite3_column_double(stmt, value_raw) * dimensionful_unit<value_type>();
            value_type e_raw = sqlite3_column_double(stmt, err_raw) * dimensionful_unit<value_type>();
    
            value_type v_nw = sqlite3_column_double(stmt, value_nw) * dimensionful_unit<value_type>();
            value_type e_nw = sqlite3_column_double(stmt, err_wiggle) * dimensionful_unit<value_type>();
            
            using container = typename DataType::container_type;
            
            data.set_raw(container(v_raw, e_raw));
            data.set_nowiggle(container(v_nw, e_nw));
          }
    
    
        template <typename DataType>
        void read_Pk_value(sqlite3_stmt* stmt, unsigned int value_raw, unsigned int value_nw, DataType& data)
          {
            using value_type = typename DataType::value_type;

            value_type v_raw = sqlite3_column_double(stmt, value_raw) * dimensionful_unit<value_type>();
            value_type v_nw = sqlite3_column_double(stmt, value_nw) * dimensionful_unit<value_type>();
    
            using container = typename DataType::container_type;
    
            data.set_raw(container(v_raw, value_type(0.0)));
            data.set_nowiggle(container(v_nw, value_type(0.0)));
          }
    
    
        void read_dd_Pk(sqlite3* db, const std::string& table, const FRW_model_token& model, const k_token& k,
                        const z_token& z, const linear_Pk_token& Pk_lin, const IR_cutoff_token& IR_cutoff,
                        const UV_cutoff_token& UV_cutoff, dd_Pk& Pk)
          {
            std::ostringstream read_stmt;
            read_stmt
              << "SELECT "
              << "Ptree_raw, err_tree_raw, P13_raw, err_13_raw, P22_raw, err_22_raw, P1loopSPT_raw, err_1loopSPT_raw, Z2_delta_raw, "
              << "Ptree_nw, err_tree_nw, P13_nw, err_13_nw, P22_nw, err_22_nw, P1loopSPT_nw, err_1loopSPT_nw, Z2_delta_nw "
              << "FROM " << table << " "
              << "WHERE mid=@mid AND zid=@zid AND kid=@kid AND Pk_id=@Pk_id AND IR_id=@IR_id AND UV_id=@UV_id;";
            
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, read_stmt.str().c_str(), read_stmt.str().length()+1, &stmt, nullptr));
            
            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), z.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), k.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), Pk_lin.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), IR_cutoff.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), UV_cutoff.get_id()));
            
            // perform read
            int result = 0;
            unsigned int count = 0;
            while((result = sqlite3_step(stmt)) != SQLITE_DONE)
              {
                if(result == SQLITE_ROW)
                  {
                    read_Pk_value(stmt, 0, 1, 9, 10, Pk.get_tree());
                    read_Pk_value(stmt, 2, 3, 11, 12, Pk.get_13());
                    read_Pk_value(stmt, 4, 5, 13, 14, Pk.get_22());
                    read_Pk_value(stmt, 6, 7, 15, 16, Pk.get_1loop_SPT());
                    read_Pk_value(stmt, 8, 17, Pk.get_Z2_delta());
    
                    ++count;
                  }
                else
                  {
                    check_stmt(db, sqlite3_clear_bindings(stmt));
                    check_stmt(db, sqlite3_finalize(stmt));
                    
                    throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_READ_PK_FAIL);
                  }
              }
            
            // clear bindings and release
            check_stmt(db, sqlite3_clear_bindings(stmt));
            check_stmt(db, sqlite3_finalize(stmt));
            
            if(count != 1) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_READ_PK_MISREAD);
          }
    
    
        void read_dd_rsd_Pk(sqlite3* db, const std::string& table, const FRW_model_token& model, const k_token& k,
                            const z_token& z, const linear_Pk_token& Pk_lin, const IR_cutoff_token& IR_cutoff,
                            const UV_cutoff_token& UV_cutoff, rsd_dd_Pk& Pk)
          {
            std::ostringstream read_stmt;
            read_stmt
              << "SELECT "
              << "Ptree_raw, err_tree_raw, P13_raw, err_13_raw, P22_raw, err_22_raw, P1loopSPT_raw, err_1loopSPT_raw, Z2_delta_raw, Z0_v_raw, Z2_v_raw, Z0_vdelta_raw, Z2_vdelta_raw, Z2_vv_raw, Z2_vvdelta_raw, Z2_vvv_raw, "
              << "Ptree_nw, err_tree_nw, P13_nw, err_13_nw, P22_nw, err_22_nw, P1loopSPT_nw, err_1loopSPT_nw, Z2_delta_nw, Z0_v_nw, Z2_v_nw, Z0_vdelta_nw, Z2_vdelta_nw, Z2_vv_nw, Z2_vvdelta_nw, Z2_vvv_nw "
              << "FROM " << table << " "
              << "WHERE mid=@mid AND zid=@zid AND kid=@kid AND Pk_id=@Pk_id AND IR_id=@IR_id AND UV_id=@UV_id;";
    
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, read_stmt.str().c_str(), read_stmt.str().length()+1, &stmt, nullptr));
    
            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), z.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), k.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), Pk_lin.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), IR_cutoff.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), UV_cutoff.get_id()));
    
            // perform read
            int result = 0;
            unsigned int count = 0;
            while((result = sqlite3_step(stmt)) != SQLITE_DONE)
              {
                if(result == SQLITE_ROW)
                  {
                    read_Pk_value(stmt, 0, 1, 16, 17, Pk.get_tree());
                    read_Pk_value(stmt, 2, 3, 18, 19, Pk.get_13());
                    read_Pk_value(stmt, 4, 5, 20, 21, Pk.get_22());
                    read_Pk_value(stmt, 6, 7, 22, 23, Pk.get_1loop_SPT());
                    read_Pk_value(stmt, 8, 24, Pk.get_Z2_delta());
                    read_Pk_value(stmt, 9, 25, Pk.get_Z0_v());
                    read_Pk_value(stmt, 10, 26, Pk.get_Z2_v());
                    read_Pk_value(stmt, 11, 27, Pk.get_Z0_vdelta());
                    read_Pk_value(stmt, 12, 28, Pk.get_Z2_vdelta());
                    read_Pk_value(stmt, 13, 29, Pk.get_Z2_vv());
                    read_Pk_value(stmt, 14, 30, Pk.get_Z2_vvdelta());
                    read_Pk_value(stmt, 15, 31, Pk.get_Z2_vvv());
            
                    ++count;
                  }
                else
                  {
                    check_stmt(db, sqlite3_clear_bindings(stmt));
                    check_stmt(db, sqlite3_finalize(stmt));
            
                    throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_READ_RSD_PK_FAIL);
                  }
              }
    
            // clear bindings and release
            check_stmt(db, sqlite3_clear_bindings(stmt));
            check_stmt(db, sqlite3_finalize(stmt));
    
            if(count != 1) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_READ_RSD_PK_MISREAD);
          }
          
      }   // namespace find_impl
    
    std::unique_ptr<oneloop_growth> find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                                         const FRW_model_token& model, const z_database& z_db)
      {
        // set up temporary table of desired z identifiers
        std::string ztab = z_table(db, mgr, policy, z_db);
        
        std::unique_ptr<oneloop_growth> payload = std::make_unique<oneloop_growth>(z_db);

        // note use of ORDER BY ... DESC which is needed to get the g, f values
        // the correct order
        // z_table() will produce a table of tokens in ascending redshift order,
        // ie. ordered into the past
        
        std::ostringstream g_select_stmt;
        g_select_stmt
          << "SELECT "
          << "g_sample.g_linear, "
          << "g_sample.A, "
          << "g_sample.B, "
          << "g_sample.D, "
          << "g_sample.E, "
          << "g_sample.F, "
          << "g_sample.G, "
          << "g_sample.J "
          << "FROM " << ztab << " "
          << "INNER JOIN (SELECT * FROM " << policy.g_factor_table() << " WHERE mid=@mid) g_sample "
          << "ON " << ztab << ".id = g_sample.zid "
          << "ORDER BY " << ztab << ".ROWID DESC;";
        
        std::ostringstream f_select_stmt;
        f_select_stmt
          << "SELECT "
          << "f_sample.f_linear, "
          << "f_sample.fA, "
          << "f_sample.fB, "
          << "f_sample.fD, "
          << "f_sample.fE, "
          << "f_sample.fF, "
          << "f_sample.fG, "
          << "f_sample.fJ "
          << "FROM " << ztab << " "
          << "INNER JOIN (SELECT * FROM " << policy.f_factor_table() << " WHERE mid=@mid) f_sample "
          << "ON " << ztab << ".id = f_sample.zid "
          << "ORDER BY " << ztab << ".ROWID DESC;";
        
        // prepare statements
        sqlite3_stmt* g_stmt;
        check_stmt(db, sqlite3_prepare_v2(db, g_select_stmt.str().c_str(), g_select_stmt.str().length()+1, &g_stmt, nullptr));
    
        sqlite3_stmt* f_stmt;
        check_stmt(db, sqlite3_prepare_v2(db, f_select_stmt.str().c_str(), f_select_stmt.str().length()+1, &f_stmt, nullptr));
    
        check_stmt(db, sqlite3_bind_int(g_stmt, sqlite3_bind_parameter_index(g_stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(f_stmt, sqlite3_bind_parameter_index(f_stmt, "@mid"), model.get_id()));
        
        int g_status = 0;
        int f_status = 0;
        
        unsigned int read_count = 0;
        
        while((g_status = sqlite3_step(g_stmt)) != SQLITE_DONE && (f_status = sqlite3_step(f_stmt)) != SQLITE_DONE)
          {
            if(g_status == SQLITE_ROW && f_status == SQLITE_ROW)
              {
                payload->push_back(sqlite3_column_double(g_stmt, 0),
                                   sqlite3_column_double(g_stmt, 1),
                                   sqlite3_column_double(g_stmt, 2),
                                   sqlite3_column_double(g_stmt, 3),
                                   sqlite3_column_double(g_stmt, 4),
                                   sqlite3_column_double(g_stmt, 5),
                                   sqlite3_column_double(g_stmt, 6),
                                   sqlite3_column_double(g_stmt, 7),
                                   sqlite3_column_double(f_stmt, 0),
                                   sqlite3_column_double(f_stmt, 1),
                                   sqlite3_column_double(f_stmt, 2),
                                   sqlite3_column_double(f_stmt, 3),
                                   sqlite3_column_double(f_stmt, 4),
                                   sqlite3_column_double(f_stmt, 5),
                                   sqlite3_column_double(f_stmt, 6),
                                   sqlite3_column_double(f_stmt, 7));
                
                ++read_count;
              }
            else
              {
                std::ostringstream msg;
                msg << ERROR_SQLITE3_FG_GROWTH_TABLE_READ_FAIL << "(" << g_status << "," << f_status << "): " << sqlite3_errmsg(db) << "]";

                check_stmt(db, sqlite3_finalize(g_stmt));
                check_stmt(db, sqlite3_finalize(f_stmt));

                throw runtime_exception(exception_type::database_error, msg.str());
              }
          }
        
        check_stmt(db, sqlite3_finalize(g_stmt));
        check_stmt(db, sqlite3_finalize(f_stmt));
    
        // drop unneeded temporary tables
        drop_temp(db, mgr, ztab);
    
        if(read_count != z_db.size()) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_FG_GROWTH_MISREAD);
        
        return std::move(payload);
      }
    
    
    std::unique_ptr<loop_integral>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
         const k_token& k, const linear_Pk_token& Pk, const IR_cutoff_token& IR_cutoff,
         const UV_cutoff_token& UV_cutoff)
      {
        delta_22_integrals delta22;
        delta_13_integrals delta13;
        rsd_22_integrals rsd22;
        rsd_13_integrals rsd13;
    
        find_impl::read_loop_kernel(db, policy.AA_table(), model, k, Pk, UV_cutoff, delta22.get_AA(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.AB_table(), model, k, Pk, UV_cutoff, delta22.get_AB(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.BB_table(), model, k, Pk, UV_cutoff, delta22.get_BB(), IR_cutoff);
    
        find_impl::read_loop_kernel(db, policy.D_table(), model, k, Pk, UV_cutoff, delta13.get_D(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.E_table(), model, k, Pk, UV_cutoff, delta13.get_E(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.F_table(), model, k, Pk, UV_cutoff, delta13.get_F(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.G_table(), model, k, Pk, UV_cutoff, delta13.get_G(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.J1_table(), model, k, Pk, UV_cutoff, delta13.get_J1(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.J2_table(), model, k, Pk, UV_cutoff, delta13.get_J2(), IR_cutoff);
    
        find_impl::read_loop_kernel(db, policy.RSD13_a_table(), model, k, Pk, UV_cutoff, rsd13.get_a(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD13_b_table(), model, k, Pk, UV_cutoff, rsd13.get_b(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD13_c_table(), model, k, Pk, UV_cutoff, rsd13.get_c(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD13_d_table(), model, k, Pk, UV_cutoff, rsd13.get_d(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD13_e_table(), model, k, Pk, UV_cutoff, rsd13.get_e(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD13_f_table(), model, k, Pk, UV_cutoff, rsd13.get_f(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD13_g_table(), model, k, Pk, UV_cutoff, rsd13.get_g(), IR_cutoff);
    
        find_impl::read_loop_kernel(db, policy.RSD22_A1_table(), model, k, Pk, UV_cutoff, rsd22.get_A1(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_A2_table(), model, k, Pk, UV_cutoff, rsd22.get_A2(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_A3_table(), model, k, Pk, UV_cutoff, rsd22.get_A3(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_A4_table(), model, k, Pk, UV_cutoff, rsd22.get_A4(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_A5_table(), model, k, Pk, UV_cutoff, rsd22.get_A5(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_B2_table(), model, k, Pk, UV_cutoff, rsd22.get_B2(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_B3_table(), model, k, Pk, UV_cutoff, rsd22.get_B3(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_B6_table(), model, k, Pk, UV_cutoff, rsd22.get_B6(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_B8_table(), model, k, Pk, UV_cutoff, rsd22.get_B8(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_B9_table(), model, k, Pk, UV_cutoff, rsd22.get_B9(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_C1_table(), model, k, Pk, UV_cutoff, rsd22.get_C1(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_C2_table(), model, k, Pk, UV_cutoff, rsd22.get_C2(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_C4_table(), model, k, Pk, UV_cutoff, rsd22.get_C4(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_D1_table(), model, k, Pk, UV_cutoff, rsd22.get_D1(), IR_cutoff);
    
        std::unique_ptr<loop_integral> payload = std::make_unique<loop_integral>(k, Pk, UV_cutoff, IR_cutoff, delta22, delta13, rsd22, rsd13);
        
        return std::move(payload);
      }
    
    
    std::unique_ptr<oneloop_Pk>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
         const k_token& k, const z_token& z, const linear_Pk_token& Pk_lin, const IR_cutoff_token& IR_cutoff,
         const UV_cutoff_token& UV_cutoff)
      {
        dd_Pk dd;
        rsd_dd_Pk rsd_mu0;
        rsd_dd_Pk rsd_mu2;
        rsd_dd_Pk rsd_mu4;
        rsd_dd_Pk rsd_mu6;
        rsd_dd_Pk rsd_mu8;
    
        find_impl::read_dd_Pk(db, policy.dd_Pk_table(), model, k, z, Pk_lin, IR_cutoff, UV_cutoff, dd);
        find_impl::read_dd_rsd_Pk(db, policy.dd_rsd_mu0_Pk_table(), model, k, z, Pk_lin, IR_cutoff, UV_cutoff, rsd_mu0);
        find_impl::read_dd_rsd_Pk(db, policy.dd_rsd_mu2_Pk_table(), model, k, z, Pk_lin, IR_cutoff, UV_cutoff, rsd_mu2);
        find_impl::read_dd_rsd_Pk(db, policy.dd_rsd_mu4_Pk_table(), model, k, z, Pk_lin, IR_cutoff, UV_cutoff, rsd_mu4);
        find_impl::read_dd_rsd_Pk(db, policy.dd_rsd_mu6_Pk_table(), model, k, z, Pk_lin, IR_cutoff, UV_cutoff, rsd_mu6);
        find_impl::read_dd_rsd_Pk(db, policy.dd_rsd_mu8_Pk_table(), model, k, z, Pk_lin, IR_cutoff, UV_cutoff, rsd_mu8);
        
        std::unique_ptr<oneloop_Pk> payload = std::make_unique<oneloop_Pk>(k, Pk_lin, IR_cutoff, UV_cutoff, z, dd, rsd_mu0, rsd_mu2, rsd_mu4, rsd_mu6, rsd_mu8);
        
        return std::move(payload);
      }
    
    
    std::unique_ptr<Matsubara_XY>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
         const linear_Pk_token& Pk, const IR_resum_token& IR_resum)
      {
        std::ostringstream read_stmt;
        read_stmt << "SELECT X, Y FROM " << policy.Matsubara_XY_table() << " WHERE mid=@mid AND Pk_id=@Pk_id AND IR_resum_id=@IR_resum_id;";
    
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, read_stmt.str().c_str(), read_stmt.str().length()+1, &stmt, nullptr));
    
        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), Pk.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_resum_id"), IR_resum.get_id()));
    
        // store value
        Mpc_units::inverse_energy2 X(0.0);
        Mpc_units::inverse_energy2 Y(0.0);
        
        // perform read
        int result = 0;
        unsigned int count = 0;
        while((result = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(result == SQLITE_ROW)
              {
                X = sqlite3_column_double(stmt, 0) * dimensionful_unit<Mpc_units::inverse_energy2>();
                Y = sqlite3_column_double(stmt, 1) * dimensionful_unit<Mpc_units::inverse_energy2>();
            
                ++count;
              }
            else
              {
                check_stmt(db, sqlite3_clear_bindings(stmt));
                check_stmt(db, sqlite3_finalize(stmt));
            
                throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_READ_MATSUBARA_XY_FAIL);
              }
          }
    
        // clear bindings and release
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));
    
        if(count != 1) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_MATSUBARA_XY_MISREAD);
        
        std::unique_ptr<Matsubara_XY> payload = std::make_unique<Matsubara_XY>(Pk, IR_resum, X, Y);
        
        return std::move(payload);
      }
    
    
    std::unique_ptr<wiggle_Pk>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const linear_Pk_token& token,
         const k_database& k_db)
      {
        // set up temporary table of desired k identifiers
        std::string ktab = k_table(db, mgr, policy, k_db);
        
        std::ostringstream read_stmt;
        read_stmt << "SELECT sample.Pk_raw, sample.Pk_nw "
                  << "FROM " << ktab << " "
                  << "INNER JOIN (SELECT * FROM " << policy.Pk_linear_table() << " WHERE Pk_id=@Pk_id) sample "
                  << "ON " << ktab << ".id = sample.kid "
                  << "ORDER BY " << ktab << ".ROWID ASC;";
    
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, read_stmt.str().c_str(), read_stmt.str().length()+1, &stmt, nullptr));
    
        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), token.get_id()));
        
        // set up databases to hold the result
        tree_Pk::database_type raw_db;
        tree_Pk_w::database_type nw_db;
    
        // perform read
        int result = 0;
        k_database::const_record_iterator t = k_db.record_cbegin();
        while((result = sqlite3_step(stmt)) != SQLITE_DONE && t != k_db.record_cend())
          {
            if(result == SQLITE_ROW)
              {
                Mpc_units::inverse_energy3 raw = sqlite3_column_double(stmt, 0) * dimensionful_unit<Mpc_units::inverse_energy3>();
                Mpc_units::inverse_energy3 nw  = sqlite3_column_double(stmt, 1) * dimensionful_unit<Mpc_units::inverse_energy3>();
                
                raw_db.add_record(*(*t), raw);
                nw_db.add_record(*(*t), nw);
                
                ++t;
              }
            else
              {
                check_stmt(db, sqlite3_clear_bindings(stmt));
                check_stmt(db, sqlite3_finalize(stmt));
            
                throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_READ_FILTERED_PK_FAIL);
              }
          }
    
        // clear bindings and release
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));
        
        // drop temporary table
        drop_temp(db, mgr, ktab);
    
        if(t != k_db.record_cend()) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_FILTERED_PK_MISREAD);
        
        std::unique_ptr<wiggle_Pk> payload = std::make_unique<wiggle_Pk>(token, nw_db, raw_db);
        
        return std::move(payload);
      }
    
  }   // namespace sqlite3_operations
