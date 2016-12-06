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
                              const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff, KernelType& kernel)
          {
            std::ostringstream read_stmt;
            read_stmt << "SELECT value, regions, evals, err, time FROM " << table << " WHERE mid=@mid AND kid=@kid AND UV_id=@UV_id AND IR_id=@IR_id;";
            
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, read_stmt.str().c_str(), read_stmt.str().length()+1, &stmt, nullptr));
            
            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), k.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), UV_cutoff.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), IR_cutoff.get_id()));
            
            // perform read
            int result = 0;
            unsigned int count = 0;
            while((result = sqlite3_step(stmt)) != SQLITE_DONE)
              {
                if(result == SQLITE_ROW)
                  {
                    kernel.value = sqlite3_column_double(stmt, 0) * dimensionful_unit<typename KernelType::value_type>();
                    kernel.regions = sqlite3_column_int(stmt, 1);
                    kernel.evaluations = sqlite3_column_int(stmt, 2);
                    kernel.error = sqlite3_column_double(stmt, 3) * dimensionful_unit<typename KernelType::value_type>();
                    kernel.time = sqlite3_column_int64(stmt, 4);
                    
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
            
            if(count > 1) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_LOOP_MOMENTUM_MISREAD);
          }
        
        
        template <typename DataType>
        void read_Pk_value(sqlite3_stmt* stmt, unsigned int value, unsigned int err, DataType& data)
          {
            data.value = sqlite3_column_double(stmt, value) * dimensionful_unit<typename DataType::value_type>();
            data.error = sqlite3_column_double(stmt, err) * dimensionful_unit<typename DataType::value_type>();
          }
    
    
        template <typename DataType>
        void read_Pk_value(sqlite3_stmt* stmt, unsigned int value, DataType& data)
          {
            data.value = sqlite3_column_double(stmt, value) * dimensionful_unit<typename DataType::value_type>();
            data.error = typename DataType::value_type(0.0);
          }
        
        
        void read_dd_Pk(sqlite3* db, const std::string& table, const FRW_model_token& model,
                        const k_token& k, const z_token& z,
                        const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff, dd_Pk& Pk)
          {
            std::ostringstream read_stmt;
            read_stmt
              << "SELECT "
              << "Ptree, err_tree, P13, err_13, P22, err_22, P1loopSPT, err_1loopSPT, Z2_delta "
              << "FROM " << table << " "
              << "WHERE mid=@mid AND zid=@zid AND kid=@kid AND IR_id=@IR_id AND UV_id=@UV_id;";
            
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, read_stmt.str().c_str(), read_stmt.str().length()+1, &stmt, nullptr));
            
            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), z.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), k.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), IR_cutoff.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), UV_cutoff.get_id()));
            
            // perform read
            int result = 0;
            unsigned int count = 0;
            while((result = sqlite3_step(stmt)) != SQLITE_DONE)
              {
                if(result == SQLITE_ROW)
                  {
                    read_Pk_value(stmt, 0, 1, Pk.get_tree());
                    read_Pk_value(stmt, 2, 3, Pk.get_13());
                    read_Pk_value(stmt, 4, 5, Pk.get_22());
                    read_Pk_value(stmt, 6, 7, Pk.get_1loop_SPT());
                    read_Pk_value(stmt, 8, Pk.get_Z2_delta());
                    
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
            
            if(count > 1) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_READ_PK_MISREAD);
          }
    
    
        void read_dd_rsd_Pk(sqlite3* db, const std::string& table, const FRW_model_token& model,
                            const k_token& k, const z_token& z,
                            const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff, rsd_dd_Pk& Pk)
          {
            std::ostringstream read_stmt;
            read_stmt
              << "SELECT "
              << "Ptree, err_tree, P13, err_13, P22, err_22, P1loopSPT, err_1loopSPT, "
              << "Z2_delta, Z0_v, Z2_v, Z0_vdelta, Z2_vdelta, Z2_vv, Z2_vvdelta, Z2_vvv "
              << "FROM " << table << " "
              << "WHERE mid=@mid AND zid=@zid AND kid=@kid AND IR_id=@IR_id AND UV_id=@UV_id;";
    
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, read_stmt.str().c_str(), read_stmt.str().length()+1, &stmt, nullptr));
    
            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), z.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), k.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), IR_cutoff.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), UV_cutoff.get_id()));
    
            // perform read
            int result = 0;
            unsigned int count = 0;
            while((result = sqlite3_step(stmt)) != SQLITE_DONE)
              {
                if(result == SQLITE_ROW)
                  {
                    read_Pk_value(stmt, 0, 1, Pk.get_tree());
                    read_Pk_value(stmt, 2, 3, Pk.get_13());
                    read_Pk_value(stmt, 4, 5, Pk.get_22());
                    read_Pk_value(stmt, 6, 7, Pk.get_1loop_SPT());
                    read_Pk_value(stmt, 8, Pk.get_Z2_delta());
                    read_Pk_value(stmt, 9, Pk.get_Z0_v());
                    read_Pk_value(stmt, 10, Pk.get_Z2_v());
                    read_Pk_value(stmt, 11, Pk.get_Z0_vdelta());
                    read_Pk_value(stmt, 12, Pk.get_Z2_vdelta());
                    read_Pk_value(stmt, 13, Pk.get_Z2_vv());
                    read_Pk_value(stmt, 14, Pk.get_Z2_vvdelta());
                    read_Pk_value(stmt, 15, Pk.get_Z2_vvv());
            
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
    
            if(count > 1) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_READ_RSD_PK_MISREAD);
          }
          
      }   // namespace find_impl
    
    oneloop_growth find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                        const FRW_model_token& model, const z_database& z_db)
      {
        // set up temporary table of desired z identifiers
        std::string ztab = z_table(db, mgr, policy, z_db);
        
        oneloop_growth payload(z_db);

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
                payload.push_back(sqlite3_column_double(g_stmt, 0),
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
        
        if(read_count != z_db.size()) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_FG_GROWTH_MISREAD);
        
        check_stmt(db, sqlite3_finalize(g_stmt));
        check_stmt(db, sqlite3_finalize(f_stmt));
    
        // drop unneeded temporary tables
        drop_temp(db, mgr, ztab);
        
        return std::move(payload);
      }
    
    
    loop_integral find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                       const FRW_model_token& model, const k_token& k,
                       const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff)
      {
        delta_22_integrals delta22;
        delta_13_integrals delta13;
        rsd_22_integrals rsd22;
        rsd_13_integrals rsd13;
    
        find_impl::read_loop_kernel(db, policy.AA_table(), model, k, IR_cutoff, UV_cutoff, delta22.get_AA());
        find_impl::read_loop_kernel(db, policy.AB_table(), model, k, IR_cutoff, UV_cutoff, delta22.get_AB());
        find_impl::read_loop_kernel(db, policy.BB_table(), model, k, IR_cutoff, UV_cutoff, delta22.get_BB());
    
        find_impl::read_loop_kernel(db, policy.D_table(), model, k, IR_cutoff, UV_cutoff, delta13.get_D());
        find_impl::read_loop_kernel(db, policy.E_table(), model, k, IR_cutoff, UV_cutoff, delta13.get_E());
        find_impl::read_loop_kernel(db, policy.F_table(), model, k, IR_cutoff, UV_cutoff, delta13.get_F());
        find_impl::read_loop_kernel(db, policy.G_table(), model, k, IR_cutoff, UV_cutoff, delta13.get_G());
        find_impl::read_loop_kernel(db, policy.J1_table(), model, k, IR_cutoff, UV_cutoff, delta13.get_J1());
        find_impl::read_loop_kernel(db, policy.J2_table(), model, k, IR_cutoff, UV_cutoff, delta13.get_J2());
    
        find_impl::read_loop_kernel(db, policy.RSD13_a_table(), model, k, IR_cutoff, UV_cutoff, rsd13.get_a());
        find_impl::read_loop_kernel(db, policy.RSD13_b_table(), model, k, IR_cutoff, UV_cutoff, rsd13.get_b());
        find_impl::read_loop_kernel(db, policy.RSD13_c_table(), model, k, IR_cutoff, UV_cutoff, rsd13.get_c());
        find_impl::read_loop_kernel(db, policy.RSD13_d_table(), model, k, IR_cutoff, UV_cutoff, rsd13.get_d());
        find_impl::read_loop_kernel(db, policy.RSD13_e_table(), model, k, IR_cutoff, UV_cutoff, rsd13.get_e());
        find_impl::read_loop_kernel(db, policy.RSD13_f_table(), model, k, IR_cutoff, UV_cutoff, rsd13.get_f());
        find_impl::read_loop_kernel(db, policy.RSD13_g_table(), model, k, IR_cutoff, UV_cutoff, rsd13.get_g());
    
        find_impl::read_loop_kernel(db, policy.RSD22_A1_table(), model, k, IR_cutoff, UV_cutoff, rsd22.get_A1());
        find_impl::read_loop_kernel(db, policy.RSD22_A2_table(), model, k, IR_cutoff, UV_cutoff, rsd22.get_A2());
        find_impl::read_loop_kernel(db, policy.RSD22_A3_table(), model, k, IR_cutoff, UV_cutoff, rsd22.get_A3());
        find_impl::read_loop_kernel(db, policy.RSD22_A4_table(), model, k, IR_cutoff, UV_cutoff, rsd22.get_A4());
        find_impl::read_loop_kernel(db, policy.RSD22_A5_table(), model, k, IR_cutoff, UV_cutoff, rsd22.get_A5());
        find_impl::read_loop_kernel(db, policy.RSD22_B2_table(), model, k, IR_cutoff, UV_cutoff, rsd22.get_B2());
        find_impl::read_loop_kernel(db, policy.RSD22_B3_table(), model, k, IR_cutoff, UV_cutoff, rsd22.get_B3());
        find_impl::read_loop_kernel(db, policy.RSD22_B6_table(), model, k, IR_cutoff, UV_cutoff, rsd22.get_B6());
        find_impl::read_loop_kernel(db, policy.RSD22_B8_table(), model, k, IR_cutoff, UV_cutoff, rsd22.get_B8());
        find_impl::read_loop_kernel(db, policy.RSD22_B9_table(), model, k, IR_cutoff, UV_cutoff, rsd22.get_B9());
        find_impl::read_loop_kernel(db, policy.RSD22_C1_table(), model, k, IR_cutoff, UV_cutoff, rsd22.get_C1());
        find_impl::read_loop_kernel(db, policy.RSD22_C2_table(), model, k, IR_cutoff, UV_cutoff, rsd22.get_C2());
        find_impl::read_loop_kernel(db, policy.RSD22_C4_table(), model, k, IR_cutoff, UV_cutoff, rsd22.get_C4());
        find_impl::read_loop_kernel(db, policy.RSD22_D1_table(), model, k, IR_cutoff, UV_cutoff, rsd22.get_D1());
    
        loop_integral payload(k, UV_cutoff, IR_cutoff, delta22, delta13, rsd22, rsd13);
        
        return std::move(payload);
      }
    
    
    oneloop_Pk find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                    const FRW_model_token& model, const k_token& k, const z_token& z,
                    const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff)
      {
        dd_Pk dd;
        rsd_dd_Pk rsd_mu0;
        rsd_dd_Pk rsd_mu2;
        rsd_dd_Pk rsd_mu4;
        rsd_dd_Pk rsd_mu6;
        rsd_dd_Pk rsd_mu8;
        
        find_impl::read_dd_Pk(db, policy.dd_Pk_table(), model, k, z, IR_cutoff, UV_cutoff, dd);
        find_impl::read_dd_rsd_Pk(db, policy.dd_rsd_mu0_Pk_table(), model, k, z, IR_cutoff, UV_cutoff, rsd_mu0);
        find_impl::read_dd_rsd_Pk(db, policy.dd_rsd_mu2_Pk_table(), model, k, z, IR_cutoff, UV_cutoff, rsd_mu2);
        find_impl::read_dd_rsd_Pk(db, policy.dd_rsd_mu4_Pk_table(), model, k, z, IR_cutoff, UV_cutoff, rsd_mu4);
        find_impl::read_dd_rsd_Pk(db, policy.dd_rsd_mu6_Pk_table(), model, k, z, IR_cutoff, UV_cutoff, rsd_mu6);
        find_impl::read_dd_rsd_Pk(db, policy.dd_rsd_mu8_Pk_table(), model, k, z, IR_cutoff, UV_cutoff, rsd_mu8);
        
        oneloop_Pk payload(k, IR_cutoff, UV_cutoff, z, dd, rsd_mu0, rsd_mu2, rsd_mu4, rsd_mu6, rsd_mu8);
        
        return std::move(payload);
      }
    
    
    Matsubara_A
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
         const IR_resum_token& IR_resum)
      {
        std::ostringstream read_stmt;
        read_stmt << "SELECT A FROM " << policy.Matsubara_A_table() << " WHERE mid=@mid AND IR_resum_id=@IR_resum_id;";
    
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, read_stmt.str().c_str(), read_stmt.str().length()+1, &stmt, nullptr));
    
        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_resum_id"), IR_resum.get_id()));
    
        // store value
        Mpc_units::inverse_energy2 value(0.0);
        
        // perform read
        int result = 0;
        unsigned int count = 0;
        while((result = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(result == SQLITE_ROW)
              {
                value = sqlite3_column_double(stmt, 0) * dimensionful_unit<Mpc_units::inverse_energy2>();
            
                ++count;
              }
            else
              {
                check_stmt(db, sqlite3_clear_bindings(stmt));
                check_stmt(db, sqlite3_finalize(stmt));
            
                throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_READ_MATSUBARA_A_FAIL);
              }
          }
    
        // clear bindings and release
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));
    
        if(count > 1) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_MATSUBARA_A_MISREAD);
        
        Matsubara_A payload(IR_resum, value);
        
        return std::move(payload);
      }
    
  }   // namespace sqlite3_operations
