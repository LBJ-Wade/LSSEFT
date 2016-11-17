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
        void read_loop_kernel(sqlite3* db, const std::string& table, const FRW_model_token& model,
                              const k_token& k, const UV_cutoff_token& UV_cutoff, const IR_cutoff_token& IR_cutoff,
                              KernelType& kernel)
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
          
      }   // namespace find_impl
    
    oneloop_growth find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy,
                        const FRW_model_token& model, z_database& z_db)
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
        
        std::ostringstream meta_select_stmt;
        meta_select_stmt
          << "SELECT time, steps FROM " << policy.gf_metadata_table() << " WHERE mid=@mid;";
        
        // prepare statements
        sqlite3_stmt* g_stmt;
        check_stmt(db, sqlite3_prepare_v2(db, g_select_stmt.str().c_str(), g_select_stmt.str().length()+1, &g_stmt, nullptr));
    
        sqlite3_stmt* f_stmt;
        check_stmt(db, sqlite3_prepare_v2(db, f_select_stmt.str().c_str(), f_select_stmt.str().length()+1, &f_stmt, nullptr));
        
        sqlite3_stmt* meta_stmt;
        check_stmt(db, sqlite3_prepare_v2(db, meta_select_stmt.str().c_str(), meta_select_stmt.str().length()+1, &meta_stmt, nullptr));
    
        check_stmt(db, sqlite3_bind_int(g_stmt, sqlite3_bind_parameter_index(g_stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(f_stmt, sqlite3_bind_parameter_index(f_stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(meta_stmt, sqlite3_bind_parameter_index(meta_stmt, "@mid"), model.get_id()));
        
        int status = 0;
        if((status = sqlite3_step(meta_stmt)) != SQLITE_DONE)
          {
            if(status == SQLITE_ROW)
              {
                payload.set_integration_metadata(sqlite3_column_int64(meta_stmt, 0),
                                                 sqlite3_column_int(meta_stmt, 1));
              }
            else
              {
                std::ostringstream msg;
                msg << ERROR_SQLITE3_FG_GROWTH_META_READ_FAIL << status << ": " << sqlite3_errmsg(db) << "]";
    
                check_stmt(db, sqlite3_finalize(g_stmt));
                check_stmt(db, sqlite3_finalize(f_stmt));
                check_stmt(db, sqlite3_finalize(meta_stmt));
    
                throw runtime_exception(exception_type::database_error, msg.str());
              }
          }
        
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
                check_stmt(db, sqlite3_finalize(meta_stmt));

                throw runtime_exception(exception_type::database_error, msg.str());
              }
          }
        
        if(read_count != z_db.size()) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_FG_GROWTH_MISREAD);
        
        check_stmt(db, sqlite3_finalize(g_stmt));
        check_stmt(db, sqlite3_finalize(f_stmt));
        check_stmt(db, sqlite3_finalize(meta_stmt));
    
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
        
        find_impl::read_loop_kernel(db, policy.AA_table(), model, k, UV_cutoff, IR_cutoff, delta22.get_AA());
        find_impl::read_loop_kernel(db, policy.AB_table(), model, k, UV_cutoff, IR_cutoff, delta22.get_AB());
        find_impl::read_loop_kernel(db, policy.BB_table(), model, k, UV_cutoff, IR_cutoff, delta22.get_BB());

        find_impl::read_loop_kernel(db, policy.D_table(), model, k, UV_cutoff, IR_cutoff, delta13.get_D());
        find_impl::read_loop_kernel(db, policy.E_table(), model, k, UV_cutoff, IR_cutoff, delta13.get_E());
        find_impl::read_loop_kernel(db, policy.F_table(), model, k, UV_cutoff, IR_cutoff, delta13.get_F());
        find_impl::read_loop_kernel(db, policy.G_table(), model, k, UV_cutoff, IR_cutoff, delta13.get_G());
        find_impl::read_loop_kernel(db, policy.J1_table(), model, k, UV_cutoff, IR_cutoff, delta13.get_J1());
        find_impl::read_loop_kernel(db, policy.J2_table(), model, k, UV_cutoff, IR_cutoff, delta13.get_J2());
        
        find_impl::read_loop_kernel(db, policy.RSD13_a_table(), model, k, UV_cutoff, IR_cutoff, rsd13.get_a());
        find_impl::read_loop_kernel(db, policy.RSD13_b_table(), model, k, UV_cutoff, IR_cutoff, rsd13.get_b());
        find_impl::read_loop_kernel(db, policy.RSD13_c_table(), model, k, UV_cutoff, IR_cutoff, rsd13.get_c());
        find_impl::read_loop_kernel(db, policy.RSD13_d_table(), model, k, UV_cutoff, IR_cutoff, rsd13.get_d());
        find_impl::read_loop_kernel(db, policy.RSD13_e_table(), model, k, UV_cutoff, IR_cutoff, rsd13.get_e());
        find_impl::read_loop_kernel(db, policy.RSD13_f_table(), model, k, UV_cutoff, IR_cutoff, rsd13.get_f());
        find_impl::read_loop_kernel(db, policy.RSD13_g_table(), model, k, UV_cutoff, IR_cutoff, rsd13.get_g());
    
        find_impl::read_loop_kernel(db, policy.RSD22_A1_table(), model, k, UV_cutoff, IR_cutoff, rsd22.get_A1());
        find_impl::read_loop_kernel(db, policy.RSD22_A2_table(), model, k, UV_cutoff, IR_cutoff, rsd22.get_A2());
        find_impl::read_loop_kernel(db, policy.RSD22_A3_table(), model, k, UV_cutoff, IR_cutoff, rsd22.get_A3());
        find_impl::read_loop_kernel(db, policy.RSD22_A4_table(), model, k, UV_cutoff, IR_cutoff, rsd22.get_A4());
        find_impl::read_loop_kernel(db, policy.RSD22_A5_table(), model, k, UV_cutoff, IR_cutoff, rsd22.get_A5());
        find_impl::read_loop_kernel(db, policy.RSD22_B2_table(), model, k, UV_cutoff, IR_cutoff, rsd22.get_B2());
        find_impl::read_loop_kernel(db, policy.RSD22_B3_table(), model, k, UV_cutoff, IR_cutoff, rsd22.get_B3());
        find_impl::read_loop_kernel(db, policy.RSD22_B6_table(), model, k, UV_cutoff, IR_cutoff, rsd22.get_B6());
        find_impl::read_loop_kernel(db, policy.RSD22_B8_table(), model, k, UV_cutoff, IR_cutoff, rsd22.get_B8());
        find_impl::read_loop_kernel(db, policy.RSD22_B9_table(), model, k, UV_cutoff, IR_cutoff, rsd22.get_B9());
        find_impl::read_loop_kernel(db, policy.RSD22_C1_table(), model, k, UV_cutoff, IR_cutoff, rsd22.get_C1());
        find_impl::read_loop_kernel(db, policy.RSD22_C2_table(), model, k, UV_cutoff, IR_cutoff, rsd22.get_C2());
        find_impl::read_loop_kernel(db, policy.RSD22_C4_table(), model, k, UV_cutoff, IR_cutoff, rsd22.get_C4());
        find_impl::read_loop_kernel(db, policy.RSD22_D1_table(), model, k, UV_cutoff, IR_cutoff, rsd22.get_D1());
    
        loop_integral payload(k, UV_cutoff, IR_cutoff, delta22, delta13, rsd22, rsd13);
        
        return std::move(payload);
      }
    
  }   // namespace sqlite3_operations
