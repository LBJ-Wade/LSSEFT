//
// Created by David Seery on 10/11/2016.
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
                                      const loop_integral_params_token& params, const k_token& k, const linear_Pk_token& Pk,
                                      const UV_cutoff_token& UV_cutoff, KernelType& kernel, const IR_cutoff_token& IR_cutoff)
          {
            std::ostringstream read_stmt;
            read_stmt
              << "SELECT raw_value, raw_regions, raw_evals, raw_err, raw_time, nw_value, nw_regions, nw_evals, nw_err, nw_time FROM "
              << table << " WHERE mid=@mid AND params_id=@params_id AND kid=@kid AND Pk_id=@Pk_id AND UV_id=@UV_id AND IR_id=@IR_id;";
            
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, read_stmt.str().c_str(), read_stmt.str().length()+1, &stmt, nullptr));
            
            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@params_id"), params.get_id()));
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
    
    
        void read_dd_Pk(sqlite3* db, const std::string& table, const FRW_model_token& model,
                        const growth_params_token& growth_params, const loop_integral_params_token& loop_params,
                        const k_token& k, const z_token& z, const linear_Pk_token& init_Pk,
                        const boost::optional<linear_Pk_token>& final_Pk, const IR_cutoff_token& IR_cutoff,
                        const UV_cutoff_token& UV_cutoff, dd_Pk& Pk)
          {
            std::ostringstream read_stmt;
            read_stmt
              << "SELECT "
              << "Ptree_raw, err_tree_raw, P13_raw, err_13_raw, P22_raw, err_22_raw, P1loopSPT_raw, err_1loopSPT_raw, Z2_d_raw, "
              << "Ptree_nw, err_tree_nw, P13_nw, err_13_nw, P22_nw, err_22_nw, P1loopSPT_nw, err_1loopSPT_nw, Z2_d_nw "
              << "FROM " << table << " "
              << "WHERE mid=@mid AND growth_params=@growth_params AND loop_params=@loop_params AND "
              << "zid=@zid AND kid=@kid AND init_Pk_id=@init_Pk_id "
              << "AND ((@final_Pk_id IS NULL AND final_Pk_id IS NULL) OR final_Pk_id=@final_Pk_id) AND IR_id=@IR_id AND UV_id=@UV_id;";
            
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, read_stmt.str().c_str(), read_stmt.str().length()+1, &stmt, nullptr));
            
            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@growth_params"), growth_params.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@loop_params"), loop_params.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), z.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), k.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@init_Pk_id"), init_Pk.get_id()));
            if(final_Pk)
              {
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@final_Pk_id"), final_Pk->get_id()));
              }
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
            
            if(count != 1)
              {
                std::ostringstream msg;
                msg << ERROR_SQLITE3_READ_PK_MISREAD << " (count=" << count << ")";
                throw runtime_exception(exception_type::database_error, msg.str());
              }
          }
    
    
        void read_dd_rsd_Pk(sqlite3* db, const std::string& table, const FRW_model_token& model,
                            const growth_params_token& growth_params, const loop_integral_params_token& loop_params,
                            const k_token& k, const z_token& z, const linear_Pk_token& init_Pk,
                            const boost::optional<linear_Pk_token>& final_Pk, const IR_cutoff_token& IR_cutoff,
                            const UV_cutoff_token& UV_cutoff, rsd_dd_Pk& Pk)
          {
            std::ostringstream read_stmt;
            read_stmt
              << "SELECT "
              << "Ptree_raw, err_tree_raw, P13_raw, err_13_raw, P22_raw, err_22_raw, P1loopSPT_raw, err_1loopSPT_raw, "
              << "Z2_d_raw, Z0_v_raw, Z2_v_raw, Z0_vd_raw, Z2_vd_raw, Z2_vv_A_raw, Z2_vv_B_raw, Z2_vvd_raw, Z2_vvv_raw, Z2_total_raw, "
              << "Ptree_nw, err_tree_nw, P13_nw, err_13_nw, P22_nw, err_22_nw, P1loopSPT_nw, err_1loopSPT_nw, "
              << "Z2_d_nw, Z0_v_nw, Z2_v_nw, Z0_vd_nw, Z2_vd_nw, Z2_vv_A_nw, Z2_vv_B_nw, Z2_vvd_nw, Z2_vvv_nw, Z2_total_nw "
              << "FROM " << table << " "
              << "WHERE mid=@mid AND growth_params=@growth_params AND loop_params=@loop_params "
              << "AND zid=@zid AND kid=@kid AND init_Pk_id=@init_Pk_id "
              << "AND ((@final_Pk_id IS NULL AND final_Pk_id IS NULL) OR final_Pk_id=@final_Pk_id) AND IR_id=@IR_id AND UV_id=@UV_id;";
    
            constexpr unsigned int Ptree_raw = 0;
            constexpr unsigned int err_tree_raw = 1;
            constexpr unsigned int P13_raw = 2;
            constexpr unsigned int err_13_raw = 3;
            constexpr unsigned int P22_raw = 4;
            constexpr unsigned int err_22_raw = 5;
            constexpr unsigned int P1loopSPT_raw = 6;
            constexpr unsigned int err_1loopSPT_raw = 7;
            constexpr unsigned int Z2_d_raw = 8;
            constexpr unsigned int Z0_v_raw = 9;
            constexpr unsigned int Z2_v_raw = 10;
            constexpr unsigned int Z0_vd_raw = 11;
            constexpr unsigned int Z2_vd_raw = 12;
            constexpr unsigned int Z2_vv_A_raw = 13;
            constexpr unsigned int Z2_vv_B_raw = 14;
            constexpr unsigned int Z2_vvd_raw = 15;
            constexpr unsigned int Z2_vvv_raw = 16;
            constexpr unsigned int Z2_total_raw = 17;
            constexpr unsigned int Ptree_nw = 18;
            constexpr unsigned int err_tree_nw = 19;
            constexpr unsigned int P13_nw = 20;
            constexpr unsigned int err_13_nw = 21;
            constexpr unsigned int P22_nw = 22;
            constexpr unsigned int err_22_nw = 23;
            constexpr unsigned int P1loopSPT_nw = 24;
            constexpr unsigned int err_1loopSPT_nw = 25;
            constexpr unsigned int Z2_d_nw = 26;
            constexpr unsigned int Z0_v_nw = 27;
            constexpr unsigned int Z2_v_nw = 28;
            constexpr unsigned int Z0_vd_nw = 29;
            constexpr unsigned int Z2_vd_nw = 30;
            constexpr unsigned int Z2_vv_A_nw = 31;
            constexpr unsigned int Z2_vv_B_nw = 32;
            constexpr unsigned int Z2_vvd_nw = 33;
            constexpr unsigned int Z2_vvv_nw = 34;
            constexpr unsigned int Z2_total_nw = 35;
            
            // prepare statement
            sqlite3_stmt* stmt;
            check_stmt(db, sqlite3_prepare_v2(db, read_stmt.str().c_str(), read_stmt.str().length()+1, &stmt, nullptr));
    
            // bind parameter values
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@growth_params"), growth_params.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@loop_params"), loop_params.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@zid"), z.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@kid"), k.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@init_Pk_id"), init_Pk.get_id()));
            if(final_Pk)
              {
                check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@final_Pk_id"), final_Pk->get_id()));
              }
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@IR_id"), IR_cutoff.get_id()));
            check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@UV_id"), UV_cutoff.get_id()));
    
            // perform read
            int result = 0;
            unsigned int count = 0;
            while((result = sqlite3_step(stmt)) != SQLITE_DONE)
              {
                if(result == SQLITE_ROW)
                  {
                    read_Pk_value(stmt, Ptree_raw, err_tree_raw, Ptree_nw, err_tree_nw, Pk.get_tree());
                    read_Pk_value(stmt, P13_raw, err_13_raw, P13_nw, err_13_nw, Pk.get_13());
                    read_Pk_value(stmt, P22_raw, err_22_raw, P22_nw, err_22_nw, Pk.get_22());
                    read_Pk_value(stmt, P1loopSPT_raw, err_1loopSPT_raw, P1loopSPT_nw, err_1loopSPT_nw, Pk.get_1loop_SPT());
                    read_Pk_value(stmt, Z2_d_raw, Z2_d_nw, Pk.get_Z2_delta());
                    read_Pk_value(stmt, Z0_v_raw, Z0_v_nw, Pk.get_Z0_v());
                    read_Pk_value(stmt, Z2_v_raw, Z2_v_nw, Pk.get_Z2_v());
                    read_Pk_value(stmt, Z0_vd_raw, Z0_vd_nw, Pk.get_Z0_vdelta());
                    read_Pk_value(stmt, Z2_vd_raw, Z2_vd_nw, Pk.get_Z2_vdelta());
                    read_Pk_value(stmt, Z2_vv_A_raw, Z2_vv_A_nw, Pk.get_Z2_vv_A());
                    read_Pk_value(stmt, Z2_vv_B_raw, Z2_vv_B_nw, Pk.get_Z2_vv_B());
                    read_Pk_value(stmt, Z2_vvd_raw, Z2_vvd_nw, Pk.get_Z2_vvdelta());
                    read_Pk_value(stmt, Z2_vvv_raw, Z2_vvv_nw, Pk.get_Z2_vvv());
                    read_Pk_value(stmt, Z2_total_raw, Z2_total_nw, Pk.get_Z2_total());
            
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
    
    
            if(count != 1)
              {
                std::ostringstream msg;
                msg << ERROR_SQLITE3_READ_RSD_PK_MISREAD << " (count=" << count << ")";
                throw runtime_exception(exception_type::database_error, msg.str());
              }
          }
          
      }   // namespace find_impl
    
    
    std::unique_ptr<oneloop_growth>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
             const growth_params_token& params, const z_database& z_db)
      {
        // set up temporary table of desired z identifiers
        std::string ztab = z_table(db, mgr, policy, z_db);
        
        std::unique_ptr<oneloop_growth> payload = std::make_unique<oneloop_growth>(params, z_db);

        // note use of ORDER BY ... DESC which is needed to get the D, f values
        // in the correct order
        // z_table() will produce a table of tokens in ascending redshift order,
        // ie. ordered into the past
        
        std::ostringstream D_select_stmt;
        D_select_stmt
          << "SELECT "
          << "D_sample.D_linear, "
          << "D_sample.A, "
          << "D_sample.B, "
          << "D_sample.D, "
          << "D_sample.E, "
          << "D_sample.F, "
          << "D_sample.G, "
          << "D_sample.J "
          << "FROM " << ztab << " "
          << "INNER JOIN (SELECT * FROM " << policy.D_factor_table() << " WHERE mid=@mid AND params_id=@params_id) D_sample "
          << "ON " << ztab << ".id = D_sample.zid "
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
          << "INNER JOIN (SELECT * FROM " << policy.f_factor_table() << " WHERE mid=@mid AND params_id=@params_id) f_sample "
          << "ON " << ztab << ".id = f_sample.zid "
          << "ORDER BY " << ztab << ".ROWID DESC;";
        
        constexpr unsigned int D_linear = 0;
        constexpr unsigned int D_A = 1;
        constexpr unsigned int D_B = 2;
        constexpr unsigned int D_D = 3;
        constexpr unsigned int D_E = 4;
        constexpr unsigned int D_F = 5;
        constexpr unsigned int D_G = 6;
        constexpr unsigned int D_J = 7;

        constexpr unsigned int f_linear = 0;
        constexpr unsigned int f_A = 1;
        constexpr unsigned int f_B = 2;
        constexpr unsigned int f_D = 3;
        constexpr unsigned int f_E = 4;
        constexpr unsigned int f_F = 5;
        constexpr unsigned int f_G = 6;
        constexpr unsigned int f_J = 7;
        
        // prepare statements
        sqlite3_stmt* D_stmt;
        check_stmt(db, sqlite3_prepare_v2(db, D_select_stmt.str().c_str(), D_select_stmt.str().length()+1, &D_stmt, nullptr));
    
        sqlite3_stmt* f_stmt;
        check_stmt(db, sqlite3_prepare_v2(db, f_select_stmt.str().c_str(), f_select_stmt.str().length()+1, &f_stmt, nullptr));
    
        check_stmt(db, sqlite3_bind_int(D_stmt, sqlite3_bind_parameter_index(D_stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(f_stmt, sqlite3_bind_parameter_index(f_stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(D_stmt, sqlite3_bind_parameter_index(D_stmt, "@params_id"), params.get_id()));
        check_stmt(db, sqlite3_bind_int(f_stmt, sqlite3_bind_parameter_index(f_stmt, "@params_id"), params.get_id()));
        
        int D_status = 0;
        int f_status = 0;
        
        unsigned int read_count = 0;
        
        while((D_status = sqlite3_step(D_stmt)) != SQLITE_DONE && (f_status = sqlite3_step(f_stmt)) != SQLITE_DONE)
          {
            if(D_status == SQLITE_ROW && f_status == SQLITE_ROW)
              {
                payload->push_back(sqlite3_column_double(D_stmt, D_linear),
                                   sqlite3_column_double(D_stmt, D_A),
                                   sqlite3_column_double(D_stmt, D_B),
                                   sqlite3_column_double(D_stmt, D_D),
                                   sqlite3_column_double(D_stmt, D_E),
                                   sqlite3_column_double(D_stmt, D_F),
                                   sqlite3_column_double(D_stmt, D_G),
                                   sqlite3_column_double(D_stmt, D_J),
                                   sqlite3_column_double(f_stmt, f_linear),
                                   sqlite3_column_double(f_stmt, f_A),
                                   sqlite3_column_double(f_stmt, f_B),
                                   sqlite3_column_double(f_stmt, f_D),
                                   sqlite3_column_double(f_stmt, f_E),
                                   sqlite3_column_double(f_stmt, f_F),
                                   sqlite3_column_double(f_stmt, f_G),
                                   sqlite3_column_double(f_stmt, f_J));
                
                ++read_count;
              }
            else
              {
                std::ostringstream msg;
                msg << ERROR_SQLITE3_DF_GROWTH_TABLE_READ_FAIL << "(" << D_status << "," << f_status << "): " << sqlite3_errmsg(db) << "]";

                check_stmt(db, sqlite3_finalize(D_stmt));
                check_stmt(db, sqlite3_finalize(f_stmt));

                throw runtime_exception(exception_type::database_error, msg.str());
              }
          }
        
        check_stmt(db, sqlite3_finalize(D_stmt));
        check_stmt(db, sqlite3_finalize(f_stmt));
    
        // drop unneeded temporary tables
        drop_temp(db, mgr, ztab);
    
        if(read_count != z_db.size()) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_DF_GROWTH_MISREAD);
        
        return std::move(payload);
      }
    
    
    std::unique_ptr<loop_integral>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
         const loop_integral_params_token& params, const k_token& k, const linear_Pk_token& Pk,
         const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff)
      {
        delta_22_integrals delta22;
        delta_13_integrals delta13;
        rsd_22_integrals rsd22;
        rsd_13_integrals rsd13;
    
        find_impl::read_loop_kernel(db, policy.AA_table(), model, params, k, Pk, UV_cutoff, delta22.get_AA(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.AB_table(), model, params, k, Pk, UV_cutoff, delta22.get_AB(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.BB_table(), model, params, k, Pk, UV_cutoff, delta22.get_BB(), IR_cutoff);
    
        find_impl::read_loop_kernel(db, policy.D_table(), model, params, k, Pk, UV_cutoff, delta13.get_D(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.E_table(), model, params, k, Pk, UV_cutoff, delta13.get_E(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.F_table(), model, params, k, Pk, UV_cutoff, delta13.get_F(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.G_table(), model, params, k, Pk, UV_cutoff, delta13.get_G(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.J1_table(), model, params, k, Pk, UV_cutoff, delta13.get_J1(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.J2_table(), model, params, k, Pk, UV_cutoff, delta13.get_J2(), IR_cutoff);
    
        find_impl::read_loop_kernel(db, policy.RSD13_a_table(), model, params, k, Pk, UV_cutoff, rsd13.get_a(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD13_b_table(), model, params, k, Pk, UV_cutoff, rsd13.get_b(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD13_c_table(), model, params, k, Pk, UV_cutoff, rsd13.get_c(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD13_d_table(), model, params, k, Pk, UV_cutoff, rsd13.get_d(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD13_e_table(), model, params, k, Pk, UV_cutoff, rsd13.get_e(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD13_f_table(), model, params, k, Pk, UV_cutoff, rsd13.get_f(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD13_g_table(), model, params, k, Pk, UV_cutoff, rsd13.get_g(), IR_cutoff);
    
        find_impl::read_loop_kernel(db, policy.RSD22_A1_table(), model, params, k, Pk, UV_cutoff, rsd22.get_A1(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_A2_table(), model, params, k, Pk, UV_cutoff, rsd22.get_A2(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_A3_table(), model, params, k, Pk, UV_cutoff, rsd22.get_A3(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_A4_table(), model, params, k, Pk, UV_cutoff, rsd22.get_A4(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_A5_table(), model, params, k, Pk, UV_cutoff, rsd22.get_A5(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_B2_table(), model, params, k, Pk, UV_cutoff, rsd22.get_B2(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_B3_table(), model, params, k, Pk, UV_cutoff, rsd22.get_B3(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_B6_table(), model, params, k, Pk, UV_cutoff, rsd22.get_B6(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_B8_table(), model, params, k, Pk, UV_cutoff, rsd22.get_B8(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_B9_table(), model, params, k, Pk, UV_cutoff, rsd22.get_B9(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_C1_table(), model, params, k, Pk, UV_cutoff, rsd22.get_C1(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_C2_table(), model, params, k, Pk, UV_cutoff, rsd22.get_C2(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_C4_table(), model, params, k, Pk, UV_cutoff, rsd22.get_C4(), IR_cutoff);
        find_impl::read_loop_kernel(db, policy.RSD22_D1_table(), model, params, k, Pk, UV_cutoff, rsd22.get_D1(), IR_cutoff);
    
        std::unique_ptr<loop_integral> payload = std::make_unique<loop_integral>(k, params, Pk, UV_cutoff, IR_cutoff, delta22, delta13, rsd22, rsd13);
        
        return std::move(payload);
      }
    
    
    std::unique_ptr<oneloop_Pk>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
             const growth_params_token& growth_params, const loop_integral_params_token& loop_params, const k_token& k,
             const z_token& z, const linear_Pk_token& init_Pk_lin, const boost::optional<linear_Pk_token>& final_Pk_lin,
             const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff)
      {
        dd_Pk dd;
        rsd_dd_Pk rsd_mu0;
        rsd_dd_Pk rsd_mu2;
        rsd_dd_Pk rsd_mu4;
        rsd_dd_Pk rsd_mu6;
        rsd_dd_Pk rsd_mu8;
    
        find_impl::read_dd_Pk(db, policy.dd_Pk_table(), model, growth_params, loop_params, k, z, init_Pk_lin, final_Pk_lin, IR_cutoff, UV_cutoff, dd);
        find_impl::read_dd_rsd_Pk(db, policy.dd_rsd_mu0_Pk_table(), model, growth_params, loop_params, k, z, init_Pk_lin, final_Pk_lin, IR_cutoff, UV_cutoff, rsd_mu0);
        find_impl::read_dd_rsd_Pk(db, policy.dd_rsd_mu2_Pk_table(), model, growth_params, loop_params, k, z, init_Pk_lin, final_Pk_lin, IR_cutoff, UV_cutoff, rsd_mu2);
        find_impl::read_dd_rsd_Pk(db, policy.dd_rsd_mu4_Pk_table(), model, growth_params, loop_params, k, z, init_Pk_lin, final_Pk_lin, IR_cutoff, UV_cutoff, rsd_mu4);
        find_impl::read_dd_rsd_Pk(db, policy.dd_rsd_mu6_Pk_table(), model, growth_params, loop_params, k, z, init_Pk_lin, final_Pk_lin, IR_cutoff, UV_cutoff, rsd_mu6);
        find_impl::read_dd_rsd_Pk(db, policy.dd_rsd_mu8_Pk_table(), model, growth_params, loop_params, k, z, init_Pk_lin, final_Pk_lin, IR_cutoff, UV_cutoff, rsd_mu8);
    
        std::unique_ptr<oneloop_Pk> payload = std::make_unique<oneloop_Pk>(k, growth_params, loop_params, init_Pk_lin, final_Pk_lin,
                                                                           IR_cutoff, UV_cutoff, z, dd, rsd_mu0, rsd_mu2, rsd_mu4, rsd_mu6, rsd_mu8);
        
        return std::move(payload);
      }
    
    
    std::unique_ptr<Matsubara_XY>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
         const MatsubaraXY_params_token& params, const linear_Pk_token& Pk, const IR_resum_token& IR_resum)
      {
        std::ostringstream read_stmt;
        read_stmt << "SELECT X, Y FROM " << policy.Matsubara_XY_table() << " WHERE mid=@mid AND params_id=@params_id AND Pk_id=@Pk_id AND IR_resum_id=@IR_resum_id;";
    
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, read_stmt.str().c_str(), read_stmt.str().length()+1, &stmt, nullptr));
    
        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@params_id"), params.get_id()));
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
        
        if(count != 1)
          {
            std::ostringstream msg;
            msg << ERROR_SQLITE3_MATSUBARA_XY_MISREAD << " (count=" << count << ")";
            throw runtime_exception(exception_type::database_error, msg.str());
          }
        
        std::unique_ptr<Matsubara_XY> payload = std::make_unique<Matsubara_XY>(params, Pk, IR_resum, X, Y);
        
        return std::move(payload);
      }
    
  }   // namespace sqlite3_operations
