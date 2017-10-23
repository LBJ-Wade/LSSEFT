//
// Created by David Seery on 11/08/2015.
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


#include "utilities.h"

#include "create.h"

#include "defaults.h"

namespace sqlite3_operations
  {
    
    namespace create_impl
      {
        
        void model_table(sqlite3* db, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << policy.FRW_model_table() << "("
              << "id INTEGER PRIMARY KEY, "
              << "name TEXT, "
              << "omega_m DOUBLE, "
              << "omega_cc DOUBLE, "
              << "h DOUBLE, "
              << "T_CMB DOUBLE, "
              << "Neff DOUBLE, "
              << "f_baryon DOUBLE, "
              << "z_star DOUBLE, "
              << "z_drag DOUBLE, "
              << "z_eq DOUBLE, "
              << "A_curv DOUBLE, "
              << "ns DOUBLE, "
              << "k_piv DOUBLE"
              << ");";
    
            exec(db, stmt.str());
          }
        
        
        void z_config_table(sqlite3* db, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << policy.redshift_config_table() << "("
              << "id INTEGER PRIMARY KEY, "
              << "z DOUBLE);";
    
            exec(db, stmt.str());
          }

        
        void wavenumber_config_table(sqlite3* db, const std::string& table_name, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << table_name << "("
              << "id INTEGER PRIMARY KEY, "
              << "k DOUBLE);";
            exec(db, stmt.str());
          }
        
        
        void transfer_function_table(sqlite3* db, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << policy.transfer_table() << "("
              << "mid INTEGER, "
              << "kid INTEGER, "
              << "zid INTEGER, "
              << "delta_m DOUBLE, "
              << "delta_r DOUBLE, "
              << "theta_m DOUBLE, "
              << "theta_r DOUBLE, "
              << "Phi DOUBLE, "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id), "
              << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id));";
    
            exec(db, stmt.str());
          }
        
        
        void oneloop_g_table(sqlite3* db, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << policy.D_factor_table() << "("
              << "mid INTEGER, "
              << "params_id INTEGER, "
              << "zid INTEGER, "
              << "D_linear DOUBLE, "
              << "A DOUBLE, "
              << "B DOUBLE, "
              << "D DOUBLE, "
              << "E DOUBLE, "
              << "F DOUBLE, "
              << "G DOUBLE, "
              << "J DOUBLE, "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (params_id) REFERENCES " << policy.growth_config_table() << "(id), "
              << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id));";
    
            exec(db, stmt.str());
          }
    
    
        void oneloop_f_table(sqlite3* db, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << policy.f_factor_table() << "("
              << "mid INTEGER, "
              << "params_id INTEGER, "
              << "zid INTEGER, "
              << "f_linear DOUBLE, "
              << "fA DOUBLE, "
              << "fB DOUBLE, "
              << "fD DOUBLE, "
              << "fE DOUBLE, "
              << "fF DOUBLE, "
              << "fG DOUBLE, "
              << "fJ DOUBLE, "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (params_id) REFERENCES " << policy.growth_config_table() << "(id), "
              << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id));";
    
            exec(db, stmt.str());
          }
        
        
        void oneloop_momentum_integral_table(sqlite3* db, const std::string& table_name, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << table_name << "("
              << "mid INTEGER, "
              << "params_id INTEGER, "
              << "kid INTEGER, "
              << "Pk_id INTEGER, "
              << "IR_id INTEGER, "
              << "UV_id INTEGER, "
              << "raw_value DOUBLE, "
              << "raw_regions DOUBLE, "
              << "raw_evals DOUBLE, "
              << "raw_err DOUBLE, "
              << "raw_time DOUBLE, "
              << "nw_value DOUBLE, "
              << "nw_regions DOUBLE, "
              << "nw_evals DOUBLE, "
              << "nw_err DOUBLE, "
              << "nw_time DOUBLE, "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (params_id) REFERENCES " << policy.growth_config_table() << "(id), "
              << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id), "
              << "FOREIGN KEY (params_id) REFERENCES " << policy.loop_integral_config_table() << "(id), "
              << "FOREIGN KEY (Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (IR_id) REFERENCES " << policy.IR_config_table() << "(id), "
              << "FOREIGN KEY (UV_id) REFERENCES " << policy.UV_config_table() << "(id));";
            
            exec(db, stmt.str());
          }


        void oneloop_rsd_Pk_table(sqlite3* db, const std::string& table_name, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << table_name << "("
              << "mid INTEGER, "
              << "growth_params INTEGER, "
              << "loop_params INTEGER, "
              << "zid INTEGER, "
              << "kid INTEGER, "
              << "init_Pk_id INTEGER, "
              << "final_Pk_id INTEGER, "
              << "IR_id INTEGER, "
              << "UV_id INTEGER, "
              << "Ptree_raw DOUBLE, "
              << "err_tree_raw DOUBLE, "
              << "P13_raw DOUBLE, "
              << "err_13_raw DOUBLE, "
              << "P22_raw DOUBLE, "
              << "err_22_raw DOUBLE, "
              << "P1loopSPT_raw DOUBLE, "
              << "err_1loopSPT_raw DOUBLE, "
              << "Ptree_nw DOUBLE, "
              << "err_tree_nw DOUBLE, "
              << "P13_nw DOUBLE, "
              << "err_13_nw DOUBLE, "
              << "P22_nw DOUBLE, "
              << "err_22_nw DOUBLE, "
              << "P1loopSPT_nw DOUBLE, "
              << "err_1loopSPT_nw DOUBLE, "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (growth_params) REFERENCES " << policy.growth_config_table() << "(id), "
              << "FOREIGN KEY (loop_params) REFERENCES " << policy.loop_integral_config_table() << "(id), "
              << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id), "
              << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id), "
              << "FOREIGN KEY (init_Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (final_Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (IR_id) REFERENCES " << policy.IR_config_table() << "(id), "
              << "FOREIGN KEY (UV_id) REFERENCES " << policy.UV_config_table() << "(id));";
            
            exec(db, stmt.str());
          }
        
        
        void multipole_Pk_table(sqlite3* db, const std::string& table_name, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << table_name << "("
              << "mid INTEGER, "
              << "growth_params INTEGER, "
              << "loop_params INTEGER, "
              << "XY_params INTEGER, "
              << "zid INTEGER, "
              << "kid INTEGER, "
              << "init_Pk_id INTEGER, "
              << "final_Pk_id INTEGER, "
              << "IR_cutoff_id INTEGER, "
              << "UV_cutoff_id INTEGER, "
              << "IR_resum_id INTEGER, "
              << "Ptree DOUBLE, "
              << "Ptree_err DOUBLE, "
              << "Ptree_resum DOUBLE, "
              << "Ptree_resum_err DOUBLE, "
              << "P13 DOUBLE, "
              << "P13_err DOUBLE, "
              << "P13_resum DOUBLE, "
              << "P13_resum_err DOUBLE, "
              << "P22 DOUBLE, "
              << "P22_err DOUBLE, "
              << "P22_resum DOUBLE, "
              << "P22_resum_err DOUBLE, "
              << "P1loopSPT DOUBLE, "
              << "P1loopSPT_err DOUBLE, "
              << "P1loopSPT_resum DOUBLE, "
              << "P1loopSPT_resum_err DOUBLE, "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (growth_params) REFERENCES " << policy.growth_config_table() << "(id), "
              << "FOREIGN KEY (loop_params) REFERENCES " << policy.loop_integral_config_table() << "(id), "
              << "FOREIGN KEY (XY_params) REFERENCES " << policy.MatsubaraXY_config_table() << "(id), "
              << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id), "
              << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id), "
              << "FOREIGN KEY (init_Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (final_Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (IR_cutoff_id) REFERENCES " << policy.IR_config_table() << "(id), "
              << "FOREIGN KEY (UV_cutoff_id) REFERENCES " << policy.UV_config_table() << "(id) "
              << "FOREIGN KEY (IR_resum_id) REFERENCES " << policy.IR_resum_config_table() << "(id));";
            
            exec(db, stmt.str());
          }
        
        
        void Matsubara_table(sqlite3* db, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << policy.Matsubara_XY_table() << "("
              << "mid INTEGER, "
              << "params_id INTEGER, "
              << "Pk_id INTEGER, "
              << "IR_resum_id INTEGER, "
              << "X DOUBLE, "
              << "Y DOUBLE, "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (params_id) REFERENCES " << policy.MatsubaraXY_config_table() << "(id), "
              << "FOREIGN KEY (Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (IR_resum_id) REFERENCES " << policy.IR_resum_config_table() << "(id));";
    
            exec(db, stmt.str());
          }
        
        
        void Pk_linear_config_table(sqlite3* db, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << policy.Pk_linear_config_table() << "("
              << "id INTEGER PRIMARY KEY, "
              << "mid INTEGER, "
              << "path TEXT, "
              << "md5_hash TEXT, "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id));";
            
            exec(db, stmt.str());
          }
        
        
        void Pk_linear_data_table(sqlite3* db, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << policy.Pk_linear_table() << "("
              << "Pk_id INTEGER, "
              << "params_id INTEGER, "
              << "kid INTEGER, "
              << "Pk_raw DOUBLE, "
              << "Pk_nw DOUBLE, "
              << "Pk_ref DOUBLE, "
              << "Pk_nw_err DOUBLE, "
              << "regions DOUBLE, "
              << "evaluations DOUBLE, "
              << "time DOUBLE, "
              << "FOREIGN KEY (Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (params_id) REFERENCES " << policy.filter_config_table() << "(id), "
              << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id));";
            
            exec(db, stmt.str());
          }
        
        
        void filter_params_config_table(sqlite3* db, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << policy.filter_config_table() << "("
              << "id INTEGER PRIMARY KEY, "
              << "amplitude DOUBLE, "
              << "pivot DOUBLE, "
              << "idx DOUBLE, "
              << "abserr DOUBLE, "
              << "relerr DOUBLE"
              << ");";
              
            exec(db, stmt.str());
          }
    
    
        void oneloop_params_config_table(sqlite3* db, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << policy.loop_integral_config_table() << "("
              << "id INTEGER PRIMARY KEY, "
              << "abserr_13 DOUBLE, "
              << "relerr_13 DOUBLE, "
              << "abserr_22 DOUBLE, "
              << "relerr_22 DOUBLE"
              << ");";
        
            exec(db, stmt.str());
          }
    
    
        void MatsubaraXY_params_config_table(sqlite3* db, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << policy.MatsubaraXY_config_table() << "("
              << "id INTEGER PRIMARY KEY, "
              << "abserr DOUBLE, "
              << "relerr DOUBLE, "
              << "qmin DOUBLE, "
              << "qmax DOUBLE"
              << ");";
        
            exec(db, stmt.str());
          }
    
    
        void growth_params_config_table(sqlite3* db, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << policy.growth_config_table() << "("
              << "id INTEGER PRIMARY KEY, "
              << "abserr DOUBLE, "
              << "relerr DOUBLE, "
              << "use_EdS INTEGER, "
              << "use_EdS_ics INTEGER"
              << ");";
        
            exec(db, stmt.str());
          }


        void counterterms_table(sqlite3* db, const std::string& table_name, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << table_name << "("
              << "mid INTEGER, "
              << "growth_params INTEGER, "
              << "XY_params INTEGER, "
              << "zid INTEGER, "
              << "kid INTEGER, "
              << "init_Pk_id INTEGER, "
              << "final_Pk_id INTEGER, "
              << "IR_cutoff_id INTEGER, "
              << "UV_cutoff_id INTEGER, "
              << "IR_resum_id INTEGER, "
              << "P0_k0_raw DOUBLE, "
              << "P0_k0_raw_err DOUBLE, "
              << "P0_k0_resum DOUBLE, "
              << "P0_k0_resum_err DOUBLE, "
              << "P2_k0_raw DOUBLE, "
              << "P2_k0_raw_err DOUBLE, "
              << "P2_k0_resum DOUBLE, "
              << "P2_k0_resum_err DOUBLE, "
              << "P4_k0_raw DOUBLE, "
              << "P4_k0_raw_err DOUBLE, "
              << "P4_k0_resum DOUBLE, "
              << "P4_k0_resum_err DOUBLE, "
              << "P0_k2_raw DOUBLE, "
              << "P0_k2_raw_err DOUBLE, "
              << "P0_k2_resum DOUBLE, "
              << "P0_k2_resum_err DOUBLE, "
              << "P2_k2_raw DOUBLE, "
              << "P2_k2_raw_err DOUBLE, "
              << "P2_k2_resum DOUBLE, "
              << "P2_k2_resum_err DOUBLE, "
              << "P4_k2_raw DOUBLE, "
              << "P4_k2_raw_err DOUBLE, "
              << "P4_k2_resum DOUBLE, "
              << "P4_k2_resum_err DOUBLE, "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (growth_params) REFERENCES " << policy.growth_config_table() << "(id), "
              << "FOREIGN KEY (XY_params) REFERENCES " << policy.MatsubaraXY_config_table() << "(id), "
              << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id), "
              << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id), "
              << "FOREIGN KEY (init_Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (final_Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (IR_cutoff_id) REFERENCES " << policy.IR_config_table() << "(id), "
              << "FOREIGN KEY (UV_cutoff_id) REFERENCES " << policy.UV_config_table() << "(id) "
              << "FOREIGN KEY (IR_resum_id) REFERENCES " << policy.IR_resum_config_table() << "(id));";

            exec(db, stmt.str());
          }
        
      }
    
    
    void create_tables(sqlite3* db, const sqlite3_policy& policy)
      {
        create_impl::model_table(db, policy);
        create_impl::z_config_table(db, policy);
    
        create_impl::wavenumber_config_table(db, policy.wavenumber_config_table(), policy);
        create_impl::wavenumber_config_table(db, policy.IR_config_table(), policy);
        create_impl::wavenumber_config_table(db, policy.UV_config_table(), policy);
        create_impl::wavenumber_config_table(db, policy.IR_resum_config_table(), policy);
        
        create_impl::Pk_linear_config_table(db, policy);
        
        create_impl::oneloop_params_config_table(db, policy);
        create_impl::MatsubaraXY_params_config_table(db, policy);
        create_impl::growth_params_config_table(db, policy);
        create_impl::filter_params_config_table(db, policy);
    
        create_impl::transfer_function_table(db, policy);
    
        create_impl::oneloop_g_table(db, policy);
        create_impl::oneloop_f_table(db, policy);

        create_impl::Pk_linear_data_table(db, policy);

        create_impl::Matsubara_table(db, policy);

        create_impl::counterterms_table(db, policy.counterterms_c0_table(), policy);
        create_impl::counterterms_table(db, policy.counterterms_c2_table(), policy);
        create_impl::counterterms_table(db, policy.counterterms_c4_table(), policy);

#include "autogenerated/create_stmts.cpp"
      }
    
  }   // namespace sqlite3_operations
