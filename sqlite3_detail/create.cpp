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
              << "Phi DOUBLE"
#ifdef LSSEFT_STRICT_DATABASE_CONSISTENCY
              << ", "
              << "PRIMARY KEY (zid, kid), "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id), "
              << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id));";
#else
              << ");";
#endif
    
            exec(db, stmt.str());
          }
        
        
        void oneloop_g_table(sqlite3* db, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << policy.g_factor_table() << "("
              << "mid INTEGER, "
              << "zid INTEGER, "
              << "g_linear DOUBLE, "
              << "A DOUBLE, "
              << "B DOUBLE, "
              << "D DOUBLE, "
              << "E DOUBLE, "
              << "F DOUBLE, "
              << "G DOUBLE, "
              << "J DOUBLE"
#ifdef LSSEFT_STRICT_DATABASE_CONSISTENCY
              << ", "
              << "PRIMARY KEY (mid, zid), "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id));";
#else
              << ");";
#endif
    
            exec(db, stmt.str());
          }
    
    
        void oneloop_f_table(sqlite3* db, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << policy.f_factor_table() << "("
              << "mid INTEGER, "
              << "zid INTEGER, "
              << "f_linear DOUBLE, "
              << "fA DOUBLE, "
              << "fB DOUBLE, "
              << "fD DOUBLE, "
              << "fE DOUBLE, "
              << "fF DOUBLE, "
              << "fG DOUBLE, "
              << "fJ DOUBLE"
#ifdef LSSEFT_STRICT_DATABASE_CONSISTENCY
              << ", "
              << "PRIMARY KEY (mid, zid), "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id));";
#else
              << ");";
#endif
    
            exec(db, stmt.str());
          }
        
        
        void oneloop_momentum_integral_table(sqlite3* db, const std::string& table_name, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << table_name << "("
              << "mid INTEGER, "
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
              << "nw_time DOUBLE"
#ifdef LSSEFT_STRICT_DATABASE_CONSISTENCY
              << ", "
              << "PRIMARY KEY (mid, kid, Pk_id, IR_id, UV_id), "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id), "
              << "FOREIGN KEY (Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (IR_id) REFERENCES " << policy.IR_config_table() << "(id), "
              << "FOREIGN KEY (UV_id) REFERENCES " << policy.UV_config_table() << "(id));";
#else
            << ");";
#endif
            
            exec(db, stmt.str());
          }
        
        
        void oneloop_Pk_table(sqlite3* db, const std::string& table_name, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << table_name << "("
              << "mid INTEGER, "
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
              << "Z2_d_raw DOUBLE, "
              << "Ptree_nw DOUBLE, "
              << "err_tree_nw DOUBLE, "
              << "P13_nw DOUBLE, "
              << "err_13_nw DOUBLE, "
              << "P22_nw DOUBLE, "
              << "err_22_nw DOUBLE, "
              << "P1loopSPT_nw DOUBLE, "
              << "err_1loopSPT_nw DOUBLE, "
              << "Z2_d_nw DOUBLE"
#ifdef LSSEFT_STRICT_DATABASE_CONSISTENCY
              << ", "
              << "PRIMARY KEY (mid, zid, kid, init_Pk_id, final_Pk_id, IR_id, UV_id), "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id), "
              << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id), "
              << "FOREIGN KEY (init_Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (final_Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (IR_id) REFERENCES " << policy.IR_config_table() << "(id), "
              << "FOREIGN KEY (UV_id) REFERENCES " << policy.UV_config_table() << "(id));";
#else
            << ");";
#endif
            
            exec(db, stmt.str());
          }
        
        
        void oneloop_Pk_resum_table(sqlite3* db, const std::string& table_name, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << table_name << "("
              << "mid INTEGER, "
              << "zid INTEGER, "
              << "kid INTEGER, "
              << "init_Pk_id INTEGER, "
              << "final_Pk_id INTEGER, "
              << "IR_cutoff_id INTEGER, "
              << "UV_cutoff_id INTEGER, "
              << "IR_resum_id INTEGER, "
              << "Ptree DOUBLE, "
              << "err_tree DOUBLE, "
              << "P13 DOUBLE, "
              << "err_13 DOUBLE, "
              << "P22 DOUBLE, "
              << "err_22 DOUBLE, "
              << "P1loop_SPT DOUBLE, "
              << "err_1loop_SPT DOUBLE, "
              << "Z2_d DOUBLE"
#ifdef LSSEFT_STRICT_DATABASE_CONSISTENCY
              << ", "
              << "PRIMARY KEY (mid, zid, kid, init_Pk_id, final_Pk_id, IR_cutoff_id, UV_cutoff_id, IR_resum_id), "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id), "
              << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id), "
              << "FOREIGN KEY (init_Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (final_Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (IR_cutoff_id) REFERENCES " << policy.IR_config_table() << "(id), "
              << "FOREIGN KEY (UV_cutoff_id) REFERENCES " << policy.UV_config_table() << "(id), "
              << "FOREIGN KEY (IR_resum_id) REFERENCES " << policy.IR_resum_config_table() << "(id));";
#else
            << ");";
#endif
    
            exec(db, stmt.str());
          }
        
        
        void oneloop_rsd_Pk_table(sqlite3* db, const std::string& table_name, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << table_name << "("
              << "mid INTEGER, "
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
              << "Z2_d_raw DOUBLE, "
              << "Z0_v_raw DOUBLE, "
              << "Z2_v_raw DOUBLE, "
              << "Z0_vd_raw DOUBLE, "
              << "Z2_vd_raw DOUBLE, "
              << "Z2_vv_raw DOUBLE, "
              << "Z2_vvd_raw DOUBLE, "
              << "Z2_vvv_raw DOUBLE, "
              << "Z2_total_raw DOUBLE, "
              << "Ptree_nw DOUBLE, "
              << "err_tree_nw DOUBLE, "
              << "P13_nw DOUBLE, "
              << "err_13_nw DOUBLE, "
              << "P22_nw DOUBLE, "
              << "err_22_nw DOUBLE, "
              << "P1loopSPT_nw DOUBLE, "
              << "err_1loopSPT_nw DOUBLE, "
              << "Z2_d_nw DOUBLE, "
              << "Z0_v_nw DOUBLE, "
              << "Z2_v_nw DOUBLE, "
              << "Z0_vd_nw DOUBLE, "
              << "Z2_vd_nw DOUBLE, "
              << "Z2_vv_nw DOUBLE, "
              << "Z2_vvd_nw DOUBLE, "
              << "Z2_vvv_nw DOUBLE, "
              << "Z2_total_nw DOUBLE"
#ifdef LSSEFT_STRICT_DATABASE_CONSISTENCY
              << ", "
              << "PRIMARY KEY (mid, zid, kid, init_Pk_id, final_Pk_id, IR_id, UV_id), "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id), "
              << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id), "
              << "FOREIGN KEY (init_Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (final_Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (IR_id) REFERENCES " << policy.IR_config_table() << "(id), "
              << "FOREIGN KEY (UV_id) REFERENCES " << policy.UV_config_table() << "(id));";
#else
            << ");";
#endif
            
            exec(db, stmt.str());
          }
        
        
        void multipole_Pk_table(sqlite3* db, const std::string& table_name, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << table_name << "("
              << "mid INTEGER, "
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
              << "Z2_d DOUBLE, "
              << "Z2_d_resum DOUBLE, "
              << "Z0_v DOUBLE, "
              << "Z0_v_resum DOUBLE, "
              << "Z2_v DOUBLE, "
              << "Z2_v_resum DOUBLE, "
              << "Z0_vd DOUBLE, "
              << "Z0_vd_resum DOUBLE, "
              << "Z2_vd DOUBLE, "
              << "Z2_vd_resum DOUBLE, "
              << "Z2_vv DOUBLE, "
              << "Z2_vv_resum DOUBLE, "
              << "Z2_vvd DOUBLE, "
              << "Z2_vvd_resum DOUBLE, "
              << "Z2_vvv DOUBLE, "
              << "Z2_vvv_resum DOUBLE, "
              << "Z2_mu0 DOUBLE, "
              << "Z2_mu0_resum DOUBLE, "
              << "Z2_mu2 DOUBLE, "
              << "Z2_mu2_resum DOUBLE, "
              << "Z2_mu4 DOUBLE, "
              << "Z2_mu4_resum DOUBLE, "
              << "Z2_mu6 DOUBLE, "
              << "Z2_mu6_resum DOUBLE, "
              << "Z2_mu8 DOUBLE, "
              << "Z2_mu8_resum DOUBLE "
#ifdef LSSEFT_STRICT_DATABASE_CONSISTENCY
              << ", "
              << "PRIMARY KEY (mid, zid, kid, init_Pk_id, final_Pk_id, IR_cutoff_id, UV_cutoff_id, IR_resum_id), "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id), "
              << "FOREIGN KEY (zid) REFERENCES " << policy.redshift_config_table() << "(id), "
              << "FOREIGN KEY (init_Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (final_Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (IR_cutoff_id) REFERENCES " << policy.IR_config_table() << "(id), "
              << "FOREIGN KEY (UV_cutoff_id) REFERENCES " << policy.UV_config_table() << "(id) "
              << "FOREIGN KEY (IR_resum_id) REFERENCES " << policy.IR_resum_config_table() << "(id));";
#else
              << ");";
#endif
            
            exec(db, stmt.str());
          }
        
        
        void Matsubara_table(sqlite3* db, const sqlite3_policy& policy)
          {
            std::ostringstream stmt;
            stmt
              << "CREATE TABLE " << policy.Matsubara_XY_table() << "("
              << "mid INTEGER, "
              << "Pk_id INTEGER, "
              << "IR_resum_id INTEGER, "
              << "X DOUBLE, "
              << "Y DOUBLE"
#ifdef LSSEFT_STRICT_DATABASE_CONSISTENCY
              << ", "
              << "PRIMARY KEY (mid, IR_resum_id), "
              << "FOREIGN KEY (mid) REFERENCES " << policy.FRW_model_table() << "(id), "
              << "FOREIGN KEY (Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (IR_resum_id) REFERENCES " << policy.IR_resum_config_table() << "(id));";
#else
              << ");";
#endif
    
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
              << "Pk_ref DOUBLE"
#ifdef LSSEFT_STRICT_DATABASE_CONSISTENCY
              << ", "
              << "PRIMARY KEY (Pk_id, kid), "
              << "FOREIGN KEY (Pk_id) REFERENCES " << policy.Pk_linear_config_table() << "(id), "
              << "FOREIGN KEY (params_id) REFERENCES " << policy.filter_config_table() << "(id), "
              << "FOREIGN KEY (kid) REFERENCES " << policy.wavenumber_config_table() << "(id));";
#else
              << ");";
#endif
            
            exec(db, stmt.str());
          }
        
        
        void filter_data_config_table(sqlite3* db, const sqlite3_policy& policy)
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
        
      }
    
    
    void create_tables(sqlite3* db, const sqlite3_policy& policy)
      {
        create_impl::model_table(db, policy);
        create_impl::z_config_table(db, policy);
    
        create_impl::wavenumber_config_table(db, policy.wavenumber_config_table(), policy);
        create_impl::wavenumber_config_table(db, policy.IR_config_table(), policy);
        create_impl::wavenumber_config_table(db, policy.UV_config_table(), policy);
        create_impl::wavenumber_config_table(db, policy.IR_resum_config_table(), policy);
        
        create_impl::transfer_function_table(db, policy);
        
        create_impl::oneloop_g_table(db, policy);
        create_impl::oneloop_f_table(db, policy);
        
        create_impl::Pk_linear_config_table(db, policy);
        create_impl::Pk_linear_data_table(db, policy);
        
        create_impl::filter_data_config_table(db, policy);
        
        create_impl::oneloop_momentum_integral_table(db, policy.AA_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.AB_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.BB_table(), policy);
        
        create_impl::oneloop_momentum_integral_table(db, policy.D_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.E_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.F_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.G_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.J1_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.J2_table(), policy);
        
        create_impl::oneloop_momentum_integral_table(db, policy.RSD13_a_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD13_b_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD13_c_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD13_d_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD13_e_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD13_f_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD13_g_table(), policy);
        
        create_impl::oneloop_momentum_integral_table(db, policy.RSD22_A1_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD22_A2_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD22_A3_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD22_A4_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD22_A5_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD22_B2_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD22_B3_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD22_B6_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD22_B8_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD22_B9_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD22_C1_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD22_C2_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD22_C4_table(), policy);
        create_impl::oneloop_momentum_integral_table(db, policy.RSD22_D1_table(), policy);
        
        create_impl::oneloop_Pk_table(db, policy.dd_Pk_table(), policy);
        create_impl::oneloop_Pk_resum_table(db, policy.dd_Pk_resum_table(), policy);
        
        create_impl::oneloop_rsd_Pk_table(db, policy.dd_rsd_mu0_Pk_table(), policy);
        create_impl::oneloop_rsd_Pk_table(db, policy.dd_rsd_mu2_Pk_table(), policy);
        create_impl::oneloop_rsd_Pk_table(db, policy.dd_rsd_mu4_Pk_table(), policy);
        create_impl::oneloop_rsd_Pk_table(db, policy.dd_rsd_mu6_Pk_table(), policy);
        create_impl::oneloop_rsd_Pk_table(db, policy.dd_rsd_mu8_Pk_table(), policy);
        
        create_impl::multipole_Pk_table(db, policy.P0_table(), policy);
        create_impl::multipole_Pk_table(db, policy.P2_table(), policy);
        create_impl::multipole_Pk_table(db, policy.P4_table(), policy);
    
        create_impl::Matsubara_table(db, policy);
      }
    
  }   // namespace sqlite3_operations
