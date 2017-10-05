//
// Created by David Seery on 22/03/2017.
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

#include "database/data_manager.h"


void data_manager::setup_write(transfer_work_list& work)
  {
    sqlite3_operations::write_performance_pragmas(this->handle, this->arg_cache.is_network_mode());

    sqlite3_operations::drop_index(this->handle, this->policy.transfer_table(), { "mid", "kid", "zid" });
  }


void data_manager::setup_growth_write()
  {
    sqlite3_operations::write_performance_pragmas(this->handle, this->arg_cache.is_network_mode());
    
    sqlite3_operations::drop_index(this->handle, this->policy.D_factor_table(), { "mid", "params_id", "zid" });
    sqlite3_operations::drop_index(this->handle, this->policy.f_factor_table(), { "mid", "params_id", "zid" });
  }


void data_manager::setup_write(loop_integral_work_list& work)
  {
    sqlite3_operations::write_performance_pragmas(this->handle, this->arg_cache.is_network_mode());

    sqlite3_operations::drop_index(this->handle, this->policy.AA_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.AB_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.BB_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    
    sqlite3_operations::drop_index(this->handle, this->policy.D_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.E_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.F_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.G_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.J1_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.J2_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    
    sqlite3_operations::drop_index(this->handle, this->policy.RSD13_a_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD13_b_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD13_c_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD13_d_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD13_e_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD13_f_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD13_g_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    
    sqlite3_operations::drop_index(this->handle, this->policy.RSD22_A1_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD22_A2_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD22_A3_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD22_A4_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD22_A5_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD22_B2_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD22_B3_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD22_B6_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD22_B8_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD22_B9_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD22_C1_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD22_C2_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD22_C4_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.RSD22_D1_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
  }


void data_manager::setup_write(filter_Pk_work_list& work)
  {
    sqlite3_operations::write_performance_pragmas(this->handle, this->arg_cache.is_network_mode());

    sqlite3_operations::drop_index(this->handle, this->policy.Pk_linear_table(), "kid");
  }


void data_manager::setup_write(Matsubara_XY_work_list& work)
  {
    sqlite3_operations::write_performance_pragmas(this->handle, this->arg_cache.is_network_mode());
  }


void data_manager::setup_write(one_loop_Pk_work_list& work)
  {
    sqlite3_operations::write_performance_pragmas(this->handle, this->arg_cache.is_network_mode());
    
    sqlite3_operations::drop_index(this->handle, this->policy.dd_Pk_table(),
                                   { "mid", "growth_params", "loop_params", "kid", "zid", "init_Pk_id", "final_Pk_id",
                                     "IR_id", "UV_id", });
    
    sqlite3_operations::drop_index(this->handle, this->policy.dd_rsd_mu0_Pk_table(),
                                   { "mid", "growth_params", "loop_params", "kid", "zid", "init_Pk_id", "final_Pk_id",
                                     "IR_id", "UV_id", });
    sqlite3_operations::drop_index(this->handle, this->policy.dd_rsd_mu2_Pk_table(),
                                   { "mid", "growth_params", "loop_params", "kid", "zid", "init_Pk_id", "final_Pk_id",
                                     "IR_id", "UV_id", });
    sqlite3_operations::drop_index(this->handle, this->policy.dd_rsd_mu4_Pk_table(),
                                   { "mid", "growth_params", "loop_params", "kid", "zid", "init_Pk_id", "final_Pk_id",
                                     "IR_id", "UV_id", });
    sqlite3_operations::drop_index(this->handle, this->policy.dd_rsd_mu6_Pk_table(),
                                   { "mid", "growth_params", "loop_params", "kid", "zid", "init_Pk_id", "final_Pk_id",
                                     "IR_id", "UV_id", });
    sqlite3_operations::drop_index(this->handle, this->policy.dd_rsd_mu8_Pk_table(),
                                   { "mid", "growth_params", "loop_params", "kid", "zid", "init_Pk_id", "final_Pk_id",
                                     "IR_id", "UV_id", });
  }

void data_manager::setup_write(one_loop_resum_Pk_work_list& work)
  {
    sqlite3_operations::write_performance_pragmas(this->handle, this->arg_cache.is_network_mode());

    sqlite3_operations::drop_index(this->handle, this->policy.dd_Pk_resum_table(),
                                   { "mid", "growth_params", "loop_params", "XY_params", "kid", "zid", "init_Pk_id",
                                     "final_Pk_id", "IR_cutoff_id", "UV_cutoff_id", "IR_resum_id" });
  }


void data_manager::setup_write(multipole_Pk_work_list& work)
  {
    sqlite3_operations::write_performance_pragmas(this->handle, this->arg_cache.is_network_mode());
    
    sqlite3_operations::drop_index(this->handle, this->policy.P0_table(),
                                   { "mid", "growth_params", "loop_params", "XY_params", "kid", "zid", "init_Pk_id",
                                     "final_Pk_id", "IR_cutoff_id", "UV_cutoff_id", "IR_resum_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.P2_table(),
                                   { "mid", "growth_params", "loop_params", "XY_params", "kid", "zid", "init_Pk_id",
                                     "final_Pk_id", "IR_cutoff_id", "UV_cutoff_id", "IR_resum_id" });
    sqlite3_operations::drop_index(this->handle, this->policy.P4_table(),
                                   { "mid", "growth_params", "loop_params", "XY_params", "kid", "zid", "init_Pk_id",
                                     "final_Pk_id", "IR_cutoff_id", "UV_cutoff_id", "IR_resum_id" });
  }


void data_manager::finalize_write(transfer_work_list& work)
  {
    sqlite3_operations::default_pragmas(this->handle);
    
    sqlite3_operations::create_index(this->handle, this->policy.transfer_table(), { "mid", "kid", "zid" });
    
    sqlite3_operations::analyze(this->handle);
  }


void data_manager::finalize_growth_write()
  {
    sqlite3_operations::default_pragmas(this->handle);
    
    sqlite3_operations::create_index(this->handle, this->policy.D_factor_table(), { "mid", "params_id", "zid" });
    sqlite3_operations::create_index(this->handle, this->policy.f_factor_table(), { "mid", "params_id", "zid" });
    
    sqlite3_operations::analyze(this->handle);
  }


void data_manager::finalize_write(loop_integral_work_list& work)
  {
    sqlite3_operations::default_pragmas(this->handle);
    
    sqlite3_operations::create_index(this->handle, this->policy.AA_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.AB_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.BB_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    
    sqlite3_operations::create_index(this->handle, this->policy.D_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.E_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.F_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.G_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.J1_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.J2_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    
    sqlite3_operations::create_index(this->handle, this->policy.RSD13_a_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD13_b_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD13_c_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD13_d_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD13_e_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD13_f_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD13_g_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    
    sqlite3_operations::create_index(this->handle, this->policy.RSD22_A1_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD22_A2_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD22_A3_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD22_A4_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD22_A5_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD22_B2_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD22_B3_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD22_B6_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD22_B8_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD22_B9_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD22_C1_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD22_C2_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD22_C4_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    sqlite3_operations::create_index(this->handle, this->policy.RSD22_D1_table(), { "mid", "params_id", "kid", "Pk_id", "IR_id", "UV_id" });
    
    sqlite3_operations::analyze(this->handle);
  }


void data_manager::finalize_write(filter_Pk_work_list& work)
  {
    sqlite3_operations::default_pragmas(this->handle);

    sqlite3_operations::create_index(this->handle, this->policy.Pk_linear_table(), "kid");
  }


void data_manager::finalize_write(Matsubara_XY_work_list& work)
  {
    sqlite3_operations::default_pragmas(this->handle);
  }


void data_manager::finalize_write(one_loop_Pk_work_list& work)
  {
    sqlite3_operations::default_pragmas(this->handle);
    
    sqlite3_operations::create_index(this->handle, this->policy.dd_Pk_table(),
                                     { "mid", "growth_params", "loop_params", "kid", "zid", "init_Pk_id", "final_Pk_id",
                                       "IR_id", "UV_id", });
    
    sqlite3_operations::create_index(this->handle, this->policy.dd_rsd_mu0_Pk_table(),
                                     { "mid", "growth_params", "loop_params", "kid", "zid", "init_Pk_id", "final_Pk_id",
                                       "IR_id", "UV_id", });
    sqlite3_operations::create_index(this->handle, this->policy.dd_rsd_mu2_Pk_table(),
                                     { "mid", "growth_params", "loop_params", "kid", "zid", "init_Pk_id", "final_Pk_id",
                                       "IR_id", "UV_id", });
    sqlite3_operations::create_index(this->handle, this->policy.dd_rsd_mu4_Pk_table(),
                                     { "mid", "growth_params", "loop_params", "kid", "zid", "init_Pk_id", "final_Pk_id",
                                       "IR_id", "UV_id", });
    sqlite3_operations::create_index(this->handle, this->policy.dd_rsd_mu6_Pk_table(),
                                     { "mid", "growth_params", "loop_params", "kid", "zid", "init_Pk_id", "final_Pk_id",
                                       "IR_id", "UV_id", });
    sqlite3_operations::create_index(this->handle, this->policy.dd_rsd_mu8_Pk_table(),
                                     { "mid", "growth_params", "loop_params", "kid", "zid", "init_Pk_id", "final_Pk_id",
                                       "IR_id", "UV_id", });
    
    sqlite3_operations::analyze(this->handle);
  }


void data_manager::finalize_write(one_loop_resum_Pk_work_list& work)
  {
    sqlite3_operations::default_pragmas(this->handle);
    
    sqlite3_operations::create_index(this->handle, this->policy.dd_Pk_resum_table(),
                                     { "mid", "growth_params", "loop_params", "XY_params", "kid", "zid", "init_Pk_id",
                                       "final_Pk_id", "IR_cutoff_id", "UV_cutoff_id", "IR_resum_id" });
    
    sqlite3_operations::analyze(this->handle);
  }


void data_manager::finalize_write(multipole_Pk_work_list& work)
  {
    sqlite3_operations::default_pragmas(this->handle);
    
    sqlite3_operations::create_index(this->handle, this->policy.P0_table(),
                                     { "mid", "growth_params", "loop_params", "XY_params", "kid", "zid", "init_Pk_id",
                                       "final_Pk_id", "IR_cutoff_id", "UV_cutoff_id", "IR_resum_id" });
    sqlite3_operations::create_index(this->handle, this->policy.P2_table(),
                                     { "mid", "growth_params", "loop_params", "XY_params", "kid", "zid", "init_Pk_id",
                                       "final_Pk_id", "IR_cutoff_id", "UV_cutoff_id", "IR_resum_id" });
    sqlite3_operations::create_index(this->handle, this->policy.P4_table(),
                                     { "mid", "growth_params", "loop_params", "XY_params", "kid", "zid", "init_Pk_id",
                                       "final_Pk_id", "IR_cutoff_id", "UV_cutoff_id", "IR_resum_id" });
    
    sqlite3_operations::analyze(this->handle);
  }
