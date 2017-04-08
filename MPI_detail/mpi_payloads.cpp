//
// Created by David Seery on 21/11/2015.
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


#include "mpi_payloads.h"


namespace MPI_detail
  {

    new_transfer_integration build_payload(const FRW_model& model, transfer_work_list::const_iterator& t)
      {
        return new_transfer_integration(model, *(*t), t->get_token(), t->get_z_db());
      }


    new_loop_momentum_integration build_payload(const FRW_model& model, loop_integral_work_list::const_iterator& t)
      {
        return new_loop_momentum_integration(model, *(*t), t->get_k_token(), t->get_UV_cutoff(), t->get_UV_token(),
                                             t->get_IR_cutoff(), t->get_IR_token(), t->get_tree_Pk_db(),
                                             t->get_params_token(), t->get_params());
      }


    new_filter_Pk build_payload(const FRW_model& model, filter_Pk_work_list::const_iterator& t)
      {
        return new_filter_Pk(model, *(*t), t->get_k_token(), t->get_Pk_token(), t->get_linear_Pk(), t->get_params_token(), t->get_params());
      }
    
    
    new_Matsubara_XY build_payload(const FRW_model&, Matsubara_XY_work_list::const_iterator& t)
      {
        return new_Matsubara_XY(t->get_IR_resum(), t->get_IR_resum_token(), t->get_linear_Pk(), t->get_params_token(), t->get_params());
      }
    
    
    new_one_loop_Pk build_payload(const FRW_model&, one_loop_Pk_work_list::const_iterator& t)
      {
        return new_one_loop_Pk(*(*t), t->get_gf_factors(), t->get_loop_data(), t->get_init_linear_Pk(), t->get_final_linear_Pk());
      }
    
    
    new_one_loop_resum_Pk build_payload(const FRW_model&, one_loop_resum_Pk_work_list::const_iterator& t)
      {
        return new_one_loop_resum_Pk(*(*t), t->get_Matsubara_XY(), t->get_Pk_data(), t->get_Df_data(),
                                     t->get_init_linear_Pk(), t->get_final_linear_Pk());
      }
    
    
    new_multipole_Pk build_payload(const FRW_model&, multipole_Pk_work_list::const_iterator& t)
      {
        return new_multipole_Pk(*(*t), t->get_Matsubara_XY(), t->get_Pk_data(), t->get_Df_data(),
                                t->get_init_linear_Pk(), t->get_final_linear_Pk());
      }
    
  }
