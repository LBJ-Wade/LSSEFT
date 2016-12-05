//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "mpi_payloads.h"


namespace MPI_detail
  {

    new_transfer_integration build_payload(const FRW_model& model, transfer_work_list::const_iterator& t)
      {
        return new_transfer_integration(model, *(*t), t->get_token(), t->get_z_db());
      }


    new_loop_momentum_integration build_payload(const FRW_model& model, loop_momentum_work_list::const_iterator& t)
      {
        return new_loop_momentum_integration(model, *(*t), t->get_k_token(),
                                             t->get_UV_cutoff(), t->get_UV_token(),
                                             t->get_IR_cutoff(), t->get_IR_token(),
                                             t->get_tree_Pk_db());
      }
    
    
    new_one_loop_Pk build_payload(const FRW_model&, one_loop_Pk_work_list::const_iterator& t)
      {
        return new_one_loop_Pk(*(*t), t->get_gf_factors(), t->get_loop_data(), t->get_tree_Pk_db());
      }
    
    
    new_multipole_Pk build_payload(const FRW_model&, multipole_Pk_work_list::const_iterator& t)
      {
        return new_multipole_Pk(*(*t), t->get_Matsubara_A(), t->get_Pk_data(), t->get_gf_data(), t->get_tree_Pk_db());
      }
    
    
    new_Matsubara_A build_payload(const FRW_model&, Matsubara_A_work_list::const_iterator& t)
      {
        return new_Matsubara_A(t->get_IR_resum(), t->get_IR_resum_token(), t->get_tree_Pk_db());
      }
    
    
    new_filter_Pk build_payload(const FRW_model& model, filter_Pk_work_list::const_iterator& t)
      {
        return new_filter_Pk(model, *(*t), t->get_k_token(), t->get_Pk_token(), t->get_linear_Pk());
      }
    
  }
