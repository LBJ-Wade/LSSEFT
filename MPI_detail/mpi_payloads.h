//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_MPI_PAYLOADS_H
#define LSSEFT_MPI_PAYLOADS_H


#include <list>

#include "cosmology/types.h"

#include "mpi_operations.h"


namespace MPI_detail
  {

    //! build payload for transfer-function integration
    new_transfer_integration build_payload(const FRW_model& model, transfer_work_list::const_iterator& t);
    
    //! build payload for linear Pk filter calculation
    new_filter_Pk build_payload(const FRW_model& model, filter_Pk_work_list::const_iterator& t);

    //! build payload for loop integration
    new_loop_momentum_integration build_payload(const FRW_model& model, loop_momentum_work_list::const_iterator& t);
    
    //! build payload for Matsubara X & Y coefficient calculation
    new_Matsubara_XY build_payload(const FRW_model&, Matsubara_XY_work_list::const_iterator& t);
    
    //! build payload for one-loop P(k) calculation
    new_one_loop_Pk build_payload(const FRW_model&, one_loop_Pk_work_list::const_iterator& t);
    
    //! build payload for resummed one-loop P(k) calculation
    new_one_loop_resum_Pk build_payload(const FRW_model&, one_loop_resum_Pk_work_list::const_iterator& t);
    
    //! build payload for multipole P(k) calculation
    new_multipole_Pk build_payload(const FRW_model&, multipole_Pk_work_list::const_iterator& t);


  }   // namespace MPI_detail


#endif //LSSEFT_MPI_PAYLOADS_H
