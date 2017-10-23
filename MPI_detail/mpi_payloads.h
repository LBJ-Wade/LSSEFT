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

#ifndef LSSEFT_MPI_PAYLOADS_H
#define LSSEFT_MPI_PAYLOADS_H


#include <list>

#include "cosmology/types.h"

#include "mpi_operations.h"


namespace MPI_detail
  {

    //! build payload for transfer-function integration
    new_transfer_integration build_payload(const FRW_model& model, transfer_work_list::const_iterator& t);

    //! build payload for loop integration
    new_loop_momentum_integration build_payload(const FRW_model& model, loop_integral_work_list::const_iterator& t);
    
    //! build payload for linear Pk filter calculation
    new_filter_Pk build_payload(const FRW_model& model, filter_Pk_work_list::const_iterator& t);
    
    //! build payload for Matsubara X & Y coefficient calculation
    new_Matsubara_XY build_payload(const FRW_model&, Matsubara_XY_work_list::const_iterator& t);
    
    //! build payload for one-loop P(k) calculation
    new_one_loop_Pk build_payload(const FRW_model&, one_loop_Pk_work_list::const_iterator& t);

    //! build payload for multipole P(k) calculation
    new_multipole_Pk build_payload(const FRW_model&, multipole_Pk_work_list::const_iterator& t);

    //! build payload for counterterm calculation
    new_counterterm build_payload(const FRW_model&, counterterm_work_list::const_iterator& t);
    
    
  }   // namespace MPI_detail


#endif //LSSEFT_MPI_PAYLOADS_H
