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

#ifndef LSSEFT_MPI_TRAITS_H
#define LSSEFT_MPI_TRAITS_H


#include "cosmology/types.h"

#include "mpi_operations.h"


namespace MPI_detail
  {

    template <typename WorkItem> struct work_item_traits;

    
    template <> struct work_item_traits<transfer_work_record>
      {
        typedef new_transfer_integration   outgoing_payload_type;
        typedef transfer_integration_ready incoming_payload_type;

        static constexpr unsigned int new_task_message() { return(MESSAGE_NEW_TRANSFER_TASK); }
        static constexpr unsigned int new_item_message() { return(MESSAGE_NEW_TRANSFER_INTEGRATION); }
      };
    
    
    template <> struct work_item_traits<filter_Pk_work_record>
      {
        typedef new_filter_Pk   outgoing_payload_type;
        typedef filter_Pk_ready incoming_payload_type;
    
        static constexpr unsigned int new_task_message() { return(MESSAGE_NEW_FILTER_PK_TASK); }
        static constexpr unsigned int new_item_message() { return(MESSAGE_NEW_FILTER_PK); }
      };

    
    template <> struct work_item_traits<loop_momentum_work_record>
      {
        typedef new_loop_momentum_integration   outgoing_payload_type;
        typedef loop_momentum_integration_ready incoming_payload_type;

        static constexpr unsigned int new_task_message() { return(MESSAGE_NEW_LOOP_INTEGRAL_TASK); }
        static constexpr unsigned int new_item_message() { return(MESSAGE_NEW_LOOP_INTEGRATION); }
      };
    
    
    template<> struct work_item_traits<Matsubara_XY_work_record>
      {
        typedef new_Matsubara_XY   outgoing_payload_type;
        typedef Matsubara_XY_ready incoming_payload_type;
        
        static constexpr unsigned int new_task_message() { return(MESSAGE_NEW_MATSUBARA_XY_TASK); }
        static constexpr unsigned int new_item_message() { return(MESSAGE_NEW_MATSUBARA_XY); }
      };
    
    
    template <> struct work_item_traits<one_loop_Pk_work_record>
      {
        typedef new_one_loop_Pk   outgoing_payload_type;
        typedef one_loop_Pk_ready incoming_payload_type;
        
        static constexpr unsigned int new_task_message() { return(MESSAGE_NEW_ONE_LOOP_PK_TASK); }
        static constexpr unsigned int new_item_message() { return(MESSAGE_NEW_ONE_LOOP_PK); }
      };
    
    
    template <> struct work_item_traits<one_loop_resum_Pk_work_record>
      {
        typedef new_one_loop_resum_Pk   outgoing_payload_type;
        typedef one_loop_resum_Pk_ready incoming_payload_type;
    
        static constexpr unsigned int new_task_message() { return(MESSAGE_NEW_ONE_LOOP_RESUM_PK_TASK); }
        static constexpr unsigned int new_item_message() { return(MESSAGE_NEW_ONE_LOOP_RESUM_PK); }
      };
    
    
    template <> struct work_item_traits<multipole_Pk_work_record>
      {
        typedef new_multipole_Pk   outgoing_payload_type;
        typedef multipole_Pk_ready incoming_payload_type;
        
        static constexpr unsigned int new_task_message() { return(MESSAGE_NEW_MULTIPOLE_PK_TASK); }
        static constexpr unsigned int new_item_message() { return(MESSAGE_NEW_MULTIPOLE_PK); }
      };

  }   // namespace MPI_detail


#endif //LSSEFT_MPI_TRAITS_H
