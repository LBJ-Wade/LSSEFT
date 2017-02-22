//
// Created by David Seery on 10/08/2015.
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


#include "MPI_detail/mpi_traits.h"
#include "MPI_detail/mpi_payloads.h"

#include "slave_controller.h"

#include "cosmology/transfer_integrator.h"
#include "cosmology/oneloop_momentum_integrator.h"
#include "cosmology/oneloop_Pk_calculator.h"
#include "cosmology/multipole_Pk_calculator.h"
#include "cosmology/Pk_filter.h"
#include "cosmology/Matsubara_XY_calculator.h"

#include "cosmology/types.h"


slave_controller::slave_controller(boost::mpi::environment& me, boost::mpi::communicator& mw, argument_cache& ac)
  : mpi_env(me),
    mpi_world(mw),
    arg_cache(ac),
    local_env(),
    err_handler(arg_cache, local_env)
  {
  }


void slave_controller::execute()
  {
    // wait until we receive an EXIT instruction from master
    bool finished = false;

    while(!finished)
      {
        // wait for a message from the master node
        boost::mpi::status stat = this->mpi_world.probe(MPI_detail::RANK_MASTER);

        switch(stat.tag())
          {
            case MPI_detail::MESSAGE_NEW_TRANSFER_TASK:
              {
                this->mpi_world.recv(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_NEW_TRANSFER_TASK);
                this->process_task<transfer_work_record>();
                break;
              }
            
            case MPI_detail::MESSAGE_NEW_FILTER_PK_TASK:
              {
                this->mpi_world.recv(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_NEW_FILTER_PK_TASK);
                this->process_task<filter_Pk_work_record>();
                break;
              }

            case MPI_detail::MESSAGE_NEW_LOOP_INTEGRAL_TASK:
              {
                this->mpi_world.recv(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_NEW_LOOP_INTEGRAL_TASK);
                this->process_task<loop_momentum_work_record>();
                break;
              }
    
            case MPI_detail::MESSAGE_NEW_MATSUBARA_XY_TASK:
              {
                this->mpi_world.recv(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_NEW_MATSUBARA_XY_TASK);
                this->process_task<Matsubara_XY_work_record>();
                break;
              }
            
            case MPI_detail::MESSAGE_NEW_ONE_LOOP_PK_TASK:
              {
                this->mpi_world.recv(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_NEW_ONE_LOOP_PK_TASK);
                this->process_task<one_loop_Pk_work_record>();
                break;
              }
    
            case MPI_detail::MESSAGE_NEW_ONE_LOOP_RESUM_PK_TASK:
              {
                this->mpi_world.recv(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_NEW_ONE_LOOP_RESUM_PK_TASK);
                this->process_task<one_loop_resum_Pk_work_record>();
                break;
              }
            
            case MPI_detail::MESSAGE_NEW_MULTIPOLE_PK_TASK:
              {
                this->mpi_world.recv(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_NEW_MULTIPOLE_PK_TASK);
                this->process_task<multipole_Pk_work_record>();
                break;
              }

            case MPI_detail::MESSAGE_TERMINATE:
              {
                this->mpi_world.recv(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_TERMINATE);
                finished = true;
                break;
              }

            default:
              {
                assert(false);
              }
          }
      }
  }


template <typename WorkItem>
void slave_controller::process_task()
  {
    // pass an acknowledgment to the master process
    this->mpi_world.isend(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_WORKER_READY);

    bool finished = false;
    
    while(!finished)
      {
        // wait for a message from the master node
        boost::mpi::status stat = this->mpi_world.probe(MPI_detail::RANK_MASTER);

        switch(stat.tag())
          {
            case MPI_detail::work_item_traits<WorkItem>::new_item_message():
              {
                typename MPI_detail::work_item_traits<WorkItem>::outgoing_payload_type payload;
                this->mpi_world.recv(stat.source(), MPI_detail::work_item_traits<WorkItem>::new_item_message(), payload);
                this->process_item(payload);
                break;
              }

            case MPI_detail::MESSAGE_END_OF_WORK:
              {
                this->mpi_world.recv(stat.source(), MPI_detail::MESSAGE_END_OF_WORK);
                finished = true;

                // send acknowledgement to master
                this->mpi_world.isend(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_END_OF_WORK_ACK);
                break;
              }

            default:
              {
                assert(false);
              }
          }
      }
  }


void slave_controller::process_item(MPI_detail::new_transfer_integration& payload)
  {
    const FRW_model& model = payload.get_model();
    const Mpc_units::energy& k = payload.get_k();
    const k_token& tok = payload.get_token();
    const z_database& z_db = payload.get_z_db();
    
//    std::cout << "Worker " << this->worker_number() << " processing transfer item: id = " << tok.get_id() << " for k = "
//              << k * Mpc_units::Mpc << " h/Mpc = " << k / Mpc_units::eV << " eV" << '\n';

    transfer_integrator integrator;
    transfer_function sample = integrator.integrate(model, k, tok, z_db);

    // inform master process that we have completed work on this integration
    MPI_detail::transfer_integration_ready return_payload(sample);
    boost::mpi::request ack = this->mpi_world.isend(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_WORK_PRODUCT_READY, return_payload);
    ack.wait();
  }


void slave_controller::process_item(MPI_detail::new_filter_Pk& payload)
  {
    const FRW_model& model = payload.get_model();
    const Mpc_units::energy& k = payload.get_k();
    const filterable_Pk& Pk_lin = payload.get_Pk_linear();
    
    const k_token& k_tok = payload.get_k_token();
    const linear_Pk_token& Pk_tok = payload.get_Pk_token();
    
    filtered_Pk_value sample;
    try
      {
        Pk_filter filter;
        auto out = filter(model, Pk_lin, k);
        sample = filtered_Pk_value(k_tok, Pk_tok, out.first, Pk_lin(k), out.second);
      }
    catch(runtime_exception& xe)
      {
        if(xe.get_exception_code() == exception_type::filter_failure)
          {
            std::cerr << "lsseft: " << xe.what() << '\n';
            sample = filtered_Pk_value(k_tok, Pk_tok, 0.0, Pk_lin(k), 0.0);
            sample.mark_failed();
          }
        else
          {
            throw;
          }
      }
    
    // inform master process we have finished work on this calculation
    MPI_detail::filter_Pk_ready return_payload(sample);
    boost::mpi::request ack = this->mpi_world.isend(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_WORK_PRODUCT_READY, return_payload);
    ack.wait();
  }


void slave_controller::process_item(MPI_detail::new_loop_momentum_integration& payload)
  {
    const FRW_model& model = payload.get_model();
    const Mpc_units::energy& k = payload.get_k();
    const Mpc_units::energy& UV_cutoff = payload.get_UV_cutoff();
    const Mpc_units::energy& IR_cutoff = payload.get_IR_cutoff();

    const k_token& k_tok = payload.get_k_token();
    const UV_cutoff_token& UV_tok = payload.get_UV_token();
    const IR_cutoff_token& IR_tok = payload.get_IR_token();
    const initial_filtered_Pk& Pk = payload.get_tree_power_spectrum();
    
//    std::cout << "Worker " << this->worker_number() << " processing loop integral item: k-id = " << k_tok.get_id()
//              << " for k = " << k * Mpc_units::Mpc << " h/Mpc, IR cutoff = " << IR_cutoff * Mpc_units::Mpc
//              << " h/Mpc, UV cutoff = " << UV_cutoff * Mpc_units::Mpc << " h/Mpc" << '\n';

    oneloop_momentum_integrator integrator;
    loop_integral sample = integrator.integrate(model, k, k_tok, UV_cutoff, UV_tok, IR_cutoff, IR_tok, Pk);

    // inform master process that we have completed work on this integration
    MPI_detail::loop_momentum_integration_ready return_payload(sample);
    boost::mpi::request ack = this->mpi_world.isend(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_WORK_PRODUCT_READY, return_payload);
    ack.wait();
  }


void slave_controller::process_item(MPI_detail::new_Matsubara_XY& payload)
  {
    const Mpc_units::energy& IR_resum = payload.get_IR_resum();
    const initial_filtered_Pk& Pk = payload.get_tree_power_spectrum();
    
    const IR_resum_token& IR_resum_tok = payload.get_IR_resum_token();
    
    Matsubara_XY_calculator calculator;
    Matsubara_XY item = calculator.calculate_Matsubara_XY(IR_resum, IR_resum_tok, Pk);
    
    // inform master process that the calculation is finished
    MPI_detail::Matsubara_XY_ready return_payload(item);
    boost::mpi::request ack = this->mpi_world.isend(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_WORK_PRODUCT_READY, return_payload);
    ack.wait();
  }


void slave_controller::process_item(MPI_detail::new_one_loop_Pk& payload)
  {
    const Mpc_units::energy& k = payload.get_k();
    const oneloop_growth& gf_factors = payload.get_gf_factors();
    const loop_integral& loop_data = payload.get_loop_data();
    const initial_filtered_Pk& Pk = payload.get_tree_power_spectrum();
    
    const k_token& k_tok = loop_data.get_k_token();
    const IR_cutoff_token& IR_tok = loop_data.get_IR_token();
    const UV_cutoff_token& UV_tok = loop_data.get_UV_token();
    
//    std::cout << "Worker " << this->worker_number() << " processing 1-loop P(k) for"
//              << " k-id = " << k_tok.get_id() << ", IR-id = " << IR_tok.get_id() << ", UV-id = " << UV_tok.get_id()
//              << "; " << gf_factors.size() << " redshifts to process" << '\n';
    
    oneloop_Pk_calculator calculator;
    std::list<oneloop_Pk> sample = calculator.calculate_dd(k, k_tok, IR_tok, UV_tok, gf_factors, loop_data, Pk);
    
    // inform master process that the calculation is finished
    std::list<boost::mpi::request> acks;
    for(oneloop_Pk& record : sample)
      {
        MPI_detail::one_loop_Pk_ready return_payload(record);
        acks.push_back(this->mpi_world.isend(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_WORK_PRODUCT_READY, return_payload));
      }
    boost::mpi::wait_all(acks.begin(), acks.end());
  }


void slave_controller::process_item(MPI_detail::new_one_loop_resum_Pk& payload)
  {
    const Mpc_units::energy& k = payload.get_k();
    const Matsubara_XY& XY = payload.get_Matsubara_XY();
    const oneloop_Pk& oneloop_data = payload.get_oneloop_Pk_data();
    const oneloop_growth_record& gf_data = payload.get_gf_data();
    const initial_filtered_Pk& Pk = payload.get_tree_power_spectrum();
    
    oneloop_Pk_calculator calculator;
    oneloop_resum_Pk sample = calculator.calculate_resum_dd(k, XY, oneloop_data, gf_data, Pk);
    
    // inform master process that the calculation is finished
    MPI_detail::one_loop_resum_Pk_ready return_payload(sample);
    boost::mpi::request ack = this->mpi_world.isend(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_WORK_PRODUCT_READY, return_payload);
    ack.wait();
  }


void slave_controller::process_item(MPI_detail::new_multipole_Pk& payload)
  {
    const Mpc_units::energy& k = payload.get_k();
    const Matsubara_XY& XY = payload.get_Matsubara_XY();
    const oneloop_Pk& oneloop_data = payload.get_oneloop_Pk_data();
    const oneloop_growth_record& gf_data = payload.get_gf_data();
    const initial_filtered_Pk& Pk = payload.get_tree_power_spectrum();
    
//    std::cout << "Worker " << this->worker_number() << " processing multipole P(k) for"
//              << " k-id = " << oneloop_data.get_k_token().get_id()
//              << ", z-id = " << oneloop_data.get_z_token().get_id()
//              << ", IR-id = " << oneloop_data.get_IR_token().get_id()
//              << ", UV-id = " << oneloop_data.get_UV_token().get_id()
//              << ", IR-resum-id = " << IR_resum_tok.get_id() << '\n';
    
    multipole_Pk_calculator calculator;
    multipole_Pk sample = calculator.calculate_Legendre(k, XY, oneloop_data, gf_data, Pk);
    
    // inform master process that the calculation is finished
    MPI_detail::multipole_Pk_ready return_payload(sample);
    boost::mpi::request ack = this->mpi_world.isend(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_WORK_PRODUCT_READY, return_payload);
    ack.wait();
  }
