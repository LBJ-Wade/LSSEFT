//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include "MPI_detail/mpi_traits.h"
#include "MPI_detail/mpi_payloads.h"

#include "slave_controller.h"

#include "cosmology/transfer_integrator.h"
#include "cosmology/oneloop_momentum_integrator.h"
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

            case MPI_detail::MESSAGE_NEW_LOOP_INTEGRAL_TASK:
              {
                this->mpi_world.recv(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_NEW_LOOP_INTEGRAL_TASK);
                this->process_task<loop_momentum_work_record>();
                break;
              }

            case MPI_detail::MESSAGE_TERMINATE:
              {
                this->mpi_world.recv(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_TERMINATE);
                finished = true;
                break;
              }

            default:
              assert(false);
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
              assert(false);
          }
      }
  }


void slave_controller::process_item(MPI_detail::new_transfer_integration& payload)
  {
    FRW_model model = payload.get_model();
    eV_units::energy k = payload.get_k();
    k_token tok = payload.get_token();
    std::shared_ptr<z_database> z_db = payload.get_z_db();

    std::cout << "Worker " << this->worker_number() << " processing transfer item: id = " << tok.get_id() << " for k = " << k * eV_units::Mpc << " h/Mpc = " << static_cast<double>(k) << " eV" << '\n';

    transfer_integrator integrator;
    transfer_function sample = integrator.integrate(model, k, tok, z_db);

    // inform master process that we have completed work on this integration
    MPI_detail::transfer_integration_ready return_payload(sample);
    boost::mpi::request ack = this->mpi_world.isend(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_WORK_PRODUCT_READY, return_payload);
  }


void slave_controller::process_item(MPI_detail::new_loop_momentum_integration& payload)
  {
    FRW_model model = payload.get_model();
    eV_units::energy k = payload.get_k();
    eV_units::energy UV_cutoff = payload.get_UV_cutoff();
    eV_units::energy IR_cutoff = payload.get_IR_cutoff();
    k_token k_tok = payload.get_k_token();
    UV_token UV_tok = payload.get_UV_token();
    IR_token IR_tok = payload.get_IR_token();
    std::shared_ptr<tree_power_spectrum> Pk = payload.get_tree_power_spectrum();

    std::cout << "Worker " << this->worker_number() << " processing loop integral item: id = " << k_tok.get_id() << " for k = " << k * eV_units::Mpc << " h/Mpc, IR cutoff = " << IR_cutoff * eV_units::Mpc << " h/Mpc, UV cutoff = " << UV_cutoff * eV_units::Mpc << " h/Mpc" << '\n';

    oneloop_momentum_integrator integrator;
    loop_integral sample = integrator.integrate(model, k, k_tok, UV_cutoff, UV_tok, IR_cutoff, IR_tok, Pk);

    // inform master process that we have completed work on this integration
    MPI_detail::loop_momentum_integration_ready return_payload(sample);
    boost::mpi::request ack = this->mpi_world.isend(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_WORK_PRODUCT_READY, return_payload);
  }
