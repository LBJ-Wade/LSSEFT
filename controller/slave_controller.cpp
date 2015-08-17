//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "slave_controller.h"

#include "cosmology/transfer_integrator.h"


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
              this->mpi_world.recv(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_NEW_TRANSFER_TASK);
              this->process_transfer_task();
              break;

            case MPI_detail::MESSAGE_TERMINATE:
              this->mpi_world.recv(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_TERMINATE);
              finished = true;
              break;

            default:
              assert(false);
          }
      }
  }


void slave_controller::process_transfer_task()
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
            case MPI_detail::MESSAGE_NEW_TRANSFER_INTEGRATION:
              {
                MPI_detail::new_transfer_integration payload;
                this->mpi_world.recv(stat.source(), MPI_detail::MESSAGE_NEW_TRANSFER_INTEGRATION, payload);
                this->transfer_integration(payload);
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


void slave_controller::transfer_integration(MPI_detail::new_transfer_integration& payload)
  {
    FRW_model model = payload.get_model();
    eV_units::energy k = payload.get_k();
    wavenumber_token tok = payload.get_token();
    std::shared_ptr<redshift_database> z_db = payload.get_z_db();

    std::cout << "Worker " << this->worker_number() << " beginning transfer task: id = " << tok.get_id() << " for k = " << k * eV_units::Mpc << " h/Mpc = " << static_cast<double>(k) << " eV" << '\n';

    transfer_integrator integrator;
    transfer_function sample = integrator.integrate(model, k, tok, z_db);

    // inform master process that we have completed work on this integration
    MPI_detail::transfer_integration_ready return_payload(sample);
    boost::mpi::request ack = this->mpi_world.isend(MPI_detail::RANK_MASTER, MPI_detail::MESSAGE_TRANSFER_INTEGRATION_READY, return_payload);
  }
