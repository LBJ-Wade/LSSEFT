//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include "core.h"

#include "master_controller.h"
#include "scheduler.h"

#include "cosmology/FRW_model.h"
#include "cosmology/concepts/range.h"

#include "database/data_manager.h"

#include "MPI_detail/mpi_operations.h"

#include "utilities/formatter.h"

#include "localizations/messages.h"

#include "boost/program_options.hpp"


master_controller::master_controller(boost::mpi::environment& me, boost::mpi::communicator& mw, argument_cache& ac)
  : mpi_env(me),
    mpi_world(mw),
    arg_cache(ac),
    local_env(),
    err_handler(arg_cache, local_env)
  {
  }


void master_controller::process_arguments(int argc, char* argv[])
  {
    // set up BOOST::program_options descriptors for command-line arguments
    boost::program_options::options_description generic("Generic");
    generic.add_options()
      (LSSEFT_SWITCH_HELP, LSSEFT_HELP_HELP)
      (LSSEFT_SWITCH_VERSION, LSSEFT_HELP_VERSION)
      (LSSEFT_SWITCH_NO_COLOUR, LSSEFT_HELP_NO_COLOUR);

    boost::program_options::options_description configuration("Configuration options");
    configuration.add_options()
      (LSSEFT_SWITCH_VERBOSE, LSSEFT_HELP_VERBOSE)
      (LSSEFT_SWITCH_DATABASE, boost::program_options::value<std::string>(), LSSEFT_HELP_DATABASE);

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      (LSSEFT_SWITCH_NO_COLOR, LSSEFT_HELP_NO_COLOR);

    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(configuration).add(hidden);

    boost::program_options::options_description output_options;
    output_options.add(generic).add(configuration);

    boost::program_options::variables_map option_map;
    boost::program_options::store(boost::program_options::parse_command_line(argc, argv, cmdline_options), option_map);
    boost::program_options::notify(option_map);

    bool emitted_version = false;

    if(option_map.count(LSSEFT_SWITCH_VERSION))
      {
        std::cout << LSSEFT_NAME << " " << LSSEFT_VERSION << " " << LSSEFT_COPYRIGHT << '\n';
        emitted_version = true;
      }

    if(option_map.count(LSSEFT_SWITCH_HELP))
      {
        if(!emitted_version) std::cout << LSSEFT_NAME << " " << LSSEFT_VERSION << " " << LSSEFT_COPYRIGHT << '\n';
        std::cout << output_options << '\n';
      }

    if(option_map.count(LSSEFT_SWITCH_VERBOSE_LONG)) this->arg_cache.set_verbose(true);
    if(option_map.count(LSSEFT_SWITCH_NO_COLOUR) || option_map.count(LSSEFT_SWITCH_NO_COLOR))
      this->arg_cache.set_colour_output(false);

    if(option_map.count(LSSEFT_SWITCH_DATABASE_LONG))
      {
        this->arg_cache.set_database_path(option_map[LSSEFT_SWITCH_DATABASE_LONG].as<std::string>());
      }
  }


void master_controller::execute()
  {
    if(!this->arg_cache.get_database_set())
      {
        this->err_handler.error(ERROR_NO_DATABASE);
        return;
      }

    // set up
    data_manager dmgr(this->arg_cache.get_database_path());
    FRW_model cosmology_model;

    std::unique_ptr<FRW_model_token> model = dmgr.tokenize(cosmology_model);

    // set up a list of wavenumber to sample, measured in h/Mpc
    stepping_range<eV_units::energy> wavenumber_samples(0.05, 0.3, 30, 1.0 / eV_units::Mpc, spacing_type::linear);

    // set up a list of redshifts at which to sample
    stepping_range<double> redshift_samples(0.01, 1000.0, 250, 1.0, spacing_type::logarithmic_bottom);

    // exchange these sample ranges for iterable databases
    std::unique_ptr<redshift_database> z_db = dmgr.build_db(redshift_samples);
    std::unique_ptr<wavenumber_database> k_db = dmgr.build_db(wavenumber_samples);

    // generate targets
    // for 1-loop calculation, we need the linear transfer function at the initial redshift,
    // plus the time-dependent 1-loop kernels through the subsequent evolution

    // build a work list for transfer functions
    std::unique_ptr<transfer_work_list> work_list = dmgr.build_transfer_work_list(*model, *k_db, *z_db);

    // distribute this work list among the slave processes
    this->scatter(cosmology_model, *model, *work_list, dmgr);

    // instruct slave processes to terminate
    this->terminate_workers();
  }


void master_controller::scatter(const FRW_model& model, const FRW_model_token& token, transfer_work_list& work,
                                data_manager& dmgr)
  {
    if(this->mpi_world.size() == 1) throw runtime_exception(exception_type::runtime_error, ERROR_TOO_FEW_WORKS);

    // instruct slave processes to await transfer function tasks
    std::unique_ptr<scheduler> sch = this->set_up_workers(MPI_detail::MESSAGE_NEW_TRANSFER_TASK);

    bool sent_closedown = false;
    transfer_work_list::iterator next_work_item = work.begin();

    while(!sch->all_inactive())
      {
        // check whether all work is exhausted
        if(next_work_item == work.end() && !sent_closedown)
          {
            sent_closedown = true;
            this->close_down_workers();
          }

        // check whether any workers are waiting for assignments
        if(next_work_item != work.end() && sch->is_assignable())
          {
            std::vector<unsigned int> unassigned_list = sch->make_assignment();
            std::vector<boost::mpi::request> requests;

            for(std::vector<unsigned int>::const_iterator t = unassigned_list.begin();
                next_work_item != work.end() && t != unassigned_list.end(); ++t)
              {
                // assign next work item to this worker
                MPI_detail::new_transfer_integration payload(model, *(*next_work_item), next_work_item->get_token(),
                                                             next_work_item->get_z_db());
                requests.push_back(this->mpi_world.isend(this->worker_rank(*t), MPI_detail::MESSAGE_NEW_TRANSFER_INTEGRATION, payload));

                sch->mark_assigned(*t);
                ++next_work_item;
              }

            // wait for all messages to be received
            boost::mpi::wait_all(requests.begin(), requests.end());
          }

        // check whether any messages are waiting in the queue
        boost::optional<boost::mpi::status> stat = this->mpi_world.iprobe();

        while(stat) // consume messages until no more are available
          {
            switch(stat->tag())
              {
                case MPI_detail::MESSAGE_TRANSFER_INTEGRATION_READY:
                  {
                    MPI_detail::transfer_integration_ready payload;
                    this->mpi_world.recv(stat->source(), MPI_detail::MESSAGE_TRANSFER_INTEGRATION_READY, payload);
                    sch->mark_unassigned(this->worker_number(stat->source()));
                    dmgr.store(token, payload.get_data());
                    break;
                  }

                case MPI_detail::MESSAGE_END_OF_WORK_ACK:
                  {
                    this->mpi_world.recv(stat->source(), MPI_detail::MESSAGE_END_OF_WORK_ACK);
                    sch->mark_inactive(this->worker_number(stat->source()));
                    break;
                  }

                default:
                  assert(false);
              }

            stat = this->mpi_world.iprobe();
          }
      }
  }


void master_controller::terminate_workers()
  {
    // send terminate message to all workers
    std::vector<boost::mpi::request> requests(this->mpi_world.size() - 1);
    for(unsigned int i = 0; i < this->mpi_world.size() - 1; ++i)
      {
        requests[i] = this->mpi_world.isend(this->worker_rank(i), MPI_detail::MESSAGE_TERMINATE);
      }

    // wait for all messages to be received, then exit ourselves
    boost::mpi::wait_all(requests.begin(), requests.end());
  }


std::unique_ptr<scheduler> master_controller::set_up_workers(unsigned int tag)
  {
    // send specified message to all workers
    std::vector<boost::mpi::request> requests(this->mpi_world.size() - 1);
    for(unsigned int i = 0; i < this->mpi_world.size() - 1; ++i)
      {
        requests[i] = this->mpi_world.isend(this->worker_rank(i), tag);
      }

    // wait for all messages to be received
    boost::mpi::wait_all(requests.begin(), requests.end());

    // create scheduler object
    std::unique_ptr<scheduler> sch(new scheduler(this->mpi_world.size() - 1));

    // wait for messages from workers reporting that they are correctly set up
    while(!sch->is_ready())
      {
        boost::mpi::status stat = this->mpi_world.probe();

        switch(stat.tag())
          {
            case MPI_detail::MESSAGE_WORKER_READY:
              {
                this->mpi_world.recv(stat.source(), MPI_detail::MESSAGE_WORKER_READY);
                sch->initialize_worker(this->worker_number(stat.source()));
              }
          }
      }

    return (sch);
  }


void master_controller::close_down_workers()
  {
    // send end-of-work message to all workers
    std::vector<boost::mpi::request> requests(this->mpi_world.size() - 1);
    for(unsigned int i = 0; i < this->mpi_world.size() - 1; ++i)
      {
        requests[i] = this->mpi_world.isend(this->worker_rank(i), MPI_detail::MESSAGE_END_OF_WORK);
      }

    // wait for all messages to be received
    boost::mpi::wait_all(requests.begin(), requests.end());
  }
