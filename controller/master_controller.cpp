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


#include "core.h"

#include "master_controller.h"

#include "utilities/formatter.h"

#include "localizations/messages.h"

#include "boost/program_options.hpp"
#include "boost/optional.hpp"


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
      (LSSEFT_SWITCH_DATABASE, boost::program_options::value<std::string>(), LSSEFT_HELP_DATABASE)
      (LSSEFT_SWITCH_INITIAL_POWERSPEC, boost::program_options::value<std::string>(), LSSEFT_HELP_INITIAL_POWERSPEC)
      (LSSEFT_SWITCH_FINAL_POWERSPEC, boost::program_options::value<std::string>(), LSSEFT_HELP_FINAL_POWERSPEC)
      (LSSEFT_SWITCH_EDS_MODE, LSSEFT_HELP_EDS_MODE);

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

    if(option_map.count(LSSEFT_SWITCH_INITIAL_POWERSPEC_LONG))
      {
        this->arg_cache.set_initial_powerspectrum_path(option_map[LSSEFT_SWITCH_INITIAL_POWERSPEC_LONG].as<std::string>());
      }
    
    if(option_map.count(LSSEFT_SWITCH_FINAL_POWERSPEC_LONG))
      {
        this->arg_cache.set_final_powerspectrum_path(option_map[LSSEFT_SWITCH_FINAL_POWERSPEC_LONG].as<std::string>());
      }
    
    if(option_map.count(LSSEFT_SWITCH_EDS_MODE)) this->arg_cache.set_EdS_mode(true);
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
    std::unique_ptr<scheduler> sch = std::make_unique<scheduler>(this->mpi_world.size() - 1);

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
                break;
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


void master_controller::integrate_loop_growth(const FRW_model& model, const FRW_model_token& token, z_database& z_db, data_manager& dmgr,
                                              const growth_params_token& params_tok, const growth_params& params)
  {
    dmgr.setup_growth_write();
    
    oneloop_growth_integrator integrator(params, params_tok);
    growth_integrator_data data = integrator.integrate(model, z_db);
    dmgr.store(token, *data.container);
    
    dmgr.finalize_growth_write();
  }
