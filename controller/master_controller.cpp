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


#include <MPI_detail/mpi_traits.h>
#include "core.h"

#include "master_controller.h"
#include "scheduler.h"

#include "cosmology/FRW_model.h"
#include "cosmology/oneloop_growth_integrator.h"
#include "cosmology/concepts/range.h"
#include "cosmology/concepts/power_spectrum.h"

#include "database/data_manager.h"

#include "MPI_detail/mpi_traits.h"
#include "MPI_detail/mpi_payloads.h"

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
      (LSSEFT_SWITCH_DATABASE, boost::program_options::value<std::string>(), LSSEFT_HELP_DATABASE)
      (LSSEFT_SWITCH_INITIAL_POWERSPEC, boost::program_options::value<std::string>(), LSSEFT_HELP_INITIAL_POWERSPEC)
      (LSSEFT_SWITCH_FINAL_POWERSPEC, boost::program_options::value<std::string>(), LSSEFT_HELP_FINAL_POWERSPEC);

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

    // fix the background cosmological model
    FRW_model cosmology_model;
    std::unique_ptr<FRW_model_token> model = dmgr.tokenize(cosmology_model);

    // set up a list of wavenumbers to sample for the transfer functions, measured in h/Mpc
    stepping_range<Mpc_units::energy> transfer_k_samples(0.01, 0.8, 30, 1.0 / Mpc_units::Mpc, spacing_type::linear);

    // set up a list of redshifts at which to sample the transfer functions
    stepping_range<double> hi_redshift_samples(1000.0, 1500.0, 5, 1.0, spacing_type::linear);

    // set up a list of redshifts at which to sample the late-time growth functions
    stepping_range<double> lo_redshift_samples(0.0, 50.0, 250, 1.0, spacing_type::logarithmic_bottom);

    // set up a list of UV cutoffs, measured in h/Mpc, to be used with the loop integrals
    stepping_range<Mpc_units::energy> UV_cutoffs(0.8, 0.8, 0, 1.0 / Mpc_units::Mpc, spacing_type::logarithmic_bottom);

    // set up a list of IR cutoffs, measured in h/Mpc, to be used with the loop integrals
    stepping_range<Mpc_units::energy> IR_cutoffs(1E-4, 1E-4, 0, 1.0 / Mpc_units::Mpc, spacing_type::logarithmic_bottom);

    // set up a list of k at which to compute the loop integrals
    stepping_range<Mpc_units::energy> loop_k_samples(0.005, 0.5, 200, 1.0 / Mpc_units::Mpc, spacing_type::logarithmic_bottom);
    
    // set up a list of IR resummation scales, measured in h/Mpc
    stepping_range<Mpc_units::energy> IR_resummation(0.3, 0.3, 0, 1.0 / Mpc_units::Mpc, spacing_type::linear);

    // exchange these sample ranges for iterable databases
    std::unique_ptr<z_database> hi_z_db              = dmgr.build_redshift_db(hi_redshift_samples);
    std::unique_ptr<z_database> lo_z_db              = dmgr.build_redshift_db(lo_redshift_samples);
    std::unique_ptr<k_database> transfer_k_db        = dmgr.build_k_db(transfer_k_samples);

    std::unique_ptr<UV_cutoff_database> UV_cutoff_db = dmgr.build_UV_cutoff_db(UV_cutoffs);
    std::unique_ptr<IR_cutoff_database> IR_cutoff_db = dmgr.build_IR_cutoff_db(IR_cutoffs);
    std::unique_ptr<IR_resum_database>  IR_resum_db  = dmgr.build_IR_resum_db(IR_resummation);
    std::unique_ptr<k_database>         loop_k_db    = dmgr.build_k_db(loop_k_samples);

    
    // GENERATE TARGETS

    // for 1-loop calculation, we need the linear transfer function at the initial redshift,
    // (although we can also ingest a linear power spectrum generated by some other means, eg. with CAMB; see below)
    // plus the time-dependent 1-loop kernels through the subsequent evolution

    // build a work list for transfer functions
    std::unique_ptr<transfer_work_list> transfer_work = dmgr.build_transfer_work_list(*model, *transfer_k_db, *hi_z_db);
    
    // distribute this work list among the worker processes
    if(transfer_work) this->scatter(cosmology_model, *model, *transfer_work, dmgr);

    // build a work list for linear and one-loop growth functions (and their growth rates);
    // we inherit ownership of its lifetime using std::unique_ptr<>
    std::unique_ptr<z_database> loop_growth_work = dmgr.build_oneloop_work_list(*model, *lo_z_db);

    // compute linear and one-loop growth functions, if needed; can be done on master process since there is only one integration
    if(loop_growth_work) this->integrate_oneloop(cosmology_model, *model, *loop_growth_work, dmgr);
    
    if(this->arg_cache.is_initial_powerspectrum_set())
      {
        // STEP 1 - READ IN AND FILTER LINEAR POWER SPECTRA
        
        // read in initial linear power spectrum in CAMB format, and ask the database to tokenize it
        // we manage its lifetime using std::shared_ptr<> because we want to share ownership with
        // the MPI work records and payloads
        std::shared_ptr<initial_Pk> init_Pk_lin_db = std::make_shared<initial_Pk>(this->arg_cache.get_initial_powerspectrum_path());
        std::unique_ptr<linear_Pk_token> init_Pk_lin = dmgr.tokenize(*model, *init_Pk_lin_db);
    
        // build a work list for filtering the linear power spectrum in wiggle/no-wiggle components
        std::shared_ptr<filterable_Pk> filterable_init_Pk_lin_db = make_filterable(*init_Pk_lin_db);
        auto init_filter_work = dmgr.build_filter_Pk_work_list(*init_Pk_lin, filterable_init_Pk_lin_db);
    
        // distribute this work list among the worker processes
        if(init_filter_work) this->scatter(cosmology_model, *model, *init_filter_work, dmgr);
        
        // exchange our linear power spectrum for the filtered version;
        // we manage its lifetime using std::shared_ptr<> since ownership is shared with the
        // MPI work records and payloads
        std::shared_ptr<initial_filtered_Pk> initial_Pk_wig = dmgr.build_wiggle_Pk(*init_Pk_lin, *init_Pk_lin_db);
    
        // if a final linear power spectrum has been specified, then read this in too, tokenize it,
        // and obtain a container
        std::shared_ptr<final_filtered_Pk> final_Pk_wig;
        if(this->arg_cache.is_final_powerspectrum_set())
          {
            // read in final power spectrum and tokenize it
            std::shared_ptr<final_Pk> final_Pk_lin_db = std::make_shared<final_Pk>(this->arg_cache.get_final_powerspectrum_path());
            std::unique_ptr<linear_Pk_token> final_Pk_lin = dmgr.tokenize(*model, *final_Pk_lin_db);
            
            // build a work list for filtering
            std::shared_ptr<filterable_Pk> filterable_final_Pk_lin_db = make_filterable(*final_Pk_lin_db);
            auto final_filter_work = dmgr.build_filter_Pk_work_list(*final_Pk_lin, filterable_final_Pk_lin_db);
            
            // distribute this work list among the worker processes
            if(final_filter_work) this->scatter(cosmology_model, *model, *final_filter_work, dmgr);
            
            // get filtered version of power spectrum
            final_Pk_wig = dmgr.build_wiggle_Pk(*final_Pk_lin, *final_Pk_lin_db);
          }

        
        // STEP 2 - COMPUTE LOOP INTEGRALS

        // build a work list for the loop integrals
        std::unique_ptr<loop_momentum_work_list> loop_momentum_work =
          dmgr.build_loop_momentum_work_list(*model, *loop_k_db, *IR_cutoff_db, *UV_cutoff_db, initial_Pk_wig);

        // distribute this work list among the worker processes
        if(loop_momentum_work) this->scatter(cosmology_model, *model, *loop_momentum_work, dmgr);
    
        
        // STEP 3 - COMPUTE MATSUBARA X & Y COEFFICIENTS
        
        // build a work list of the Matsubara X & Y coefficients associated with these resummation scales
        std::unique_ptr<Matsubara_XY_work_list> Matsubara_work =
          dmgr.build_Matsubara_XY_work_list(*model, *IR_resum_db, initial_Pk_wig);
    
        // distribute this work list among the worker processes
        if(Matsubara_work) this->scatter(cosmology_model, *model, *Matsubara_work, dmgr);
        
        
        // STEP 4 - COMPUTE ONE-LOOP POWER SPECTRA IN REAL AND REDSHIFT SPACE
        
        // build a work list for the individual power spectrum components
        std::unique_ptr<one_loop_Pk_work_list> Pk_work =
          dmgr.build_one_loop_Pk_work_list(*model, *lo_z_db, *loop_k_db, *IR_cutoff_db, *UV_cutoff_db, initial_Pk_wig);

        // distribute this work list among the worker processes
        if(Pk_work) this->scatter(cosmology_model, *model, *Pk_work, dmgr);
        
        // build a work list for the resummed power spectrum components
        std::unique_ptr<one_loop_resum_Pk_work_list> Pk_resum_work =
          dmgr.build_one_loop_resum_Pk_work_list(*model, *lo_z_db, *loop_k_db, *IR_cutoff_db, *UV_cutoff_db, *IR_resum_db, initial_Pk_wig);
        
        // distribute this work list among the worker processes
        if(Pk_resum_work) this->scatter(cosmology_model, *model, *Pk_resum_work, dmgr);
        
        
        // STEP 5 - COMPUTE MULTIPOLE DECOMPOSIITON OF REDSHIFT-SPACE POWER SPECTRUM
        
        // build a work list for the resummed multipole power spectra
        std::unique_ptr<multipole_Pk_work_list> multipole_Pk_work =
          dmgr.build_multipole_Pk_work_list(*model, *lo_z_db, *loop_k_db, *IR_cutoff_db, *UV_cutoff_db, *IR_resum_db, initial_Pk_wig);
    
        // distribute this work list among the worker processes
        if(multipole_Pk_work) this->scatter(cosmology_model, *model, *multipole_Pk_work, dmgr);
      }

    // instruct slave processes to terminate
    this->terminate_workers();
  }


template <typename WorkItemList>
void master_controller::scatter(const FRW_model& model, const FRW_model_token& token, WorkItemList& work, data_manager& dmgr)
  {
    using WorkItem = typename WorkItemList::value_type;
    
    boost::timer::cpu_timer timer;

    if(this->mpi_world.size() == 1) throw runtime_exception(exception_type::runtime_error, ERROR_TOO_FEW_WORKERS);

    // instruct slave processes to await transfer function tasks
    std::unique_ptr<scheduler> sch = this->set_up_workers(MPI_detail::work_item_traits<WorkItem>::new_task_message());

    bool sent_closedown = false;
    auto next_work_item = work.cbegin();

    while(!sch->all_inactive())
      {
        // check whether all work is exhausted
        if(next_work_item == work.cend() && !sent_closedown)
          {
            sent_closedown = true;
            this->close_down_workers();
          }

        // check whether any workers are waiting for assignments
        if(next_work_item != work.cend() && sch->is_assignable())
          {
            std::vector<unsigned int> unassigned_list = sch->make_assignment();
            std::vector<boost::mpi::request> requests;

            for(std::vector<unsigned int>::const_iterator t = unassigned_list.begin();
                next_work_item != work.cend() && t != unassigned_list.end(); ++t)
              {
                // assign next work item to this worker
                requests.push_back(this->mpi_world.isend(this->worker_rank(*t),
                                                         MPI_detail::work_item_traits<WorkItem>::new_item_message(),
                                                         MPI_detail::build_payload(model, next_work_item)));

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
                case MPI_detail::MESSAGE_WORK_PRODUCT_READY:
                  {
                    this->store_payload<WorkItem>(token, stat->source(), dmgr);
                    sch->mark_unassigned(this->worker_number(stat->source()));
                    break;
                  }

                case MPI_detail::MESSAGE_END_OF_WORK_ACK:
                  {
                    this->mpi_world.recv(stat->source(), MPI_detail::MESSAGE_END_OF_WORK_ACK);
                    sch->mark_inactive(this->worker_number(stat->source()));
                    break;
                  }

                default:
                  {
                    assert(false);
                  }
              }

            stat = this->mpi_world.iprobe();
          }
      }

    timer.stop();
    std::cout << "lsseft: completed work in time " << format_time(timer.elapsed().wall) << '\n';
  }


template <typename WorkItem>
void master_controller::store_payload(const FRW_model_token& token, unsigned int source, data_manager& dmgr)
  {
    typename MPI_detail::work_item_traits<WorkItem>::incoming_payload_type payload;

    this->mpi_world.recv(source, MPI_detail::MESSAGE_WORK_PRODUCT_READY, payload);
    dmgr.store(token, payload.get_data());
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


void master_controller::integrate_oneloop(const FRW_model& model, const FRW_model_token& token, z_database& z_db,
                                          data_manager& dmgr)
  {
    oneloop_growth_integrator integrator;

    growth_integrator_data data = integrator.integrate(model, z_db);

    dmgr.store(token, *data.container);
  }
