//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "master_controller.h"

#include "core.h"
#include "localizations/messages.h"

#include "cosmology/FRW_model.h"
#include "cosmology/concepts/range.h"

#include "database/data_manager.h"

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
      (LSSEFT_SWITCH_HELP,      LSSEFT_HELP_HELP)
      (LSSEFT_SWITCH_VERSION,   LSSEFT_HELP_VERSION)
      (LSSEFT_SWITCH_NO_COLOUR, LSSEFT_HELP_NO_COLOUR);

    boost::program_options::options_description configuration("Configuration options");
    configuration.add_options()
      (LSSEFT_SWITCH_VERBOSE,                                                   LSSEFT_HELP_VERBOSE)
      (LSSEFT_SWITCH_DATABASE,  boost::program_options::value< std::string >(), LSSEFT_HELP_DATABASE);

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

    if(option_map.count(LSSEFT_SWITCH_VERBOSE_LONG))                                          this->arg_cache.set_verbose(true);
    if(option_map.count(LSSEFT_SWITCH_NO_COLOUR) || option_map.count(LSSEFT_SWITCH_NO_COLOR)) this->arg_cache.set_colour_output(false);

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
    stepping_range<eV_units::energy> wavenumber_samples(0.05, 0.3, 30, 1.0/eV_units::Mpc, spacing_type::linear);

    // set up a list of redshifts at which to sample
    stepping_range<double> redshift_samples(0.01, 1000.0, 250, 1.0, spacing_type::logarithmic_bottom);

    // exchange these sample ranges for iterable databases
    std::unique_ptr<redshift_database>   z_db = dmgr.build_db(redshift_samples);
    std::unique_ptr<wavenumber_database> k_db = dmgr.build_db(wavenumber_samples);

    // generate targets
    // for 1-loop calculation, we need the linear transfer function at the initial redshift,
    // plus the time-dependent 1-loop kernels through the subsequent evolution

    // build a work list for transfer functions
    std::unique_ptr<transfer_work_list> work_list = dmgr.build_transfer_work_list(*model, *k_db, *z_db);

    // distribute this work list among the slave processes
    this->scatter(*work_list);
  }


void master_controller::scatter(transfer_work_list& work)
  {

  }
