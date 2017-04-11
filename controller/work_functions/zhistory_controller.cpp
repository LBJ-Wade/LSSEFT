//
// Created by David Seery on 11/04/2017.
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

#include "controller/master_controller.h"
#include "cosmology/models/Planck_defaults.h"


void master_controller::execute()
  {
    if(!this->arg_cache.is_database_set())
      {
        this->err_handler.error(ERROR_NO_DATABASE);
        return;
      }
    
    // set up
    data_manager dmgr(this->arg_cache.get_database_path(), this->err_handler, this->arg_cache);
    
    // fix the background cosmological model
    // here, that's taken to have parameters matching the MDR1 simulation
    FRW_model cosmology_model;
    std::unique_ptr<FRW_model_token> model = dmgr.tokenize(cosmology_model);
    
    // set up a list of redshifts at which to sample the late-time growth functions
    stepping_range<double> lo_redshift_samples(0.001, 50.0, 500, 1.0, spacing_type::logarithmic_bottom);
    
    std::unique_ptr<z_database> lo_z_db = dmgr.build_redshift_db(lo_redshift_samples);
    
    // SET UP PARAMETERS
    
    // set up parameters for growth function
    // allow specification of full or EdS growth function on command line, but always use EdS ics
    growth_params Df_params(this->arg_cache.use_EdS());
    std::unique_ptr<growth_params_token> growth_tok = dmgr.tokenize(Df_params);
    
    
    // GENERATE TARGETS
    
    std::unique_ptr<z_database> loop_growth_work = dmgr.build_loop_growth_work_list(*model, *lo_z_db, *growth_tok);
    
    // compute linear and one-loop growth functions, if needed; can be done on master process since there is only one integration
    if(loop_growth_work) this->integrate_loop_growth(cosmology_model, *model, *loop_growth_work, dmgr, *growth_tok, Df_params);
    
    // instruct slave processes to terminate
    this->terminate_workers();
  }
