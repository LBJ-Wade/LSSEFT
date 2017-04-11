//
// Created by David Seery on 02/04/2017.
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
#include "cosmology/models/MDR1_sim.h"


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
    FRW_model cosmology_model(MDR1::name, MDR1::omega_m, MDR1::omega_cc, MDR1::h, MDR1::T_CMB, MDR1::Neff,
                              MDR1::f_baryon, MDR1::z_star, MDR1::z_drag, MDR1::z_eq, MDR1::Acurv, MDR1::ns, MDR1::kpiv);
    std::unique_ptr<FRW_model_token> model = dmgr.tokenize(cosmology_model);
    
    // set up a list of wavenumbers to sample for the transfer functions, measured in h/Mpc
    // stepping_range<Mpc_units::energy> transfer_k_samples(0.01, 0.8, 30, 1.0 / Mpc_units::Mpc, spacing_type::linear);
    
    // set up a list of redshifts at which to sample the transfer functions
    // stepping_range<double> hi_redshift_samples(1000.0, 1500.0, 5, 1.0, spacing_type::linear);
    
    // set up a list of redshifts at which to sample the late-time growth functions; we need the z=50 point
    // to define where the integration starts
    stepping_range<double> z50(50.0, 50.0, 0, 1.0, spacing_type::linear);
    stepping_range<double> z0(0.0, 0.0, 0, 1.0, spacing_type::linear);
    stepping_range<double> z025(0.25, 0.25, 0, 1.0, spacing_type::linear);
    stepping_range<double> z05(0.5, 0.5, 0, 1.0, spacing_type::linear);
    stepping_range<double> z075(0.75, 0.75, 0, 1.0, spacing_type::linear);
    stepping_range<double> z1(1.0, 1.0, 0, 1.0, spacing_type::linear);
    auto lo_redshift_samples = z0 + z025 + z05 + z075 + z1 + z50;
    
    // set up a list of UV cutoffs, measured in h/Mpc, to be used with the loop integrals
    stepping_range<Mpc_units::energy> UV_cutoffs(1.4, 1.4, 0, 1.0 / Mpc_units::Mpc, spacing_type::logarithmic_bottom);
    
    // set up a list of IR cutoffs, measured in h/Mpc, to be used with the loop integrals
    stepping_range<Mpc_units::energy> IR_cutoffs(1E-4, 1E-4, 0, 1.0 / Mpc_units::Mpc, spacing_type::logarithmic_bottom);
    
    // set up a list of k at which to compute the loop integrals
    stepping_range<Mpc_units::energy> loop_k_samples(0.005, 1.0, 500, 1.0 / Mpc_units::Mpc, spacing_type::logarithmic_bottom);
    
    // set up a list of IR resummation scales, measured in h/Mpc
    stepping_range<Mpc_units::energy> IR_resummation(1.4, 1.4, 0, 1.0 / Mpc_units::Mpc, spacing_type::linear);
    
    // exchange these sample ranges for iterable databases
    // std::unique_ptr<z_database> hi_z_db              = dmgr.build_redshift_db(hi_redshift_samples);
    std::unique_ptr<z_database> lo_z_db              = dmgr.build_redshift_db(lo_redshift_samples);
    // std::unique_ptr<k_database> transfer_k_db        = dmgr.build_k_db(transfer_k_samples);
    
    std::unique_ptr<UV_cutoff_database> UV_cutoff_db = dmgr.build_UV_cutoff_db(UV_cutoffs);
    std::unique_ptr<IR_cutoff_database> IR_cutoff_db = dmgr.build_IR_cutoff_db(IR_cutoffs);
    std::unique_ptr<IR_resum_database>  IR_resum_db  = dmgr.build_IR_resum_db(IR_resummation);
    std::unique_ptr<k_database>         loop_k_db    = dmgr.build_k_db(loop_k_samples);
    
    
    // SET UP PARAMETERS
    
    // set up parameters for filter
    Pk_filter_params filter_params(0.25, 0.05 / Mpc_units::Mpc, 0.02);
    std::unique_ptr<filter_params_token> filter_tok = dmgr.tokenize(filter_params);
    
    // set up parameters for growth function;
    // allow specification of full or EdS growth function on command line, but always use EdS ics
    growth_params Df_params(this->arg_cache.use_EdS());
    std::unique_ptr<growth_params_token> growth_tok = dmgr.tokenize(Df_params);
    
    // set up parameters for Matsubara X&Y integral
    // defauts to qmin = 10 Mpc/h, qmax = 300 Mpc/h
    MatsubaraXY_params XY_params;
    std::unique_ptr<MatsubaraXY_params_token> XY_tok = dmgr.tokenize(XY_params);
    
    // set up parameters for loop integral
    loop_integral_params loop_params;
    std::unique_ptr<loop_integral_params_token> loop_tok = dmgr.tokenize(loop_params);
    
    
    // GENERATE TARGETS
    
    // for 1-loop calculation, we need the linear transfer function at the initial redshift,
    // (although we can also ingest a linear power spectrum generated by some other means, eg. with CAMB; see below)
    // plus the time-dependent 1-loop kernels through the subsequent evolution
    
    // build a work list for transfer functions
    // std::unique_ptr<transfer_work_list> transfer_work = dmgr.build_transfer_work_list(*model, *transfer_k_db, *hi_z_db);
    
    // distribute this work list among the worker processes
    // if(transfer_work) this->scatter(cosmology_model, *model, *transfer_work, dmgr);
    
    // build a work list for linear and one-loop growth functions (and their growth rates);
    // we inherit ownership of its lifetime using std::unique_ptr<>
    std::unique_ptr<z_database> loop_growth_work = dmgr.build_loop_growth_work_list(*model, *lo_z_db, *growth_tok);
    
    // compute linear and one-loop growth functions, if needed; can be done on master process since there is only one integration
    if(loop_growth_work) this->integrate_loop_growth(cosmology_model, *model, *loop_growth_work, dmgr, *growth_tok, Df_params);
    
    if(this->arg_cache.is_initial_powerspectrum_set())
      {
        // STEP 1 - READ IN AND FILTER LINEAR POWER SPECTRA
        
        // read in initial linear power spectrum in CAMB format, and ask the database to tokenize it
        // we manage its lifetime using std::shared_ptr<> because we want to share ownership with
        // the MPI work records and payloads
        std::shared_ptr<initial_Pk> init_Pk_lin_db = std::make_shared<initial_Pk>(this->arg_cache.get_initial_powerspectrum_path());
        std::unique_ptr<linear_Pk_token> init_Pk_tok = dmgr.tokenize(*model, *init_Pk_lin_db);
        
        // build a work list for filtering the linear power spectrum in wiggle/no-wiggle components
        std::shared_ptr<filterable_Pk> filterable_init_Pk_lin_db = make_filterable(*init_Pk_lin_db);
        auto init_filter_work = dmgr.build_filter_Pk_work_list(*init_Pk_tok, filterable_init_Pk_lin_db, *filter_tok, filter_params);
        
        // distribute this work list among the worker processes
        if(init_filter_work) this->scatter(cosmology_model, *model, *init_filter_work, dmgr);
        
        // exchange our linear power spectrum for the filtered version;
        // we manage its lifetime using std::shared_ptr<> since ownership is shared with the
        // MPI work records and payloads
        std::shared_ptr<initial_filtered_Pk> init_Pk_filt = dmgr.build_wiggle_Pk(*init_Pk_tok, *init_Pk_lin_db);
        
        // if a final linear power spectrum has been specified, then read this in too, tokenize it,
        // and obtain a container
        std::unique_ptr<linear_Pk_token> final_Pk_tok;
        std::shared_ptr<final_filtered_Pk> final_Pk_filt;
        if(this->arg_cache.is_final_powerspectrum_set())
          {
            // read in final power spectrum and tokenize it
            std::shared_ptr<final_Pk> final_Pk_lin_db = std::make_shared<final_Pk>(this->arg_cache.get_final_powerspectrum_path());
            final_Pk_tok = dmgr.tokenize(*model, *final_Pk_lin_db);
            
            // build a work list for filtering
            std::shared_ptr<filterable_Pk> filterable_final_Pk_lin_db = make_filterable(*final_Pk_lin_db);
            auto final_filter_work = dmgr.build_filter_Pk_work_list(*final_Pk_tok, filterable_final_Pk_lin_db, *filter_tok, filter_params);
            
            // distribute this work list among the worker processes
            if(final_filter_work) this->scatter(cosmology_model, *model, *final_filter_work, dmgr);
            
            // get filtered version of power spectrum
            final_Pk_filt = dmgr.build_wiggle_Pk(*final_Pk_tok, *final_Pk_lin_db);
            
            // rescale final power spectrum
            dmgr.rescale_final_Pk(*model, *growth_tok, *final_Pk_filt, *lo_z_db);
          }
        
        
        // STEP 2 - COMPUTE MATSUBARA X & Y COEFFICIENTS
        
        // build a work list of the Matsubara X & Y coefficients associated with these resummation scales
        std::unique_ptr<Matsubara_XY_work_list> Matsubara_work =
          dmgr.build_Matsubara_XY_work_list(*model, *IR_resum_db, init_Pk_filt, *XY_tok, XY_params);
        
        // distribute this work list among the worker processes
        if(Matsubara_work) this->scatter(cosmology_model, *model, *Matsubara_work, dmgr);
        
        
        // STEP 3 - COMPUTE LOOP INTEGRALS
        
        // build a work list for the loop integrals
        std::unique_ptr<loop_integral_work_list> loop_momentum_work =
          dmgr.build_loop_momentum_work_list(*model, *loop_k_db, *IR_cutoff_db, *UV_cutoff_db, init_Pk_filt, *loop_tok, loop_params);
        
        // distribute this work list among the worker processes
        if(loop_momentum_work) this->scatter(cosmology_model, *model, *loop_momentum_work, dmgr);
        
        
        // STEP 4 - COMPUTE ONE-LOOP POWER SPECTRA IN REAL AND REDSHIFT SPACE
        
        // build a work list for the individual power spectrum components
        std::unique_ptr<one_loop_Pk_work_list> Pk_work =
          dmgr.build_one_loop_Pk_work_list(*model, *growth_tok, *loop_tok, *lo_z_db, *loop_k_db,
                                           *IR_cutoff_db, *UV_cutoff_db, init_Pk_filt, final_Pk_filt);
        
        // distribute this work list among the worker processes
        if(Pk_work) this->scatter(cosmology_model, *model, *Pk_work, dmgr);
        
        // build a work list for the resummed power spectrum components
        std::unique_ptr<one_loop_resum_Pk_work_list> Pk_resum_work =
          dmgr.build_one_loop_resum_Pk_work_list(*model, *growth_tok, *loop_tok, *XY_tok, *lo_z_db,
                                                 *loop_k_db, *IR_cutoff_db, *UV_cutoff_db, *IR_resum_db, init_Pk_filt,
                                                 final_Pk_filt);
        
        // distribute this work list among the worker processes
        if(Pk_resum_work) this->scatter(cosmology_model, *model, *Pk_resum_work, dmgr);
        
        
        // STEP 5 - COMPUTE MULTIPOLE DECOMPOSIITON OF REDSHIFT-SPACE POWER SPECTRUM
        
        // build a work list for the resummed multipole power spectra
        std::unique_ptr<multipole_Pk_work_list> multipole_Pk_work =
          dmgr.build_multipole_Pk_work_list(*model, *growth_tok, *loop_tok, *XY_tok, *lo_z_db,
                                            *loop_k_db, *IR_cutoff_db, *UV_cutoff_db, *IR_resum_db, init_Pk_filt,
                                            final_Pk_filt);
        
        // distribute this work list among the worker processes
        if(multipole_Pk_work) this->scatter(cosmology_model, *model, *multipole_Pk_work, dmgr);
      }
    
    // instruct slave processes to terminate
    this->terminate_workers();
  }
