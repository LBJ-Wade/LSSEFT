//
// Created by David Seery on 14/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_ONE_LOOP_PK_CALCULATOR_H
#define LSSEFT_ONE_LOOP_PK_CALCULATOR_H


#include <list>

#include "FRW_model.h"
#include "cosmology/concepts/oneloop_Pk.h"
#include "concepts/oneloop_growth.h"
#include "concepts/loop_integral.h"
#include "cosmology/concepts/power_spectrum.h"

#include "database/tokens.h"

#include "defaults.h"


class oneloop_Pk_calculator
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor is default
    oneloop_Pk_calculator() = default;
    
    //! destructor is default
    ~oneloop_Pk_calculator() = default;
    
    
    // INTERFACE
    
  public:
    
    //! calculate power spectra
    std::list<oneloop_Pk>
    calculate(const Mpc_units::energy& k, const k_token& k_tok, const IR_cutoff_token& IR_tok, const UV_cutoff_token& UV_tok,
              const oneloop_growth& gf_factors, const loop_integral& loop_data, const tree_Pk& Ptree);
    
  private:
    
    //! compute delta-delta power spectrum and counterterm
    dd_Pk compute_dd(const Mpc_units::energy& k, const oneloop_growth_record& val, const loop_integral& loop_data,
                     const tree_Pk& Ptree);
    
    //! compute mu^0 coefficient of rsd delta-delta power spectrum and counterterms
    rsd_dd_Pk compute_rsd_dd_mu0(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                 const loop_integral& loop_data,
                                 const tree_Pk& Ptree);
    
    //! compute mu^2 coefficient of rsd delta-delta power spectrum and counterterms
    rsd_dd_Pk compute_rsd_dd_mu2(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                 const loop_integral& loop_data,
                                 const tree_Pk& Ptree);
    
    //! compute mu^4 coefficient of rsd delta-delta power spectrum and counterterms
    rsd_dd_Pk compute_rsd_dd_mu4(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                 const loop_integral& loop_data,
                                 const tree_Pk& Ptree);
    
    //! compute mu^6 coefficient of rsd delta-delta power spectrum and counterterms
    rsd_dd_Pk compute_rsd_dd_mu6(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                 const loop_integral& loop_data,
                                 const tree_Pk& Ptree);
    
    //! compute mu^8 coefficient of rsd delta-delta power spectrum and counterterms
    rsd_dd_Pk compute_rsd_dd_mu8(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                 const loop_integral& loop_data,
                                 const tree_Pk& Ptree);
    
  };


#endif //LSSEFT_ONE_LOOP_PK_CALCULATOR_H
