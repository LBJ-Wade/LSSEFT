//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "oneloop_momentum_integrator.h"


loop_integral oneloop_momentum_integrator::integrate(const FRW_model& model, const eV_units::energy& k,
                                                     const k_token& k_tok, const eV_units::energy& UV_cutoff,
                                                     const UV_token& UV_tok, const eV_units::energy& IR_cutoff,
                                                     const IR_token& IR_tok, std::shared_ptr<tree_power_spectrum>& Pk)
  {
    return loop_integral(k, k_tok, UV_cutoff, UV_tok, IR_cutoff, IR_tok, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  }
