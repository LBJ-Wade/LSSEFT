//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "loop_integral.h"


loop_integral::loop_integral(const Mpc_units::energy& _k, const k_token& kt, const Mpc_units::energy& UV,
                             const UV_token& UVt, const Mpc_units::energy& IR, const IR_token& IRt,
                             bool f, const inverse_energy3_kernel& _AA, const inverse_energy3_kernel& _AB, const inverse_energy3_kernel& _BB,
                             const dimless_kernel& _D, const dimless_kernel& _E, const dimless_kernel& _F, const dimless_kernel& _G)
  : k(_k),
    k_tok(kt),
    UV_cutoff(UV),
    UV_tok(UVt),
    IR_cutoff(IR),
    IR_tok(IRt),
    fail(f),
    AA(_AA),
    AB(_AB),
    BB(_BB),
    D(_D),
    E(_E),
    F(_F),
    G(_G)
  {
  }


void loop_integral::set_integration_metadata(boost::timer::nanosecond_type t)
  {
    this->integration_time = t;
  }
