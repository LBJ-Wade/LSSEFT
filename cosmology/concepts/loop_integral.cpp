//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "loop_integral.h"


loop_integral::loop_integral(const eV_units::energy& _k, const k_token& kt, const eV_units::energy& UV,
                             const UV_token& UVt, const eV_units::energy& IR, const IR_token& IRt,
                             double _A, double _B, double _D, double _E, double _F, double _G)
  : k(_k),
    k_tok(kt),
    UV_cutoff(UV),
    UV_tok(UVt),
    IR_cutoff(IR),
    IR_tok(IRt),
    A(_A),
    B(_B),
    D(_D),
    E(_E),
    F(_F),
    G(_G)
  {
  }


void loop_integral::set_integration_metadata(boost::timer::nanosecond_type t, size_t s)
  {
    this->integration_time = t;
    this->steps = s;
  }
