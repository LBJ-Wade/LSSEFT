//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "loop_integral.h"


loop_integral::loop_integral(const eV_units::energy& _k, const k_token& kt, const eV_units::energy& UV,
                             const UV_token& UVt, const eV_units::energy& IR, const IR_token& IRt, double v)
  : k(_k),
    k_tok(kt),
    UV_cutoff(UV),
    UV_tok(UVt),
    IR_cutoff(IR),
    IR_tok(IRt),
    value(v)
  {
  }


void loop_integral::set_integration_metadata(boost::timer::nanosecond_type t, size_t s)
  {
    this->integration_time = t;
    this->steps = s;
  }
