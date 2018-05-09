//
// Created by David Seery on 21/11/2015.
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

#include "loop_integral.h"


loop_integral::loop_integral(const k_token& kt, const loop_integral_params_token& pt, const linear_Pk_token& Pkt,
                             const UV_cutoff_token& UVt, const IR_cutoff_token& IRt, const kernels& ker_)
  : k(kt),
    params(pt),
    Pk_lin(Pkt),
    UV_cutoff(UVt),
    IR_cutoff(IRt),
    ker(ker_)
  {
  }


loop_integral::loop_integral()
  : k(0),
    params(0),
    Pk_lin(0),
    UV_cutoff(0),
    IR_cutoff(0),
    ker()
  {
  }

