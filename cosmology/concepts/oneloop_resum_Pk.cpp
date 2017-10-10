//
// Created by David Seery on 09/12/2016.
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

#include "oneloop_resum_Pk.h"


oneloop_resum_Pk::oneloop_resum_Pk(const k_token& kt, const growth_params_token& gp,
                                   const loop_integral_params_token& lp, const MatsubaraXY_params_token& XYp,
                                   const linear_Pk_token& Pkt_i, const boost::optional<linear_Pk_token>& Pkt_f,
                                   const IR_cutoff_token& IRt, const UV_cutoff_token& UVt, const z_token& zt,
                                   const IR_resum_token& IRrt, const resum_dd_Pk& Pkr)
  : k(kt),
    growth_params(gp),
    loop_params(lp),
    XY_params(XYp),
    init_Pk(Pkt_i),
    final_Pk(Pkt_f),
    IR_cutoff(IRt),
    UV_cutoff(UVt),
    z(zt),
    IR_resum(IRrt),
    Pk_resum(Pkr)
  {
  }


oneloop_resum_Pk::oneloop_resum_Pk()
  : k(0.0),
    growth_params(0),
    loop_params(0),
    XY_params(0),
    init_Pk(0),
    final_Pk(boost::none),
    IR_cutoff(0),
    UV_cutoff(0),
    z(0),
    IR_resum(0),
    Pk_resum()
  {
  }