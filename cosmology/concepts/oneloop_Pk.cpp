//
// Created by David Seery on 14/11/2016.
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

#include "oneloop_Pk.h"


oneloop_Pk::oneloop_Pk(const k_token& kt, const growth_params_token& gt, const loop_integral_params_token& lt,
                       const linear_Pk_token& Pkt_i, const boost::optional<linear_Pk_token>& Pkt_f,
                       const IR_cutoff_token& IRt, const UV_cutoff_token& UVt, const z_token& zt,
                       const rsd_dd_Pk& _rsd_mu0, const rsd_dd_Pk& _rsd_mu2, const rsd_dd_Pk& _rsd_mu4,
                       const rsd_dd_Pk& _rsd_mu6, const rsd_dd_Pk& _rsd_mu8)
  : k(kt),
    loop_params(lt),
    growth_params(gt),
    init_Pk(Pkt_i),
    final_Pk(Pkt_f),
    UV_cutoff(UVt),
    IR_cutoff(IRt),
    z(zt),
    rsd_dd_mu0(_rsd_mu0),
    rsd_dd_mu2(_rsd_mu2),
    rsd_dd_mu4(_rsd_mu4),
    rsd_dd_mu6(_rsd_mu6),
    rsd_dd_mu8(_rsd_mu8)
  {
  }


oneloop_Pk::oneloop_Pk()
  : k(0),
    loop_params(0),
    growth_params(0),
    init_Pk(0),
    final_Pk(boost::none),
    UV_cutoff(0),
    IR_cutoff(0),
    z(0),
    rsd_dd_mu0(),
    rsd_dd_mu2(),
    rsd_dd_mu4(),
    rsd_dd_mu6(),
    rsd_dd_mu8()
  {
  }


rsd_dd_Pk::rsd_dd_Pk(const Pk_value& _Pt, const Pk_value& _P13, const Pk_value& _P22)
  : Ptree(_Pt),
    P13(_P13),
    P22(_P22),
    P1loopSPT(_Pt + _P13 + _P22)
  {
  }


rsd_dd_Pk::rsd_dd_Pk()
  : Ptree(),
    P13(),
    P22(),
    P1loopSPT()
  {
  }
