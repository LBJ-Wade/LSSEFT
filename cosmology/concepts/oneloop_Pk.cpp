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


oneloop_Pk::oneloop_Pk(const k_token& kt, const linear_Pk_token& Pkt_i, const boost::optional<linear_Pk_token>& Pkt_f,
                       const IR_cutoff_token& IRt, const UV_cutoff_token& UVt, const z_token& zt,
                       const dd_Pk& _dd, const rsd_dd_Pk& _rsd_mu0, const rsd_dd_Pk& _rsd_mu2, const rsd_dd_Pk& _rsd_mu4,
                       const rsd_dd_Pk& _rsd_mu6, const rsd_dd_Pk& _rsd_mu8)
  : k(kt),
    init_Pk(Pkt_i),
    final_Pk(Pkt_f),
    UV_cutoff(UVt),
    IR_cutoff(IRt),
    z(zt),
    dd(_dd),
    rsd_dd_mu0(_rsd_mu0),
    rsd_dd_mu2(_rsd_mu2),
    rsd_dd_mu4(_rsd_mu4),
    rsd_dd_mu6(_rsd_mu6),
    rsd_dd_mu8(_rsd_mu8)
  {
  }


oneloop_Pk::oneloop_Pk()
  : k(0),
    init_Pk(0),
    final_Pk(boost::none),
    UV_cutoff(0),
    IR_cutoff(0),
    z(0),
    dd(),
    rsd_dd_mu0(),
    rsd_dd_mu2(),
    rsd_dd_mu4(),
    rsd_dd_mu6(),
    rsd_dd_mu8()
  {
  }


rsd_dd_Pk::rsd_dd_Pk(const Pk_value& _Pt, const Pk_value& _P13, const Pk_value& _P22,
                     const k2_Pk_value& _Z2d, const Pk_value& _Z0v, const k2_Pk_value& _Z2v,
                     const Pk_value& _Z0vd, const k2_Pk_value& _Z2vd,
                     const k2_Pk_value& _Z2vv, const k2_Pk_value& _Z2vvd, const k2_Pk_value& _Z2vvv)
  : Ptree(_Pt),
    P13(_P13),
    P22(_P22),
    P1loopSPT(_Pt + _P13 + _P22),
    Z2_delta(_Z2d),
    Z0_v(_Z0v),
    Z2_v(_Z2v),
    Z0_vdelta(_Z0vd),
    Z2_vdelta(_Z2vd),
    Z2_vv(_Z2vv),
    Z2_vvdelta(_Z2vvd),
    Z2_vvv(_Z2vvv)
  {
  }


rsd_dd_Pk::rsd_dd_Pk()
  : Ptree(),
    P13(),
    P22(),
    P1loopSPT(),
    Z2_delta(),
    Z0_v(),
    Z2_v(),
    Z0_vdelta(),
    Z2_vdelta(),
    Z2_vv(),
    Z2_vvdelta(),
    Z2_vvv()
  {
  }
