//
// Created by David Seery on 18/11/2016.
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

#include "multipole_Pk.h"


Pk_ell::Pk_ell(const Mpc_units::inverse_energy3& _Pt, const Mpc_units::inverse_energy3& _Pt_resum,
               const Mpc_units::inverse_energy3& _P13, const Mpc_units::inverse_energy3& _P13_resum,
               const Mpc_units::inverse_energy3& _P22, const Mpc_units::inverse_energy3& _P22_resum,
               const Mpc_units::inverse_energy3& _PSPT, const Mpc_units::inverse_energy3& _PSPT_resum,
               const Mpc_units::inverse_energy& _Z2d, const Mpc_units::inverse_energy3& _Z0v,
               const Mpc_units::inverse_energy& _Z2v, const Mpc_units::inverse_energy3& _Z0vd,
               const Mpc_units::inverse_energy& _Z2vd, const Mpc_units::inverse_energy& _Z2vv,
               const Mpc_units::inverse_energy& _Z2vvd, const Mpc_units::inverse_energy& _Z2vvv)
  : Ptree(_Pt),
    Ptree_resum(_Pt_resum),
    P13(_P13),
    P13_resum(_P13_resum),
    P22(_P22),
    P22_resum(_P22_resum),
    PSPT(_PSPT),
    PSPT_resum(_PSPT_resum),
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


Pk_ell::Pk_ell()
  : Ptree(0.0),
    Ptree_resum(0.0),
    P13(0.0),
    P13_resum(0.0),
    P22(0.0),
    P22_resum(0.0),
    PSPT(0.0),
    PSPT_resum(0.0),
    Z2_delta(0.0),
    Z0_v(0.0),
    Z2_v(0.0),
    Z0_vdelta(0.0),
    Z2_vdelta(0.0),
    Z2_vv(0.0),
    Z2_vvdelta(0.0),
    Z2_vvv(0.0)
  {
  }


multipole_Pk::multipole_Pk(const k_token& kt, const linear_Pk_token& Pkt, const IR_cutoff_token& IRt,
                           const UV_cutoff_token& UVt, const z_token& zt, const IR_resum_token& IRrt, const Pk_ell& _P0,
                           const Pk_ell& _P2, const Pk_ell& _P4)
  : k(kt),
    Pk(Pkt),
    IR_cutoff(IRt),
    UV_cutoff(UVt),
    z(zt),
    IR_resum(IRrt),
    P0(_P0),
    P2(_P2),
    P4(_P4)
  {
  }


multipole_Pk::multipole_Pk()
  : k(0),
    Pk(0),
    IR_cutoff(0),
    UV_cutoff(0),
    z(0),
    IR_resum(0),
    P0(),
    P2(),
    P4()
  {
  }
