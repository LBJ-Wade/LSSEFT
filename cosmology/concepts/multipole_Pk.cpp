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


Pk_ell::Pk_ell(const Pk_resum& _tree, const Pk_resum& _P13, const Pk_resum& _P22, const Pk_resum& _PSPT)
  : Ptree(std::move(_tree)),
    P13(std::move(_P13)),
    P22(std::move(_P22)),
    PSPT(std::move(_PSPT))
  {
  }


multipole_Pk::multipole_Pk(const k_token kt, const growth_params_token& gp, const loop_integral_params_token& lp,
                           const MatsubaraXY_params_token& XYp, const linear_Pk_token Pkt_i,
                           const boost::optional<linear_Pk_token> Pkt_f, const IR_cutoff_token IRt, const UV_cutoff_token UVt,
                           const z_token zt, const IR_resum_token IRrt, const Pk_ell _P0, const Pk_ell _P2, const Pk_ell _P4)
  : k(std::move(kt)),
    growth_params(gp),
    loop_params(lp),
    XY_params(XYp),
    init_Pk(std::move(Pkt_i)),
    final_Pk(std::move(Pkt_f)),
    IR_cutoff(std::move(IRt)),
    UV_cutoff(std::move(UVt)),
    z(std::move(zt)),
    IR_resum(std::move(IRrt)),
    P0(std::move(_P0)),
    P2(std::move(_P2)),
    P4(std::move(_P4))
  {
  }


multipole_Pk::multipole_Pk()
  : k(0),
    growth_params(0),
    loop_params(0),
    XY_params(0),
    init_Pk(0),
    final_Pk(boost::none),
    IR_cutoff(0),
    UV_cutoff(0),
    z(0),
    IR_resum(0),
    P0(Pk_resum(), Pk_resum(), Pk_resum(), Pk_resum()),
    P2(Pk_resum(), Pk_resum(), Pk_resum(), Pk_resum()),
    P4(Pk_resum(), Pk_resum(), Pk_resum(), Pk_resum())
  {
  }
