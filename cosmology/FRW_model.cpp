//
// Created by David Seery on 11/08/2015.
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

#include <iostream>

#include "FRW_model.h"


FRW_model::FRW_model(std::string nm, double om, double occ, double h_, Mpc_units::energy tc, double ne, double fb, double zs,
                     double zd, double ze, double Ac, double n, Mpc_units::energy kp)
  : name(nm),
    omega_m(om),
    omega_cc(occ),
    h(h_),
    T_CMB(tc),
    Neff(ne),
    f_baryon(fb),
    z_star(zs),
    z_drag(zd),
    z_eq(ze),
    A_curv(Ac),
    ns(n),
    k_piv(kp)
  {
  }
