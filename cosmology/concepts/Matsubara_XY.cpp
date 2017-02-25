//
// Created by David Seery on 21/11/2016.
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

#include "Matsubara_XY.h"


Matsubara_XY::Matsubara_XY(const linear_Pk_token& Pkt, const IR_resum_token& IRt, const Mpc_units::inverse_energy2& _X,
                           const Mpc_units::inverse_energy2& _Y)
  : IR_resum_tok(IRt),
    Pk_lin(Pkt),
    X(_X),
    Y(_Y)
  {
  }


Matsubara_XY::Matsubara_XY()
  : IR_resum_tok(0),
    Pk_lin(0),
    X(0.0),
    Y(0.0)
  {
  }
