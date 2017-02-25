//
// Created by David Seery on 05/12/2016.
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

#include "filtered_Pk_value.h"


filtered_Pk_value::filtered_Pk_value(const k_token& kt, const linear_Pk_token& Pt, Mpc_units::inverse_energy3 _Pk_nw,
                         Mpc_units::inverse_energy3 _Pk_raw, Mpc_units::inverse_energy3 _Pk_ref)
  : fail(false),
    k_tok(kt),
    Pk_tok(Pt),
    Pk_nw(std::move(_Pk_nw)),
    Pk_raw(std::move(_Pk_raw)),
    Pk_ref(std::move(_Pk_ref))
  {
  }


filtered_Pk_value::filtered_Pk_value()
  : fail(false),
    k_tok(0),
    Pk_tok(0),
    Pk_nw(0.0),
    Pk_raw(0.0),
    Pk_ref(0.0)
  {
  }
