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

#include "wiggle.h"


wiggle_Pk::wiggle_Pk(const linear_Pk_token& t, const tree_Pk_w::database_type& w, const tree_Pk::database_type& r)
  : tok(t),
    nowiggle(w),
    raw(r)
  {
  }


bool wiggle_Pk::is_valid(const Mpc_units::energy& k) const
  {
    return this->nowiggle.is_valid(k, 0, 0) && this->raw.is_valid(k, 0, 0);
  }
