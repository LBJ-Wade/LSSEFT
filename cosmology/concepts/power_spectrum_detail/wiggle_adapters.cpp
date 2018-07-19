//
// Created by David Seery on 06/12/2016.
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

#include "wiggle_adapters.h"


wiggle_Pk_raw_adapter::wiggle_Pk_raw_adapter(const initial_filtered_Pk& w, const Mpc_units::energy& klo_, const Mpc_units::energy& khi_)
  : container(w),
    klo(klo_),
    khi(khi_)
  {
  }


Mpc_units::inverse_energy3 wiggle_Pk_raw_adapter::operator()(const Mpc_units::energy& k) const
  {
    if(k < this->klo || k > this->khi) return 0.0;
    return this->container.Pk_raw(k);
  }


wiggle_Pk_wiggle_adapter::wiggle_Pk_wiggle_adapter(const initial_filtered_Pk& w, const Mpc_units::energy& klo_, const Mpc_units::energy& khi_)
  : container(w),
    klo(klo_),
    khi(khi_)
  {
  }


Mpc_units::inverse_energy3 wiggle_Pk_wiggle_adapter::operator()(const Mpc_units::energy& k) const
  {
    if(k < this->klo || k > this->khi) return 0.0;
    return this->container.Pk_wiggle(k);
  }


wiggle_Pk_nowiggle_adapter::wiggle_Pk_nowiggle_adapter(const initial_filtered_Pk& w, const Mpc_units::energy& klo_, const Mpc_units::energy& khi_)
  : container(w),
    klo(klo_),
    khi(khi_)
  {
  }


Mpc_units::inverse_energy3 wiggle_Pk_nowiggle_adapter::operator()(const Mpc_units::energy& k) const
  {
    if(k < this->klo || k > this->khi) return 0.0;
    return this->container.Pk_nowiggle(k);
  }
