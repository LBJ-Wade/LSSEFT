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

#ifndef LSSEFT_POWER_SPECTRUM_TYPES_H
#define LSSEFT_POWER_SPECTRUM_TYPES_H


#include "generic.h"
#include "units/Mpc_units.h"

#include "boost/filesystem.hpp"


namespace Pk_database_impl
  {
    
    template <int m>
    struct PowerSpectrumTag
      {
        enum { SpectrumClass=m };
      };
    
    typedef PowerSpectrumTag<0> TreeTag;
    typedef PowerSpectrumTag<1> WiggleTag;
    typedef PowerSpectrumTag<2> NoWiggleTag;
    typedef PowerSpectrumTag<3> ApproxTag;
    
  }


// convenience type for tree-level power spectrum
typedef generic_Pk< Pk_database_impl::TreeTag, Mpc_units::inverse_energy3, true > tree_Pk;

// convenience types for wiggle/no-wiggle tree-level power spectra
typedef generic_Pk< Pk_database_impl::WiggleTag, Mpc_units::inverse_energy3, true > tree_Pk_w;
typedef generic_Pk< Pk_database_impl::NoWiggleTag, Mpc_units::inverse_energy3, true > tree_Pk_nw;

// approximation to a power spectrum from a fitting function
typedef generic_Pk< Pk_database_impl::ApproxTag, Mpc_units::inverse_energy3, true > approx_Pk;


#endif //LSSEFT_POWER_SPECTRUM_TYPES_H
