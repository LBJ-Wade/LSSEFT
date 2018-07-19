//
// Created by David Seery on 19/07/2018.
// --@@
// Copyright (c) 2018 University of Sussex. All rights reserved.
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

#ifndef LSSEFT_GENERIC_POWER_SPECTRUM_H
#define LSSEFT_GENERIC_POWER_SPECTRUM_H


#include "units/Mpc_units.h"


// generic class representing a evaluable power spectrum
template <typename Dimension>
class generic_Pk
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor is default
    generic_Pk() = default;

    //! destructor is default
    ~generic_Pk() = default;


    // INTERFACE

  public:

    //! evaluate spline
    virtual Dimension operator()(const Mpc_units::energy& k) const = 0;

  };


#endif //LSSEFT_GENERIC_POWER_SPECTRUM_H
