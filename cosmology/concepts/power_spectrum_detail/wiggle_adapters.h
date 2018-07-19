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

#ifndef LSSEFT_POWER_SPECTRUM_SPLINE_H
#define LSSEFT_POWER_SPECTRUM_SPLINE_H


#include "units/Mpc_units.h"

#include "generic.h"
#include "wiggle.h"


// define adapter classes for wiggle power spectrum
class wiggle_Pk_raw_adapter: public generic_Pk<Mpc_units::inverse_energy3>
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor captures wiggle_Pk container
    wiggle_Pk_raw_adapter(const initial_filtered_Pk& w, const Mpc_units::energy& klo_, const Mpc_units::energy& khi_);
    
    //! destructor is default
    ~wiggle_Pk_raw_adapter() = default;
    
    
    // INTERFACE
    
  public:
    
    //! evaluate spline
    Mpc_units::inverse_energy3 operator()(const Mpc_units::energy& k) const override final;
    
    
    // INTERNAL DATA
    
  private:
    
    //! capture wiggle_Pk container
    const initial_filtered_Pk& container;

    //! UV cutoff
    Mpc_units::energy khi;

    //! IR cutoff
    Mpc_units::energy klo;
    
  };


class wiggle_Pk_wiggle_adapter: public generic_Pk<Mpc_units::inverse_energy3>
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! constructor captures wiggle_Pk container
    wiggle_Pk_wiggle_adapter(const initial_filtered_Pk& w, const Mpc_units::energy& klo_, const Mpc_units::energy& khi_);
    
    //! destructor is default
    ~wiggle_Pk_wiggle_adapter() = default;
    
    
    // INTERFACE
  
  public:
    
    //! evaluate spline
    Mpc_units::inverse_energy3 operator()(const Mpc_units::energy& k) const override final;
    
    
    // INTERNAL DATA
  
  private:
    
    //! capture wiggle_Pk container
    const initial_filtered_Pk& container;

    //! UV cutoff
    Mpc_units::energy khi;

    //! IR cutoff
    Mpc_units::energy klo;

  };


class wiggle_Pk_nowiggle_adapter: public generic_Pk<Mpc_units::inverse_energy3>
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! constructor captures wiggle_Pk container
    wiggle_Pk_nowiggle_adapter(const initial_filtered_Pk& w, const Mpc_units::energy& klo_, const Mpc_units::energy& khi_);
    
    //! destructor is default
    ~wiggle_Pk_nowiggle_adapter() = default;
    
    
    // INTERFACE
  
  public:
    
    //! evaluate spline
    Mpc_units::inverse_energy3 operator()(const Mpc_units::energy& k) const override final;
    
    
    // INTERNAL DATA
  
  private:
    
    //! capture wiggle_Pk container
    const initial_filtered_Pk& container;

    //! UV cutoff
    Mpc_units::energy khi;

    //! IR cutoff
    Mpc_units::energy klo;

  };


#endif //LSSEFT_POWER_SPECTRUM_SPLINE_H
