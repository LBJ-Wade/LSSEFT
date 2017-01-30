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

#include "wiggle.h"


// generic class representing a splined power spectrum
class spline_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor is default
    spline_Pk() = default;
    
    //! destructor is default
    ~spline_Pk() = default;
    
    
    // INTERFACE
    
  public:
    
    //! evaluate spline
    virtual Mpc_units::inverse_energy3 operator()(const Mpc_units::energy& k) const = 0;
    
  };


// define adapter classes for wiggle power spectrum
class wiggle_Pk_raw_adapter: public spline_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor captures wiggle_Pk container
    wiggle_Pk_raw_adapter(const wiggle_Pk& w);
    
    //! destructor is default
    ~wiggle_Pk_raw_adapter() = default;
    
    
    // INTERFACE
    
  public:
    
    //! evaluate spline
    Mpc_units::inverse_energy3 operator()(const Mpc_units::energy& k) const override final;
    
    
    // INTERNAL DATA
    
  private:
    
    //! capture wiggle_Pk container
    const wiggle_Pk& container;
    
  };


class wiggle_Pk_wiggle_adapter: public spline_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! constructor captures wiggle_Pk container
    wiggle_Pk_wiggle_adapter(const wiggle_Pk& w);
    
    //! destructor is default
    ~wiggle_Pk_wiggle_adapter() = default;
    
    
    // INTERFACE
  
  public:
    
    //! evaluate spline
    Mpc_units::inverse_energy3 operator()(const Mpc_units::energy& k) const override final;
    
    
    // INTERNAL DATA
  
  private:
    
    //! capture wiggle_Pk container
    const wiggle_Pk& container;
    
  };


class wiggle_Pk_nowiggle_adapter: public spline_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! constructor captures wiggle_Pk container
    wiggle_Pk_nowiggle_adapter(const wiggle_Pk& w);
    
    //! destructor is default
    ~wiggle_Pk_nowiggle_adapter() = default;
    
    
    // INTERFACE
  
  public:
    
    //! evaluate spline
    Mpc_units::inverse_energy3 operator()(const Mpc_units::energy& k) const override final;
    
    
    // INTERNAL DATA
  
  private:
    
    //! capture wiggle_Pk container
    const wiggle_Pk& container;
    
  };


#endif //LSSEFT_POWER_SPECTRUM_SPLINE_H
