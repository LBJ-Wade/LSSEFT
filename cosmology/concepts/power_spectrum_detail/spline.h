//
// Created by David Seery on 06/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
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
