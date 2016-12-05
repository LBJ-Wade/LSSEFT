//
// Created by David Seery on 05/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
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
    
  }


// convenience type for tree-level power spectrum
typedef generic_Pk< Pk_database_impl::TreeTag, Mpc_units::inverse_energy3, true > tree_Pk;

// convenience types for wiggle/no-wiggle tree-level power spectra
typedef generic_Pk< Pk_database_impl::WiggleTag, Mpc_units::inverse_energy3, true > wiggle_Pk;
typedef generic_Pk< Pk_database_impl::NoWiggleTag, Mpc_units::inverse_energy3, true > nowiggle_Pk;


#endif //LSSEFT_POWER_SPECTRUM_TYPES_H
