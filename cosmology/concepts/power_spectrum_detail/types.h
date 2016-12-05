//
// Created by David Seery on 05/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_POWER_SPECTRUM_TYPES_H
#define LSSEFT_POWER_SPECTRUM_TYPES_H


#include "generic.h"
#include "units/Mpc_units.h"

#include "boost/filesystem.hpp"


namespace power_spectrum_database_impl
  {
    
    template <int m>
    struct PowerSpectrumTag
      {
        enum { SpectrumClass=m };
      };
    
    typedef PowerSpectrumTag<0> TreeTag;
    
  }


// convenience type for tree-level power spectrum
typedef generic_power_spectrum< power_spectrum_database_impl::TreeTag, Mpc_units::inverse_energy3, true > tree_power_spectrum;


#endif //LSSEFT_POWER_SPECTRUM_TYPES_H
