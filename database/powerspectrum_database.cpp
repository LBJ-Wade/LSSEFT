//
// Created by David Seery on 04/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include <limits>

#include "powerspectrum_database.h"


powerspectrum_database::powerspectrum_database()
  : k_min(std::numeric_limits<double>::max()),
    k_max(std::numeric_limits<double>::min())
  {
  }


void powerspectrum_database::add_record(const eV_units::energy& k, const eV_units::inverse_energy3& Pk)
  {
    std::pair<database_type::iterator, bool> emplaced_value = this->database.emplace(k, Pk_record(k, Pk));
    assert(emplaced_value.second);

    if(k > this->k_max) this->k_max = k;
    if(k < this->k_min) this->k_min = k;
  }
