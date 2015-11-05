//
// Created by David Seery on 04/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "powerspectrum_database.h"


void powerspectrum_database::add_record(const eV_units::energy& k, double Pk)
  {
    std::pair<database_type::iterator, bool> emplaced_value = this->database.emplace(k, Pk_record(k, Pk));
    assert(emplaced_value.second);
  }
