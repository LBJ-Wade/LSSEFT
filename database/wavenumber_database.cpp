//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include <assert.h>

#include "wavenumber_database.h"


wavenumber_record::wavenumber_record(const eV_units::energy& _k, const wavenumber_token& tok)
  : k(_k),
    token(tok)
  {
  }


void wavenumber_database::add_record(const eV_units::energy& k, const wavenumber_token& tok)
  {
    std::pair<database_type::iterator, bool> emplaced_value = this->database.emplace(tok.get_id(), wavenumber_record(k, tok));
    assert(emplaced_value.second);
  }
