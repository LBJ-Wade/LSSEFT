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
    std::pair<database_type::iterator, bool> emplaced_value = this->database.emplace(static_cast<double>(k), wavenumber_record(k, tok));
    assert(emplaced_value.second);

    this->key_index.emplace(tok.get_id(), emplaced_value.first);
  }


wavenumber_database::record_iterator wavenumber_database::lookup(wavenumber_token tok)
  {
    key_index_type::iterator t = this->key_index.find(tok.get_id());
    database_type::iterator dt = t->second;

    return record_iterator(dt);
  }


wavenumber_database::const_record_iterator wavenumber_database::lookup(wavenumber_token tok) const
  {
    key_index_type::const_iterator t = this->key_index.find(tok.get_id());
    database_type::const_iterator dt = t->second;

    return const_record_iterator(dt);
  }


void wavenumber_database::rebuild_key_index()
  {
    this->key_index.clear();

    for(database_type::iterator t = this->database.begin(); t != this->database.end(); ++t)
      {
        this->key_index.emplace(t->second.get_token().get_id(), t);
      }
  }
