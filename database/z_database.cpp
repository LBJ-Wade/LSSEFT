//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include <assert.h>

#include "z_database.h"


void z_database::add_record(double z, const z_token& tok)
  {
    std::pair<database_type::iterator, bool> emplaced_value = this->database.emplace(z, z_record(z, tok));
    assert(emplaced_value.second);

    this->key_index.emplace(tok.get_id(), emplaced_value.first);
  }


z_database::record_iterator z_database::lookup(z_token tok)
  {
    key_index_type::iterator t = this->key_index.find(tok.get_id());          // find has logarithmic complexity
    database_type::iterator dt = t->second;

    return record_iterator(dt);
  }


z_database::const_record_iterator z_database::lookup(z_token tok) const
  {
    key_index_type::const_iterator t = this->key_index.find(tok.get_id());    // find has logarithmic complexity
    database_type::const_iterator dt = t->second;

    return const_record_iterator(dt);
  }


void z_database::rebuild_key_index()
  {
    this->key_index.clear();

    for(database_type::iterator t = this->database.begin(); t != this->database.end(); ++t)
      {
        this->key_index.emplace(t->second.get_token().get_id(), t);
      }
  }
