//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include <assert.h>

#include "redshift_database.h"


redshift_record::redshift_record(double _z, const redshift_token& tok)
  : z(_z),
    token(tok)
  {
  }


void redshift_database::add_record(double z, const redshift_token& tok)
  {
    std::pair<database_type::iterator, bool> emplaced_value = this->database.emplace(z, redshift_record(z, tok));
    assert(emplaced_value.second);

    this->key_index.emplace(tok.get_id(), emplaced_value.first);
  }


redshift_database::record_iterator redshift_database::lookup(redshift_token tok)
  {
    key_index_type::iterator t = this->key_index.find(tok.get_id());          // find has logarithmic complexity
    database_type::iterator dt = t->second;

    return record_iterator(dt);
  }


redshift_database::const_record_iterator redshift_database::lookup(redshift_token tok) const
  {
    key_index_type::const_iterator t = this->key_index.find(tok.get_id());    // find has logarithmic complexity
    database_type::const_iterator dt = t->second;

    return const_record_iterator(dt);
  }


void redshift_database::rebuild_key_index()
  {
    this->key_index.clear();

    for(database_type::iterator t = this->database.begin(); t != this->database.end(); ++t)
      {
        this->key_index.emplace(t->second.get_token().get_id(), t);
      }
  }
