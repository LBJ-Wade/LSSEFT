//
// Created by David Seery on 12/08/2015.
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
