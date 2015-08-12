//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include <assert.h>

#include "redshift_database.h"


redshift_record::redshift_record(double _z, std::shared_ptr<redshift_token> tok)
  : z(_z),
    token(tok)
  {
  }


void redshift_database::add_record(double z, std::shared_ptr<redshift_token> tok)
  {
    std::pair<database_type::iterator, bool> emplaced_value = this->database.emplace(tok->get_id(), redshift_record(z, tok));
    assert(emplaced_value.second);
  }
