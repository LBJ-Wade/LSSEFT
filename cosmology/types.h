//
// Created by David Seery on 13/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_TYPES_H
#define LSSEFT_TYPES_H


#include <map>

#include "database/tokens.h"
#include "database/redshift_database.h"


//! work record for a transfer function calculation
class transfer_work_record
  {

  public:

    // CONSTRUCTOR, DESTRUCTOR

    //! constructor
    transfer_work_record(const std::shared_ptr<wavenumber_token>& k, const std::shared_ptr<redshift_database>& z)
      : k_token(k),
        z_db(z)
      {
      }

    //! destructor is defaults
    ~transfer_work_record() = default;


    // INTERFACE

  public:

    //! dereference to get k_token
    const std::shared_ptr<wavenumber_token>& operator*() { return(this->k_token); }

    //! get redshift database
    const std::shared_ptr<redshift_database>& get_z_db() { return(this->z_db); }


    // INTERNAL DATA

  private:

    //! wavenumber token
    std::shared_ptr<wavenumber_token> k_token;

    //! redshift database
    std::shared_ptr<redshift_database> z_db;

  };

//! list of work for transfer function calculation
typedef std::map< unsigned int, transfer_work_record> transfer_work_list;


#endif //LSSEFT_TYPES_H
