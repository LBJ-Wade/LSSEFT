//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_WAVENUMBER_H
#define LSSEFT_SQLITE3_WAVENUMBER_H


#include <memory>

#include "database/transaction_manager.h"
#include "sqlite3_policy.h"

#include "units/eV_units.h"

#include "boost/optional.hpp"

#include "sqlite3.h"


namespace sqlite3_operations
  {

    //! lookup id for a wavenumber
    boost::optional<unsigned int> lookup_wavenumber(sqlite3* db, transaction_manager& mgr,
                                                    const eV_units::energy& k, const sqlite3_policy& policy,
                                                    double tol);

    //! insert a wavenumber
    unsigned int insert_wavenumber(sqlite3* db, transaction_manager& mgr,
                                   const eV_units::energy& k, const sqlite3_policy& policy);

  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_WAVENUMBER_H
