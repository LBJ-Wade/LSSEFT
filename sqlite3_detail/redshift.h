//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_REDSHIFT_H
#define LSSEFT_SQLITE3_REDSHIFT_H


#include <memory>

#include "database/transaction_manager.h"
#include "sqlite3_policy.h"

#include "boost/optional.hpp"

#include "sqlite3.h"


namespace sqlite3_operations
  {

    //! lookup id for a redshift
    boost::optional<unsigned int> lookup_redshift(sqlite3* db, std::shared_ptr<transaction_manager>& mgr, double z,
                                                  const sqlite3_policy& policy, double tol);

    //! insert a redshift
    unsigned int insert_redshift(sqlite3* db, std::shared_ptr<transaction_manager>& mgr, double z,
                                 const sqlite3_policy& policy);

  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_REDSHIFT_H
