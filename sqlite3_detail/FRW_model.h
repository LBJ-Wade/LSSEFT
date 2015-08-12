//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_FRW_MODEL_H
#define LSSEFT_SQLITE3_FRW_MODEL_H


#include <memory>

#include "cosmology/FRW_model.h"

#include "database/transaction_manager.h"
#include "sqlite3_policy.h"

#include "boost/optional.hpp"

#include "sqlite3.h"


namespace sqlite3_operations
  {

    //! lookup id for an FRW model
    boost::optional<unsigned int> lookup_FRW_model(sqlite3* db, std::shared_ptr<transaction_manager>& mgr, const FRW_model& obj, const sqlite3_policy& policy, double tol);

    //! insert an FRW model
    unsigned int insert_FRW_model(sqlite3* db, std::shared_ptr<transaction_manager>& mgr, const FRW_model& obj, const sqlite3_policy& policy);

  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_FRW_MODEL_H
