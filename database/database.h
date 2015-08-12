//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_DATABASE_H
#define LSSEFT_DATABASE_H


#include <memory>
#include <sqlite3_detail/sqlite3_policy.h>

#include "tokens.h"
#include "transaction_manager.h"

#include "cosmology/FRW_model.h"

#include "boost/filesystem/operations.hpp"

#include "sqlite3.h"


class database
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor opens SQLite3 handle, and creates tables if this is a new database
    database(const boost::filesystem::path& c);

    //! destructor closes SQLite3 handle
    ~database();


    // TOKENS
    // tokens are the basic unit of currency used when interacting with the database

  public:

    FRW_model_token tokenize_FRW_model(const FRW_model& obj);


    // TRANSACTIONS

  protected:

    //! open a transaction; throws an exception if a transaction is already held open
    std::shared_ptr<transaction_manager> open_transaction();

    //! begin a new transaction on the database
    void begin_transaction();

    //! commit a transaction on the database
    void commit_transaction();

    //! rollback a transaction on the database
    void rollback_transaction();

    //! release a transaction
    void release_transaction();


    // LOOKUP OR INSERT RECORDS

  protected:

    //! lookup or insert a new FRW model
    unsigned int lookup_or_insert_model(std::shared_ptr<transaction_manager>& mgr, const FRW_model& obj);


    // INTERNAL DATA

  private:

    //! path to data container
    boost::filesystem::path container;

    //! SQLite3 handler for container
    sqlite3* handle;


    // TRANSACTIONS

    //! current transaction manager, if one exists
    std::weak_ptr<transaction_manager> current_transaction;


    // SQLite3 policies

    //! sqlite3_policy object
    sqlite3_policy policy;


    // SEARCH TOLERANCES

    //! tolerance to use when searching for FRW model parameters
    double FRW_model_tol;

    //! tolerance to use when searching for time configurations
    double t_tol;

    //! tolerance to use when searching for wavenumber configurations
    double k_tol;

  };


#endif //LSSEFT_DATABASE_H
