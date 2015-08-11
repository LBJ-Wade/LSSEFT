//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_TRANSACTION_MANAGER_H
#define LSSEFT_TRANSACTION_MANAGER_H


#include <functional>


class transaction_manager
  {

    // TYPEDEFS FOR HANDLER CALLBACKS

  public:

    typedef std::function<void()> open_handler;
    typedef std::function<void()> commit_handler;
    typedef std::function<void()> rollback_handler;
    typedef std::function<void()> release_handler;


    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    transaction_manager(open_handler& o, commit_handler& c, rollback_handler& r, release_handler& l);

    //! destructor
    ~transaction_manager();

    // allow moving
    transaction_manager(transaction_manager&& obj) = default;

    // disable copying (transactions are atomic and indivisible)
    transaction_manager(const transaction_manager& obj) = delete;

    // disable assignment (transaction are atomic and indivisible)
    transaction_manager& operator=(const transaction_manager& obj) = delete;


    // TRANSACTION MANAGEMENT

  public:

    //! commit this transaction
    void commit();

    //! rollback this transaction
    void rollback();


    // INTERNAL DATA

  private:

    // CALLBACKS

    //! open callback
    open_handler do_open;

    //! commit callback
    commit_handler do_commit;

    //! rollback callback
    rollback_handler do_rollback;

    //! release callback
    release_handler do_release;


    // STATUS

    //! has this transaction manager been rolled back?
    bool rolled_back;

    //! has this transaction manager been committed?
    bool committed;

  };


#endif //LSSEFT_TRANSACTION_MANAGER_H
