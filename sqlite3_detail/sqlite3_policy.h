//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_SQLITE3_POLICY_H
#define LSSEFT_SQLITE3_POLICY_H


#include <string>


class sqlite3_policy
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    sqlite3_policy();

    //! destructor is default
    ~sqlite3_policy() = default;


    // GET TABLE NAMES

    //! FRW model table
    const std::string& FRW_model_table() const { return(this->FRW_model); }

    //! time configuration table
    const std::string& redshift_config_table() const { return(this->redshift_config); }

    //! wavenumber configuration table
    const std::string& wavenumber_config_table() const { return(this->wavenumber_config); }

    //! transfer function table
    const std::string& transfer_table() const { return(this->transfer); }

    //! 1-loop kernels table
    const std::string& oneloop_table() const { return(this->oneloop); }

    //! temporary table
    const std::string& temp_table() const { return(this->temp); }


    // INTERNAL DATA

  private:

    // TABLE NAMES

    //! FRW model table
    std::string FRW_model;

    //! redshift configuration table
    std::string redshift_config;

    //! wavenumber configuration table
    std::string wavenumber_config;

    //! transfer function table
    std::string transfer;

    //! 1-loop kernels table
    std::string oneloop;

    //! temporary table name
    std::string temp;

  };


#endif //LSSEFT_SQLITE3_POLICY_H
