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

    //! delta_m transfer function table
    const std::string& delta_m_table() const { return(this->delta_m); }

    //! delta_r transfer function table
    const std::string& delta_r_table() const { return(this->delta_r); }

    //! theta_m transfer function table
    const std::string& theta_m_table() const { return(this->theta_m); }

    //! theta_r transfer function table
    const std::string& theta_r_table() const { return(this->theta_r); }

    //! Phi transfer function table
    const std::string Phi_table() const { return(this->Phi); }

    //! 1-loop kernels table
    const std::string& oneloop_table() const { return(this->oneloop); }


    // INTERNAL DATA

  private:

    // TABLE NAMES

    //! FRW model table
    std::string FRW_model;

    //! redshift configuration table
    std::string redshift_config;

    //! wavenumber configuration table
    std::string wavenumber_config;

    //! delta_m transfer function table
    std::string delta_m;

    //! delta_r transfer function table
    std::string delta_r;

    //! theta_m transfer function table
    std::string theta_m;

    //! theta_r transfer function table
    std::string theta_r;

    //! Phi transfer function table
    std::string Phi;

    //! 1-loop kernels table
    std::string oneloop;

  };


#endif //LSSEFT_SQLITE3_POLICY_H
