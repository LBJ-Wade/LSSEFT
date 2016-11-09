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

    //! IR configuration table
    const std::string& IR_config_table() const { return(this->IR_config); }

    //! UV configuration table
    const std::string& UV_config_table() const { return(this->UV_config); }

    //! transfer function table
    const std::string& transfer_table() const { return(this->transfer); }

    //! 1-loop growth factors table
    const std::string& oneloop_table() const { return(this->oneloop); }

    //! 1-loop AA integral table
    const std::string& AA_table() const { return(this->AA); }
    
    //! 1-loop AA integral table
    const std::string& AB_table() const { return(this->AB); }
    
    //! 1-loop AA integral table
    const std::string& BB_table() const { return(this->BB); }
    
    //! 1-loop AA integral table
    const std::string& D_table() const { return(this->D); }
    
    //! 1-loop AA integral table
    const std::string& E_table() const { return(this->E); }
    
    //! 1-loop AA integral table
    const std::string& F_table() const { return(this->F); }
    
    //! 1-loop AA integral table
    const std::string& G_table() const { return(this->G); }
    
    //! 1-loop AA integral table
    const std::string& J1_table() const { return(this->J1); }
    
    //! 1-loop AA integral table
    const std::string& J2_table() const { return(this->J2); }

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

    //! IR cutoff configuration table
    std::string IR_config;

    //! UV cutoff configuration table
    std::string UV_config;

    //! transfer function table
    std::string transfer;

    //! 1-loop growth factors table
    std::string oneloop;

    //! 1-loop momentum integrals - AA
    std::string AA;
    
    //! 1-loop momentum integrals - AB
    std::string AB;
    
    //! 1-loop momentum integrals - BB
    std::string BB;
    
    //! 1-loop momentum integrals - D
    std::string D;
    
    //! 1-loop momentum integrals - E
    std::string E;
    
    //! 1-loop momentum integrals - F
    std::string F;
    
    //! 1-loop momentum integrals - G
    std::string G;
    
    //! 1-loop momentum integrals - J1
    std::string J1;
    
    //! 1-loop momentum integrals - J2
    std::string J2;

    //! temporary table name
    std::string temp;

  };


#endif //LSSEFT_SQLITE3_POLICY_H
