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
    const std::string& growth_factor_table() const { return(this->growth_factor); }
    
    //! 1-loop growth rate table
    const std::string& growth_rate_table() const { return(this->growth_rate); }

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
    const std::string FRW_model;

    //! redshift configuration table
    const std::string redshift_config;

    //! wavenumber configuration table
    const std::string wavenumber_config;

    //! IR cutoff configuration table
    const std::string IR_config;

    //! UV cutoff configuration table
    const std::string UV_config;

    //! transfer function table
    const std::string transfer;

    //! 1-loop growth factors table
    const std::string growth_factor;
    
    //! 1-loop growth rate table
    const std::string growth_rate;

    //! 1-loop momentum integrals - AA
    const std::string AA;
    
    //! 1-loop momentum integrals - AB
    const std::string AB;
    
    //! 1-loop momentum integrals - BB
    const std::string BB;
    
    //! 1-loop momentum integrals - D
    const std::string D;
    
    //! 1-loop momentum integrals - E
    const std::string E;
    
    //! 1-loop momentum integrals - F
    const std::string F;
    
    //! 1-loop momentum integrals - G
    const std::string G;
    
    //! 1-loop momentum integrals - J1
    const std::string J1;
    
    //! 1-loop momentum integrals - J2
    const std::string J2;

    //! temporary table name
    const std::string temp;

  };


#endif //LSSEFT_SQLITE3_POLICY_H
