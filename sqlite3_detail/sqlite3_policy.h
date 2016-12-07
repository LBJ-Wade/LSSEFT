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

    //! IR-cutoff configuration table
    const std::string& IR_config_table() const { return(this->IR_config); }

    //! UV-cutoff configuration table
    const std::string& UV_config_table() const { return(this->UV_config); }
    
    //! IR-resummation configuration table
    const std::string& IR_resum_config_table() const { return(this->IR_resum_config); }
    
    //! linear power spectrum configuration table
    const std::string& Pk_linear_config_table() const { return(this->Pk_linear_config); }
    
    //! linear power spectrum table
    const std::string& Pk_linear_table() const { return(this->Pk_linear); }
    
    //! transfer function table
    const std::string& transfer_table() const { return(this->transfer); }

    //! 1-loop growth g-factor table
    const std::string& g_factor_table() const { return(this->growth_g_factor); }
    
    //! 1-loop growth f-factor table
    const std::string& f_factor_table() const { return(this->growth_f_factor); }

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
    
    //! 1-loop RSD13_a integral table
    const std::string& RSD13_a_table() const { return this->RSD13_a; }
    
    //! 1-loop RSD13_b integral table
    const std::string& RSD13_b_table() const { return this->RSD13_b; }
    
    //! 1-loop RSD13_c integral table
    const std::string& RSD13_c_table() const { return this->RSD13_c; }
    
    //! 1-loop RSD13_d integral table
    const std::string& RSD13_d_table() const { return this->RSD13_d; }
    
    //! 1-loop RSD13_e integral table
    const std::string& RSD13_e_table() const { return this->RSD13_e; }
    
    //! 1-loop RSD13_f integral table
    const std::string& RSD13_f_table() const { return this->RSD13_f; }
    
    //! 1-loop RSD13_g integral table
    const std::string& RSD13_g_table() const { return this->RSD13_g; }
    
    //! 1-loop RSD22_A1 integral table
    const std::string& RSD22_A1_table() const { return this->RSD22_A1; }
    
    //! 1-loop RSD22_A2 integral table
    const std::string& RSD22_A2_table() const { return this->RSD22_A2; }
    
    //! 1-loop RSD22_A3 integral table
    const std::string& RSD22_A3_table() const { return this->RSD22_A3; }
    
    //! 1-loop RSD22_A3 integral table
    const std::string& RSD22_A4_table() const { return this->RSD22_A4; }
    
    //! 1-loop RSD22_A5 integral table
    const std::string& RSD22_A5_table() const { return this->RSD22_A5; }
    
    //! 1-loop RSD22_B2 integral table
    const std::string& RSD22_B2_table() const { return this->RSD22_B2; }
    
    //! 1-loop RSD22_B3 integral table
    const std::string& RSD22_B3_table() const { return this->RSD22_B3; }
    
    //! 1-loop RSD22_B6 integral table
    const std::string& RSD22_B6_table() const { return this->RSD22_B6; }
    
    //! 1-loop RSD22_B8 integral table
    const std::string& RSD22_B8_table() const { return this->RSD22_B8; }
    
    //! 1-loop RSD22_B9 integral table
    const std::string& RSD22_B9_table() const { return this->RSD22_B9; }
    
    //! 1-loop RSD22_C1 integral table
    const std::string& RSD22_C1_table() const { return this->RSD22_C1; }
    
    //! 1-loop RSD22_C2 integral table
    const std::string& RSD22_C2_table() const { return this->RSD22_C2; }
    
    //! 1-loop RSD22_C4 integral table
    const std::string& RSD22_C4_table() const { return this->RSD22_C4; }
    
    //! 1-loop RSD22_D1 integral table
    const std::string& RSD22_D1_table() const { return this->RSD22_D1; }
    
    //! 1-loop delta-delta power spectrum table
    const std::string& dd_Pk_table() const { return(this->dd_Pk); }
    
    //! 1-loop delta-delta rsd power spectrum mu^0 table
    const std::string& dd_rsd_mu0_Pk_table() const { return(this->dd_rsd_mu0_Pk); }
    
    //! 1-loop delta-delta rsd power spectrum mu^2 table
    const std::string& dd_rsd_mu2_Pk_table() const { return(this->dd_rsd_mu2_Pk); }
    
    //! 1-loop delta-delta rsd power spectrum mu^4 table
    const std::string& dd_rsd_mu4_Pk_table() const { return(this->dd_rsd_mu4_Pk); }
    
    //! 1-loop delta-delta rsd power spectrum mu^6 table
    const std::string& dd_rsd_mu6_Pk_table() const { return(this->dd_rsd_mu6_Pk); }
    
    //! 1-loop delta-delta rsd power spectrum mu^8 table
    const std::string& dd_rsd_mu8_Pk_table() const { return(this->dd_rsd_mu8_Pk); }
    
    //! 1-loop multipole P0 table
    const std::string& P0_table() const { return this->P0; }
    
    //! 1-loop multipole P2 table
    const std::string& P2_table() const { return this->P2; }
    
    //! 1-loop multipole P4 table
    const std::string& P4_table() const { return this->P4; }
    
    //! Matsubara-XY table
    const std::string& Matsubara_XY_table() const { return this->Matsubara_XY; }

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
    
    //! IR-resummation configuration table
    const std::string IR_resum_config;
    
    //! linear power spectrum P(k) configuration table
    const std::string Pk_linear_config;
    
    //! linear power spectrum data table
    const std::string Pk_linear;

    //! transfer function table
    const std::string transfer;

    //! 1-loop growth g-factor table
    const std::string growth_g_factor;
    
    //! 1-loop growth f-factor table
    const std::string growth_f_factor;
    
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
    
    //! 1-loop momentum integrals - RSD13_a
    const std::string RSD13_a;
    
    //! 1-loop momentum integrals - RSD13_b
    const std::string RSD13_b;
    
    //! 1-loop momentum integrals - RSD13_c
    const std::string RSD13_c;
    
    //! 1-loop momentum integrals - RSD13_d
    const std::string RSD13_d;
    
    //! 1-loop momentum integrals - RSD13_e
    const std::string RSD13_e;
    
    //! 1-loop momentum integrals - RSD13_f
    const std::string RSD13_f;
    
    //! 1-loop momentum integrals - RSD13_g
    const std::string RSD13_g;
    
    //! 1-loop momentum integrals - RSD22_A1
    const std::string RSD22_A1;
    
    //! 1-loop momentum integrals - RSD22_A2
    const std::string RSD22_A2;
    
    //! 1-loop momentum integrals - RSD22_A3
    const std::string RSD22_A3;
    
    //! 1-loop momentum integrals - RSD22_A4
    const std::string RSD22_A4;
    
    //! 1-loop momentum integrals - RSD22_A5
    const std::string RSD22_A5;
    
    //! 1-loop momentum integrals - RSD22_B2
    const std::string RSD22_B2;
    
    //! 1-loop momentum integrals - RSD22_B3
    const std::string RSD22_B3;
    
    //! 1-loop momentum integrals - RSD22_B6
    const std::string RSD22_B6;
    
    //! 1-loop momentum integrals - RSD22_B8
    const std::string RSD22_B8;
    
    //! 1-loop momentum integrals - RSD22_B9
    const std::string RSD22_B9;
    
    //! 1-loop momentum integrals - RSD22_C1
    const std::string RSD22_C1;
    
    //! 1-loop momentum integrals - RSD22_C2
    const std::string RSD22_C2;
    
    //! 1-loop momentum integrals - RSD22_C4
    const std::string RSD22_C4;
    
    //! 1-loop momentum integrals - RSD22_D1
    const std::string RSD22_D1;
    
    //! 1-loop delta delta power spectrum table
    const std::string dd_Pk;
    
    //! 1-loop delta_s delta_s power spectrum table, mu^0 coefficient
    const std::string dd_rsd_mu0_Pk;
    
    //! 1-loop delta_s delta_s power spectrum table, mu^2 coefficient
    const std::string dd_rsd_mu2_Pk;
    
    //! 1-loop delta_s delta_s power spectrum table, mu^4 coefficient
    const std::string dd_rsd_mu4_Pk;
    
    //! 1-loop delta_s delta_s power spectrum table, mu^6 coefficient
    const std::string dd_rsd_mu6_Pk;
    
    //! 1-loop delta_s delta_s power spectrum table, mu^8 coefficient
    const std::string dd_rsd_mu8_Pk;
    
    //! 1-loop multipole P0 table
    const std::string P0;
    
    //! 1-loop multipole P2 table
    const std::string P2;
    
    //! 1-loop multipole P4 table
    const std::string P4;
    
    //! Matsubara X & Y coefficients
    const std::string Matsubara_XY;

    //! temporary table name
    const std::string temp;

  };


#endif //LSSEFT_SQLITE3_POLICY_H
