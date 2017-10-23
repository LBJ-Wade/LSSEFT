//
// Created by David Seery on 12/08/2015.
// --@@ // Copyright (c) 2017 University of Sussex. All rights reserved.
//
// This file is part of the Sussex Effective Field Theory for
// Large-Scale Structure platform (LSSEFT).
//
// LSSEFT is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// LSSEFT is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with LSSEFT.  If not, see <http://www.gnu.org/licenses/>.
//
// @license: GPL-2
// @contributor: David Seery <D.Seery@sussex.ac.uk>
// --@@
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
    
    //! filtering parameter configurations
    const std::string& filter_config_table() const { return(this->filter_config); }
    
    //! oneloop parameter configurations
    const std::string& loop_integral_config_table() const { return(this->loop_integral_config); }
    
    //! MatsubaraXY parameter configurations
    const std::string& MatsubaraXY_config_table() const { return(this->MatsubaraXY_config); }
    
    //! growth parameter configurations
    const std::string& growth_config_table() const { return(this->growth_config); }
    
    //! transfer function table
    const std::string& transfer_table() const { return(this->transfer); }

    //! 1-loop growth g-factor table
    const std::string& D_factor_table() const { return(this->growth_D_factor); }
    
    //! 1-loop growth f-factor table
    const std::string& f_factor_table() const { return(this->growth_f_factor); }

    //! Matsubara-XY table
    const std::string& Matsubara_XY_table() const { return this->Matsubara_XY; }

    //! counterterms table - P0
    const std::string& counterterms_c0_table() const { return this->counterterms_c0; }

    //! counterterms table - P2
    const std::string& counterterms_c2_table() const { return this->counterterms_c2; }

    //! counterterms table - P4
    const std::string& counterterms_c4_table() const { return this->counterterms_c4; }

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
    
    //! filtering parameter configurations
    const std::string filter_config;
    
    //! one-loop parameter configurations
    const std::string loop_integral_config;
    
    //! MatsubaraXY parameter configurations
    const std::string MatsubaraXY_config;
    
    //! growth function parameter configurations
    const std::string growth_config;

    //! transfer function table
    const std::string transfer;

    //! 1-loop growth g-factor table
    const std::string growth_D_factor;
    
    //! 1-loop growth f-factor table
    const std::string growth_f_factor;

    //! Matsubara X & Y coefficients
    const std::string Matsubara_XY;

    //! counterterms table - P0
    const std::string counterterms_c0;

    //! counterterms table - P2
    const std::string counterterms_c2;

    //! counterterms table - P4
    const std::string counterterms_c4;

    //! temporary table name
    const std::string temp;

  };


#endif //LSSEFT_SQLITE3_POLICY_H
