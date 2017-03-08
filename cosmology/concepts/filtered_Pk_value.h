//
// Created by David Seery on 05/12/2016.
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

#ifndef LSSEFT_FILTERED_PK_H
#define LSSEFT_FILTERED_PK_H


#include "database/tokens.h"
#include "units/Mpc_units.h"

#include "boost/serialization/serialization.hpp"

class filtered_Pk_value
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    filtered_Pk_value(const k_token& kt, const linear_Pk_token& Pt, Mpc_units::inverse_energy3 _Pk_nw,
                      Mpc_units::inverse_energy3 _Pk_raw, Mpc_units::inverse_energy3 _Pk_ref);
    
    //! empty constructor, used only for receiving MPI payloads
    filtered_Pk_value();
    
    //! destructor is default
    ~filtered_Pk_value() = default;
    
    
    // INTERFACE
  
  public:
    
    //! get failure state
    bool get_fail() const { return this->fail; }
    
    //! set failed flag
    void mark_failed() { this->fail = true; }
    
    
    //! get wavenumber token
    const k_token& get_k_token() const { return this->k_tok; }
    
    //! get power spectrum token
    const linear_Pk_token& get_Pk_token() const { return this->Pk_tok; }
    
    //! get no-wiggle power spectrum
    const Mpc_units::inverse_energy3 get_Pk_nowiggle() const { return this->Pk_nw; }
    
    //! get raw power spectrum
    const Mpc_units::inverse_energy3 get_Pk_raw() const { return this->Pk_raw; }
    
    //! get reference power spectrum
    const Mpc_units::inverse_energy3 get_Pk_ref() const { return this->Pk_ref; }
    
    
    // INTERNAL DATA
    
  private:
    
    //! failure state
    bool fail;
    
    
    // CONFIGURATION DATA
    
    //! wavenumber token
    k_token k_tok;
    
    //! power spectrum token
    linear_Pk_token Pk_tok;
    
    
    // PAYLOAD DATA
    
    //! no-wiggle power spectrum
    Mpc_units::inverse_energy3 Pk_nw;
    
    //! raw power spectrum
    Mpc_units::inverse_energy3 Pk_raw;
    
    //! reference power spectrum
    Mpc_units::inverse_energy3 Pk_ref;
 
 
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & fail;
        ar & k_tok;
        ar & Pk_tok;
        ar & Pk_nw;
        ar & Pk_raw;
        ar & Pk_ref;
      }
    
  };


#endif //LSSEFT_FILTERED_PK_H
