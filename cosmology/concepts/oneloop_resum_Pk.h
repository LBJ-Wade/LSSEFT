//
// Created by David Seery on 09/12/2016.
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

#ifndef LSSEFT_ONE_LOOP_RESUM_PK_H
#define LSSEFT_ONE_LOOP_RESUM_PK_H


#include "oneloop_Pk.h"


// define convenience types for resummed versions of the basic power spectrum (~ 1/k^3) and k^2 * basic power spectrum (~ 1/k)
// these differ from Pk_value and k2_Pk_value by dropping the need to keep raw & nowiggle information
typedef Pk_value_group<Mpc_units::inverse_energy3> resum_Pk_value;
typedef Pk_value_group<Mpc_units::inverse_energy>  k2_resum_Pk_value;

// define a container class for the power spectrum data
typedef generic_dd_Pk<resum_Pk_value, k2_resum_Pk_value> resum_dd_Pk;


class oneloop_resum_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    oneloop_resum_Pk(const k_token& kt, const linear_Pk_token& Pkt, const IR_cutoff_token& IRt,
                     const UV_cutoff_token& UVt, const z_token& zt, const IR_resum_token& IRrt,
                     const resum_dd_Pk& Pkr);
    
    //! empty constructor, used when receiving an MPI payload
    oneloop_resum_Pk();
    
    //! destructor is default
    ~oneloop_resum_Pk() = default;
    
    
    // INTERFACE
  
  public:
    
    //! get wavenumber token
    const k_token& get_k_token() const { return this->k; }
    
    //! get power spectrum token
    const linear_Pk_token& get_Pk_token() const { return this->Pk_lin; }
    
    //! get UV cutoff token
    const UV_cutoff_token& get_UV_cutoff_token() const { return this->UV_cutoff; }
    
    //! get IR cutoff token
    const IR_cutoff_token& get_IR_cutoff_token() const { return this->IR_cutoff; }
    
    //! get z token
    const z_token& get_z_token() const { return this->z; }
    
    //! get IR resummation token
    const IR_resum_token& get_IR_resum_token() const { return this->IR_resum; }
    
    //! get resummed power spectrum
    const resum_dd_Pk& get_Pk_resum() const { return this->Pk_resum; }
    
    
    // INTERNAL DATA
  
  private:
    
    // CONFIGURATION DATA
    
    //! wavenumber token
    k_token k;
    
    //! power spectrum token
    linear_Pk_token Pk_lin;
    
    //! UV cutoff token
    UV_cutoff_token UV_cutoff;
    
    //! IR cutoff token
    IR_cutoff_token IR_cutoff;
    
    //! redshift token
    z_token z;
    
    //! IR resummation scale token
    IR_resum_token IR_resum;
    
    
    // VALUES
    
    //! resummed power spectrum
    resum_dd_Pk Pk_resum;
    
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & k;
        ar & Pk_lin;
        ar & UV_cutoff;
        ar & IR_cutoff;
        ar & z;
        ar & IR_resum;
        ar & Pk_resum;
      }
    
  };


#endif //LSSEFT_ONE_LOOP_RESUM_PK_H
