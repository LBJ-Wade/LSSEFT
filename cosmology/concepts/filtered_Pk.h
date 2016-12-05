//
// Created by David Seery on 05/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_FILTERED_PK_H
#define LSSEFT_FILTERED_PK_H


#include "database/tokens.h"
#include "units/Mpc_units.h"

#include "boost/serialization/serialization.hpp"

class filtered_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    filtered_Pk(const k_token& kt, const linear_Pk_token& Pt, Mpc_units::inverse_energy3 _Pk_w, Mpc_units::inverse_energy3 _Pk_raw);
    
    //! empty constructor, used only for receiving MPI payloads
    filtered_Pk();
    
    //! destructor is default
    ~filtered_Pk() = default;
    
    
    // INTERFACE
  
  public:
    
    //! get wavenumber token
    const k_token& get_k_token() const { return this->k_tok; }
    
    //! get power spectrum token
    const linear_Pk_token& get_Pk_token() const { return this->Pk_tok; }
    
    //! get no-wiggle power spectrum
    const Mpc_units::inverse_energy3 get_Pk_w() const { return this->Pk_w; }
    
    //! get raw power spectrum
    const Mpc_units::inverse_energy3 get_Pk_raw() const { return this->Pk_raw; }
    
    
    // INTERNAL DATA
    
  private:
    
    // CONFIGURATION DATA
    
    //! wavenumber token
    k_token k_tok;
    
    //! power spectrum token
    linear_Pk_token Pk_tok;
    
    
    // PAYLOAD DATA
    
    //! no-wiggle power spectrum
    Mpc_units::inverse_energy3 Pk_w;
    
    //! raw power spectrum
    Mpc_units::inverse_energy3 Pk_raw;
 
 
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & k_tok;
        ar & Pk_tok;
        ar & Pk_w;
        ar & Pk_raw;
      }
    
  };


#endif //LSSEFT_FILTERED_PK_H
