//
// Created by David Seery on 21/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_MATSUBARA_A_H
#define LSSEFT_MATSUBARA_A_H


#include "units/Mpc_units.h"
#include "database/tokens.h"


class Matsubara_A
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor
    Matsubara_A(const IR_resum_token& _IRt, const Mpc_units::inverse_energy2& _A);
    
    //! empty constructor, used for receiving MPI payloads
    Matsubara_A();
    
    //! destructor is default
    ~Matsubara_A() = default;
    
    
    // INTERFACE
  
  public:

    //! overload * operator to return value
    const Mpc_units::inverse_energy2 operator*() const { return this->A; }
    
    //! allow automatic conversion to Mpc_units::inverse_energy2
    operator Mpc_units::inverse_energy2() const { return this->A; }
    
    //! accessor for value
    const Mpc_units::inverse_energy2 get_value() const { return this->A; }

    //! accessor for database token
    const IR_resum_token& get_token() const { return this->IR_resum_tok; }
    
    
    
    // INTERNAL DATA
  
  private:
    
    //! token
    IR_resum_token IR_resum_tok;
    
    //! value
    Mpc_units::inverse_energy2 A;
    
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & IR_resum_tok;
        ar & A;
      }
    
    
  };


#endif //LSSEFT_MATSUBARA_A_H
