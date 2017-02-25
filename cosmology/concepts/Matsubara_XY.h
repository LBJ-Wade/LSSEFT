//
// Created by David Seery on 21/11/2016.
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

#ifndef LSSEFT_MATSUBARA_XY_H
#define LSSEFT_MATSUBARA_XY_H


#include "units/Mpc_units.h"
#include "database/tokens.h"


class Matsubara_XY
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor
    Matsubara_XY(const linear_Pk_token& Pkt, const IR_resum_token& IRt, const Mpc_units::inverse_energy2& _X,
                 const Mpc_units::inverse_energy2& _Y);
    
    //! empty constructor, used for receiving MPI payloads
    Matsubara_XY();
    
    //! destructor is default
    ~Matsubara_XY() = default;
    
    
    // OPERATOR OVERLOADS
  
  public:

    //! overload * operator to return X + Y value
    const Mpc_units::inverse_energy2 operator*() const { return this->X + this->Y; }
    
    //! allow automatic conversion to Mpc_units::inverse_energy2
    operator Mpc_units::inverse_energy2() const { return *(*this); }
    
    
    // ACCESSORS
    
  public:
    
    //! accessor: X coefficient
    const Mpc_units::inverse_energy2 get_X() const { return this->X; }
    
    //! accessor: Y coefficient
    const Mpc_units::inverse_energy2 get_Y() const { return this->Y; }

    //! accessor for database token
    const IR_resum_token& get_IR_resum_token() const { return this->IR_resum_tok; }
    
    //! accessor for power spectrum token
    const linear_Pk_token& get_Pk_token() const { return this->Pk_lin; }
    
    
    // INTERNAL DATA
  
  private:
    
    //! IR resummation scale token
    IR_resum_token IR_resum_tok;
    
    //! linear power spectrum token
    linear_Pk_token Pk_lin;
    
    //! X
    Mpc_units::inverse_energy2 X;
    
    //! Y
    Mpc_units::inverse_energy2 Y;
    
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & IR_resum_tok;
        ar & Pk_lin;
        ar & X;
        ar & Y;
      }
    
    
  };


#endif //LSSEFT_MATSUBARA_XY_H
