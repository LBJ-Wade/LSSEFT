//
// Created by David Seery on 08/03/2017.
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

#ifndef LSSEFT_PK_RESUM_VALUE_H
#define LSSEFT_PK_RESUM_VALUE_H


#include "Pk_value.h"


template <typename ValueType>
class resum_Pk_component
  {
    
    // TYPES
  
  public:
    
    typedef ValueType value_type;
    typedef Pk_value_group<ValueType> element_type;
    
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! value constructor
    resum_Pk_component(const element_type _raw, const element_type _resum)
      : raw(std::move(_raw)),
        resum(std::move(_resum))
      {
      }
    
    //! empty constructor
    resum_Pk_component()
      : raw(),
        resum()
      {
      }
    
    //! destructor is default
    ~resum_Pk_component() = default;
    
    
    // ACCESSORS
  
  public:
    
    //! get raw value
    const element_type& get_raw() const { return this->raw; }
    
    //! get resummed value
    const element_type& get_resum() const { return this->resum; }
    
    //! set raw value
    resum_Pk_component<ValueType>& set_raw(Pk_value_group<ValueType> r) { this->raw = std::move(r); return *this; }
    
    //! set resummed value
    resum_Pk_component<ValueType>& set_resum(Pk_value_group<ValueType> r) { this->resum = std::move(r); return *this; }
    
    
    // INTERNAL DATA
  
  private:
    
    //! raw value
    element_type raw;
    
    //! resummed value
    element_type resum;
    
    
    // enable boost::serialization support
    friend class boost::serialization::access;
    
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & raw;
        ar & resum;
      }
    
  };


typedef resum_Pk_component<Mpc_units::inverse_energy> k2_Pk_resum;
typedef resum_Pk_component<Mpc_units::inverse_energy3> Pk_resum;


//! addition
template <typename ValueType>
resum_Pk_component<ValueType> operator+(const resum_Pk_component<ValueType>& a, const resum_Pk_component<ValueType>& b)
  {
    return resum_Pk_component<ValueType>{a.get_raw() + b.get_raw(), a.get_resum() + b.get_resum()};
  }


//! subtraction
template <typename ValueType>
resum_Pk_component<ValueType> operator-(const resum_Pk_component<ValueType>& a, const resum_Pk_component<ValueType>& b)
  {
    return resum_Pk_component<ValueType>{a.get_raw() - b.get_raw(), a.get_resum() - b.get_resum()};
  }


//! multiplication by number
template <typename ValueType>
resum_Pk_component<ValueType> operator*(double a, const resum_Pk_component<ValueType>& b)
  {
    return resum_Pk_component<ValueType>{a * b.get_raw(), a * b.get_resum()};
  }


template <typename ValueType>
resum_Pk_component<ValueType> operator*(const resum_Pk_component<ValueType>& a, double b)
  {
    return b*a;
  }


#endif //LSSEFT_PK_RESUM_VALUE_H
