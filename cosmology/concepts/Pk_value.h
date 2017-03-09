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

#ifndef LSSEFT_PK_VALUE_H
#define LSSEFT_PK_VALUE_H


#include "loop_integral.h"
#include "units/Mpc_units.h"

#include "boost/serialization/serialization.hpp"


template <typename ValueType>
class Pk_value_group
  {
  
  public:
    
    typedef ValueType value_type;
    typedef ValueType error_type;
    
    //! value constructor
    Pk_value_group(value_type v, value_type e=value_type(0.0))
      : value(std::move(v)),
        error(std::move(e))
      {
      }
    
    //! empty constructor
    Pk_value_group()
      : value(value_type(0.0)),
        error(value_type(0.0))
      {
      }
    
    
    // ACCESSORS
  
  public:
    
    //! get value
    const value_type& get_value() const { return this->value; }
    
    //! get error
    const value_type& get_error() const { return this->error; }
    
    
    // INTERNAL DATA
  
  private:
    
    //! value for this Pk component
    value_type value;
    
    //! erorr estimate for this Pk component
    value_type error;
  
  
  private:
    
    // enable boost::serialization support
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & value;
        ar & error;
      }
    
  };


// convenience typedefs for representing power spectra and counterterms
typedef Pk_value_group<Mpc_units::inverse_energy3> Pk_element;
typedef Pk_value_group<Mpc_units::inverse_energy>  k2_Pk_element;


//! overload + and -
template <typename ValueType>
Pk_value_group<ValueType> operator+(const Pk_value_group<ValueType>& A, const Pk_value_group<ValueType>& B)
  {
    auto value = A.get_value() + B.get_value();
    
    double A_error = static_cast<double>(A.get_error());
    double B_error = static_cast<double>(B.get_error());
    double error = std::sqrt(A_error * A_error + B_error * B_error);
    
    return Pk_value_group<ValueType>(value, ValueType(error));
  }

template <typename ValueType>
Pk_value_group<ValueType> operator-(const Pk_value_group<ValueType>& A, const Pk_value_group<ValueType>& B)
  {
    auto value = A.get_value() - B.get_value();
    
    double A_error = static_cast<double>(A.get_error());
    double B_error = static_cast<double>(B.get_error());
    double error = std::sqrt(A_error * A_error + B_error * B_error);
    
    return Pk_value_group<ValueType>(value, ValueType(error));
  }

template <typename ValueType>
Pk_value_group<ValueType> operator+(const Pk_value_group<ValueType>& A, const ValueType& B)
  {
    auto value = A.get_value() + B;
    auto error = A.get_error();
    
    return Pk_value_group<ValueType>(value, error);
  }

template <typename ValueType>
Pk_value_group<ValueType> operator+(const ValueType& A, const Pk_value_group<ValueType>& B)
  {
    auto value = B.get_value() + A;
    auto error = B.get_error();
    
    return Pk_value_group<ValueType>(value, error);
  }

template <typename ValueType>
Pk_value_group<ValueType> operator-(const Pk_value_group<ValueType>& A, const ValueType& B)
  {
    auto value = A.get_value() - B;
    auto error = A.get_error();
    
    return Pk_value_group<ValueType>(value, error);
  }

template <typename ValueType>
Pk_value_group<ValueType> operator-(const ValueType& A, const Pk_value_group<ValueType>& B)
  {
    auto value = A - B.get_value();
    auto error = B.get_error();
    
    return Pk_value_group<ValueType>(value, error);
  }


//! overload scalar multiplication
template <typename PkType, typename ScalarType, typename std::enable_if_t< !std::is_integral<ScalarType>::value >* = nullptr>
auto operator*(const ScalarType& S, const Pk_value_group<PkType>& P) -> Pk_value_group<decltype(S*P.get_value())>
  {
    auto value = S * P.get_value();
    auto error = std::abs(S) * P.get_error();
    
    return Pk_value_group<decltype(value)>(value, error);
  }

template <typename PkType, typename ScalarType, typename std::enable_if_t< !std::is_integral<ScalarType>::value >* = nullptr>
auto operator*(const Pk_value_group<PkType>& P, const ScalarType& S) -> Pk_value_group<decltype(S*P.get_value())>
  {
    auto value = S * P.get_value();
    auto error = std::abs(S) * P.get_error();
    
    return Pk_value_group<decltype(value)>(value, error);
  }


//! overload scalar division
template <typename PkType, typename ScalarType, typename std::enable_if_t< !std::is_integral<ScalarType>::value >* = nullptr>
auto operator/(const Pk_value_group<PkType>& P, const ScalarType& S) -> Pk_value_group<decltype(P.get_value()/S)>
  {
    auto value = P.get_value() / S;
    auto error = P.get_error() / std::abs(S);
    
    return Pk_value_group<decltype(value)>(value, error);
  }


// AUTOMATIC CONVERSION OF LOOP INTEGRAL KERNELS TO PK VALUES


// addition and subtraction of loop integral results gives a Pk_value_group
template <typename ValueType>
Pk_value_group<ValueType> operator+(const loop_integral_result<ValueType>& A, const loop_integral_result<ValueType>& B)
  {
    auto value = A.value + B.value;
    double error = std::sqrt(static_cast<double>(A.error)*static_cast<double>(A.error) + static_cast<double>(B.error)*static_cast<double>(B.error));
    
    return Pk_value_group<ValueType>(value, ValueType(error));
  }

template <typename ValueType>
Pk_value_group<ValueType> operator-(const loop_integral_result<ValueType>& A, const loop_integral_result<ValueType>& B)
  {
    auto value = A.value - B.value;
    double error = std::sqrt(static_cast<double>(A.error)*static_cast<double>(A.error) + static_cast<double>(B.error)*static_cast<double>(B.error));
    
    return Pk_value_group<ValueType>(value, ValueType(error));
  }


// addition and subtraction of loop integral results with Pk_value_group<>s gives Pk_value_group<>
template <typename ValueType>
Pk_value_group<ValueType> operator+(const Pk_value_group<ValueType>& A, const loop_integral_result<ValueType>& B)
  {
    auto value = A.get_value() + B.value;
    
    double A_error = static_cast<double>(A.get_error());
    double B_error = static_cast<double>(B.error);
    double error = std::sqrt(A_error * A_error + B_error * B_error);
    
    return Pk_value_group<ValueType>(value, ValueType(error));
  }


template <typename ValueType>
Pk_value_group<ValueType> operator+(const loop_integral_result<ValueType>& A, const Pk_value_group<ValueType>& B)
  {
    auto value = A.value + B.get_value();
    
    double A_error = static_cast<double>(A.error);
    double B_error = static_cast<double>(B.get_error());
    double error = std::sqrt(A_error * A_error + B_error * B_error);
    
    return Pk_value_group<ValueType>(value, ValueType(error));
  }


template <typename ValueType>
Pk_value_group<ValueType> operator-(const Pk_value_group<ValueType>& A, const loop_integral_result<ValueType>& B)
  {
    auto value = A.get_value() - B.value;
    
    double A_error = static_cast<double>(A.get_error());
    double B_error = static_cast<double>(B.error);
    double error = std::sqrt(A_error * A_error + B_error * B_error);
    
    return Pk_value_group<ValueType>(value, ValueType(error));
  }


template <typename ValueType>
Pk_value_group<ValueType> operator-(const loop_integral_result<ValueType>& A, const Pk_value_group<ValueType>& B)
  {
    auto value = A.value - B.get_value();
    
    double A_error = static_cast<double>(A.error);
    double B_error = static_cast<double>(B.get_error());
    double error = std::sqrt(A_error * A_error + B_error * B_error);
    
    return Pk_value_group<ValueType>(value, ValueType(error));
  }


// allow multiplication of loop integral results by scalar-type objects
template <typename ScalarType, typename LoopType, typename std::enable_if_t< !std::is_integral<ScalarType>::value >* = nullptr>
auto operator*(const ScalarType& S, const loop_integral_result<LoopType>& L) -> Pk_value_group<decltype(S*L.value)>
  {
    return Pk_value_group<decltype(S*L.value)>(S * L.value, std::abs(S) * L.error);
  };


template <typename ScalarType, typename LoopType, typename std::enable_if_t< !std::is_integral<ScalarType>::value >* = nullptr>
auto operator*(const loop_integral_result<LoopType>& L, const ScalarType& S) -> Pk_value_group<decltype(S*L.value)>
  {
    return Pk_value_group<decltype(S*L.value)>(S * L.value, std::abs(S) * L.error);
  };


// multiplication of loop integral results by Pk_value_group<> gives a Pk_value_group<>
template <typename PkType, typename LoopType>
auto operator*(const Pk_value_group<PkType> P, const loop_integral_result<LoopType>& L) -> Pk_value_group<decltype(P.get_value()*L.value)>
  {
    auto value = P.get_value() * L.value;
    
    double relP_err = P.get_error() / P.get_value();
    double relL_err = L.error / L.value;
    
    if(!std::isfinite(relP_err)) relP_err = 0.0;
    if(!std::isfinite(relL_err)) relL_err = 0.0;
    
    double quadrature = std::sqrt(relP_err * relP_err + relL_err * relL_err);
    
    return Pk_value_group<decltype(value)>(value, std::abs(value) * quadrature);
  };


template <typename PkType, typename LoopType>
auto operator*(const loop_integral_result<LoopType>& L, const Pk_value_group<PkType> P) -> Pk_value_group<decltype(P.get_value()*L.value)>
  {
    auto value = P.get_value() * L.value;
    
    double relP_err = P.get_error() / P.get_value();
    double relL_err = L.error / L.value;
    
    if(!std::isfinite(relP_err)) relP_err = 0.0;
    if(!std::isfinite(relL_err)) relL_err = 0.0;
    
    double quadrature = std::sqrt(relP_err * relP_err + relL_err * relL_err);
    
    return Pk_value_group<decltype(value)>(value, std::abs(value) * quadrature);
  };


// allow multiplication of Pk_value_group<> objects by other Pk_value_group<> objects
template <typename PkTypeA, typename PkTypeB>
auto operator*(const Pk_value_group<PkTypeA> PA, const Pk_value_group<PkTypeB>& PB) -> Pk_value_group<decltype(PA.get_value()*PB.get_value())>
  {
    auto value = PA.get_value() * PB.get_value();
    
    double relPA_err = PA.get_error() / PA.get_value();
    double relPB_err = PB.get_error() / PB.get_value();
    
    if(!std::isfinite(relPA_err)) relPA_err = 0.0;
    if(!std::isfinite(relPB_err)) relPB_err = 0.0;
    
    double quadrature = std::sqrt(relPA_err * relPA_err + relPB_err * relPB_err);
    
    return Pk_value_group<decltype(value)>(value, std::abs(value) * quadrature);
  };


#endif //LSSEFT_PK_VALUE_H
