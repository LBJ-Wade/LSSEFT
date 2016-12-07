//
// Created by David Seery on 14/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_ONE_LOOP_PK_H
#define LSSEFT_ONE_LOOP_PK_H


#include "loop_integral.h"

#include "database/tokens.h"
#include "units/Mpc_units.h"

#include "boost/timer/timer.hpp"
#include "boost/serialization/serialization.hpp"


template <typename ValueType>
class Pk_value_group
  {
    
  public:
    
    typedef ValueType value_type;
    
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


//! overload + and -
template <typename ValueType>
Pk_value_group<ValueType> operator+(const Pk_value_group<ValueType>& A, const Pk_value_group<ValueType>& B)
  {
    auto value = A.get_value() + B.get_value();
    
    double A_error = static_cast<double>(A.get_error());
    double B_error = static_cast<double>(B.get_error());
    double error = std::sqrt(A_error*A_error + B_error*B_error);
    
    return Pk_value_group<ValueType>(value, ValueType(error));
  }

template <typename ValueType>
Pk_value_group<ValueType> operator-(const Pk_value_group<ValueType>& A, const Pk_value_group<ValueType>& B)
  {
    auto value = A.get_value() - B.get_value();
    
    double A_error = static_cast<double>(A.get_error());
    double B_error = static_cast<double>(B.get_error());
    double error = std::sqrt(A_error*A_error + B_error*B_error);
    
    return Pk_value_group<ValueType>(value, ValueType(error));
  }


//! overload scalar multiplication
template <typename PkType, typename ScalarType>
auto operator*(const ScalarType& S, const Pk_value_group<PkType>& P) -> Pk_value_group<decltype(S*P.get_value())>
  {
    auto value = S * P.get_value();
    auto error = std::abs(S) * P.get_error();
    
    return Pk_value_group<decltype(value)>(value, error);
  }

template <typename PkType, typename ScalarType>
auto operator*(const Pk_value_group<PkType>& P, const ScalarType& S) -> Pk_value_group<decltype(S*P.get_value())>
  {
    auto value = S * P.get_value();
    auto error = std::abs(S) * P.get_error();
    
    return Pk_value_group<decltype(value)>(value, error);
  }


//! overload scalar division
template <typename PkType, typename ScalarType>
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
    double error = std::sqrt(A_error*A_error + B_error*B_error);
    
    return Pk_value_group<ValueType>(value, error);
  }

template <typename ValueType>
Pk_value_group<ValueType> operator+(const loop_integral_result<ValueType>& A, const Pk_value_group<ValueType>& B)
  {
    auto value = A.value + B.get_value();
    
    double A_error = static_cast<double>(A.error);
    double B_error = static_cast<double>(B.get_error());
    double error = std::sqrt(A_error*A_error + B_error*B_error);
    
    return Pk_value_group<ValueType>(value, error);
  }

template <typename ValueType>
Pk_value_group<ValueType> operator-(const Pk_value_group<ValueType>& A, const loop_integral_result<ValueType>& B)
  {
    auto value = A.get_value() - B.value;
    
    double A_error = static_cast<double>(A.get_error());
    double B_error = static_cast<double>(B.error);
    double error = std::sqrt(A_error*A_error + B_error*B_error);
    
    return Pk_value_group<ValueType>(value, error);
  }

template <typename ValueType>
Pk_value_group<ValueType> operator-(const loop_integral_result<ValueType>& A, const Pk_value_group<ValueType>& B)
  {
    auto value = A.value - B.get_value();
    
    double A_error = static_cast<double>(A.error);
    double B_error = static_cast<double>(B.get_error());
    double error = std::sqrt(A_error*A_error + B_error*B_error);
    
    return Pk_value_group<ValueType>(value, error);
  }


// multiplication of a loop integral result by a scalar converts to a Pk_value_group<>
template <typename ValueType>
Pk_value_group<ValueType> operator*(double A, const loop_integral_result<ValueType>& B)
  {
    return Pk_value_group<ValueType>(A * B.value, std::abs(A) * B.error);
  }

template <typename ValueType>
Pk_value_group<ValueType> operator*(const loop_integral_result<ValueType>& A, double B)
  {
    return Pk_value_group<ValueType>(B * A.value, std::abs(B) * A.error);
  }


// allow multiplication of loop integral results by scalar-type objects
template <typename ScalarType, typename LoopType>
auto operator*(const ScalarType& S, const loop_integral_result<LoopType>& L) -> Pk_value_group<decltype(S*L.value)>
{
    return Pk_value_group<decltype(S*L.value)>(S * L.value, std::abs(S) * L.error);
  };

template <typename ScalarType, typename LoopType>
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
    double relPB_err = PB.get_value() / PB.get_value();
    
    if(!std::isfinite(relPA_err)) relPA_err = 0.0;
    if(!std::isfinite(relPB_err)) relPB_err = 0.0;
    
    double quadrature = std::sqrt(relPA_err * relPA_err + relPB_err * relPB_err);
    
    return Pk_value_group<decltype(value)>(value, std::abs(value) * quadrature);
  };


template <typename ValueType>
class dimensionful_Pk_component
  {
    
  public:
    
    typedef ValueType value_type;
    typedef Pk_value_group<ValueType> container_type;
    
    //! value constructor; error is zero if not specified which allows
    //! assignment-on-construction to a fixed number
    dimensionful_Pk_component(Pk_value_group<ValueType> r, Pk_value_group<ValueType> w)
      : raw(std::move(r)),
        wiggle(std::move(w))
      {
      }
    
    //! empty constructor
    dimensionful_Pk_component()
      : raw(),
        wiggle()
      {
      }
    
    
    // CONVERSIONS
    
  public:
    
//    //! allow direct assignment from a value_type object, in which case there is zero error
//    dimensionful_Pk_component<ValueType>& operator=(value_type v)
//      {
//        this->value = std::move(v);
//        this->error = value_type(0.0);
//        return *this;
//      }


    // ACCESSOR
    
  public:
    
    //! get raw data
    const Pk_value_group<ValueType>& get_raw() const { return this->raw; }
    
    //! get wiggle data
    const Pk_value_group<ValueType>& get_wiggle() const { return this->wiggle; }
    
    //! get no-wiggle data
    const Pk_value_group<ValueType> get_nowiggle() const { return this->raw - this->wiggle; }
    
    //! set raw data
    void set_raw(Pk_value_group<ValueType> r) { this->raw = r; }
    
    //! set wiggle data
    void set_wiggle(Pk_value_group<ValueType> w) { this->wiggle = w; }
    
    
    // DATA
    
  private:
    
    //! raw data
    Pk_value_group<ValueType> raw;
    
    //! wiggle data
    Pk_value_group<ValueType> wiggle;
  
  
  private:
    
    // enable boost::serialization support
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & raw;
        ar & wiggle;
      }
    
  };


// define convenience types for basic power spectrum (~ 1/k^3) and k^2 * basic power spectrum (~ 1/k)
typedef dimensionful_Pk_component<Mpc_units::inverse_energy3> Pk_value;
typedef dimensionful_Pk_component<Mpc_units::inverse_energy>  k2_Pk_value;


//! overload + and - so that power spectrum and counterterm values can be added
template <typename ValueType>
dimensionful_Pk_component<ValueType> operator+(const dimensionful_Pk_component<ValueType>& A, const dimensionful_Pk_component<ValueType>& B)
  {
    return dimensionful_Pk_component<ValueType>(A.get_raw() + B.get_raw(), A.get_wiggle() + B.get_wiggle());
  }

template <typename ValueType>
dimensionful_Pk_component<ValueType> operator-(const dimensionful_Pk_component<ValueType>& A, const dimensionful_Pk_component<ValueType>& B)
  {
    return dimensionful_Pk_component<ValueType>(A.get_raw() - B.get_raw(), A.get_wiggle() - B.get_wiggle());
  }


//! overload scalar multiplication
template <typename ScalarType, typename PkType>
auto operator*(const ScalarType& S, const dimensionful_Pk_component<PkType>& P) -> dimensionful_Pk_component<decltype(S*P.get_raw().get_value())>
  {
    return dimensionful_Pk_component<decltype(S*P.get_raw().get_value())>(S * P.get_raw(), S * P.get_wiggle());
  }

template <typename ScalarType, typename PkType>
auto operator*(const dimensionful_Pk_component<PkType>& P, const ScalarType& S) -> dimensionful_Pk_component<decltype(S*P.get_raw().get_value())>
  {
    return dimensionful_Pk_component<decltype(S*P.get_raw().get_value())>(S * P.get_raw(), S * P.get_wiggle());
  }


//! overload scalar division
template <typename ScalarType, typename PkType>
auto operator/(const dimensionful_Pk_component<PkType>& P, const ScalarType& S) -> dimensionful_Pk_component<decltype(P.get_raw().get_value()/S)>
  {
    return dimensionful_Pk_component<decltype(P.get_raw().get_value()/S)>(P.get_raw() / S, P.get_wiggle() / S);
  }



// AUTOMATIC CONVERSION OF LOOP INTEGRAL KERNELS TO PK COMPONENTS


// multiplication of any loop integral output by a scalar converts to a dimensionful_Pk_component<>
template <typename ScalarType, typename LoopType>
auto operator*(const ScalarType& S, const loop_integral_output<LoopType>& P) -> dimensionful_Pk_component<decltype(S*P.get_raw().value)>
  {
    return dimensionful_Pk_component<decltype(S*P.get_raw().value)>(S * P.get_raw(), S * P.get_wiggle());
  }

template <typename ScalarType, typename LoopType>
auto operator*(const loop_integral_output<LoopType>& P, const ScalarType& S) -> dimensionful_Pk_component<decltype(S*P.get_raw().value)>
  {
    return dimensionful_Pk_component<decltype(S*P.get_raw().value)>(S * P.get_raw(), S * P.get_wiggle());
  }


// multiplication of loop integral output by dimensionful_Pk_component<>
template <typename LoopType, typename PkType>
auto operator*(const dimensionful_Pk_component<PkType>& P, const loop_integral_output<LoopType>& L) -> dimensionful_Pk_component<decltype(P.get_raw().get_value()*L.get_raw().value)>
  {
    return dimensionful_Pk_component<decltype(P.get_raw().get_value()*L.get_raw().value)>(P.get_raw() * L.get_raw(), P.get_wiggle() * L.get_wiggle());
  }

template <typename LoopType, typename PkType>
auto operator*(const loop_integral_output<LoopType>& L, const dimensionful_Pk_component<PkType>& P) -> dimensionful_Pk_component<decltype(P.get_raw().get_value()*L.get_raw().value)>
  {
    return dimensionful_Pk_component<decltype(P.get_raw().get_value()*L.get_raw().value)>(P.get_raw() * L.get_raw(), P.get_wiggle() * L.get_wiggle());
  }

  
// addition and subtraction of loop integral outputs gives dimensionful_Pk_component<>
template <typename ValueType>
dimensionful_Pk_component<ValueType> operator+(const loop_integral_output<ValueType>& A, const loop_integral_output<ValueType>& B)
  {
    return dimensionful_Pk_component<ValueType>(A.get_raw() + B.get_raw(), A.get_wiggle() + B.get_wiggle());
  }

template <typename ValueType>
dimensionful_Pk_component<ValueType> operator-(const loop_integral_output<ValueType>& A, const loop_integral_output<ValueType>& B)
  {
    return dimensionful_Pk_component<ValueType>(A.get_raw() - B.get_raw(), A.get_wiggle() - B.get_wiggle());
  }


// addition and subtraction of loop integral outputs with dimensionful_Pk_component<>s gives dimensionful_Pk_component<>
template <typename ValueType>
dimensionful_Pk_component<ValueType> operator+(const dimensionful_Pk_component<ValueType>& A, const loop_integral_output<ValueType>& B)
  {
    return dimensionful_Pk_component<ValueType>(A.get_raw() + B.get_raw(), A.get_wiggle + B.get_wiggle());
  }

template <typename ValueType>
dimensionful_Pk_component<ValueType> operator+(const loop_integral_output<ValueType>& A, const dimensionful_Pk_component<ValueType>& B)
  {
    return dimensionful_Pk_component<ValueType>(A.get_raw() + B.get_raw(), A.get_wiggle + B.get_wiggle());
  }

template <typename ValueType>
dimensionful_Pk_component<ValueType> operator-(const dimensionful_Pk_component<ValueType>& A, const loop_integral_output<ValueType>& B)
  {
    return dimensionful_Pk_component<ValueType>(A.get_raw() - B.get_raw(), A.get_wiggle - B.get_wiggle());
  }

template <typename ValueType>
dimensionful_Pk_component<ValueType> operator-(const loop_integral_output<ValueType>& A, const dimensionful_Pk_component<ValueType>& B)
  {
    return dimensionful_Pk_component<ValueType>(A.get_raw() - B.get_raw(), A.get_wiggle - B.get_wiggle());
  }


// multiplication of dimensionful_Pk_component<>s
template <typename PkTypeA, typename PkTypeB>
auto operator*(const dimensionful_Pk_component<PkTypeA> PA, const dimensionful_Pk_component<PkTypeB>& PB) -> dimensionful_Pk_component<decltype(PA.get_raw().get_value()*PB.get_raw().get_value())>
  {
    return dimensionful_Pk_component<decltype(PA.get_raw().get_value()*PB.get_raw().get_value())>(PA.get_raw() * PB.get_raw(), PA.get_wiggle() * PB.get_wiggle());
  }


class dd_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    dd_Pk(const Pk_value& _Pt, const Pk_value& _P13, const Pk_value& _P22, const k2_Pk_value& _Z2d);
    
    //! empty constructor for use when overwriting with MPI payloads
    dd_Pk();
    
    //! destructor is default
    ~dd_Pk() = default;
    
    
    // INTERFACE
    
  public:
    
    //! get tree value
    Pk_value& get_tree() { return this->Ptree; }
    const Pk_value& get_tree() const { return this->Ptree; }
    
    //! get 13 value
    Pk_value& get_13() { return this->P13; }
    const Pk_value& get_13() const { return this->P13; }
    
    //! get 22 value
    Pk_value& get_22() { return this->P22; }
    const Pk_value& get_22() const { return this->P22; }
    
    //! get total SPT power spectrum = 13 + 22
    Pk_value& get_1loop_SPT() { return this->P1loopSPT; }
    const Pk_value& get_1loop_SPT() const { return this->P1loopSPT; }
    
    
    // COUNTERTERMS
    
    //! get EFT counterterm
    k2_Pk_value& get_Z2_delta() { return this->Z2_delta; }
    const k2_Pk_value& get_Z2_delta() const { return this->Z2_delta; }
    
    
    // INTERNAL DATA
    
  private:
    
    //! tree power spectrum
    Pk_value Ptree;
    
    //! 13 terms
    Pk_value P13;
    
    //! 22 terms
    Pk_value P22;
    
    //! total 1-loop SPT value
    Pk_value P1loopSPT;
    
    //! coefficient of the counterterm Z2_delta
    k2_Pk_value Z2_delta;
    
    
    // enable boost::serialization support
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & Ptree;
        ar & P13;
        ar & P22;
        ar & P1loopSPT;
        ar & Z2_delta;
      }
    
  };


class rsd_dd_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! value constructor
    rsd_dd_Pk(const Pk_value& _Pt, const Pk_value& _P13, const Pk_value& _P22,
              const k2_Pk_value& _Z2d, const Pk_value& _Z0v, const k2_Pk_value& _Z2v,
              const Pk_value& _Z0vd, const k2_Pk_value& _Z2vd,
              const k2_Pk_value& _Z2vv, const k2_Pk_value& _Z2vvd, const k2_Pk_value& _Z2vvv);
    
    //! empty constructor for use when overwriting with MPI payloads
    rsd_dd_Pk();
    
    //! destructor is default
    ~rsd_dd_Pk() = default;
    
    
    // INTERFACE
  
  public:
    
    //! get tree value
    Pk_value& get_tree() { return this->Ptree; }
    const Pk_value& get_tree() const { return this->Ptree; }
    
    //! get 13 value
    Pk_value& get_13() { return this->P13; }
    const Pk_value& get_13() const { return this->P13; }
    
    //! get 22 value
    Pk_value& get_22() { return this->P22; }
    const Pk_value& get_22() const { return this->P22; }
    
    //! get total 13 + 22
    Pk_value& get_1loop_SPT() { return this->P1loopSPT; }
    const Pk_value& get_1loop_SPT() const { return this->P1loopSPT; }
    
    
    // COUNTERTERMS
    
    //! get Z2_delta counterterm
    k2_Pk_value& get_Z2_delta() { return this->Z2_delta; }
    const k2_Pk_value& get_Z2_delta() const { return this->Z2_delta; }
    
    //! get Z0_v counterterm
    Pk_value& get_Z0_v() { return this->Z0_v; }
    const Pk_value& get_Z0_v() const { return this->Z0_v; }
    
    //! get Z2_v counterterm
    k2_Pk_value& get_Z2_v() { return this->Z2_v; }
    const k2_Pk_value& get_Z2_v() const { return this->Z2_v; }
    
    //! get Z0_vdelta counterterm
    Pk_value& get_Z0_vdelta() { return this->Z0_vdelta; }
    const Pk_value& get_Z0_vdelta() const { return this->Z0_vdelta; }
    
    //! get Z2_vdelta counterterm
    k2_Pk_value& get_Z2_vdelta() { return this->Z2_vdelta; }
    const k2_Pk_value& get_Z2_vdelta() const { return this->Z2_vdelta; }
    
    //! get Z2_vv counterterm
    k2_Pk_value& get_Z2_vv() { return this->Z2_vv; }
    const k2_Pk_value& get_Z2_vv() const { return this->Z2_vv; }
    
    //! get Z2_vvdelta counterterm
    k2_Pk_value& get_Z2_vvdelta() { return this->Z2_vvdelta; }
    const k2_Pk_value& get_Z2_vvdelta() const { return this->Z2_vvdelta; }
    
    //! get Z2_vvv counterterm
    k2_Pk_value& get_Z2_vvv() { return this->Z2_vvv; }
    const k2_Pk_value& get_Z2_vvv() const { return this->Z2_vvv; }
    
    
    // INTERNAL DATA
  
  private:
    
    //! tree power spectrum
    Pk_value Ptree;
    
    //! 13 terms
    Pk_value P13;
    
    //! 22 terms
    Pk_value P22;
    
    //! total 1-loop SPT value
    Pk_value P1loopSPT;
    
    //! coefficient of the counterterm Z2_delta
    k2_Pk_value Z2_delta;
    
    //! coefficient of the counterterm Z0_v
    Pk_value Z0_v;
    
    //! coefficient of the counterterm Z2_v
    k2_Pk_value Z2_v;
    
    //! coefficient of the counterterm Z0_vd
    Pk_value Z0_vdelta;
    
    //! coefficient of the counterterm Z2_vd
    k2_Pk_value Z2_vdelta;
    
    //! coefficient of the counterterm Z2_vv
    k2_Pk_value Z2_vv;
    
    //! coefficient of the counterterm Z2_vvdelta
    k2_Pk_value Z2_vvdelta;
    
    //! coefficient of the counterterm Z2_vvv
    k2_Pk_value Z2_vvv;
    
    
    // enable boost::serialization support
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & Ptree;
        ar & P13;
        ar & P22;
        ar & P1loopSPT;
        ar & Z2_delta;
        ar & Z0_v;
        ar & Z2_v;
        ar & Z0_vdelta;
        ar & Z2_vdelta;
        ar & Z2_vv;
        ar & Z2_vvdelta;
        ar & Z2_vvv;
      }
    
  };


class oneloop_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    oneloop_Pk(const k_token& kt, const linear_Pk_token& Pkt, const IR_cutoff_token& IRt,
               const UV_cutoff_token& UVt, const z_token& zt,
               const dd_Pk& _dd, const rsd_dd_Pk& _rsd_mu0, const rsd_dd_Pk& _rsd_mu2, const rsd_dd_Pk& _rsd_mu4,
               const rsd_dd_Pk& _rsd_mu6, const rsd_dd_Pk& _rsd_mu8);
    
    //! empty constructor, used when receiving an MPI payload
    oneloop_Pk();
    
    //! destructor is default
    ~oneloop_Pk() = default;
    
    
    // INTERFACE
    
  public:
    
    //! get wavenumber token
    const k_token& get_k_token() const { return this->k; }
    
    //! get power spectrum token
    const linear_Pk_token& get_Pk_token() const { return this->Pk_lin; }
    
    //! get UV cutoff token
    const UV_cutoff_token& get_UV_token() const { return this->UV_cutoff; }
    
    //! get IR cutoff token
    const IR_cutoff_token& get_IR_token() const { return this->IR_cutoff; }
    
    //! get z token
    const z_token& get_z_token() const { return this->z; }
    
    
    //! get delta-delta power spectrum
    const dd_Pk& get_dd() const { return this->dd; }
    
    //! get delta-delta RSD power spectrum mu^0 coefficient
    const rsd_dd_Pk& get_dd_rsd_mu0() const { return this->rsd_dd_mu0; }
    
    //! get delta-delta RSD power spectrum mu^2 coefficient
    const rsd_dd_Pk& get_dd_rsd_mu2() const { return this->rsd_dd_mu2; }
    
    //! get delta-delta RSD power spectrum mu^4 coefficient
    const rsd_dd_Pk& get_dd_rsd_mu4() const { return this->rsd_dd_mu4; }
    
    //! get delta-delta RSD power spectrum mu^6 coefficient
    const rsd_dd_Pk& get_dd_rsd_mu6() const { return this->rsd_dd_mu6; }
    
    //! get delta-delta RSD power spectrum mu^8 coefficient
    const rsd_dd_Pk& get_dd_rsd_mu8() const { return this->rsd_dd_mu8; }
    
    
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
    
    
    // VALUES
    
    //! delta-delta power spectrum
    dd_Pk dd;
    
    //! mu^0 term in delta_s-delta_s power spectrum
    rsd_dd_Pk rsd_dd_mu0;
    
    //! mu^2 term in delta_s-delta_s power spectrum
    rsd_dd_Pk rsd_dd_mu2;
    
    //! mu^4 term in delta_s-delta_s power spectrum
    rsd_dd_Pk rsd_dd_mu4;
    
    //! mu^6 term in delta_s-delta_s power spectrum
    rsd_dd_Pk rsd_dd_mu6;
    
    //! mu^8 term in delta_s-delta_s power spectrum
    rsd_dd_Pk rsd_dd_mu8;
    
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & k;
        ar & UV_cutoff;
        ar & IR_cutoff;
        ar & z;
        ar & dd;
        ar & rsd_dd_mu0;
        ar & rsd_dd_mu2;
        ar & rsd_dd_mu4;
        ar & rsd_dd_mu6;
        ar & rsd_dd_mu8;
      }
    
    
  };


#endif //LSSEFT_ONE_LOOP_PK_H
