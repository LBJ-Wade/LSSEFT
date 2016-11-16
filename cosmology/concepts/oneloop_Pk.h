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
class dimensionful_Pk_component
  {
    
  public:
    
    typedef ValueType value_type;
    
    //! value constructor; error is zero if not specified which allows
    //! assignment-on-construction to a fixed number
    dimensionful_Pk_component(value_type v, value_type e=value_type(0.0))
      : value(std::move(v)),
        error(std::move(e))
      {
      }
    
    //! empty constructor
    dimensionful_Pk_component()
      : value(value_type(0.0)),
        error(value_type(0.0))
      {
      }
    
    
    // CONVERSIONS
    
  public:
    
    //! allow direct assignment from a value_type object, in which case there is zero error
    dimensionful_Pk_component<ValueType>& operator=(value_type v)
      {
        this->value = std::move(v);
        this->error = value_type(0.0);
        return *this;
      }
    
    
    // DATA
    
  public:
    
    value_type value;
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


typedef dimensionful_Pk_component<Mpc_units::inverse_energy3> Pk_value;
typedef dimensionful_Pk_component<Mpc_units::inverse_energy>  k2_Pk_value;


//! overload + and - so that power spectrum and counterterm values can be added
template <typename ValueType>
dimensionful_Pk_component<ValueType> operator+(const dimensionful_Pk_component<ValueType>& A, const dimensionful_Pk_component<ValueType>& B)
  {
    dimensionful_Pk_component<ValueType> res;
    
    res.value = A.value + B.value;
    res.error = std::sqrt(static_cast<double>(A.error)*static_cast<double>(A.error) + static_cast<double>(B.error)*static_cast<double>(B.error));
    
    return std::move(res);
  }

template <typename ValueType>
dimensionful_Pk_component<ValueType> operator-(const dimensionful_Pk_component<ValueType>& A, const dimensionful_Pk_component<ValueType>& B)
  {
    dimensionful_Pk_component<ValueType> res;
    
    res.value = A.value - B.value;
    res.error = std::sqrt(static_cast<double>(A.error)*static_cast<double>(A.error) + static_cast<double>(B.error)*static_cast<double>(B.error));
    
    return std::move(res);
  }


//! overload scalar multiplication
template <typename ValueType>
dimensionful_Pk_component<ValueType> operator*(double A, const dimensionful_Pk_component<ValueType>& B)
  {
    dimensionful_Pk_component<ValueType> res;
    
    res.value = A * B.value;
    res.error = std::abs(A) * B.error;
    
    return std::move(res);
  }

template <typename ValueType>
dimensionful_Pk_component<ValueType> operator*(const dimensionful_Pk_component<ValueType>& A, double B)
  {
    dimensionful_Pk_component<ValueType> res;
    
    res.value = B * A.value;
    res.error = std::abs(B) * A.error;
    
    return std::move(res);
  }

//! overload scalar division
template <typename ValueType>
dimensionful_Pk_component<ValueType> operator/(const dimensionful_Pk_component<ValueType>& A, double B)
  {
    dimensionful_Pk_component<ValueType> res;
    
    res.value = A.value / B;
    res.error = A.error / std::abs(B);
    
    return std::move(res);
  }


// overload multiplication of a loop integral result to give a Pk component
template <typename ValueType>
dimensionful_Pk_component<ValueType> operator*(double A, const loop_integral_output<ValueType>& B)
  {
    dimensionful_Pk_component<ValueType> res;
    
    res.value = A * B.value;
    res.error = std::abs(A) * B.error;
    
    return std::move(res);
  }

template <typename ValueType>
dimensionful_Pk_component<ValueType> operator*(const loop_integral_output<ValueType>& A, double B)
  {
    dimensionful_Pk_component<ValueType> res;
    
    res.value = B * A.value;
    res.error = std::abs(B) * A.error;
    
    return std::move(res);
  }

inline dimensionful_Pk_component<Mpc_units::inverse_energy3>
operator*(const Mpc_units::inverse_energy3& A, const loop_integral_output<double>& B)
  {
    dimensionful_Pk_component<Mpc_units::inverse_energy3> res;
    
    res.value = A * B.value;
    res.error = std::abs(A) * B.error;
    
    return std::move(res);
  }

inline dimensionful_Pk_component<Mpc_units::inverse_energy3>
operator*(const loop_integral_output<double>& A, const Mpc_units::inverse_energy3& B)
  {
    dimensionful_Pk_component<Mpc_units::inverse_energy3> res;
    
    res.value = B * A.value;
    res.error = std::abs(B) * A.error;
    
    return std::move(res);
  }


// overload addition/subtraction of loop integral results to give a Pk component
template <typename ValueType>
dimensionful_Pk_component<ValueType> operator+(const loop_integral_output<ValueType>& A, const loop_integral_output<ValueType>& B)
  {
    dimensionful_Pk_component<ValueType> res;
    
    res.value = A.value + B.value;
    res.error = std::sqrt(static_cast<double>(A.error)*static_cast<double>(A.error) + static_cast<double>(B.error)*static_cast<double>(B.error));
    
    return std::move(res);
  }

template <typename ValueType>
dimensionful_Pk_component<ValueType> operator-(const loop_integral_output<ValueType>& A, const loop_integral_output<ValueType>& B)
  {
    dimensionful_Pk_component<ValueType> res;
    
    res.value = A.value - B.value;
    res.error = std::sqrt(static_cast<double>(A.error)*static_cast<double>(A.error) + static_cast<double>(B.error)*static_cast<double>(B.error));
    
    return std::move(res);
  }


// overload addition of loop integral and Pk components to give another Pk component
template <typename ValueType>
dimensionful_Pk_component<ValueType> operator+(const dimensionful_Pk_component<ValueType>& A, const loop_integral_output<ValueType>& B)
  {
    dimensionful_Pk_component<ValueType> res;
    
    res.value = A.value + B.value;
    res.error = std::sqrt(static_cast<double>(A.error)*static_cast<double>(A.error) + static_cast<double>(B.error)*static_cast<double>(B.error));
    
    return std::move(res);
  }

template <typename ValueType>
dimensionful_Pk_component<ValueType> operator+(const loop_integral_output<ValueType>& A, const dimensionful_Pk_component<ValueType>& B)
  {
    dimensionful_Pk_component<ValueType> res;
    
    res.value = A.value + B.value;
    res.error = std::sqrt(static_cast<double>(A.error)*static_cast<double>(A.error) + static_cast<double>(B.error)*static_cast<double>(B.error));
    
    return std::move(res);
  }


// overload subtraction of loop integral and Pk components to give another Pk component
template <typename ValueType>
dimensionful_Pk_component<ValueType> operator-(const dimensionful_Pk_component<ValueType>& A, const loop_integral_output<ValueType>& B)
  {
    dimensionful_Pk_component<ValueType> res;
    
    res.value = A.value - B.value;
    res.error = std::sqrt(static_cast<double>(A.error)*static_cast<double>(A.error) + static_cast<double>(B.error)*static_cast<double>(B.error));
    
    return std::move(res);
  }

template <typename ValueType>
dimensionful_Pk_component<ValueType> operator-(const loop_integral_output<ValueType>& A, const dimensionful_Pk_component<ValueType>& B)
  {
    dimensionful_Pk_component<ValueType> res;
    
    res.value = A.value - B.value;
    res.error = std::sqrt(static_cast<double>(A.error)*static_cast<double>(A.error) + static_cast<double>(B.error)*static_cast<double>(B.error));
    
    return std::move(res);
  }


// overload multiplication of dimensionless Pk components with 1/Mpc^3 dimensionful quantities
dimensionful_Pk_component<Mpc_units::inverse_energy3>
inline operator*(const Mpc_units::inverse_energy3& A, const dimensionful_Pk_component<double>& B)
  {
    dimensionful_Pk_component<Mpc_units::inverse_energy3> res;
    
    res.value = A * B.value;
    res.error = std::abs(A) * B.error;
    
    return std::move(res);
  }

dimensionful_Pk_component<Mpc_units::inverse_energy3>
inline operator*(const dimensionful_Pk_component<double>& A, const Mpc_units::inverse_energy3& B)
  {
    dimensionful_Pk_component<Mpc_units::inverse_energy3> res;
    
    res.value = B * A.value;
    res.error = std::abs(B) * A.error;
    
    return std::move(res);
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
      }
    
  };


class oneloop_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    oneloop_Pk(const k_token& kt, const UV_token& UVt, const IR_token& IRt, const z_token& zt,
               const dd_Pk& _dd, const rsd_dd_Pk& _rsd_mu0, const rsd_dd_Pk& _rsd_mu2,
               const rsd_dd_Pk& _rsd_mu4, const rsd_dd_Pk& _rsd_mu6, const rsd_dd_Pk& _rsd_mu8);
    
    //! empty constructor, used when receiving an MPI payload
    oneloop_Pk();
    
    //! destructor is default
    ~oneloop_Pk() = default;
    
    
    // INTERFACE
    
  public:
    
    //! get wavenumber token
    const k_token& get_k_token() const { return this->k; }
    
    //! get UV cutoff token
    const UV_token& get_UV_token() const { return this->UV_cutoff; }
    
    //! get IR cutoff token
    const IR_token& get_IR_token() const { return this->IR_cutoff; }
    
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
    
    //! UV cutoff token
    UV_token UV_cutoff;
    
    //! IR cutoff token
    IR_token IR_cutoff;
    
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
