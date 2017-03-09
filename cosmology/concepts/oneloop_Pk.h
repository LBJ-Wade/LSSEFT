//
// Created by David Seery on 14/11/2016.
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

#ifndef LSSEFT_ONE_LOOP_PK_H
#define LSSEFT_ONE_LOOP_PK_H


#include "Pk_value.h"

#include "database/tokens.h"
#include "units/Mpc_units.h"

#include "power_spectrum.h"

#include "boost/timer/timer.hpp"
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/optional.hpp"


template <typename ValueType>
class raw_wiggle_Pk_component
  {
    
    // TYPES
    
  public:
    
    typedef ValueType value_type;
    typedef Pk_value_group<ValueType> container_type;
    
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor; error is zero if not specified which allows
    //! assignment-on-construction to a fixed number
    raw_wiggle_Pk_component(container_type r, container_type nw)
      : raw(std::move(r)),
        nowiggle(std::move(nw))
      {
      }
    
    //! empty constructor
    raw_wiggle_Pk_component()
      : raw(),
        nowiggle()
      {
      }

    
    // ACCESSORS
    
  public:
    
    //! get raw data
    const container_type& get_raw() const { return this->raw; }
    
    //! get wiggle data
    const container_type get_wiggle() const { return this->raw - this->nowiggle; }
    
    //! get no-wiggle data
    const container_type& get_nowiggle() const { return this->nowiggle; }
    
    //! set raw data
    raw_wiggle_Pk_component<ValueType>& set_raw(Pk_value_group<ValueType> r) { this->raw = std::move(r); return *this; }
    
    //! set wiggle data
    raw_wiggle_Pk_component<ValueType>& set_nowiggle(Pk_value_group<ValueType> w) { this->nowiggle = std::move(w); return *this; }
    
    
    // DATA
    
  private:
    
    //! raw data
    container_type raw;
    
    //! wiggle data
    container_type nowiggle;
  
  
  private:
    
    // enable boost::serialization support
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & raw;
        ar & nowiggle;
      }
    
  };


// define convenience types for basic power spectrum (~ 1/k^3) and k^2 * basic power spectrum (~ 1/k)
typedef raw_wiggle_Pk_component<Mpc_units::inverse_energy3> Pk_value;
typedef raw_wiggle_Pk_component<Mpc_units::inverse_energy>  k2_Pk_value;


// build Pk_value from a initial_filtered_Pk
template <typename PkContainer>
Pk_value build_Pk_value(const Mpc_units::energy& k, const PkContainer& Pk)
  {
    Pk_value_group<Mpc_units::inverse_energy3> Ptree_raw(Pk.Pk_raw(k));
    Pk_value_group<Mpc_units::inverse_energy3> Ptree_nw(Pk.Pk_nowiggle(k));

    return Pk_value(Ptree_raw, Ptree_nw);
  }


//! overload + and - so that power spectrum and counterterm values can be added
template <typename ValueType>
raw_wiggle_Pk_component<ValueType> operator+(const raw_wiggle_Pk_component<ValueType>& A, const raw_wiggle_Pk_component<ValueType>& B)
  {
    return raw_wiggle_Pk_component<ValueType>(A.get_raw() + B.get_raw(), A.get_nowiggle() + B.get_nowiggle());
  }

template <typename ValueType>
raw_wiggle_Pk_component<ValueType> operator-(const raw_wiggle_Pk_component<ValueType>& A, const raw_wiggle_Pk_component<ValueType>& B)
  {
    return raw_wiggle_Pk_component<ValueType>(A.get_raw() - B.get_raw(), A.get_nowiggle() - B.get_nowiggle());
  }


//! overload scalar multiplication
template <typename ScalarType, typename PkType, typename std::enable_if_t< !std::is_integral<ScalarType>::value >* = nullptr>
auto operator*(const ScalarType& S, const raw_wiggle_Pk_component<PkType>& P) -> raw_wiggle_Pk_component<decltype(S*P.get_raw().get_value())>
  {
    return raw_wiggle_Pk_component<decltype(S*P.get_raw().get_value())>(S * P.get_raw(), S * P.get_nowiggle());
  }

template <typename ScalarType, typename PkType, typename std::enable_if_t< !std::is_integral<ScalarType>::value >* = nullptr>
auto operator*(const raw_wiggle_Pk_component<PkType>& P, const ScalarType& S) -> raw_wiggle_Pk_component<decltype(S*P.get_raw().get_value())>
  {
    return raw_wiggle_Pk_component<decltype(S*P.get_raw().get_value())>(S * P.get_raw(), S * P.get_nowiggle());
  }


//! overload scalar division
template <typename ScalarType, typename PkType, typename std::enable_if_t< !std::is_integral<ScalarType>::value >* = nullptr>
auto operator/(const raw_wiggle_Pk_component<PkType>& P, const ScalarType& S) -> raw_wiggle_Pk_component<decltype(P.get_raw().get_value()/S)>
  {
    return raw_wiggle_Pk_component<decltype(P.get_raw().get_value()/S)>(P.get_raw() / S, P.get_nowiggle() / S);
  }



// AUTOMATIC CONVERSION OF LOOP INTEGRAL KERNELS TO PK COMPONENTS


// multiplication of any loop integral output by a scalar converts to a raw_wiggle_Pk_component<>
template <typename ScalarType, typename LoopType, typename std::enable_if_t< !std::is_integral<ScalarType>::value >* = nullptr>
auto operator*(const ScalarType& S, const loop_integral_output<LoopType>& P) -> raw_wiggle_Pk_component<decltype(S*P.get_raw().value)>
  {
    return raw_wiggle_Pk_component<decltype(S*P.get_raw().value)>(S * P.get_raw(), S * P.get_nowiggle());
  }

template <typename ScalarType, typename LoopType, typename std::enable_if_t< !std::is_integral<ScalarType>::value >* = nullptr>
auto operator*(const loop_integral_output<LoopType>& P, const ScalarType& S) -> raw_wiggle_Pk_component<decltype(S*P.get_raw().value)>
  {
    return raw_wiggle_Pk_component<decltype(S*P.get_raw().value)>(S * P.get_raw(), S * P.get_nowiggle());
  }


// multiplication of loop integral output by raw_wiggle_Pk_component<>
template <typename LoopType, typename PkType>
auto operator*(const raw_wiggle_Pk_component<PkType>& P, const loop_integral_output<LoopType>& L) -> raw_wiggle_Pk_component<decltype(P.get_raw().get_value()*L.get_raw().value)>
  {
    return raw_wiggle_Pk_component<decltype(P.get_raw().get_value()*L.get_raw().value)>(P.get_raw() * L.get_raw(), P.get_nowiggle() * L.get_nowiggle());
  }

template <typename LoopType, typename PkType>
auto operator*(const loop_integral_output<LoopType>& L, const raw_wiggle_Pk_component<PkType>& P) -> raw_wiggle_Pk_component<decltype(P.get_raw().get_value()*L.get_raw().value)>
  {
    return raw_wiggle_Pk_component<decltype(P.get_raw().get_value()*L.get_raw().value)>(P.get_raw() * L.get_raw(), P.get_nowiggle() * L.get_nowiggle());
  }

  
// addition and subtraction of loop integral outputs gives raw_wiggle_Pk_component<>
template <typename ValueType>
raw_wiggle_Pk_component<ValueType> operator+(const loop_integral_output<ValueType>& A, const loop_integral_output<ValueType>& B)
  {
    return raw_wiggle_Pk_component<ValueType>(A.get_raw() + B.get_raw(), A.get_nowiggle() + B.get_nowiggle());
  }

template <typename ValueType>
raw_wiggle_Pk_component<ValueType> operator-(const loop_integral_output<ValueType>& A, const loop_integral_output<ValueType>& B)
  {
    return raw_wiggle_Pk_component<ValueType>(A.get_raw() - B.get_raw(), A.get_nowiggle() - B.get_nowiggle());
  }


// addition and subtraction of loop integral outputs with raw_wiggle_Pk_component<>s gives raw_wiggle_Pk_component<>
template <typename ValueType>
raw_wiggle_Pk_component<ValueType> operator+(const raw_wiggle_Pk_component<ValueType>& A, const loop_integral_output<ValueType>& B)
  {
    return raw_wiggle_Pk_component<ValueType>(A.get_raw() + B.get_raw(), A.get_nowiggle() + B.get_nowiggle());
  }

template <typename ValueType>
raw_wiggle_Pk_component<ValueType> operator+(const loop_integral_output<ValueType>& A, const raw_wiggle_Pk_component<ValueType>& B)
  {
    return raw_wiggle_Pk_component<ValueType>(A.get_raw() + B.get_raw(), A.get_nowiggle() + B.get_nowiggle());
  }

template <typename ValueType>
raw_wiggle_Pk_component<ValueType> operator-(const raw_wiggle_Pk_component<ValueType>& A, const loop_integral_output<ValueType>& B)
  {
    return raw_wiggle_Pk_component<ValueType>(A.get_raw() - B.get_raw(), A.get_nowiggle() - B.get_nowiggle());
  }

template <typename ValueType>
raw_wiggle_Pk_component<ValueType> operator-(const loop_integral_output<ValueType>& A, const raw_wiggle_Pk_component<ValueType>& B)
  {
    return raw_wiggle_Pk_component<ValueType>(A.get_raw() - B.get_raw(), A.get_nowiggle() - B.get_nowiggle());
  }


// multiplication of raw_wiggle_Pk_component<>s
template <typename PkTypeA, typename PkTypeB>
auto operator*(const raw_wiggle_Pk_component<PkTypeA> PA, const raw_wiggle_Pk_component<PkTypeB>& PB) -> raw_wiggle_Pk_component<decltype(PA.get_raw().get_value()*PB.get_raw().get_value())>
  {
    return raw_wiggle_Pk_component<decltype(PA.get_raw().get_value()*PB.get_raw().get_value())>(PA.get_raw() * PB.get_raw(), PA.get_nowiggle() * PB.get_nowiggle());
  }


// template class for a density-density 1-loop power spectrum (templated so it can be re-used for the resummed version, too)
template <typename PkValueType, typename k2PkValueType>
class generic_dd_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    generic_dd_Pk(const PkValueType& _Pt, const PkValueType& _P13, const PkValueType& _P22, const k2PkValueType& _Z2d);
    
    //! value constructor
    generic_dd_Pk(const PkValueType& _Pt, const PkValueType& _P13, const PkValueType& _P22, const PkValueType& _P1loopSPT, const k2PkValueType& _Z2d);
    
    //! empty constructor for use when overwriting with MPI payloads
    generic_dd_Pk();
    
    //! destructor is default
    ~generic_dd_Pk() = default;
    
    
    // INTERFACE
    
  public:
    
    //! get tree value
    PkValueType& get_tree() { return this->Ptree; }
    const PkValueType& get_tree() const { return this->Ptree; }
    
    //! get 13 value
    PkValueType& get_13() { return this->P13; }
    const PkValueType& get_13() const { return this->P13; }
    
    //! get 22 value
    PkValueType& get_22() { return this->P22; }
    const PkValueType& get_22() const { return this->P22; }
    
    //! get total SPT power spectrum = 13 + 22
    PkValueType& get_1loop_SPT() { return this->P1loopSPT; }
    const PkValueType& get_1loop_SPT() const { return this->P1loopSPT; }
    
    
    // COUNTERTERMS
    
    //! get EFT counterterm
    k2PkValueType& get_Z2_delta() { return this->Z2_delta; }
    const k2PkValueType& get_Z2_delta() const { return this->Z2_delta; }
    
    
    // INTERNAL DATA
    
  private:
    
    //! tree power spectrum
    PkValueType Ptree;
    
    //! 13 terms
    PkValueType P13;
    
    //! 22 terms
    PkValueType P22;
    
    //! total 1-loop SPT value
    PkValueType P1loopSPT;
    
    //! coefficient of the counterterm Z2_delta
    k2PkValueType Z2_delta;
    
    
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


template <typename PkValueType, typename k2PkValueType>
generic_dd_Pk<PkValueType, k2PkValueType>::generic_dd_Pk(const PkValueType& _Pt, const PkValueType& _P13, const PkValueType& _P22, const k2PkValueType& _Z2d)
  : Ptree(_Pt),
    P13(_P13),
    P22(_P22),
    P1loopSPT(_Pt + _P13 + _P22),
    Z2_delta(_Z2d)
  {
  }


template <typename PkValueType, typename k2PkValueType>
generic_dd_Pk<PkValueType, k2PkValueType>::generic_dd_Pk(const PkValueType& _Pt, const PkValueType& _P13,
                                                         const PkValueType& _P22, const PkValueType& _P1loopSPT,
                                                         const k2PkValueType& _Z2d)
  : Ptree(_Pt),
    P13(_P13),
    P22(_P22),
    P1loopSPT(_P1loopSPT),
    Z2_delta(_Z2d)
  {
  }


template <typename PkValueType, typename k2PkValueType>
generic_dd_Pk<PkValueType, k2PkValueType>::generic_dd_Pk()
  : Ptree(),
    P13(),
    P22(),
    P1loopSPT(),
    Z2_delta()
  {
  }


// the general one-loop delta-delta power spectrum uses raw/nowiggle parts
typedef generic_dd_Pk<Pk_value, k2_Pk_value> dd_Pk;


class rsd_dd_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! value constructor
    rsd_dd_Pk(const Pk_value& _Pt, const Pk_value& _P13, const Pk_value& _P22,
              const k2_Pk_value& _Z2d, const Pk_value& _Z0v, const k2_Pk_value& _Z2v,
              const Pk_value& _Z0vd, const k2_Pk_value& _Z2vd,
              const k2_Pk_value& _Z2vv, const k2_Pk_value& _Z2vvd, const k2_Pk_value& _Z2vvv,
              const k2_Pk_value& _Z2_total);
    
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
    
    //± get total Z2_total counterterm
    k2_Pk_value& get_Z2_total() { return this->Z2_total; }
    const k2_Pk_value& get_Z2_total() const { return this->Z2_total; }
    
    
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
    
    //! coefficient of the total counterterm for this power of mu
    k2_Pk_value Z2_total;
    
    
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
        ar & Z2_total;
      }
    
  };


class oneloop_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    oneloop_Pk(const k_token& kt, const linear_Pk_token& Pkt_i, const boost::optional<linear_Pk_token>& Pkt_f,
               const IR_cutoff_token& IRt, const UV_cutoff_token& UVt, const z_token& zt,
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
    
    //! get initial power spectrum token
    const linear_Pk_token& get_init_Pk_token() const { return this->init_Pk; }
    
    //! get final power spectrum token, if provided
    const boost::optional<linear_Pk_token>& get_final_Pk_token() const { return this->final_Pk; }
    
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
    
    //! initial power spectrum token
    linear_Pk_token init_Pk;
    
    //! final power spectrum token
    boost::optional<linear_Pk_token> final_Pk;
    
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
        ar & init_Pk;
        ar & final_Pk;
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
