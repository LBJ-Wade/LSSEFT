//
// Created by David Seery on 21/11/2015.
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

#ifndef LSSEFT_LOOP_INTEGRAL_H
#define LSSEFT_LOOP_INTEGRAL_H


#include "database/tokens.h"
#include "units/Mpc_units.h"

#include "boost/timer/timer.hpp"
#include "boost/serialization/serialization.hpp"


template <typename ValueType>
class loop_integral_result
  {

  public:
    
    typedef ValueType value_type;
    
    //! constructor initializes zero values, which should be overwritten later
    loop_integral_result()
      : value(value_type(0.0)),
        error(value_type(0.0)),
        regions(0),
        evaluations(0),
        time(0)
      {
      }
    
    //! destructor is default
    ~loop_integral_result() = default;

    
    // DATA
  
  public:
    
    value_type                    value;
    value_type                    error;
    
    unsigned int                  regions;
    unsigned int                  evaluations;
    boost::timer::nanosecond_type time;
  
  private:
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & value;
        ar & regions;
        ar & evaluations;
        ar & error;
        ar & time;
      }

  };


template <typename ValueType>
class loop_integral_output
  {
  
  public:
    
    typedef ValueType value_type;

    //! constructor initializes result values to zero
    //! these should be overwritten later
    loop_integral_output()
      : raw(),
        nowiggle()
      {
      }
    
    //! destructor is default
    ~loop_integral_output() = default;
    
    
    // ACCESSORS
    
  public:
    
    //! get raw value
    loop_integral_result<ValueType>& get_raw() { return this->raw; }
    const loop_integral_result<ValueType>& get_raw() const { return this->raw; }
    
    //! get no-wiggle value
    loop_integral_result<ValueType>& get_nowiggle() { return this->nowiggle; }
    const loop_integral_result<ValueType>& get_nowiggle() const { return this->nowiggle; }
    
    
    // INTERNAL DATA
    
  private:
    
    //! raw result
    loop_integral_result<ValueType> raw;
    
    //! wiggle result
    loop_integral_result<ValueType> nowiggle;
    
  
  private:
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & raw;
        ar & nowiggle;
      }
    
  };


typedef loop_integral_output<Mpc_units::inverse_energy3> inverse_energy3_integral;
typedef loop_integral_output<double>                     dimless_integral;


template <typename DimensionfulType>
inline DimensionfulType dimensionful_unit();


template <>
inline double dimensionful_unit<double>()
  {
    return 1.0;
  }

template <>
inline Mpc_units::inverse_energy dimensionful_unit<Mpc_units::inverse_energy>()
  {
    return Mpc_units::Mpc;
  }

template <>
inline Mpc_units::inverse_energy2 dimensionful_unit<Mpc_units::inverse_energy2>()
  {
    return Mpc_units::Mpc2;
  }

template <>
inline Mpc_units::inverse_energy3 dimensionful_unit<Mpc_units::inverse_energy3>()
  {
    return Mpc_units::Mpc3;
  }


class delta_22_integrals
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! value constructor
    delta_22_integrals(const inverse_energy3_integral& _AA, const inverse_energy3_integral& _AB,
                       const inverse_energy3_integral& _BB);
    
    //! empty constructor for use when overwriting with MPI payloads
    delta_22_integrals();
    
    //! destructor is default
    ~delta_22_integrals() = default;
    
    
    // INTERFACE
  
  public:
    
    //! get failure state
    bool get_fail() const { return this->fail; }
    
    //! set failed flag
    void mark_failed() { this->fail = true; }
    
    
    //! get AA-value
    inverse_energy3_integral& get_AA() { return (this->AA); }
    const inverse_energy3_integral& get_AA() const { return (this->AA); }
    
    //! get AB-value
    inverse_energy3_integral& get_AB() { return (this->AB); }
    const inverse_energy3_integral& get_AB() const { return (this->AB); }
    
    //! get BB-value
    inverse_energy3_integral& get_BB() { return (this->BB); }
    const inverse_energy3_integral& get_BB() const { return (this->BB); }
    
    
    // INTERNAL DATA
  
  private:
    
    //! failure state
    bool fail;
    
    //! AA-type integral P_22
    inverse_energy3_integral AA;
    
    //! AB-type integral P_22
    inverse_energy3_integral AB;
    
    //! BB-type integral P_22
    inverse_energy3_integral BB;
    
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & fail;
        ar & AA;
        ar & AB;
        ar & BB;
      }
    
  };


class delta_13_integrals
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    delta_13_integrals(const dimless_integral& _D, const dimless_integral& _E, const dimless_integral& _F, const dimless_integral& _G,
                       const dimless_integral& _J1, const dimless_integral& _J2);
    
    //! empty constructor for use when overwriting with MPI payloads
    delta_13_integrals();
    
    //! destructor
    ~delta_13_integrals() = default;
    
    
    // INTERFACE
    
  public:
    
    //! get failure state
    bool get_fail() const { return this->fail; }
    
    //! set failed flag
    void mark_failed() { this->fail = true; }
    
    
    //! get D-value
    dimless_integral& get_D() { return(this->D); }
    const dimless_integral& get_D() const { return(this->D); }
    
    //! get E-value
    dimless_integral& get_E() { return(this->E); }
    const dimless_integral& get_E() const { return(this->E); }
    
    //! get F-value
    dimless_integral& get_F() { return(this->F); }
    const dimless_integral& get_F() const { return(this->F); }
    
    //! get G-value
    dimless_integral& get_G() { return(this->G); }
    const dimless_integral& get_G() const { return(this->G); }
    
    //! get J1-value
    dimless_integral& get_J1() { return(this->J1); }
    const dimless_integral& get_J1() const { return(this->J1); }
    
    //! get J1-value
    dimless_integral& get_J2() { return(this->J2); }
    const dimless_integral& get_J2() const { return(this->J2); }

    
    // INTERNAL DATA
    
  private:
    
    //! failure state
    bool fail;
    
    //! D-type integral P_13
    dimless_integral D;
    
    //! E-type integral P_13
    dimless_integral E;
    
    //! F-type integral P_13
    dimless_integral F;
    
    //! G-type integral P_13
    dimless_integral G;
    
    //! J1-type kernel P_13
    dimless_integral J1;
    
    //! J2-type kernel P_13
    dimless_integral J2;
    
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & fail;
        ar & D;
        ar & E;
        ar & F;
        ar & G;
        ar & J1;
        ar & J2;
      }
    
  };


class rsd_22_integrals
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    rsd_22_integrals(const inverse_energy3_integral& _A1, const inverse_energy3_integral& _A2,
                     const inverse_energy3_integral& _A3, const inverse_energy3_integral& _A4,
                     const inverse_energy3_integral& _A5, const inverse_energy3_integral& _B2,
                     const inverse_energy3_integral& _B3, const inverse_energy3_integral& _B6,
                     const inverse_energy3_integral& _B8, const inverse_energy3_integral& _B9,
                     const inverse_energy3_integral& _C1, const inverse_energy3_integral& _C2,
                     const inverse_energy3_integral& _C4, const inverse_energy3_integral& _D1);
    
    //! empty constructor for use when overwriting with MPI payloads
    rsd_22_integrals();
    
    //! destructor is default
    ~rsd_22_integrals() = default;
    
    
    // INTERFACE
    
  public:
    
    //! get failure state
    bool get_fail() const { return this->fail; }
    
    //! set failed flag
    void mark_failed() { this->fail = true; }
    
    
    //! get A1-value
    inverse_energy3_integral& get_A1() { return(this->A1); }
    const inverse_energy3_integral& get_A1() const { return(this->A1); }
    
    //! get A2-value
    inverse_energy3_integral& get_A2() { return(this->A2); }
    const inverse_energy3_integral& get_A2() const { return(this->A2); }
    
    //! get A3-value
    inverse_energy3_integral& get_A3() { return(this->A3); }
    const inverse_energy3_integral& get_A3() const { return(this->A3); }
    
    //! get A4-value
    inverse_energy3_integral& get_A4() { return(this->A4); }
    const inverse_energy3_integral& get_A4() const { return(this->A4); }
    
    //! get A5-value
    inverse_energy3_integral& get_A5() { return(this->A5); }
    const inverse_energy3_integral& get_A5() const { return(this->A5); }
    
    //! get B2-value
    inverse_energy3_integral& get_B2() { return(this->B2); }
    const inverse_energy3_integral& get_B2() const { return(this->B2); }
    
    //! get B3-value
    inverse_energy3_integral& get_B3() { return(this->B3); }
    const inverse_energy3_integral& get_B3() const { return(this->B3); }
    
    //! get B6-value
    inverse_energy3_integral& get_B6() { return(this->B6); }
    const inverse_energy3_integral& get_B6() const { return(this->B6); }
    
    //! get B8-value
    inverse_energy3_integral& get_B8() { return(this->B8); }
    const inverse_energy3_integral& get_B8() const { return(this->B8); }
    
    //! get B9-value
    inverse_energy3_integral& get_B9() { return(this->B9); }
    const inverse_energy3_integral& get_B9() const { return(this->B9); }
    
    //! get C1-value
    inverse_energy3_integral& get_C1() { return(this->C1); }
    const inverse_energy3_integral& get_C1() const { return(this->C1); }
    
    //! get C2-value
    inverse_energy3_integral& get_C2() { return(this->C2); }
    const inverse_energy3_integral& get_C2() const { return(this->C2); }
    
    //! get C4-value
    inverse_energy3_integral& get_C4() { return(this->C4); }
    const inverse_energy3_integral& get_C4() const { return(this->C4); }
    
    //! get D1-value
    inverse_energy3_integral& get_D1() { return(this->D1); }
    const inverse_energy3_integral& get_D1() const { return(this->D1); }
    
    
    // INTERNAL DATA
    
    //! failure state
    bool fail;
    
    
    //! A1-type integral P_22
    inverse_energy3_integral A1;
    
    //! A2-type integral P_22
    inverse_energy3_integral A2;
    
    //! A3-type integral P_22
    inverse_energy3_integral A3;
    
    //! A4-type integral P_22
    inverse_energy3_integral A4;
    
    //! A5-type integral P_22
    inverse_energy3_integral A5;
    
    //! B2-type integral P_22
    inverse_energy3_integral B2;
    
    //! B3-type integral P_22
    inverse_energy3_integral B3;
    
    //! B6-type integral P_22
    inverse_energy3_integral B6;
    
    //! B8-type integral P_22
    inverse_energy3_integral B8;
    
    //! B9-type integral P_22
    inverse_energy3_integral B9;
    
    //! C1-type integral P_22
    inverse_energy3_integral C1;
    
    //! C2-type integral P_22
    inverse_energy3_integral C2;
    
    //! C4-type integral P_22
    inverse_energy3_integral C4;
    
    //! D1-type integral P_22
    inverse_energy3_integral D1;

    
  private:
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & fail;
        ar & A1;
        ar & A2;
        ar & A3;
        ar & A4;
        ar & A5;
        ar & B2;
        ar & B3;
        ar & B6;
        ar & B8;
        ar & B9;
        ar & C1;
        ar & C2;
        ar & C4;
        ar & D1;
      }
    
  };


class rsd_13_integrals
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! value constructor
    rsd_13_integrals(const dimless_integral& _a, const dimless_integral& _b, const dimless_integral& _c,
                     const dimless_integral& _d, const dimless_integral& _e, const dimless_integral& _f,
                     const dimless_integral& _g);
    
    //! empty constructor for use when overwriting with MPI payloads
    rsd_13_integrals();
    
    //! destructor is default
    ~rsd_13_integrals() = default;
    
    
    // INTERFACE
  
  public:
    
    //! get failure state
    bool get_fail() const { return this->fail; }
    
    //! set failed flag
    void mark_failed() { this->fail = true; }
    
    
    //! get a-value
    dimless_integral& get_a() { return(this->a); }
    const dimless_integral& get_a() const { return(this->a); }
    
    //! get b-value
    dimless_integral& get_b() { return(this->b); }
    const dimless_integral& get_b() const { return(this->b); }
    
    //! get c-value
    dimless_integral& get_c() { return(this->c); }
    const dimless_integral& get_c() const { return(this->c); }
    
    //! get d-value
    dimless_integral& get_d() { return(this->d); }
    const dimless_integral& get_d() const { return(this->d); }
    
    //! get e-value
    dimless_integral& get_e() { return(this->e); }
    const dimless_integral& get_e() const { return(this->e); }
    
    //! get f-value
    dimless_integral& get_f() { return(this->f); }
    const dimless_integral& get_f() const { return(this->f); }
    
    //! get g-value
    dimless_integral& get_g() { return(this->g); }
    const dimless_integral& get_g() const { return(this->g); }
    
    
    // INTERNAL DATA
    
  private:
    
    //! failure state
    bool fail;
    
    
    //! a-type integral P_13
    dimless_integral a;
    
    //! b-type integral P_13
    dimless_integral b;
    
    //! c-type integral P_13
    dimless_integral c;
    
    //! d-type integral P_13
    dimless_integral d;
    
    //! e-type integral P_13
    dimless_integral e;
    
    //! f-type integral P_13
    dimless_integral f;
    
    //! g-type integral P_13
    dimless_integral g;
  
  
  private:
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & fail;
        ar & a;
        ar & b;
        ar & c;
        ar & d;
        ar & e;
        ar & f;
        ar & g;
      }
    
  };


class loop_integral
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! value constructor
    loop_integral(const k_token& kt, const loop_integral_params_token& pt, const linear_Pk_token& Pt,
                  const UV_cutoff_token& UVt, const IR_cutoff_token& IRt, const delta_22_integrals& d22,
                  const delta_13_integrals& d13, const rsd_22_integrals& r22, const rsd_13_integrals& r13);
    
    //! empty constructor, used for constructing an empty container to be overwritten by an MPI payload
    loop_integral();

    //! destructor is default
    ~loop_integral() = default;


    // INTERFACE

  public:
    
    //! get parameters token
    const loop_integral_params_token& get_params_token() const { return this->params; }

    //! get wavenumber token
    const k_token& get_k_token() const { return this->k; }
    
    //! get linear power spectrum token
    const linear_Pk_token& get_Pk_token() const { return this->Pk_lin; }

    //! get UV cutoff token
    const UV_cutoff_token& get_UV_token() const { return this->UV_cutoff; }

    //! get IR cutoff token
    const IR_cutoff_token& get_IR_token() const { return this->IR_cutoff; }
    
    
    
    //! get delta-22 container
    const delta_22_integrals& get_delta22() const { return this->delta22; }
    
    //! get delta-13 container
    const delta_13_integrals& get_delta13() const { return this->delta13; }
    
    //! get rsd-22 container
    const rsd_22_integrals& get_rsd22() const { return this->rsd22; }
    
    //! get rsd-13 container
    const rsd_13_integrals& get_rsd13() const { return this->rsd13; }


    // INTERNAL DATA

  private:

    // CONFIGURATION DATA
    
    //! parameters token
    loop_integral_params_token params;

    //! wavenumber token
    k_token k;
    
    //! linear power spectrum token
    linear_Pk_token Pk_lin;

    //! UV cutoff token
    UV_cutoff_token UV_cutoff;

    //! IR cutoff token
    IR_cutoff_token IR_cutoff;


    // VALUES
    
    //! delta-13 integrals
    delta_13_integrals delta13;
    
    //! delta-22 integrals
    delta_22_integrals delta22;
    
    //! rsd-13 integrals
    rsd_13_integrals rsd13;
    
    //! rsd-22 integrals
    rsd_22_integrals rsd22;


    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & k;
        ar & UV_cutoff;
        ar & IR_cutoff;
        ar & delta13;
        ar & delta22;
        ar & rsd13;
        ar & rsd22;
      }

  };


#endif //LSSEFT_LOOP_INTEGRAL_H
