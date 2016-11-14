//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_LOOP_INTEGRAL_H
#define LSSEFT_LOOP_INTEGRAL_H


#include "database/tokens.h"
#include "units/Mpc_units.h"

#include "boost/timer/timer.hpp"
#include "boost/serialization/serialization.hpp"


class inverse_energy3_kernel
  {

  public:
    
    typedef Mpc_units::inverse_energy3 value_type;
    
    inverse_energy3_kernel()
      : value(value_type(0.0)),
        regions(0),
        evaluations(0),
        error(0.0),
        time(0)
      {
      }

    value_type                    value;
    unsigned int                  regions;
    unsigned int                  evaluations;
    double                        error;
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


class dimless_kernel
  {

  public:

    typedef double value_type;
    
    dimless_kernel()
      : value(0.0),
        regions(0),
        evaluations(0),
        error(0.0),
        time(0)
      {
      }

    value_type                    value;
    unsigned int                  regions;
    unsigned int                  evaluations;
    double                        error;
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


template <typename DimensionfulType>
inline DimensionfulType dimensionful_unit();


template <>
inline double dimensionful_unit<double>()
  {
    return 1.0;
  }

template <>
inline Mpc_units::inverse_energy3 dimensionful_unit<Mpc_units::inverse_energy3>()
  {
    return Mpc_units::Mpc3;
  }

template <>
inline Mpc_units::inverse_energy dimensionful_unit<Mpc_units::inverse_energy>()
  {
    return Mpc_units::Mpc;
  }


class delta_22_integrals
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! value constructor
    delta_22_integrals(const inverse_energy3_kernel& _AA, const inverse_energy3_kernel& _AB,
                       const inverse_energy3_kernel& _BB);
    
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
    inverse_energy3_kernel& get_AA() { return (this->AA); }
    const inverse_energy3_kernel& get_AA() const { return (this->AA); }
    
    //! get AB-value
    inverse_energy3_kernel& get_AB() { return (this->AB); }
    const inverse_energy3_kernel& get_AB() const { return (this->AB); }
    
    //! get BB-value
    inverse_energy3_kernel& get_BB() { return (this->BB); }
    const inverse_energy3_kernel& get_BB() const { return (this->BB); }
    
    
    // INTERNAL DATA
  
  private:
    
    //! failure state
    bool fail;
    
    //! AA-type integral P_22
    inverse_energy3_kernel AA;
    
    //! AB-type integral P_22
    inverse_energy3_kernel AB;
    
    //! BB-type integral P_22
    inverse_energy3_kernel BB;
    
    
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
    delta_13_integrals(const dimless_kernel& _D, const dimless_kernel& _E, const dimless_kernel& _F, const dimless_kernel& _G,
                       const dimless_kernel& _J1, const dimless_kernel& _J2);
    
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
    dimless_kernel& get_D() { return(this->D); }
    const dimless_kernel& get_D() const { return(this->D); }
    
    //! get E-value
    dimless_kernel& get_E() { return(this->E); }
    const dimless_kernel& get_E() const { return(this->E); }
    
    //! get F-value
    dimless_kernel& get_F() { return(this->F); }
    const dimless_kernel& get_F() const { return(this->F); }
    
    //! get G-value
    dimless_kernel& get_G() { return(this->G); }
    const dimless_kernel& get_G() const { return(this->G); }
    
    //! get J1-value
    dimless_kernel& get_J1() { return(this->J1); }
    const dimless_kernel& get_J1() const { return(this->J1); }
    
    //! get J1-value
    dimless_kernel& get_J2() { return(this->J2); }
    const dimless_kernel& get_J2() const { return(this->J2); }

    
    // INTERNAL DATA
    
  private:
    
    //! failure state
    bool fail;
    
    //! D-type integral P_13
    dimless_kernel D;
    
    //! E-type integral P_13
    dimless_kernel E;
    
    //! F-type integral P_13
    dimless_kernel F;
    
    //! G-type integral P_13
    dimless_kernel G;
    
    //! J1-type kernel P_13
    dimless_kernel J1;
    
    //! J2-type kernel P_13
    dimless_kernel J2;
    
    
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
//    rsd_22_integrals() = default;
    
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
    
    
    // INTERNAL DATA
    
    //! failure state
    bool fail;
  
    
  private:
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & fail;
      }
    
  };


class rsd_13_integrals
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! value constructor
//    rsd_13_integrals() = default;
    
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
    
    
    // INTERNAL DATA
    
    //! failure state
    bool fail;
  
    
  private:
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & fail;
      }
    
  };


class loop_integral
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! value constructor
    loop_integral(const k_token& kt, const UV_token& UVt, const IR_token& IRt,
                  const delta_22_integrals& d22, const delta_13_integrals& d13,
                  const rsd_22_integrals& r22, const rsd_13_integrals& r13);
    
    //! empty constructor, used for constructing an empty container to be overwritten by an MPI payload
    loop_integral();

    //! destructor is default
    ~loop_integral() = default;


    // INTERFACE

  public:

    //! get wavenumber token
    const k_token& get_k_token() const { return this->k; }

    //! get UV cutoff token
    const UV_token& get_UV_token() const { return this->UV_cutoff; }

    //! get IR cutoff token
    const IR_token& get_IR_token() const { return this->IR_cutoff; }
    
    
    
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

    //! wavenumber token
    k_token k;

    //! UV cutoff token
    UV_token UV_cutoff;

    //! IR cutoff token
    IR_token IR_cutoff;


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
