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

    value_type value = value_type(0.0);
    unsigned int regions;
    unsigned int evaluations;
    double       error;

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
      }

  };


class dimless_kernel
  {

  public:

    typedef double value_type;

    value_type   value;
    unsigned int regions;
    unsigned int evaluations;
    double       error;

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
      }

  };


class loop_integral
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    loop_integral(const Mpc_units::energy& _k, const k_token& kt,
                  const Mpc_units::energy& UV, const UV_token& UVt,
                  const Mpc_units::energy& IR, const IR_token& IRt,
                  bool f, const inverse_energy3_kernel& _A, const inverse_energy3_kernel& _B,
                  const dimless_kernel& _D, const dimless_kernel& _E, const dimless_kernel& _F, const dimless_kernel& _G);

    //! destructor is default
    ~loop_integral() = default;


    // INTERFACE

  public:

    //! get failure state
    bool get_fail() const { return(this->fail); }

    //! get wavenumber token
    const k_token& get_k_token() const { return(this->k_tok); }

    //! get UV cutoff token
    const UV_token& get_UV_token() const { return(this->UV_tok); }

    //! get IR cutoff token
    const IR_token& get_IR_token() const { return(this->IR_tok); }

    //! get A-value
    const inverse_energy3_kernel& get_A() const { return(this->A); }

    //! get B-value
    const inverse_energy3_kernel& get_B() const { return(this->B); }

    //! get D-value
    const dimless_kernel& get_D() const { return(this->D); }

    //! get E-value
    const dimless_kernel& get_E() const { return(this->E); }

    //! get F-value
    const dimless_kernel& get_F() const { return(this->F); }

    //! get G-value
    const dimless_kernel& get_G() const { return(this->G); }


    // METADATA

  public:

    //! store integration time
    void set_integration_metadata(boost::timer::nanosecond_type t);

    //! get integration time
    boost::timer::nanosecond_type get_integration_time() const { return(this->integration_time); }


    // INTERNAL DATA

  private:

    // CONFIGURATION DATA

    //! failure state
    bool fail;

    //! wavenumber
    Mpc_units::energy k;

    //! wavenumber token
    k_token k_tok;

    //! UV cutoff
    Mpc_units::energy UV_cutoff;

    //! UV cutoff token
    UV_token UV_tok;

    //! IR cutoff
    Mpc_units::energy IR_cutoff;

    //! IR cutoff token
    IR_token IR_tok;


    // VALUE

    //! A-type integral P_13
    inverse_energy3_kernel A;

    //! B-type integral P_13
    inverse_energy3_kernel B;

    //! D-type integral P_13
    dimless_kernel D;

    //! E-type integral P_13
    dimless_kernel E;

    //! F-type integral P_13
    dimless_kernel F;

    //! G-type integral P_13
    dimless_kernel G;


    // METADATA

    //! time taken to perform integration
    boost::timer::nanosecond_type integration_time;


    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & k;
        ar & k_tok;
        ar & UV_cutoff;
        ar & UV_tok;
        ar & IR_cutoff;
        ar & IR_tok;
        ar & A;
        ar & B;
        ar & D;
        ar & E;
        ar & F;
        ar & G;
        ar & integration_time;
      }

  };


#endif //LSSEFT_LOOP_INTEGRAL_H
