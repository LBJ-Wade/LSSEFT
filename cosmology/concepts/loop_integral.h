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


class loop_integral
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    loop_integral(const Mpc_units::energy& _k, const k_token& kt,
                  const Mpc_units::energy& UV, const UV_token& UVt,
                  const Mpc_units::energy& IR, const IR_token& IRt,
                  const Mpc_units::inverse_energy3& _A, const Mpc_units::inverse_energy3& _B,
                  double _D, double _E, double _F, double _G);

    //! destructor is default
    ~loop_integral() = default;


    // INTERFACE

  public:

    //! get wavenumber token
    const k_token& get_k_token() const { return(this->k_tok); }

    //! get UV cutoff token
    const UV_token& get_UV_token() const { return(this->UV_tok); }

    //! get IR cutoff token
    const IR_token& get_IR_token() const { return(this->IR_tok); }

    //! get A-value
    const Mpc_units::inverse_energy3& get_A() const { return(this->A); }

    //! get B-value
    const Mpc_units::inverse_energy3& get_B() const { return(this->B); }

    //! get D-value
    double get_D() const { return(this->D); }

    //! get E-value
    double get_E() const { return(this->E); }

    //! get F-value
    double get_F() const { return(this->F); }

    //! get G-value
    double get_G() const { return(this->G); }


    // METADATA

  public:

    //! store integration time
    void set_integration_metadata(boost::timer::nanosecond_type t, size_t s);

    //! get integration time
    boost::timer::nanosecond_type get_integration_time() const { return(this->integration_time); }

    //! get number of steps used by integrator
    size_t get_integration_steps() const { return(this->steps); }


    // INTERNAL DATA

  private:

    // CONFIGURATION DATA

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
    Mpc_units::inverse_energy3 A;

    //! B-type integral P_13
    Mpc_units::inverse_energy3 B;

    //! D-type integral P_13
    double D;

    //! E-type integral P_13
    double E;

    //! F-type integral P_13
    double F;

    //! G-type integral P_13
    double G;


    // METADATA

    //! time taken to perform integration
    boost::timer::nanosecond_type integration_time;

    //! number of steps used by integrator
    size_t steps;


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
        ar & steps;
      }

  };


#endif //LSSEFT_LOOP_INTEGRAL_H
