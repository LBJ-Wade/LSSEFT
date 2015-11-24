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

    value_type                    value = value_type(0.0);
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


class loop_integral
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    loop_integral(const Mpc_units::energy& _k, const k_token& kt,
                  const Mpc_units::energy& UV, const UV_token& UVt,
                  const Mpc_units::energy& IR, const IR_token& IRt,
                  bool f, const inverse_energy3_kernel& _AA, const inverse_energy3_kernel& _AB, const inverse_energy3_kernel& _BB,
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

    //! get AA-value
    const inverse_energy3_kernel& get_AA() const { return(this->AA); }

    //! get AB-value
    const inverse_energy3_kernel& get_AB() const { return(this->AB); }

    //! get BB-value
    const inverse_energy3_kernel& get_BB() const { return(this->BB); }

    //! get D-value
    const dimless_kernel& get_D() const { return(this->D); }

    //! get E-value
    const dimless_kernel& get_E() const { return(this->E); }

    //! get F-value
    const dimless_kernel& get_F() const { return(this->F); }

    //! get G-value
    const dimless_kernel& get_G() const { return(this->G); }


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

    //! AA-type integral P_22
    inverse_energy3_kernel AA;

    //! AB-type integral P_22
    inverse_energy3_kernel AB;

    //! BB-type integral P_22
    inverse_energy3_kernel BB;

    //! D-type integral P_13
    dimless_kernel D;

    //! E-type integral P_13
    dimless_kernel E;

    //! F-type integral P_13
    dimless_kernel F;

    //! G-type integral P_13
    dimless_kernel G;


    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & fail;
        ar & k;
        ar & k_tok;
        ar & UV_cutoff;
        ar & UV_tok;
        ar & IR_cutoff;
        ar & IR_tok;
        ar & AA;
        ar & AB;
        ar & BB;
        ar & D;
        ar & E;
        ar & F;
        ar & G;
      }

  };


#endif //LSSEFT_LOOP_INTEGRAL_H
