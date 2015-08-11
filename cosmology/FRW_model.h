//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_FRW_MODEL_H
#define LSSEFT_FRW_MODEL_H


#include "Planck_defaults.h"


class FRW_model
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    //! cosmological parameters refer to their values today; if none are provided, then we use
    //! Planck values as sensible defaults
    FRW_model(double om= Planck2013::omega_m, double occ= Planck2013::omega_cc, double h_= Planck2013::h, double tc= Planck2013::T_CMB.val);

    //! destructor is default
    ~FRW_model() = default;


    // INTERROGATE COSMOLOGICAL PARAMETERS

  public:

    //! get value of omega_m
    double get_omega_m() const { return(this->omega_m); }

    //! get value of omega_cc
    double get_omega_cc() const { return(this->omega_cc); }

    //! get value of h
    double get_h() const { return(this->h); }

    //! get value of T_CMB
    double get_T_CMB() const { return(this->T_CMB); }


    // INTERNAL DATA

  private:

    //! omega_m, measured today
    double omega_m;

    //! omega_Lambda or omega_cc, measured today
    double omega_cc;

    //! reduced Hubble constant, ie. quoted in units of 100 km/s/Mpc, measured today
    double h;

    //! CMB temperature T_CMB, measured today
    double T_CMB;

  };


#endif //LSSEFT_FRW_MODEL_H