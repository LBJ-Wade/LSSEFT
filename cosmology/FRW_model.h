//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_FRW_MODEL_H
#define LSSEFT_FRW_MODEL_H


#include "Planck_defaults.h"

#include "boost/serialization/serialization.hpp"


class FRW_model
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    //! cosmological parameters refer to their values today; if none are provided, then we use
    //! Planck values as sensible defaults
    FRW_model(double om = Planck2013::omega_m, double occ = Planck2013::omega_cc, double h_ = Planck2013::h,
              Mpc_units::energy tc = Planck2013::T_CMB, double ne = Planck2013::Neff,
              double fb = Planck2013::f_baryon, double zs = Planck2013::z_star, double zd = Planck2013::z_drag,
              double ze = Planck2013::z_eq, double Ac = Planck2013::Acurv, double n = Planck2013::ns,
              Mpc_units::energy kp = Planck2013::kpiv);

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
    Mpc_units::energy get_T_CMB() const { return(this->T_CMB); }

    //! get value of Neff
    double get_Neff() const { return(this->Neff); }
    
    //! get value of f_baryon
    double get_f_baryon() const { return this->f_baryon; }
    
    //! get zstar
    double get_z_star() const { return this->z_star; }
    
    //! get zdrag
    double get_z_drag() const { return this->z_drag; }
    
    //! get zeq
    double get_z_eq() const { return this->z_eq; }
    
    //! get value of A_curv
    double get_A_curv() const { return this->A_curv; }
    
    //! get value of n_s
    double get_ns() const { return this->ns; }
    
    //! get value of k_piv
    Mpc_units::energy get_k_piv() const { return this->k_piv; }


    // INTERNAL DATA

  private:

    //! omega_m, measured today
    double omega_m;

    //! omega_Lambda or omega_cc, measured today
    double omega_cc;

    //! reduced Hubble constant, ie. quoted in units of 100 km/s/Mpc, measured today
    double h;

    //! CMB temperature T_CMB, measured today
    Mpc_units::energy T_CMB;

    //! number of effective relativistic degrees of freedom
    double Neff;
    
    //! baryon fraction
    double f_baryon;
    
    //! z*
    double z_star;
    
    //! drag epoch
    double z_drag;
    
    //! matter-radiation equality
    double z_eq;
    
    //! amplitude of curvature perturbations
    double A_curv;
    
    //! spectral index of curvature perturbations
    double ns;
    
    //! pivot scale
    Mpc_units::energy k_piv;


    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & omega_m;
        ar & omega_cc;
        ar & h;
        ar & T_CMB;
        ar & Neff;
        ar & f_baryon;
        ar & z_star;
        ar & z_drag;
        ar & z_eq;
        ar & A_curv;
        ar & ns;
        ar & k_piv;
      }

  };


#endif //LSSEFT_FRW_MODEL_H
