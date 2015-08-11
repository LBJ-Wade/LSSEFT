//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_PLANCK_DEFAULTS_H
#define LSSEFT_PLANCK_DEFAULTS_H


#include "units/eV_units.h"


// Planck 2013+WP 'best fit' parameters from Table 2 on p.11 of arXiv:1303.5076
namespace Planck2013
  {

    constexpr double           omega_cc = 0.6817;
    constexpr double           omega_m  = 1.0 - omega_cc;
    constexpr double           h        = 0.6704;
    constexpr double           f_baryon = 0.15;   // Planck+WP best fit values is 0.154706
    constexpr eV_units::energy H0       = 100 * h * eV_units::Kilometre / (eV_units::Second * eV_units::Mpc);
    constexpr eV_units::energy T_CMB    = 2.725 * eV_units::Kelvin;

    // fluctuation two-point function
    constexpr double Acurv  = 2.2150E-9;
    constexpr double ns     = 0.9619;
    constexpr double kpiv   = 0.05; // measured in 1/Mpc
    constexpr double sigma8 = 0.8347;

    // CMB-related redshifts
    constexpr double z_star = 1090.48;
    constexpr double z_drag = 1059.25;
    constexpr double z_eq   = 3403;

  }   // namespace Planck


#endif //LSSEFT_PLANCK_DEFAULTS_H
