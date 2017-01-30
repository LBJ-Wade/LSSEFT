//
// Created by David Seery on 11/08/2015.
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

#ifndef LSSEFT_PLANCK_DEFAULTS_H
#define LSSEFT_PLANCK_DEFAULTS_H


#include "units/Mpc_units.h"


// Planck 2013+WP 'best fit' parameters from Table 2 on p.11 of arXiv:1303.5076
namespace Planck2013
  {

    constexpr double            omega_cc = 0.6817;
    constexpr double            omega_m  = 1.0 - omega_cc;
    constexpr double            h        = 0.6704;
    constexpr double            f_baryon = 0.15401;
    constexpr Mpc_units::energy H0       = 100 * h * Mpc_units::Kilometre / (Mpc_units::Second * Mpc_units::Mpc);
    constexpr Mpc_units::energy T_CMB    = 2.7255 * Mpc_units::Kelvin;
    constexpr double            Neff     = 3.046;   // Standard Model value

    // fluctuation two-point function
    constexpr double Acurv           = 2.2150E-9;
    constexpr double ns              = 0.9619;
    constexpr Mpc_units::energy kpiv = 0.05 / Mpc_units::Mpc;
    constexpr double sigma8          = 0.8347;

    // CMB-related redshifts
    constexpr double z_star = 1090.48;
    constexpr double z_drag = 1059.25;
    constexpr double z_eq   = 3403;

  }   // namespace Planck2013


// Planck 2015 TT+TE+EE+lowP+lensing+ext parameters from Table 4 on p.31 of arXiv:1502.01589
namespace Planck2015
  {
    
    constexpr double            omega_cc = 0.6911;
    constexpr double            omega_m  = 1.0 - omega_cc;
    constexpr double            h        = 0.6774;
    constexpr double            f_baryon = 0.156821;
    constexpr Mpc_units::energy H0       = 100 * h * Mpc_units::Kilometre / (Mpc_units::Second * Mpc_units::Mpc);
    constexpr Mpc_units::energy T_CMB    = 2.7255 * Mpc_units::Kelvin;
    constexpr double            Neff     = 3.046;   // Standard Model value
    
    // fluctuation two-point function
    constexpr double Acurv           = 2.1420E-9;
    constexpr double ns              = 0.9667;
    constexpr Mpc_units::energy kpiv = 0.05 / Mpc_units::Mpc;
    constexpr double sigma8          = 0.8159;
    
    // CMB-related redshifts
    constexpr double z_star = 1089.90;
    constexpr double z_drag = 1059.68;
    constexpr double z_eq   = 3371;
    
  }


#endif //LSSEFT_PLANCK_DEFAULTS_H
