//
// Created by David Seery on 13/03/2017.
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

#ifndef LSSEFT_MDR1_SIM_H
#define LSSEFT_MDR1_SIM_H


#include "units/Mpc_units.h"


// Parameters for the WizCOLA realization suite
namespace WizCOLA
  {
    
    const     std::string       name     = "Matches WizCOLA realizations";

    constexpr double            omega_cc = 0.727;
    constexpr double            omega_m  = 0.273;
    constexpr double            h        = 0.705;
    constexpr double            f_baryon = 0.0456/omega_m;
    constexpr Mpc_units::energy H0       = 100 * h * Mpc_units::Kilometre / (Mpc_units::Second * Mpc_units::Mpc);
    constexpr Mpc_units::energy T_CMB    = 2.7255 * Mpc_units::Kelvin;
    constexpr double            Neff     = 3.046;   // Standard Model value
    
    // fluctuation two-point function
    constexpr double Acurv           = 2.183e-9;
    constexpr double ns              = 0.961;
    constexpr Mpc_units::energy kpiv = 0.05 / Mpc_units::Mpc;
    constexpr double sigma8          = 0.812;
    
    // CMB-related redshifts (extracted from CAMB)
    constexpr double z_star = 1088.02;
    constexpr double z_drag = 1060.39;
    constexpr double z_eq   = 3161.49;
    
  }


#endif //LSSEFT_MDR1_SIM_H
