//
// Created by David Seery on 10/08/2015.
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

#ifndef LSSEFT_DEFAULTS_H
#define LSSEFT_DEFAULTS_H


#include "units/Mpc_units.h"


// default terminal environment variable
#define LSSEFT_TERM_ENV "TERM"

// default filesystem search path
#define LSSEFT_SHELL_PATH_ENV "PATH"

// Python executable name
#define LSSEFT_PYTHON_EXECUTABLE "python"

// default Python location
#define LSSEFT_DEFAULT_PYTHON_PATH "/usr/local/python"

// database search tolerances

constexpr double LSSEFT_DEFAULT_FRW_MODEL_PARAMETER_TOLERANCE       = 1E-5;
constexpr double LSSEFT_DEFAULT_REDSHIFT_CONFIGURATION_TOLERANCE    = 1E-5;
constexpr double LSSEFT_DEFAULT_WAVENUMBER_CONFIGURATION_TOLERANCE  = 1E-10;
constexpr double LSSEFT_DEFAULT_FILTER_CONFIGURATION_TOLERANCE      = 1E-5;
constexpr double LSSEFT_DEFAULT_ONELOOP_CONFIGURATION_TOLERANCE     = 1E-5;
constexpr double LSSEFT_DEFAULT_MATSUBARAXY_CONFIGURATION_TOLERANCE = 1E-5;
constexpr double LSSEFT_DEFAULT_GROWTH_CONFIGURATION_TOLERANCE      = 1E-5;

// default absolute and relative errors during integration
constexpr double LSSEFT_DEFAULT_ODE_ABS_ERR                         = (1E-12);
constexpr double LSSEFT_DEFAULT_ODE_REL_ERR                         = (1E-6);

constexpr double LSSEFT_DEFAULT_INTEGRAL_ABS_ERR_13                 = (1E-3);
constexpr double LSSEFT_DEFAULT_INTEGRAL_REL_ERR_13                 = (1E-3);

constexpr double LSSEFT_DEFAULT_INTEGRAL_ABS_ERR_22                 = (1E-3);
constexpr double LSSEFT_DEFAULT_INTEGRAL_REL_ERR_22                 = (1E-3);

constexpr double LSSEFT_DEFAULT_FILTER_PK_ABS_ERR                   = (1E-3);
constexpr double LSSEFT_DEFAULT_FILTER_PK_REL_ERR                   = (1E-3);

// scale 0.25 for lambda suggested by Vlah, Seljak, Chu & Feng p.23 arXiv:1509.02120
// they also suggested it should grow slightly at large k, so we choose an index of 0.04 by default
constexpr double LSSEFT_DEFAULT_FILTER_PK_AMPLITUDE                 = (0.25);
constexpr double LSSEFT_DEFAULT_FILTER_PK_PIVOT                     = (0.07);
constexpr double LSSEFT_DEFAULT_FILTER_PK_INDEX                     = (0.04);

constexpr Mpc_units::inverse_energy LSSEFT_DEFAULT_RESUM_QMIN       = 10 * Mpc_units::Mpc;
constexpr Mpc_units::inverse_energy LSSEFT_DEFAULT_RESUM_QMAX       = 300 * Mpc_units::Mpc;

// cross-over scale from series expansion of RSD mapping to exp + erf representation
constexpr double LSSEFT_SERIES_CROSSOVER = 0.15;

#endif //LSSEFT_DEFAULTS_H
