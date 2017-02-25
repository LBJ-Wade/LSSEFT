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


// default Python location

#define LSSEFT_DEFAULT_PYTHON_PATH "/usr/local/python"

// enforce strict database consistency?
#define LSSEFT_STRICT_DATABASE_CONSISTENCY

// database search tolerances

constexpr double LSSEFT_DEFAULT_FRW_MODEL_PARAMETER_TOLERANCE      = 1E-5;
constexpr double LSSEFT_DEFAULT_REDSHIFT_CONFIGURATION_TOLERANCE   = 1E-5;
constexpr double LSSEFT_DEFAULT_WAVENUMBER_CONFIGURATION_TOLERANCE = 1E-10;

// default absolute and relative errors during integration
constexpr double LSSEFT_DEFAULT_ODE_ABS_ERR                        = (1E-12);
constexpr double LSSEFT_DEFAULT_ODE_REL_ERR                        = (1E-6);

constexpr double LSSEFT_DEFAULT_INTEGRAL_ABS_ERR_13                = (1E-10);
constexpr double LSSEFT_DEFAULT_INTEGRAL_REL_ERR_13                = (1E-8);

constexpr double LSSEFT_DEFAULT_INTEGRAL_ABS_ERR_22                = (1E-10);
constexpr double LSSEFT_DEFAULT_INTEGRAL_REL_ERR_22                = (1E-6);

constexpr double LSSEFT_DEFAULT_FILTER_PK_ABS_ERR                  = (1E-10);
constexpr double LSSEFT_DEFAULT_FILTER_PK_REL_ERR                  = (1E-6);


#endif //LSSEFT_DEFAULTS_H
