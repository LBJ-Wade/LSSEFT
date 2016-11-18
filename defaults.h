//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
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


#endif //LSSEFT_DEFAULTS_H
