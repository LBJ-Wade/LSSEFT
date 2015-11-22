//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_DEFAULTS_H
#define LSSEFT_DEFAULTS_H


#define LSSEFT_DEFAULT_PYTHON_PATH "/usr/local/python"

// database search tolerances
#define LSSEFT_DEFAULT_FRW_MODEL_PARAMETER_TOLERANCE      (1E-5)
#define LSSEFT_DEFAULT_REDSHIFT_CONFIGURATION_TOLERANCE   (1E-5)
#define LSSEFT_DEFAULT_WAVENUMBER_CONFIGURATION_TOLERANCE (1E-10)

// default absolute and relative errors during integration
#define LSSEFT_DEFAULT_ODE_ABS_ERR      (1E-12)
#define LSSEFT_DEFAULT_ODE_REL_ERR      (1E-6)

#define LSSEFT_DEFAULT_INTEGRAL_ABS_ERR (1E-10)
#define LSSEFT_DEFAULT_INTEGRAL_REL_ERR (1E-3)

#define LSSEFT_DEFAULT_CUHRE_KEY        (11)


#endif //LSSEFT_DEFAULTS_H
