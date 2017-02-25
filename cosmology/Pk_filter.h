//
// Created by David Seery on 07/12/2016.
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

#ifndef LSSEFT_PK_FILTER_H
#define LSSEFT_PK_FILTER_H


#include "FRW_model.h"
#include "cosmology/concepts/power_spectrum.h"

#include "database/tokens.h"

#include "defaults.h"

#include "cuba.h"


class Pk_filter
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! constructor is default
    Pk_filter(double r = LSSEFT_DEFAULT_FILTER_PK_REL_ERR,
              double a = LSSEFT_DEFAULT_FILTER_PK_ABS_ERR)
      : rel_err(std::abs(r)),
        abs_err(std::abs(a))
      {
      }
    
    //! destructor is default
    ~Pk_filter() = default;
    
    
    // INTERFACE
    
  public:
    
    //! filter a linear power spectrum; returns estimate of the no-wiggle component and
    //! the reference power spectrum for the same scale
    std::pair< Mpc_units::inverse_energy3, Mpc_units::inverse_energy3 >
    operator()(const FRW_model& model, const filterable_Pk& Pk_lin, const Mpc_units::energy& k);
    
    
    // INTERNAL API
    
  private:
    
    //! apply filter to a given integrand
    double integrate(const double slog_min, const double slog_max, const double klog, const double lambda,
                     const filterable_Pk& Pk_lin, const approx_Pk& Papprox, integrand_t integrand);

    //! compute Eisenstein & Hu approximation to the power spectrum
    std::unique_ptr<approx_Pk> eisenstein_hu(const FRW_model& model, const filterable_Pk& Pk_lin);
    
    
    // INTERNAL DATA
    
  private:
    
    //! relative tolerance
    double rel_err;
    
    //! absolute tolerance
    double abs_err;
    
  };


#endif //LSSEFT_PK_FILTER_H
