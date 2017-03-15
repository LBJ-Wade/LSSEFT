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

#include "boost/timer/timer.hpp"
#include "boost/serialization/serialization.hpp"

#include "cuba.h"


class Pk_filter_params
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! constructor
    Pk_filter_params(double A = LSSEFT_DEFAULT_FILTER_PK_AMPLITUDE,
                     Mpc_units::energy p = LSSEFT_DEFAULT_FILTER_PK_PIVOT,
                     double n = LSSEFT_DEFAULT_FILTER_PK_INDEX,
                     double r = LSSEFT_DEFAULT_FILTER_PK_REL_ERR,
                     double a = LSSEFT_DEFAULT_FILTER_PK_ABS_ERR)
      : amplitude(A),
        pivot(p),
        index(n),
        rel_err(r),
        abs_err(a)
      {
      }
    
    
    //! destructor is default
    ~Pk_filter_params() = default;
    
    
    // INTERFACE
  
  public:
    
    //! get amplitude
    double get_amplitude() const { return this->amplitude; }
    
    //! get pivot
    const Mpc_units::energy& get_pivot() const { return this->pivot; }
    
    //! get index
    double get_index() const { return this->index; }
    
    //! get relative error
    double get_relerr() const { return this->rel_err; }
    
    //! get absolute error
    double get_abserr() const { return this->abs_err; }
    
    
    // INTERNAL DATA
  
  private:
    
    //! amplitude of window function
    double amplitude;
    
    //! pivot of window function
    Mpc_units::energy pivot;
    
    //! spectral index of window function
    double index;
    
    //! relative error used during filtering
    double rel_err;
    
    //! absolute error used during filtering
    double abs_err;
    
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & amplitude;
        ar & pivot;
        ar & index;
        ar & rel_err;
        ar & abs_err;
      };
    
  };


template <typename ValueType>
class filter_result
  {
  
  public:
    
    typedef ValueType value_type;
    
    //! constructor initializes zero values
    filter_result()
      : value(value_type(0.0)),
        error(value_type(0.0)),
        regions(0),
        evaluations(0),
        time(0)
      {
      }
    
    //! destructor is default
    ~filter_result() = default;
    
    
    // DATA
    
  public:
    
    value_type                    value;
    value_type                    error;
    
    unsigned int                  regions;
    unsigned int                  evaluations;
    boost::timer::nanosecond_type time;
  
  private:
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & value;
        ar & regions;
        ar & evaluations;
        ar & error;
        ar & time;
      }
    
  };


typedef filter_result<Mpc_units::inverse_energy3> Pk_filter_result;


class Pk_filter
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! constructor accepts Pk_filter_data package
    Pk_filter(const Pk_filter_params& p)
      : params(p)
      {
      }
    
    //! destructor is default
    ~Pk_filter() = default;
    
    
    // INTERFACE
    
  public:
    
    //! filter a linear power spectrum; returns estimate of the no-wiggle component and
    //! the reference power spectrum for the same scale
    std::pair< Pk_filter_result, Mpc_units::inverse_energy3 >
    operator()(const FRW_model& model, const filterable_Pk& Pk_lin, const Mpc_units::energy& k);
    
    
    // INTERNAL API
    
  private:
    
    //! apply filter to a given integrand
    template <typename ResultType>
    bool integrate(const double slog_min, const double slog_max, const double klog, const double lambda,
                   const filterable_Pk& Pk_lin, const approx_Pk& Papprox, integrand_t integrand,
                   ResultType& result);

    //! compute Eisenstein & Hu approximation to the power spectrum
    std::unique_ptr<approx_Pk> eisenstein_hu(const FRW_model& model, const filterable_Pk& Pk_lin);
    
    
    // INTERNAL DATA
    
  private:
    
    //! filtering parameters
    const Pk_filter_params params;
    
  };


#endif //LSSEFT_PK_FILTER_H
