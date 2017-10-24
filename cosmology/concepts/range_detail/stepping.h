//
// Created by David Seery on 12/08/2015.
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

#ifndef LSSEFT_STEPPING_RANGE_H
#define LSSEFT_STEPPING_RANGE_H


#include <assert.h>
#include <algorithm>
#include <vector>

#include "range.h"

#include "exceptions.h"
#include "localizations/messages.h"


enum class spacing_type
  {
    linear,
    logarithmic_bottom,
    logarithmic_top
  };

template <typename Value>
class stepping_range: public range<Value>
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor; defaults to linear spacing if not specified
    stepping_range(double lo, double hi, size_t st, const Value& u, spacing_type tp= spacing_type::linear);

    //! destructor is defaults
    virtual ~stepping_range() = default;


    // INTERFACE -- implements a 'range' concept

  public:

    //! get number of elements
    virtual size_t size() override                    { return(this->elements.size()); }

    //! get grid of elements
    virtual const std::vector<Value>& grid() override { return(this->elements); }

    //! overload subscript operator
    virtual Value operator[](unsigned int d) override;


    // INTERNAL API

  private:

    //! reset grid of elements and recompute it
    void populate_grid();

    //! populate grid for subcase of linear spacing
    void populate_linear_grid();

    //! populate grid for subcase of logarithmic bottom spacing
    void populate_log_bottom_grid();

    //! populate grid for subcase of logarithmic top spacing
    void populate_log_top_grid();


    // CLONE -- implements a 'range' concept

    //! clone self (note covariant return type)
    virtual stepping_range<Value>* clone() const override { return new stepping_range<Value>(dynamic_cast<const stepping_range<Value>&>(*this)); }


    // INTERNAL DATA

  private:

    //! minimum value
    double min;

    //! maximum value
    double max;

    //! number of steps
    size_t steps;

    //! unit in which values are measured
    const Value unit;

    //! spacing time
    spacing_type spacing;

    //! grid of values
    std::vector<Value> elements;

  };


template <typename Value>
stepping_range<Value>::stepping_range(double lo, double hi, size_t st, const Value& u, spacing_type tp)
  : min(lo),
    max(hi),
    steps(st),
    spacing(tp),
    unit(u)
  {
    this->populate_grid();
  }


template <typename Value>
Value stepping_range<Value>::operator[](unsigned int d)
  {
    assert(d < this->elements.size());

    if(d < this->elements.size()) return(this->elements[d]);
    else                          throw std::out_of_range(ERROR_RANGE_OUT_OF_RANGE);

  }


template <typename Value>
void stepping_range<Value>::populate_grid()
  {
    this->elements.clear();
    this->elements.reserve(this->steps+1);

    switch(this->spacing)
      {
        case spacing_type::linear:
          {
            this->populate_linear_grid();
            break;
          }

        case spacing_type::logarithmic_bottom:
          {
            this->populate_log_bottom_grid();
            break;
          }

        case spacing_type::logarithmic_top:
          {
            this->populate_log_top_grid();
            break;
          }
      }

    // sort grid into order and remove duplicates
    std::sort(this->elements.begin(), this->elements.end());

    auto last = std::unique(this->elements.begin(), this->elements.end(), range_impl::DuplicateRemovalPredicate<Value>(1E-10));
    this->elements.erase(last, this->elements.end());
  }


template <typename Value>
void stepping_range<Value>::populate_linear_grid()
  {
    if(this->steps == 0)
      {
        this->elements.push_back(this->min * this->unit);
      }
    else
      {
        for(unsigned int i = 0; i <= this->steps; ++i)
          {
            double v = this->min + (static_cast<double>(i)/this->steps)*(this->max-this->min);
            if(!std::isnan(v)) this->elements.push_back(v * this->unit);
          }
      }
  }


template <typename Value>
void stepping_range<Value>::populate_log_bottom_grid()
  {
    if(this->steps == 0)
      {
        this->elements.push_back(this->min * this->unit);
      }
    else
      {
        // if max and min are both positive, perform log-spacing as we expect
        if(this->max > 0.0 && this->min > 0.0)
          {
            for(unsigned int i = 0; i <= this->steps; ++i)
              {
                double v = this->min * std::pow(this->max/this->min, static_cast<double>(i)/this->steps);
                if(!std::isnan(v)) this->elements.push_back(v * this->unit);
              }
          }
        else
          // otherwise we shift the range to begin at 1, perform the log-spacing, and then
          // reverse the shift
          {
            double shifted_max = this->max - (this->min-1.0);
            for(unsigned int i = 0; i <= this->steps; ++i)
              {
                double v = this->min-1.0 + std::pow(shifted_max, static_cast<double>(i)/this->steps);
                if(!std::isnan(v)) this->elements.push_back(v * this->unit);
              }
          }
      }
  }


template <typename Value>
void stepping_range<Value>::populate_log_top_grid()
  {
    if(this->steps == 0)
      {
        this->elements.push_back(this->min * this->unit);
      }
    else
      {
        if(this->max > 0.0 && this->min > 0.0)
          {
            for(unsigned int i = 0; i <= this->steps; ++i)
              {
                double v = this->max + this->min - this->min * std::pow(this->max/this->min, static_cast<double>(i)/this->steps);
                if(!std::isnan(v)) this->elements.push_back(v * this->unit);
              }
          }
        else
          {
            double shifted_max = this->max - (this->min-1.0);
            for(unsigned int i = 0; i <= this->steps; ++i)
              {
                double v = this->max+1.0 - std::pow(shifted_max, static_cast<double>(i)/this->steps);
                if(!std::isnan(v)) this->elements.push_back(v * this->unit);
              }
          }
        // the result is out-of-order, but it will be sorted later
      }

  }


#endif //LSSEFT_STEPPING_RANGE_H
