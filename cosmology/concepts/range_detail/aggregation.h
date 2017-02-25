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

#ifndef LSSEFT_AGGREGATION_RANGE_H
#define LSSEFT_AGGREGATION_RANGE_H


#include <memory>
#include <algorithm>
#include <list>
#include <vector>
#include <assert.h>

#include "range.h"


// forward-declare addition overloads which will be friended in the class declaration
template <typename Value> class aggregation_range;

template <typename Value>
aggregation_range<Value> operator+(const aggregation_range<Value>& lhs, const range<Value>& rhs);

template <typename Value>
aggregation_range<Value> operator+(const range<Value>& lhs, const range<Value>& rhs);


template <typename Value>
class aggregation_range: public range<Value>
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! construct an empty aggregation range
    aggregation_range();

    //! construct an aggregation range from a single subrange
    aggregation_range(const range<Value>& a);

    //! construct an aggregation range from a pair of subranges
    aggregation_range(const range<Value>& a, const range<Value>& b);


    // OVERLOAD ARITHMETIC OPERATIONS TO ALLOW AGGREGATE RANGES TO BE COMBINED

  public:

    //! assignment
    aggregation_range<Value>& operator=(const range<Value>& rhs) { this->subrange_list.clear(); this->add_subrange(rhs); return(*this); }

    //! compound addition
    aggregation_range<Value>& operator+=(const range<Value>& rhs) { this->add_subrange(rhs); return(*this); }

    //! addition
    friend aggregation_range<Value> operator+ <>(const aggregation_range<Value>& lhs, const range<Value>& rhs);
    friend aggregation_range<Value> operator+ <>(const range<Value>& lhs, const range<Value>& rhs);


    // INTERFACE -- implements a 'range' concept

  public:

    //! get number of elements
    virtual size_t size() override                    { if(this->dirty) this->populate_grid(); return(this->elements.size()); }

    //! get grid of elements
    virtual const std::vector<Value>& grid() override { if(this->dirty) this->populate_grid(); return(this->elements); }

    //! overload subscript operator
    virtual Value operator[](unsigned int d) override;


    // ADD RANGES TO THE AGGREGATION

  public:

    //! add a subrange
    void add_subrange(const range<Value>& s);


    // INTERNAL API

  protected:

    //! clear grid and repopulate it from subranges
    void populate_grid();


    // CLONE -- implements a 'range' concept

    //! clone self (note covariant return type)
    virtual aggregation_range<Value>* clone() const override { return new aggregation_range<Value>(dynamic_cast<const aggregation_range<Value>&>(*this)); }


    // INTERNAL DATA

  protected:

    //! dirty flag - if set, the grid requires recalculation
    //! we use this approach to enable lazy evaluation; the grid isn't constructed until it is needed
    bool dirty;

    //! grid of values
    std::vector<Value> elements;

    //! set up type alias for subrange list; use std::unique_ptr<> to manage lifetime of instances
    typedef std::list< std::unique_ptr< range<Value> > > subrange_list_type;

    //! list of subranges
    subrange_list_type subrange_list;

  };


template <typename Value>
aggregation_range<Value> operator+(const aggregation_range<Value>& lhs, const range<Value>& rhs)
  {
    return(aggregation_range<Value>(lhs) += rhs);
  }


template <typename Value>
aggregation_range<Value> operator+(const range<Value>& lhs, const range<Value>& rhs)
  {
    return(aggregation_range<Value>(lhs, rhs));
  }


template <typename Value>
aggregation_range<Value>::aggregation_range()
  : dirty(true)
  {
  }


template <typename Value>
aggregation_range<Value>::aggregation_range(const range<Value>& a)
  : aggregation_range<Value>()    // delegate to empty constructor
  {
    this->subrange_list.emplace_back(a.clone());
  }


template <typename Value>
aggregation_range<Value>::aggregation_range(const range<Value>& a, const range<Value>& b)
  : aggregation_range<Value>()    // delegate to empty constructor
  {
    this->subrange_list.emplace_back(a.clone());
    this->subrange_list.emplace_back(b.clone());
  }


template <typename Value>
void aggregation_range<Value>::add_subrange(const range<Value>& s)
  {
    this->subrange_list.emplace_back(s.clone());
    this->dirty = true;
  }


template <typename Value>
Value aggregation_range<Value>::operator[](unsigned int d)
  {
    if(this->dirty) this->populate_grid();

    assert(d < this->elements.size());

    if(d < this->elements.size()) return(this->elements[d]);
    else                          throw std::out_of_range(ERROR_RANGE_OUT_OF_RANGE);
  }


template <typename Value>
void aggregation_range<Value>::populate_grid()
  {
    this->elements.clear();

    // splice subranges together
    for(typename subrange_list_type::const_iterator t = this->subrange_list.begin(); t != this->subrange_list.end(); ++t)
      {
        std::vector<Value> temp = (*t)->grid();
        this->elements.reserve(this->elements.size() + temp.size());
        this->elements.insert(this->elements.end(), temp.begin(), temp.end());
      }

    // sort resulting list into order and remove duplicates
    std::sort(this->elements.begin(), this->elements.end());

    auto last = std::unique(this->elements.begin(), this->elements.end(), range_impl::DuplicateRemovalPredicate<Value>(1E-10));
    this->elements.erase(last, this->elements.end());

    // flag grid as clean
    this->dirty = false;
  }


#endif //LSSEFT_AGGREGATION_RANGE_H
