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

#ifndef LSSEFT_ABSTRACT_RANGE_H
#define LSSEFT_ABSTRACT_RANGE_H


#include <vector>

#include "common.h"


template <typename Value>
class range
  {

    // DESTRUCTOR

  public:

    //! destructor is default; declare virtual
    virtual ~range() = default;


    // INTERFACE
    // note that these functions cannot be marked const if we wish to retain the option of lazy
    // evaluation for aggregation ranges

  public:

    //! get number of elements in the range
    virtual size_t size() = 0;

    //! get grid of elements
    virtual const std::vector<Value>& grid() = 0;

    //! subscripting operator
    virtual Value operator[](unsigned int element) = 0;


    // CLONE

  public:

    //! clone self (note covariant return type)
    virtual range<Value>* clone() const = 0;

  };


#endif //LSSEFT_ABSTRACT_RANGE_H
