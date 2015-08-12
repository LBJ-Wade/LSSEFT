//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
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
