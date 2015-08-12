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

    // INTERRFACE

  public:

    //! get number of elements in the range
    virtual size_t size() const = 0;

    //! get grid of elements
    virtual const std::vector<Value>& grid() const = 0;

    //! subscripting operator
    virtual Value operator[](unsigned int element) const = 0;

  };


#endif //LSSEFT_ABSTRACT_RANGE_H
