//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_RANGE_COMMON_H
#define LSSEFT_RANGE_COMMON_H


#include <cmath>


namespace range_impl
  {

    template <typename Value>
    class DuplicateRemovalPredicate
      {

        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! constructor
        DuplicateRemovalPredicate(double t);

        //! destructor is default
        ~DuplicateRemovalPredicate() = default;


        // OPERATOR

      public:

        bool operator()(Value& a, Value& b)
          {
            return(std::abs((static_cast<double>(a)-static_cast<double>(b))/static_cast<double>(a))<tol);
          }


        // INTERNAL DATA

      protected:

        double tol;

      };


    template <typename Value>
    DuplicateRemovalPredicate<Value>::DuplicateRemovalPredicate(double t)
      : tol(t)
      {
      }

  }


#endif //LSSEFT_RANGE_COMMON_H
