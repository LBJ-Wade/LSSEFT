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
