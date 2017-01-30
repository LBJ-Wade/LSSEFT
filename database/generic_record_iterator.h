//
// Created by David Seery on 13/08/2015.
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

#ifndef LSSEFT_GENERIC_RECORD_ITERATOR_H
#define LSSEFT_GENERIC_RECORD_ITERATOR_H


#include <type_traits>
#include <iterator>


namespace configuration_database
  {

    //! record-valued iterator: when dereferenced, yields the underlying database record
    template <typename Iterator, typename ConstIterator, typename RecordType, bool is_const_iterator=true>
    class generic_record_iterator: public std::iterator< std::bidirectional_iterator_tag, RecordType >
      {

      private:

        //! set up a type alias for a reference to a database record;
        //! we need different versions for the const and non-const types
        typedef typename std::conditional< is_const_iterator, const RecordType&, RecordType& >::type reference_type;

        //! set up type alias for a pointer to a database record
        typedef typename std::conditional< is_const_iterator, const RecordType*, RecordType* >::type pointer_type;

        //! set up type alias for the underlying iterator
        typedef typename std::conditional< is_const_iterator, ConstIterator, Iterator >::type raw_iterator_type;


      public:

        //! set up a type alias for a record
        typedef typename std::conditional< is_const_iterator, const RecordType, RecordType >::type record_type;


        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! default constructor; iterator points to nothing when it is created
        generic_record_iterator() = default;

        //! value constructor; iterator points to given record in the underlying database
        generic_record_iterator(raw_iterator_type i)
          : iter(i)
          {
          }

        //! copy constructor; allows implicit conversion from a regular iterator to a const iterator
        generic_record_iterator(const generic_record_iterator<Iterator, ConstIterator, RecordType, false>& obj)
          : iter(obj.iter)
          {
          }


        // OVERLOAD COMPARISON OPERATORS

      public:

        //! equality comparison
        bool operator==(const generic_record_iterator& obj) const
          {
            return(this->iter == obj.iter);
          }

        //! inequality comparison
        bool operator!=(const generic_record_iterator& obj) const
          {
            return(this->iter != obj.iter);
          }


        // DEREFERENCING

      public:

        //! dereference iterator to get record which it points to
        reference_type operator*()
          {
            return(this->iter->second);
          }

        //! member access operator into record which iterator points to
        pointer_type operator->()
          {
            return(&this->iter->second);
          }

        //! dereference iterator to get record which it points to -- const version
        const reference_type operator*() const
          {
            return(this->iter->second);
          }

        //! member access operator into record which iterator points to -- const version
        const pointer_type operator->() const
          {
            return(&this->iter->second);
          }


        // INCREMENT, DECREMENT

      public:

        //! prefix decrement
        generic_record_iterator& operator--()
          {
            --this->iter;
            return(*this);
          }

        //! postfix decrement
        generic_record_iterator& operator--(int)
          {
            const generic_record_iterator old(*this);
            --this->iter;
            return(old);
          }

        //! prefix increment
        generic_record_iterator& operator++()
          {
            ++this->iter;
            return(*this);
          }

        //! postfix increment
        generic_record_iterator& operator++(int)
          {
            const generic_record_iterator old(*this);
            ++this->iter;
            return(old);
          }

        // make the const generic_record_iterator a friend of the non-const generic_record_iterator,
        // so the copy constructor can access its private member variables during implicit conversion
        friend class generic_record_iterator<Iterator, ConstIterator, RecordType, true>;



        // INTERNAL DATA

      private:

        //! current value of iterator
        raw_iterator_type iter;

      };

  }   // namespace configuration_database


#endif //LSSEFT_GENERIC_RECORD_ITERATOR_H
