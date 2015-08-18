//
// Created by David Seery on 13/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_GENERIC_VALUE_ITERATOR_H
#define LSSEFT_GENERIC_VALUE_ITERATOR_H


#include <type_traits>
#include <iterator>


namespace configuration_database
  {

    //! record-valued iterator: when dereferenced, yields the value stored in the underlying database record
    template <typename Iterator, typename ConstIterator, typename ValueType, bool is_const_iterator=true>
    class generic_value_iterator : public std::iterator< std::bidirectional_iterator_tag, ValueType >
      {

      private:

        //! set up a type alias for a reference to a database record;
        //! we need different versions for the const and non-const types
        typedef typename std::conditional< is_const_iterator, const ValueType&, ValueType& >::type reference_type;

        //! set up type alias for a pointer to a database record
        typedef typename std::conditional< is_const_iterator, const ValueType*, ValueType* >::type pointer_type;

        //! set up type alias for the underlying iterator into DatabaseType
        typedef typename std::conditional< is_const_iterator, ConstIterator, Iterator >::type raw_iterator_type;


      public:

        //! set up a type alias for a record
        typedef typename std::conditional< is_const_iterator, const ValueType, ValueType >::type value_type;


        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! default constructor; iterator points to nothing when it is created
        generic_value_iterator() = default;

        //! value constructor; iterator points to given record in the underlying database
        generic_value_iterator(raw_iterator_type i)
          : iter(i)
          {
          }

        //! copy constructor; allows implicit conversion from a regular iterator to a const iterator
        generic_value_iterator(const generic_value_iterator<Iterator, ConstIterator, ValueType, false>& obj)
          : iter(obj.iter)
          {
          }

        //! equality comparison
        bool operator==(const generic_value_iterator& obj) const
          {
            return(this->iter == obj.iter);
          }

        //! inequality comparison
        bool operator!=(const generic_value_iterator& obj) const
          {
            return(this->iter != obj.iter);
          }

        //! dereference iterator to get value held by record which iterator points to
        value_type operator*()
          {
            return(*this->iter->second);
          }

        //! member access operator into value held by record which iterator points to
        pointer_type operator->()
          {
            return(&(*this->iter->second));
          }

        //! dereference iterator to get value held by record which iterator points to -- const version
        const value_type operator*() const
          {
            return(*this->iter->second);
          }

        //! member access operator into value held by record which iterator points to -- const version
        const pointer_type operator->() const
          {
            return(&(*this->iter->second));
          }

        //! prefix decrement
        generic_value_iterator& operator--()
          {
            --this->iter;
            return(*this);
          }

        //! postfix decrement
        generic_value_iterator& operator--(int)
          {
            const generic_value_iterator old(*this);
            --this->iter;
            return(old);
          }

        //! prefix increment
        generic_value_iterator& operator++()
          {
            ++this->iter;
            return(*this);
          }

        //! postfix increment
        generic_value_iterator& operator++(int)
          {
            const generic_value_iterator old(*this);
            ++this->iter;
            return(old);
          }

        // make the const generic_value_iterator a friend of the non-const generic_value_iterator,
        // so the copy constructor can access its private member variables during implicit conversion
        friend class generic_value_iterator<Iterator, ConstIterator, ValueType, true>;



        // INTERNAL DATA

      private:

        //! current value of iterator
        raw_iterator_type iter;

      };

  }   // namespace configuration_database


#endif //LSSEFT_GENERIC_VALUE_ITERATOR_H
