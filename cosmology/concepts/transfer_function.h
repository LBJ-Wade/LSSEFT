//
// Created by David Seery on 17/08/2015.
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

#ifndef LSSEFT_TRANSFER_FUNCTION_H
#define LSSEFT_TRANSFER_FUNCTION_H


#include <memory>
#include <vector>

#include "database/tokens.h"
#include "database/z_database.h"
#include "units/Mpc_units.h"

#include "boost/timer/timer.hpp"
#include "boost/serialization/serialization.hpp"


struct transfer_record
  {
    double delta_m;
    double theta_m;
    double delta_r;
    double theta_r;
    double Phi;
  };


typedef std::pair< const z_token&, transfer_record > transfer_value;


namespace transfer_function_impl
  {

    template <typename RecordIterator, typename ConstRecordIterator, typename ValueIterator, typename ConstValueIterator, bool is_const_iterator=true>
    class generic_tokenized_database_iterator: public std::iterator< std::bidirectional_iterator_tag, transfer_value >
      {

      private:

        //! type alias for raw record iterator
        typedef typename std::conditional< is_const_iterator, ConstRecordIterator, RecordIterator >::type record_iterator;

        //! type alias for raw value iterator
        typedef typename std::conditional< is_const_iterator, ConstValueIterator, ValueIterator >::type value_iterator;


        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! default constructor: points to nothing when it is constructed
        generic_tokenized_database_iterator() = default;

        //! value constructor; points to a given element in the transfer function sample
        generic_tokenized_database_iterator(record_iterator r, value_iterator dm, value_iterator tm,
                                            value_iterator dr, value_iterator tr, value_iterator P)
          : record_iter(r),
            delta_m_iter(dm),
            theta_m_iter(tm),
            delta_r_iter(dr),
            theta_r_iter(tr),
            Phi_iter(P)
          {
          }

        //! copy constructor; allows implicit conversion from a regular iterator to a const iterator
        generic_tokenized_database_iterator(const generic_tokenized_database_iterator<RecordIterator, ConstRecordIterator, ValueIterator, ConstValueIterator, false>& obj)
          : record_iter(obj.record_iter),
            delta_m_iter(obj.delta_m_iter),
            theta_m_iter(obj.theta_m_iter),
            delta_r_iter(obj.delta_r_iter),
            theta_r_iter(obj.theta_r_iter),
            Phi_iter(obj.Phi_iter)
          {
          }


        // OVERLOAD COMPARISON OPERATORS

      public:

        //! equality comparison
        bool operator==(const generic_tokenized_database_iterator& obj) const
          {
            // all should be in step, so need only compare one of them
            return(this->record_iter == obj.record_iter);
          }

        //! inequality comparison
        bool operator!=(const generic_tokenized_database_iterator& obj) const
          {
            // all iterators should be in step, so need only compare one of them
            return(this->record_iter != obj.record_iter);
          }


        // DEREFERENCING

      public:

        //! dereference iterator to get value
        transfer_value operator*() const
          {
            transfer_record rec;
            rec.delta_m = *this->delta_m_iter;
            rec.theta_m = *this->theta_m_iter;
            rec.delta_r = *this->delta_r_iter;
            rec.theta_r = *this->theta_r_iter;
            rec.Phi     = *this->Phi_iter;

            return transfer_value(this->record_iter->get_token(), rec);
          }


        // INCREMENT, DECREMENT

      public:

        //! prefix decrement
        generic_tokenized_database_iterator& operator--()
          {
            --this->delta_m_iter;
            --this->theta_m_iter;
            --this->delta_r_iter;
            --this->theta_r_iter;
            --this->Phi_iter;
            --this->record_iter;
            return(*this);
          }

        //! postfix decrement
        generic_tokenized_database_iterator& operator--(int)
          {
            const generic_tokenized_database_iterator old(*this);
            --(*this);
            return(old);
          }

        //! prefix increment
        generic_tokenized_database_iterator& operator++()
          {
            ++this->delta_m_iter;
            ++this->theta_m_iter;
            ++this->delta_r_iter;
            ++this->theta_r_iter;
            ++this->Phi_iter;
            ++this->record_iter;
            return(*this);
          }

        //! postfix increment
        generic_tokenized_database_iterator& operator++(int)
          {
            const generic_tokenized_database_iterator old(*this);
            ++(*this);
            return(old);
          }

        // make the const version a friend of the non-const version,
        // so the copy constructor can access its private members during implicit conversion
        friend class generic_tokenized_database_iterator<RecordIterator, ConstRecordIterator, ValueIterator, ConstValueIterator, true>;


        // INTERNAL DATA

      private:

        //! iterator into redshift database
        record_iterator record_iter;

        //! iterator into delta_m sample
        value_iterator delta_m_iter;

        //! iterator into theta_m sample
        value_iterator theta_m_iter;

        //! iterator into delta_r sample
        value_iterator delta_r_iter;

        //! iterator into theta_r sample
        value_iterator theta_r_iter;

        //! iterator into Phi sample
        value_iterator Phi_iter;

      };

  }   // namespace transfer_function_impl;


class transfer_function
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    transfer_function(const Mpc_units::energy& _k, const k_token& t, std::shared_ptr<z_database> z);

    //! destructor is default
    ~transfer_function() = default;


    // ITERATORS

  public:

    //! type alias for non-const iterator
    typedef transfer_function_impl::generic_tokenized_database_iterator<z_database::reverse_record_iterator, z_database::const_reverse_record_iterator,
                                                                        std::vector<double>::iterator, std::vector<double>::const_iterator, false> iterator;

    //! type alias for const iterator
    typedef transfer_function_impl::generic_tokenized_database_iterator<z_database::reverse_record_iterator, z_database::const_reverse_record_iterator,
                                                                        std::vector<double>::iterator, std::vector<double>::const_iterator,true> const_iterator;


    iterator begin()
      {
        return(iterator(this->z_db->record_rbegin(), this->delta_m->begin(), this->theta_m->begin(), this->delta_r->begin(), this->theta_r->begin(), this->Phi->begin()));
      }

    iterator end()
      {
        return(iterator(this->z_db->record_rend(), this->delta_m->end(), this->theta_m->end(), this->delta_r->end(), this->theta_r->end(), this->Phi->end()));
      }

    const_iterator begin() const
      {
        return(const_iterator(this->z_db->record_crbegin(), this->delta_m->cbegin(), this->theta_m->cbegin(), this->delta_r->cbegin(), this->theta_r->cbegin(), this->Phi->cbegin()));
      }

    const_iterator end() const
      {
        return(const_iterator(this->z_db->record_crend(), this->delta_m->cend(), this->theta_m->cend(), this->delta_r->cend(), this->theta_r->cend(), this->Phi->cend()));
      }

    const_iterator cbegin() const
      {
        return (const_iterator(this->z_db->record_crbegin(), this->delta_m->cbegin(), this->theta_m->cbegin(), this->delta_r->cbegin(), this->theta_r->cbegin(), this->Phi->cbegin()));
      }

    const_iterator cend() const
      {
        return(const_iterator(this->z_db->record_crend(), this->delta_m->cend(), this->theta_m->cend(), this->delta_r->cend(), this->theta_r->cend(), this->Phi->cend()));
      }


    // INTERFACE

  public:

    //! store components of the transfer functions
    void push_back(double delta_m, double delta_r, double theta_m, double theta_r, double Phi);

    //! get wavenumber token
    const k_token& get_k_token() const { return(this->token); }


    // METADATA

  public:

    //! store integration time
    void set_integration_metadata(boost::timer::nanosecond_type t, size_t s);

    //! get integration time
    boost::timer::nanosecond_type get_integration_time() const { return(this->integration_time); }

    //! get number of steps used by integrator
    size_t get_integration_steps() const { return(this->steps); }


    // INTERNAL DATA

  private:

    // CONFIGURATION DATA

    //! wavenumber
    Mpc_units::energy k;

    //! wavenumber token
    k_token token;

    //! redshift database; managed using a std::shared_ptr<>
    //! to avoid unnecessary duplication expense
    std::shared_ptr<z_database> z_db;


    // TRANSFER FUNCTIONS

    // these are managed using std::shared_ptr<>s to avoid
    // expensive duplication

    //! delta_m transfer function
    std::shared_ptr< std::vector<double> > delta_m;

    //! theta_m transfer function
    std::shared_ptr< std::vector<double> > theta_m;

    //! delta_r transfer function
    std::shared_ptr< std::vector<double> > delta_r;

    //! theta_r transfer function
    std::shared_ptr< std::vector<double> > theta_r;

    //! Phi transfer function
    std::shared_ptr< std::vector<double> > Phi;


    // METADATA

    //! time taken to perform integration
    boost::timer::nanosecond_type integration_time;

    //! number of steps used by integrator
    size_t steps;


    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & k;
        ar & token;
        ar & z_db;
        ar & delta_m;
        ar & theta_m;
        ar & delta_r;
        ar & theta_r;
        ar & Phi;
        ar & integration_time;
        ar & steps;
      }

  };


#endif //LSSEFT_TRANSFER_FUNCTION_H
