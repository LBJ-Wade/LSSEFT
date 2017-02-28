//
// Created by David Seery on 05/12/2016.
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

#ifndef LSSEFT_POWER_SPECTRUM_WIGGLE_H
#define LSSEFT_POWER_SPECTRUM_WIGGLE_H


#include "types.h"
#include "database/tokens.h"


template <typename Tag>
class generic_wiggle_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor -- populate with existing data
    generic_wiggle_Pk(const linear_Pk_token& t, const tree_Pk_w::database_type& nw, const tree_Pk::database_type& r);
    
    //! destructor is default
    ~generic_wiggle_Pk() = default;
    
    
    // DATABASE SERVICES
    
  public:
    
    //! get underlying power spectrum database for wiggle Pk
    const tree_Pk_w::database_type& get_nowiggle_db() const { return this->nowiggle.get_db(); }
    
    //! get underlying power spectrum database for raw Pk
    const tree_Pk::database_type& get_raw_db() const { return this->raw.get_db(); }
    
    //! ask spline to determine whether a given k-mode is acceptable for evaluation
    bool is_valid(const Mpc_units::energy& k) const;
    
    //! get smallest k-value we can evaluate
    Mpc_units::energy get_min_k() const;
    
    //! get largest k-value we can evalute
    Mpc_units::energy get_max_k() const;

    
    // TOKEN MANAGEMENT
    
  public:
    
    //! get token for linear power spectrum
    const linear_Pk_token& get_token() const { return this->tok; }
    
    
    // RESCALING
    
  public:
    
    //! set rescaling factor
    generic_wiggle_Pk<Tag>& set_rescaling(double f=1.0)
      {
        this->nowiggle.set_rescaling(f);
        this->raw.set_rescaling(f);
        return *this;
      }
    
    
    // EVALUATION
    
  public:
    
    //! evaluate spline for wiggle Pk
    Mpc_units::inverse_energy3 Pk_wiggle(const Mpc_units::energy& k) const { return this->raw(k) - this->nowiggle(k); }
    
    //! evaluate spline for raw Pk
    Mpc_units::inverse_energy3 Pk_raw(const Mpc_units::energy& k) const { return this->raw(k); }
    
    //! evaluate spline for no-wiggle Pk
    Mpc_units::inverse_energy3 Pk_nowiggle(const Mpc_units::energy& k) const { return this->nowiggle(k); }
    
    
    // INTERNAL DATA
    
  private:
    
    //! token for corresponding linear power spectrum
    linear_Pk_token tok;
    
    //! wiggle Pk container
    tree_Pk_w nowiggle;
    
    //! raw Pk container
    tree_Pk raw;
    
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & tok;
        ar & nowiggle;
        ar & raw;
      }
    
  };


template <typename Tag>
generic_wiggle_Pk<Tag>::generic_wiggle_Pk(const linear_Pk_token& t, const tree_Pk_w::database_type& w, const tree_Pk::database_type& r)
  : tok(t),
    nowiggle(w),
    raw(r)
  {
  }


template <typename Tag>
bool generic_wiggle_Pk<Tag>::is_valid(const Mpc_units::energy& k) const
  {
    return this->nowiggle.is_valid(k) && this->raw.is_valid(k);
  }


template <typename Tag>
Mpc_units::energy generic_wiggle_Pk<Tag>::get_min_k() const
  {
    return std::min(this->nowiggle.get_min_k(), this->raw.get_min_k());
  }


template <typename Tag>
Mpc_units::energy generic_wiggle_Pk<Tag>::get_max_k() const
  {
    return std::max(this->nowiggle.get_min_k(), this->raw.get_min_k());
  }


namespace boost
  {
    
    namespace serialization
      {
        
        template <typename Archive, typename Tag>
        inline void save_construct_data(Archive& ar, const generic_wiggle_Pk<Tag>* t, const unsigned int file_version)
          {
          }
    
    
        template <typename Archive, typename Tag>
        inline void load_construct_data(Archive& ar, generic_wiggle_Pk<Tag>* t, const unsigned int file_version)
          {
            // construct an object with empty databases; these will be overwritten by the standard
            // deserialization
            linear_Pk_token tok(0);
            tree_Pk_w::database_type nowiggle;
            tree_Pk::database_type raw;
            ::new(t) generic_wiggle_Pk<Tag>(tok, nowiggle, raw);
          }
        
      }   // namespace serialization
    
  }   // namespace boost


namespace generic_wiggle_Pk_impl
  {
    
    template <int m>
    struct WigglePowerSpectrumTag
      {
        enum { WigglePowerSpectrumClass = m };
      };
    
    typedef WigglePowerSpectrumTag<0> InitialTag;
    typedef WigglePowerSpectrumTag<1> FinalTag;
    
  }   // namespace generic_wiggle_Pk_impl


// convenience types for initial and final power spectrum
typedef generic_wiggle_Pk< generic_wiggle_Pk_impl::InitialTag > initial_filtered_Pk;
typedef generic_wiggle_Pk< generic_wiggle_Pk_impl::FinalTag > final_filtered_Pk;


#endif //LSSEFT_POWER_SPECTRUM_WIGGLE_H
