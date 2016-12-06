//
// Created by David Seery on 05/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_POWER_SPECTRUM_WIGGLE_H
#define LSSEFT_POWER_SPECTRUM_WIGGLE_H


#include "types.h"


class wiggle_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor -- populate with existing data
    wiggle_Pk(const tree_Pk_w::database_type& w, const tree_Pk::database_type& r);
    
    //! destructor is default
    ~wiggle_Pk() = default;
    
    
    // DATABASE SERVICES
    
  public:
    
    //! get underlying power spectrum database for wiggle Pk
    const tree_Pk_w::database_type& get_wiggle_db() const { return this->wiggle.get_db(); }
    
    //! get underlying power spectrum database for raw Pk
    const tree_Pk::database_type& get_raw_db() const { return this->raw.get_db(); }
    
    //! ask spline to determine whether a given k-mode is acceptable for evaluation
    bool is_valid(const Mpc_units::energy& k) const { return this->wiggle.is_valid(k) && this->raw.is_valid(k); }
    
    
    // EVALUATION
    
  public:
    
    //! evaluate spline for wiggle Pk
    Mpc_units::inverse_energy3 Pk_wiggle(const Mpc_units::energy& k) const { return this->wiggle(k); }
    
    //! evaluate spline for raw Pk
    Mpc_units::inverse_energy3 Pk_raw(const Mpc_units::energy& k) const { return this->raw(k); }
    
    //! evluate spline for no-wiggle Pk
    Mpc_units::inverse_energy3 Pk_nowiggle(const Mpc_units::energy& k) const { return this->raw(k) - this->wiggle(k); }
    
    
    // INTERNAL DATA
    
  private:
    
    //! wiggle Pk container
    tree_Pk_w wiggle;
    
    //! raw Pk container
    tree_Pk raw;
    
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
      }
    
  };


namespace boost
  {
    
    namespace serialization
      {
        
        template <typename Archive>
        inline void save_construct_data(Archive& ar, const wiggle_Pk* t, const unsigned int file_version)
          {
            ar << t->get_wiggle_db();
            ar << t->get_raw_db();
          }
    
    
        template <typename Archive>
        inline void load_construct_data(Archive& ar, wiggle_Pk* t, const unsigned int file_version)
          {
            tree_Pk_w::database_type wiggle;
            tree_Pk::database_type raw;

            ar >> wiggle;
            ar >> raw;
            
            ::new(t) wiggle_Pk(wiggle, raw);
          }
        
      }   // namespace serialization
    
  }   // namespace boost


#endif //LSSEFT_POWER_SPECTRUM_WIGGLE_H
