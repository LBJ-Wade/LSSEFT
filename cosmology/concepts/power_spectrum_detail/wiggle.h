//
// Created by David Seery on 05/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_POWER_SPECTRUM_WIGGLE_H
#define LSSEFT_POWER_SPECTRUM_WIGGLE_H


#include "types.h"
#include "database/tokens.h"


class wiggle_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor -- populate with existing data
    wiggle_Pk(const linear_Pk_token& t, const tree_Pk_w::database_type& nw, const tree_Pk::database_type& r);
    
    //! destructor is default
    ~wiggle_Pk() = default;
    
    
    // DATABASE SERVICES
    
  public:
    
    //! get underlying power spectrum database for wiggle Pk
    const tree_Pk_w::database_type& get_nowiggle_db() const { return this->nowiggle.get_db(); }
    
    //! get underlying power spectrum database for raw Pk
    const tree_Pk::database_type& get_raw_db() const { return this->raw.get_db(); }
    
    //! ask spline to determine whether a given k-mode is acceptable for evaluation
    bool is_valid(const Mpc_units::energy& k) const;
    
    
    // TOKEN MANAGEMENT
    
  public:
    
    //! get token for linear power spectrum
    const linear_Pk_token& get_token() const { return this->tok; }
    
    
    // EVALUATION
    
  public:
    
    //! evaluate spline for wiggle Pk
    Mpc_units::inverse_energy3 Pk_wiggle(const Mpc_units::energy& k) const { return this->raw(k) - this->nowiggle(k); }
    
    //! evaluate spline for raw Pk
    Mpc_units::inverse_energy3 Pk_raw(const Mpc_units::energy& k) const { return this->raw(k); }
    
    //! evluate spline for no-wiggle Pk
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
      }
    
  };


namespace boost
  {
    
    namespace serialization
      {
        
        template <typename Archive>
        inline void save_construct_data(Archive& ar, const wiggle_Pk* t, const unsigned int file_version)
          {
            ar << t->get_token();
            ar << t->get_nowiggle_db();
            ar << t->get_raw_db();
          }
    
    
        template <typename Archive>
        inline void load_construct_data(Archive& ar, wiggle_Pk* t, const unsigned int file_version)
          {
            linear_Pk_token tok(0);
            tree_Pk_w::database_type nowiggle;
            tree_Pk::database_type raw;

            ar >> tok;
            ar >> nowiggle;
            ar >> raw;
            
            ::new(t) wiggle_Pk(tok, nowiggle, raw);
          }
        
      }   // namespace serialization
    
  }   // namespace boost


#endif //LSSEFT_POWER_SPECTRUM_WIGGLE_H
