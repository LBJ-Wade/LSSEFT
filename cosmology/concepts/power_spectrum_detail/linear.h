//
// Created by David Seery on 05/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_POWER_SPECTRUM_LINEAR_H
#define LSSEFT_POWER_SPECTRUM_LINEAR_H


#include <string>

#include "types.h"

#include "boost/filesystem.hpp"


// place serialization overloads before linear_Pk class definition to avoid compiler problems with overload resolution
// for templated friend functions (should disappear as compilers mature)
class linear_Pk;

namespace boost
  {
    
    namespace serialization
      {

        template <typename Archive>
        void save_construct_data(Archive& ar, const linear_Pk* t, const unsigned int file_version);

        template <typename Archive>
        inline void load_construct_data(Archive& ar, linear_Pk* t, const unsigned int file_version);
        
      }
    
  }


// container class for linear power spectrum, read in from a file
class linear_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! constructor -- read in from a file in CAMB format
    linear_Pk(const boost::filesystem::path& p);
    
    //! constructor -- populate with existing data
    linear_Pk(const std::string& p, const tree_Pk::database_type& d, const std::string& h);
    
    //! destructor is default
    ~linear_Pk() = default;
    
    
    // DATABASE SERVICES
  
  public:
    
    //! get underlying power spectrum database
    const tree_Pk::database_type& get_db() const { return this->container.get_db(); }
    
    //! ask spline to determine whether a given k-mode is acceptable for evaluation
    bool is_valid(const Mpc_units::energy& k) const { return this->container.is_valid(k); }
    
    
    // EVALUATION
    
  public:
    
    //! evaluate spline
    Mpc_units::inverse_energy3 operator()(const Mpc_units::energy& k) const { return this->container(k); }
    
    
    // METADATA
    
  public:
    
    //! get file path
    const boost::filesystem::path& get_path() const { return this->path; }
    
    //! get MD5 hash for file
    const std::string& get_MD5_hash() const { return this->md5_hash; }
    
    
    // INTERNAL API
    
  private:
    
    //! evaluate MD5 hash for a file
    std::string hash(const boost::filesystem::path& p);
    
    
    // INTERNAL DATA
  
  private:
    
    // declare boost::serialization::save_construct_data() to be friend, so it can archive the tree_Pk container
    template <typename Archive>
    friend void boost::serialization::save_construct_data(Archive& ar, const linear_Pk* t, const unsigned int file_version);
    
    //! fully-qualified path to file; should come before container in declaration so that it is constructed first
    boost::filesystem::path path;
    
    //! tree power spectrum container
    tree_Pk container;
    
    //! MD5-hash for file
    std::string md5_hash;


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
        inline void save_construct_data(Archive& ar, const linear_Pk* t, const unsigned int file_version)
          {
            ar << t->get_path().string();
            ar << t->container.get_db();
            ar << t->get_MD5_hash();
          }
        
        
        template <typename Archive>
        inline void load_construct_data(Archive& ar, linear_Pk* t, const unsigned int file_version)
          {
            std::string path;
            std::string hash;
            tree_Pk::database_type db;
            
            ar >> path;
            ar >> db;
            ar >> hash;
            
            ::new(t) linear_Pk(path, db, hash);
          }
        
      }   // namespace serialization
    
  }   // namespace boost


#endif //LSSEFT_POWER_SPECTRUM_LINEAR_H
