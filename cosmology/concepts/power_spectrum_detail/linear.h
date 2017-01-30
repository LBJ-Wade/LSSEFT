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

#ifndef LSSEFT_POWER_SPECTRUM_LINEAR_H
#define LSSEFT_POWER_SPECTRUM_LINEAR_H


#include <string>

#include "types.h"

#include "boost/filesystem.hpp"


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
    bool is_valid(const Mpc_units::energy& k,
                  double bottom_clearance = SPLINE_PK_DEFAULT_BOTTOM_CLEARANCE,
                  double top_clearance = SPLINE_PK_DEFAULT_TOP_CLEARANCE) const;
    
    
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
            ar << t->get_db();
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
