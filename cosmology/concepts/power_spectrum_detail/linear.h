//
// Created by David Seery on 05/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_POWER_SPECTRUM_LINEAR_H
#define LSSEFT_POWER_SPECTRUM_LINEAR_H


#include <string>

#include "types.h"

#include "boost/filesystem.hpp"


// container class for linear power spectrum, read in from a file
class linear_power_spectrum
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! constructor -- read in from a file in CAMB format
    linear_power_spectrum(const boost::filesystem::path& p);
    
    //! destructor is default
    ~linear_power_spectrum() = default;
    
    
    // INTERFACE
  
  public:
    
    //! evaluate spline
    Mpc_units::inverse_energy3 operator()(const Mpc_units::energy& k) const { return this->container(k); }
    
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
    tree_power_spectrum container;
    
    //! MD5-hash for file
    std::string md5_hash;
    
  };


#endif //LSSEFT_POWER_SPECTRUM_LINEAR_H
