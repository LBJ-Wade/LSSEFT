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


#include <fstream>
#include <string>

#include "types.h"
#include "wiggle.h"

#include "boost/filesystem.hpp"
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/split_member.hpp"

#include "openssl/md5.h"



namespace generic_linear_Pk_impl
  {
    
    template <int m>
    struct LinearPowerSpectrumTag
      {
        enum { LinearPowerSpectrumClass=m };
      };
    
    typedef LinearPowerSpectrumTag<0> InitialTag;
    typedef LinearPowerSpectrumTag<1> FinalTag;
    
    typedef LinearPowerSpectrumTag<99> FilterableTag;
    
  }


// container class for linear power spectrum, read in from a file;
// the tag can be used to distinguish power spectra for different uses (eg. initial and final versions)
// so that they cannot be inadvertently mixed up

template <typename Tag, typename FilterPartnerType>
class generic_linear_Pk
  {
    
    // TYPE DEFINITIONS
    
  public:

    // filtered_Pk_type isn't used in generic_linear_Pk (which has no internal notion of filtering;
    // it's just a container for a power spectrum), but we store it here
    // so the information about the filtered partner type is just carried along with generic_linear_Pk
    typedef FilterPartnerType filtered_Pk_type;
    
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! constructor -- read in from a file in CAMB format
    generic_linear_Pk(const boost::filesystem::path& p);
    
    //! constructor -- populate with existing data
    generic_linear_Pk(const std::string& p, const tree_Pk::database_type& d, const std::string& h);
    
    //! constructor -- populate with existing data
    generic_linear_Pk(const boost::filesystem::path& p, const tree_Pk::database_type& d, const std::string& h);
    
    //! destructor is default
    ~generic_linear_Pk() = default;
    
    
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
    void save(Archive& ar, unsigned int version) const
      {
        ar << path.string();
        ar << container;
        ar << md5_hash;
      }
    
    template <typename Archive>
    void load(Archive& ar, unsigned int version)
      {
        std::string path_string;
        ar >> path_string;
        ar >> container;
        ar >> md5_hash;
        
        path = path_string;
      }
    
    BOOST_SERIALIZATION_SPLIT_MEMBER()
    
  };


template <typename Tag, typename FilterPartnerType>
generic_linear_Pk<Tag, FilterPartnerType>::generic_linear_Pk(const boost::filesystem::path& p)
  : path(p.is_absolute() ? p : boost::filesystem::absolute(p)),
    container(path)    // need to be sure path is initialized before container
  {
    this->md5_hash = this->hash(p);
  }


template <typename Tag, typename FilterPartnerType>
generic_linear_Pk<Tag, FilterPartnerType>::generic_linear_Pk(const std::string& p, const tree_Pk::database_type& d, const std::string& h)
  : path(p),
    container(d),
    md5_hash(h)
  {
  }


template <typename Tag, typename FilterPartnerType>
generic_linear_Pk<Tag, FilterPartnerType>::generic_linear_Pk(const boost::filesystem::path& p, const tree_Pk::database_type& d, const std::string& h)
  : path(p),
    container(d),
    md5_hash(h)
  {
  }


template <typename Tag, typename FilterPartnerType>
std::string generic_linear_Pk<Tag, FilterPartnerType>::hash(const boost::filesystem::path& p)
  {
    unsigned char result[MD5_DIGEST_LENGTH];
    
    std::ifstream in(p.string().c_str(), std::ios_base::in);
    
    if(!in)
      {
        in.close();
        return std::string();   // return empty hash
      }
    
    // read in file in 1k chunks
    constexpr unsigned int BUFFER_SIZE = 64*1024;
    std::unique_ptr<char> buffer(new char[BUFFER_SIZE]);
    
    MD5_CTX md5_context;
    MD5_Init(&md5_context);
    
    while(!in.eof())
      {
        in.read(buffer.get(), BUFFER_SIZE);
        MD5_Update(&md5_context, buffer.get(), in.gcount());
      }
    
    MD5_Final(result, &md5_context);
    
    in.close();
    
    std::ostringstream hash;
    for(unsigned char i : result)
    {
        hash << std::setw(2) << std::hex << static_cast<int>(i);
      }
    
    return hash.str();
  }


template <typename Tag, typename FilterPartnerType>
bool generic_linear_Pk<Tag, FilterPartnerType>::is_valid(const Mpc_units::energy& k, double bottom_clearance, double top_clearance) const
  {
    return this->container.is_valid(k, bottom_clearance, top_clearance);
  }


namespace boost
  {
    
    namespace serialization
      {
        
        template <typename Archive, typename Tag, typename FilterPartnerType>
        inline void save_construct_data(Archive& ar, const generic_linear_Pk<Tag, FilterPartnerType>* t, const unsigned int file_version)
          {
          }
        
        
        template <typename Archive, typename Tag, typename FilterPartnerType>
        inline void load_construct_data(Archive& ar, generic_linear_Pk<Tag, FilterPartnerType>* t, const unsigned int file_version)
          {
            // create an empty object with null contents; these contents will later be overwritte
            // by standard deserialization
            std::string path;
            std::string hash;
            tree_Pk::database_type db;
            ::new(t) generic_linear_Pk<Tag, FilterPartnerType>(path, db, hash);
          }
        
      }   // namespace serialization
    
  }   // namespace boost


// convenience type for initial and final linear power spectrum
typedef generic_linear_Pk< generic_linear_Pk_impl::InitialTag, initial_filtered_Pk > initial_Pk;
typedef generic_linear_Pk< generic_linear_Pk_impl::FinalTag, final_filtered_Pk > final_Pk;

// 'filterable_Pk' is an anonymous container type used for messaging over MPI
// initial_Pk and final_Pk are all convertible to filterable_Pk
typedef generic_linear_Pk< generic_linear_Pk_impl::FilterableTag, void > filterable_Pk;


// allow conversion of generic type to filterable type;
// a filterable power spectrum isn't used in calculations, only by the MPI backend to actually
// perform filtering
template <typename Tag, typename FilterPartnerType>
std::unique_ptr<filterable_Pk> make_filterable(const generic_linear_Pk<Tag, FilterPartnerType>& Pk)
  {
    return std::make_unique<filterable_Pk>(Pk.get_path(), Pk.get_db(), Pk.get_MD5_hash());
  };


#endif //LSSEFT_POWER_SPECTRUM_LINEAR_H
