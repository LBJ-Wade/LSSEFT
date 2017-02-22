//
// Created by David Seery on 11/08/2015.
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

#ifndef LSSEFT_TOKENS_H
#define LSSEFT_TOKENS_H


#include <functional>

#include "sqlite3_detail/sqlite3_policy.h"

#include "boost/serialization/serialization.hpp"
#include "boost/serialization/base_object.hpp"
#include "boost/serialization/assume_abstract.hpp"


// forward declare functions which will be friended in class declarations

class generic_token;

bool operator<(const generic_token& a, const generic_token& b);
bool operator==(const generic_token& a, const generic_token& b);


template <typename Token>
const std::string& tokenization_table(const sqlite3_policy& policy);


//! base class: a generic token
class generic_token
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    generic_token(unsigned int i);

    //! destructor is default
    virtual ~generic_token() = default;


    // INTERFACE

  public:

    //! return id associated with this token
    virtual unsigned int get_id() const { return(this->id); }


    // OVERLOAD COMPARISON OPERATORS

  public:

    //! overload ordering
    friend bool operator<(const generic_token& a, const generic_token& b);

    //! overload equality
    friend bool operator==(const generic_token& a, const generic_token& b);


    // INTERNAL DATA

  private:

    unsigned int id;


    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & id;
      }

  };


BOOST_SERIALIZATION_ASSUME_ABSTRACT(generic_token)


namespace boost
  {
    
    // generic_token has no default constructor, so we have to specialize load/store methods
    
    namespace serialization
      {

        template <typename Archive>
        inline void save_construct_data(Archive& ar, const generic_token* t, const unsigned int file_version)
          {
          }
        
        
        template <typename Archive>
        inline void load_construct_data(Archive& ar, generic_token* t, const unsigned int file_version)
          {
            // invoke in-place constructor with blank id; will be overwritten later during deserialization
            ::new(t) generic_token(0);
          }
        
      }   // namespace serialization
    
  }   // namespace boost


namespace std
  {
    
    template<>
    class hash<generic_token>
      {
        
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! constructor is default
        hash() = default;
        
        //! destructor is default
        ~hash() = default;
        
        
        // IMPLEMENTATION
        
      public:
        
        //! hash function
        size_t operator()(const generic_token& tok) const
          {
            std::hash<unsigned int> hasher;

            return hasher(tok.get_id());
          }
        
      };
    
  }


//! token representing an FRW model
class FRW_model_token: public generic_token
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    FRW_model_token(unsigned int i);

    //! destructor is default
    virtual ~FRW_model_token() = default;


  private:

    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & boost::serialization::base_object<generic_token>(*this);
      }

  };


//! token representing a redshift
class z_token: public generic_token
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    z_token(unsigned int i);

    //! destructor is default
    virtual ~z_token() = default;


  private:

    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & boost::serialization::base_object<generic_token>(*this);
      }

  };


// specialize tokenization for this
template <>
inline const std::string& tokenization_table<z_token>(const sqlite3_policy& policy)
  {
    return(policy.redshift_config_table());
  }


//! token representing a k-value at which we sample the power spectrum
class k_token: public generic_token
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    k_token(unsigned int i);

    //! destructor is default
    virtual ~k_token() = default;


  private:

    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & boost::serialization::base_object<generic_token>(*this);
      }

  };


// specialize tokenization for this
template <>
inline const std::string& tokenization_table<k_token>(const sqlite3_policy& policy)
  {
    return(policy.wavenumber_config_table());
  }



//! token representing a k-value corresponding to an IR cutoff
class IR_cutoff_token: public generic_token
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    IR_cutoff_token(unsigned int i);

    //! destructor is default
    virtual ~IR_cutoff_token() = default;


  private:

    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & boost::serialization::base_object<generic_token>(*this);
      }

  };


// specialize tokenization for this
template <>
inline const std::string& tokenization_table<IR_cutoff_token>(const sqlite3_policy& policy)
  {
    return(policy.IR_config_table());
  }


//! token representing a k-value corresponding to a UV cutoff
class UV_cutoff_token: public generic_token
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    UV_cutoff_token(unsigned int i);

    //! destructor is default
    virtual ~UV_cutoff_token() = default;


  private:

    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & boost::serialization::base_object<generic_token>(*this);
      }

  };


// specialize tokenization for this
template <>
inline const std::string& tokenization_table<UV_cutoff_token>(const sqlite3_policy& policy)
  {
    return(policy.UV_config_table());
  }


//! token representing a k-value associated with an IR resummation scale
class IR_resum_token: public generic_token
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor
    IR_resum_token(unsigned int i);
    
    //! destructor is default
    virtual ~IR_resum_token() = default;
    
    
  private:
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & boost::serialization::base_object<generic_token>(*this);
      }
    
  };


// specialize tokenization for IR-resummation-scale tokens
template <>
inline const std::string& tokenization_table<IR_resum_token>(const sqlite3_policy& policy)
  {
    return(policy.IR_resum_config_table());
  }


//! token representing a linear power spectrum P(k)
class linear_Pk_token: public generic_token
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! constructor
    linear_Pk_token(unsigned int i);
    
    //! destructor is default
    virtual ~linear_Pk_token() = default;
  
  
  private:
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & boost::serialization::base_object<generic_token>(*this);
      }
    
  };


// specialize tokenization for linear P(k) tokens
template <>
inline const std::string& tokenization_table<linear_Pk_token>(const sqlite3_policy& policy)
  {
    return(policy.Pk_linear_config_table());
  }


namespace boost
  {
    
    // generic_token has no default constructor, so we have to specialize load/store methods
    
    namespace serialization
      {
        
        template <typename Archive>
        inline void save_construct_data(Archive& ar, const FRW_model_token* t, const unsigned int file_version)
          {
          }
        
        
        template <typename Archive>
        inline void load_construct_data(Archive& ar, FRW_model_token* t, const unsigned int file_version)
          {
            // invoke in-place constructor with blank id; will be overwritten later during deserialization
            ::new(t) FRW_model_token(0);
          }
    
    
        template <typename Archive>
        inline void save_construct_data(Archive& ar, const z_token* t, const unsigned int file_version)
          {
          }
    
    
        template <typename Archive>
        inline void load_construct_data(Archive& ar, z_token* t, const unsigned int file_version)
          {
            // invoke in-place constructor with blank id; will be overwritten later during deserialization
            ::new(t) z_token(0);
          }
    
    
        template <typename Archive>
        inline void save_construct_data(Archive& ar, const k_token* t, const unsigned int file_version)
          {
          }
    
    
        template <typename Archive>
        inline void load_construct_data(Archive& ar, k_token* t, const unsigned int file_version)
          {
            // invoke in-place constructor with blank id; will be overwritten later during deserialization
            ::new(t) k_token(0);
          }
    
    
        template <typename Archive>
        inline void save_construct_data(Archive& ar, const IR_cutoff_token* t, const unsigned int file_version)
          {
          }
    
    
        template <typename Archive>
        inline void load_construct_data(Archive& ar, IR_cutoff_token* t, const unsigned int file_version)
          {
            // invoke in-place constructor with blank id; will be overwritten later during deserialization
            ::new(t) IR_cutoff_token(0);
          }
    
    
        template <typename Archive>
        inline void save_construct_data(Archive& ar, const UV_cutoff_token* t, const unsigned int file_version)
          {
          }
    
    
        template <typename Archive>
        inline void load_construct_data(Archive& ar, UV_cutoff_token* t, const unsigned int file_version)
          {
            // invoke in-place constructor with blank id; will be overwritten later during deserialization
            ::new(t) UV_cutoff_token(0);
          }
    
    
        template <typename Archive>
        inline void save_construct_data(Archive& ar, const IR_resum_token* t, const unsigned int file_version)
          {
          }
    
    
        template <typename Archive>
        inline void load_construct_data(Archive& ar, IR_resum_token* t, const unsigned int file_version)
          {
            // invoke in-place constructor with blank id; will be overwritten later during deserialization
            ::new(t) IR_resum_token(0);
          }
    
    
        template <typename Archive>
        inline void save_construct_data(Archive& ar, const linear_Pk_token* t, const unsigned int file_version)
          {
          }
    
    
        template <typename Archive>
        inline void load_construct_data(Archive& ar, linear_Pk_token* t, const unsigned int file_version)
          {
            // invoke in-place constructor with blank id; will be overwritten later during deserialization
            ::new(t) linear_Pk_token(0);
          }
    
    
      }   // namespace serialization
    
  }   // namespace boost


namespace std
  {
    
    template<>
    class hash<FRW_model_token>
      {
        
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! constructor is default
        hash() = default;
        
        //! destrucotr is default
        ~hash() = default;
        
        
        // IMPLEMENTATION
        
      public:
        
        //! hash function
        size_t operator()(const FRW_model_token& tok)
          {
            std::hash<generic_token> hasher;
            return hasher(tok);
          }
        
      };
    
    
    template<>
    class hash<z_token>
      {
        
        // CONSTRUCTOR, DESTRUCTOR
      
      public:
        
        //! constructor is default
        hash() = default;
        
        //! destrucotr is default
        ~hash() = default;
        
        
        // IMPLEMENTATION
      
      public:
        
        //! hash function
        size_t operator()(const z_token& tok) const
          {
            std::hash<generic_token> hasher;
            return hasher(tok);
          }
        
      };
    
    
    template<>
    class hash<k_token>
      {
        
        // CONSTRUCTOR, DESTRUCTOR
      
      public:
        
        //! constructor is default
        hash() = default;
        
        //! destrucotr is default
        ~hash() = default;
        
        
        // IMPLEMENTATION
      
      public:
        
        //! hash function
        size_t operator()(const k_token& tok) const
          {
            std::hash<generic_token> hasher;
            return hasher(tok);
          }
        
      };
    
    
    template<>
    class hash<IR_cutoff_token>
      {
        
        // CONSTRUCTOR, DESTRUCTOR
      
      public:
        
        //! constructor is default
        hash() = default;
        
        //! destrucotr is default
        ~hash() = default;
        
        
        // IMPLEMENTATION
      
      public:
        
        //! hash function
        size_t operator()(const IR_cutoff_token& tok) const
          {
            std::hash<generic_token> hasher;
            return hasher(tok);
          }
        
      };
    
    
    template<>
    class hash<UV_cutoff_token>
      {
        
        // CONSTRUCTOR, DESTRUCTOR
      
      public:
        
        //! constructor is default
        hash() = default;
        
        //! destrucotr is default
        ~hash() = default;
        
        
        // IMPLEMENTATION
      
      public:
        
        //! hash function
        size_t operator()(const UV_cutoff_token& tok) const
          {
            std::hash<generic_token> hasher;
            return hasher(tok);
          }
        
      };
    
    
    template<>
    class hash<IR_resum_token>
      {
        
        // CONSTRUCTOR, DESTRUCTOR
      
      public:
        
        //! constructor is default
        hash() = default;
        
        //! destrucotr is default
        ~hash() = default;
        
        
        // IMPLEMENTATION
      
      public:
        
        //! hash function
        size_t operator()(const IR_resum_token& tok) const
          {
            std::hash<generic_token> hasher;
            return hasher(tok);
          }
        
      };
    
    
    template<>
    class hash<linear_Pk_token>
      {
        
        // CONSTRUCTOR, DESTRUCTOR
      
      public:
        
        //! constructor is default
        hash() = default;
        
        //! destrucotr is default
        ~hash() = default;
        
        
        // IMPLEMENTATION
      
      public:
        
        //! hash function
        size_t operator()(const linear_Pk_token& tok) const
          {
            std::hash<generic_token> hasher;
            return hasher(tok);
          }
        
      };
    
  }   // namespace std



#endif //LSSEFT_TOKENS_H
