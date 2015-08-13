//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_TOKENS_H
#define LSSEFT_TOKENS_H


#include "boost/serialization/serialization.hpp"
#include "boost/serialization/base_object.hpp"
#include "boost/serialization/assume_abstract.hpp"


// forward declare functions which will be friended in class declarations

class generic_token;

bool operator<(const generic_token& a, const generic_token& b);
bool operator==(const generic_token& a, const generic_token& b);


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
class redshift_token: public generic_token
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    redshift_token(unsigned int i);

    //! destructor is default
    virtual ~redshift_token() = default;


  private:


    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & boost::serialization::base_object<generic_token>(*this);
      }

  };


//! token representing a wavenumber
class wavenumber_token: public generic_token
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    wavenumber_token(unsigned int i);

    //! destructor is default
    virtual ~wavenumber_token() = default;


  private:


    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & boost::serialization::base_object<generic_token>(*this);
      }

  };



#endif //LSSEFT_TOKENS_H
