//
// Created by David Seery on 31/07/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_ARGUMENT_CACHE_H
#define LSSEFT_ARGUMENT_CACHE_H


#include <string>


#include "boost/filesystem/operations.hpp"

#include "boost/serialization/string.hpp"
#include "boost/serialization/list.hpp"


class argument_cache
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor sets default values for all options
    argument_cache();

    //! destructor is default
    ~argument_cache() = default;


    // INTERFACE -- CONFIGURATION OPTIONS

    //! get verbose output setting
    bool get_verbose() const { return(this->verbose); }

    //! set verbose output
    void set_verbose(bool a) { this->verbose = a; }

    //! get colourized output setting
    bool get_colour_output() const { return(this->colour_output); }

    //! set colourized output
    void set_colour_output(bool a) { this->verbose = a; }


    // INTERFACE -- DATABASE

    //! get database path
    const boost::filesystem::path& get_database_path() const { return(this->database); }

    //! set database path
    void set_database_path(const std::string& p) { this->database = p; }

    //! determine if a database has been set
    bool get_database_set() const { return(!this->database.empty()); }


    // INTERNAL DATA

  private:

    //! generate verbose output?
    bool verbose;

    //! generate colourized output?
    bool colour_output;

    //! database path
    boost::filesystem::path database;

    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & verbose;
        ar & colour_output;
        ar & database;
      }

  };


#endif //LSSEFT_ARGUMENT_CACHE_H
