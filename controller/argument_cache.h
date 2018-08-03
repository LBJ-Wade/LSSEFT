//
// Created by David Seery on 31/07/2015.
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

#ifndef LSSEFT_ARGUMENT_CACHE_H
#define LSSEFT_ARGUMENT_CACHE_H


#include <string>


#include "boost/filesystem/operations.hpp"

#include "boost/serialization/serialization.hpp"
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
    
  public:

    //! get verbose output setting
    bool get_verbose() const { return(this->verbose); }

    //! set verbose output
    void set_verbose(bool a) { this->verbose = a; }

    //! get colourized output setting
    bool get_colour_output() const { return(this->colour_output); }

    //! set colourized output
    void set_colour_output(bool a) { this->verbose = a; }


    // INTERFACE -- DATABASE
    
  public:

    //! get database path
    const boost::filesystem::path& get_database_path() const { return this->database; }

    //! set database path
    void set_database_path(const std::string& p) { this->database = p; }

    //! determine if a database has been set
    bool is_database_set() const { return !this->database.empty(); }
    
    //! is network mode enables
    bool is_network_mode() const { return this->network_mode; }
    
    //! set network mode
    void set_network_mode(bool m) { this->network_mode = m; }
    
    
    // INTERFACE -- INITIAL AND FINAL POWER SPECTRA
    
  public:

    //! get initial linear power spectrum path
    const boost::filesystem::path& get_initial_powerspectrum_path() const { return(this->init_linear_Pk); }

    //! set initial linear power spectrum path
    void set_initial_powerspectrum_path(const std::string& p) { this->init_linear_Pk = p; }

    //! determine whether an initial linear power spectrum has been set
    bool is_initial_powerspectrum_set() const { return(!this->init_linear_Pk.empty()); }
    
    //! get final linear power spectrum path
    const boost::filesystem::path& get_final_powerspectrum_path() const { return(this->final_linear_Pk); }
    
    //! set final linear linear power spectrum path
    void set_final_powerspectrum_path(const std::string& p) { this->final_linear_Pk = p; }
    
    //! determine whether a final linear power spectrum has been set
    bool is_final_powerspectrum_set() const { return(!this->final_linear_Pk.empty()); }


    // INTERFACE -- PARAMETERS FILE

  public:

    //! set parameters file path
    void set_parameter_file_path(const std::string& p) { this->parameter_file = p; }

    //! determine whether a parameters file path has been set
    bool is_parameter_file_set() const { return !this->parameter_file.empty(); }

    //! get parameter file path
    const boost::filesystem::path& get_parameter_file_path() const { return this->parameter_file; }


    // INTERFACE -- KMODES FILE

  public:

    //! set kmodes file path
    void set_kmodes_file_path(const std::string& p) { this->kmodes_file = p; }

    //! determine whether a kmodes file path has been set
    bool is_kmodes_file_set() const { return !this->kmodes_file.empty(); }

    //! get kmodes file path
    const boost::filesystem::path& get_kmodes_file_path() const { return this->kmodes_file; }

    
    // INTERFACE -- CALCULATION PARAMETERS
    
  public:
    
    //! query whether we are using EdS approximations to the growth functions
    bool use_EdS() const { return this->EdS_mode; }
    
    //! set EdS mode
    void set_EdS_mode(bool m) { this->EdS_mode = m; }


    // INTERNAL DATA

  private:

    //! generate verbose output?
    bool verbose;

    //! generate colourized output?
    bool colour_output;
    
    //! use Einstein-de Sitter approximations to growth functions?
    bool EdS_mode;
    
    //! should we use network mode, ie. disable write-ahead log?
    bool network_mode;

    //! database path
    boost::filesystem::path database;

    //! path for initial linear power spectrum
    boost::filesystem::path init_linear_Pk;
    
    //! path for final linear power spectrum
    boost::filesystem::path final_linear_Pk;

    //! path for parameters file
    boost::filesystem::path parameter_file;

    //! path for kmodes file
    boost::filesystem::path kmodes_file;

    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & verbose;
        ar & colour_output;
        ar & EdS_mode;
        ar & network_mode;
        ar & database;
        ar & init_linear_Pk;
        ar & final_linear_Pk;
        ar & parameter_file;
        ar & kmodes_file;
      }

  };


#endif //LSSEFT_ARGUMENT_CACHE_H
