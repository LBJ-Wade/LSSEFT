//
// Created by David Seery on 29/10/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_TREE_POWER_SPECTRUM_H
#define LSSEFT_TREE_POWER_SPECTRUM_H

#include <memory>

#include "database/powerspectrum_database.h"

#include "boost/filesystem/operations.hpp"
#include "boost/serialization/serialization.hpp"

#include "datatable.h"
#include "bsplineapproximant.h"


class tree_power_spectrum
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor -- read in from a file in CAMB format
    tree_power_spectrum(const boost::filesystem::path& p);

    //! constructor -- from directly-supplied database
    tree_power_spectrum(const powerspectrum_database& db);

    //! destructor is default
    ~tree_power_spectrum() = default;


    // INTERFACE

  public:

    //! get power spectrum database
    const powerspectrum_database& get_db() const { return(this->database); }

    //! evaluate spline
    Mpc_units::inverse_energy3 operator()(const Mpc_units::energy& k) const;


    // INTERNAL API

  private:

    //! ingest CAMB-format powerspectrum file
    void ingest_CAMB(const boost::filesystem::path& p);

    //! recalculate spline approximant
    void recalculate_spline();


    // INTERNAL DATA

  private:

    //! power spectrum
    powerspectrum_database database;

    //! splines representing power spectrum
    std::unique_ptr<SPLINTER::DataTable> table;

    std::unique_ptr<SPLINTER::BSplineApproximant> spline;


    // enable boost::serialization support and hence automated packing for transmission over MPI
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
        inline void save_construct_data(Archive& ar, const tree_power_spectrum* t, const unsigned int file_version)
          {
            ar << t->get_db();
          }


        template <typename Archive>
        inline void load_construct_data(Archive& ar, tree_power_spectrum* t, const unsigned int file_version)
          {
            powerspectrum_database db;
            ar >> db;

            ::new(t) tree_power_spectrum(db);
          }

      }   // namespace serialization

  }   // namespace boost


#endif //LSSEFT_TREE_POWER_SPECTRUM_H
