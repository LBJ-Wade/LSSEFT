//
// Created by David Seery on 29/10/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_TREE_POWER_SPECTRUM_H
#define LSSEFT_TREE_POWER_SPECTRUM_H

#include <memory>
#include <fstream>

#include "database/powerspectrum_database.h"

#include "exceptions.h"
#include "localizations/messages.h"

#include "boost/filesystem/operations.hpp"
#include "boost/serialization/serialization.hpp"

#include "SPLINTER/datatable.h"
#include "SPLINTER/bspline.h"
#include "SPLINTER/bsplinebuilder.h"


namespace power_spectrum_impl
  {
    
    template <int m>
    struct PowerSpectrumTag
      {
        enum { SpectrumClass=m };
      };
    
    typedef PowerSpectrumTag<0> TreeTag;
    
    
    template<typename Dimension> struct DimensionTraits;
    
    template <>
    struct DimensionTraits<Mpc_units::inverse_energy3>
      {
        constexpr Mpc_units::inverse_energy3 unit() const { return Mpc_units::Mpc3; }
      };
    
    template <>
    struct DimensionTraits<Mpc_units::inverse_energy>
      {
        constexpr Mpc_units::inverse_energy unit() const { return Mpc_units::Mpc; }
      };
    
    
  }   // namespace power_spectrum_impl


template <typename Tag, typename Dimension, bool protect=false>
class power_spectrum
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor -- read in from a file in CAMB format
    power_spectrum(const boost::filesystem::path& p);

    //! constructor -- from directly-supplied database
    power_spectrum(const powerspectrum_database<Dimension>& db);

    //! destructor is default
    ~power_spectrum() = default;


    // INTERFACE

  public:

    //! get power spectrum database
    const powerspectrum_database<Dimension>& get_db() const { return(this->database); }

    //! evaluate spline
    Dimension operator()(const Mpc_units::energy& k) const;


    // INTERNAL API

  private:

    //! ingest CAMB-format powerspectrum file
    void ingest_CAMB(const boost::filesystem::path& p);

    //! recalculate spline approximant
    void recalculate_spline();


    // INTERNAL DATA

  private:

    //! power spectrum
    powerspectrum_database<Dimension> database;

    //! splines representing power spectrum
    std::unique_ptr<SPLINTER::DataTable> table;

    std::unique_ptr<SPLINTER::BSpline> spline;


    // enable boost::serialization support and hence automated packing for transmission over MPI
    friend class boost::serialization::access;


    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
      }

  };


template <typename Tag, typename Dimension, bool protect>
power_spectrum<Tag, Dimension, protect>::power_spectrum(const boost::filesystem::path& p)
  {
    if(!boost::filesystem::exists(p))
      {
        std::ostringstream msg;
        msg << ERROR_POWERSPECTRUM_FILE_NOT_THERE_A << " " << p << " " << ERROR_POWERSPECTRUM_FILE_NOT_THERE_B;
        throw runtime_exception(exception_type::runtime_error, msg.str());
      }
    
    this->ingest_CAMB(p);
    this->recalculate_spline();
  }


template <typename Tag, typename Dimension, bool protect>
power_spectrum<Tag, Dimension, protect>::power_spectrum(const powerspectrum_database<Dimension>& db)
  : database(db)
  {
    this->recalculate_spline();
  }


template <typename Tag, typename Dimension, bool protect>
void power_spectrum<Tag, Dimension, protect>::ingest_CAMB(const boost::filesystem::path& p)
  {
    std::ifstream in;
    in.open(p.string());
    
    if(!in.good())
      {
        std::ostringstream msg;
        msg << ERROR_POWERSPECTRUM_FILE_NOT_READABLE_A << " " << p << " " << ERROR_POWERSPECTRUM_FILE_NOT_READABLE_B;
        throw runtime_exception(exception_type::runtime_error, msg.str());
      }
    
    for(std::string line; std::getline(in, line); )
      {
        std::stringstream line_stream(line);
        
        double _k, _Pk;
        line_stream >> _k >> _Pk;
        
        Mpc_units::energy k = _k / Mpc_units::Mpc;
        Mpc_units::inverse_energy3 Pk = _Pk * power_spectrum_impl::DimensionTraits<Dimension>().unit();
        this->database.add_record(k, Pk);
      }
    
    in.close();
  }


template <typename Tag, typename Dimension, bool protect>
void power_spectrum<Tag, Dimension, protect>::recalculate_spline()
  {
    this->table.release();
    this->spline.release();
    
    this->table = std::make_unique<SPLINTER::DataTable>();
    
    for(typename powerspectrum_database<Dimension>::const_record_iterator t = this->database.record_begin(); t != this->database.record_end(); ++t)
      {
        this->table->addSample(t->get_wavenumber() * Mpc_units::Mpc, t->get_Pk() / power_spectrum_impl::DimensionTraits<Dimension>().unit());
      }
    
    try
      {
        this->spline = std::make_unique<SPLINTER::BSpline>(SPLINTER::BSpline::Builder(*this->table).degree(3).build());
      }
    catch(SPLINTER::Exception& xe)
      {
        throw runtime_exception(exception_type::spline_error, xe.what());
      }
  }


template <typename Tag, typename Dimension, bool protect>
Dimension power_spectrum<Tag, Dimension, protect>::operator()(const Mpc_units::energy& k) const
  {
    constexpr double TEN_PERCENT_UPPER_CLEARANCE = 0.9;
    constexpr double TEN_PERCENT_LOWER_CLEARANCE = 1.1;
    
    if(protect)
      {
        if(k > TEN_PERCENT_UPPER_CLEARANCE*this->database.get_k_max())
          {
            std::ostringstream msg;
            msg << ERROR_POWERSPECTRUM_SPLINE_TOO_BIG << " (k = " << k * Mpc_units::Mpc << " h/Mpc, k_max = " << this->database.get_k_max() * Mpc_units::Mpc << " h/Mpc)";
            throw std::overflow_error(msg.str());
          }
        if(k < TEN_PERCENT_LOWER_CLEARANCE*this->database.get_k_min())
          {
            std::ostringstream msg;
            msg << ERROR_POWERSPECTRUM_SPLINE_TOO_SMALL << " (k = " << k * Mpc_units::Mpc << " h/Mpc, k_min = " << this->database.get_k_min() * Mpc_units::Mpc << " h/Mpc)";
            throw std::overflow_error(msg.str());
          }
      }
    
    SPLINTER::DenseVector x(1);
    x(0) = k * Mpc_units::Mpc;
    
    return(this->spline->eval(x) * power_spectrum_impl::DimensionTraits<Dimension>().unit());
  }


namespace boost
  {

    namespace serialization
      {

        template <typename Archive, typename Tag, typename Dimension, bool protect>
        inline void save_construct_data(Archive& ar, const power_spectrum<Tag, Dimension, protect>* t, const unsigned int file_version)
          {
            ar << t->get_db();
          }


        template <typename Archive, typename Tag, typename Dimension, bool protect>
        inline void load_construct_data(Archive& ar, power_spectrum<Tag, Dimension, protect>* t, const unsigned int file_version)
          {
            powerspectrum_database<Dimension> db;
            ar >> db;

            ::new(t) power_spectrum<Tag, Dimension, protect>(db);
          }

      }   // namespace serialization

  }   // namespace boost


// set up convenience type for tree-level power spectrum
typedef power_spectrum< power_spectrum_impl::TreeTag, Mpc_units::inverse_energy3, true > tree_power_spectrum;


#endif //LSSEFT_TREE_POWER_SPECTRUM_H
