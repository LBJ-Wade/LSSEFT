//
// Created by David Seery on 29/10/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_GENERIC_POWER_SPECTRUM_H
#define LSSEFT_GENERIC_POWER_SPECTRUM_H

#include <memory>
#include <fstream>

#include "database/Pk_database.h"

#include "exceptions.h"
#include "localizations/messages.h"

#include "boost/filesystem/operations.hpp"
#include "boost/serialization/serialization.hpp"

#include "SPLINTER/datatable.h"
#include "SPLINTER/bspline.h"
#include "SPLINTER/bsplinebuilder.h"


constexpr double SPLINE_PK_DEFAULT_TOP_CLEARANCE = 0.9;
constexpr double SPLINE_PK_DEFAULT_BOTTOM_CLEARANCE = 1.1;


template <typename Tag, typename Dimension, bool protect=false>
class generic_Pk
  {
    
    // TYPEDEFS
    
  public:
    
    typedef Pk_database<Dimension> database_type;

    
    // CONSTRUCTOR, DESTRUCTOR

  public:
    
    //! constructor -- from directly-supplied database
    generic_Pk(const Pk_database<Dimension>& db);

    //! destructor is default
    ~generic_Pk() = default;


    // DATABASE SERVICES

  public:

    //! get power spectrum database
    const Pk_database<Dimension>& get_db() const { return(this->database); }
    
    //! ask whether it is valid to evaluate the spline at a given k-value
    bool is_valid(const Mpc_units::energy& k,
                  double bottom_clearance = SPLINE_PK_DEFAULT_BOTTOM_CLEARANCE,
                  double top_clearance = SPLINE_PK_DEFAULT_TOP_CLEARANCE) const;
    
    
    // EVALUATION
    
  public:

    //! evaluate spline, using k-range protection if enabled
    template <bool P=protect, typename std::enable_if<P>::type* = nullptr>
    Dimension operator()(const Mpc_units::energy& k) const;

    template <bool P=protect, typename std::enable_if<!P>::type* = nullptr>
    Dimension operator()(const Mpc_units::energy& k) const;
    
  private:
    
    Dimension evaluate(const Mpc_units::energy& k) const;


    // INTERNAL API

  private:

    //! recalculate spline approximant
    void recalculate_spline();


    // INTERNAL DATA

  private:

    //! power spectrum
    Pk_database<Dimension> database;

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
generic_Pk<Tag, Dimension, protect>::generic_Pk(const Pk_database<Dimension>& db)
  : database(db)
  {
    this->recalculate_spline();
  }


template <typename Tag, typename Dimension, bool protect>
void generic_Pk<Tag, Dimension, protect>::recalculate_spline()
  {
    this->table.release();
    this->spline.release();
    
    this->table = std::make_unique<SPLINTER::DataTable>();
    
    for(typename Pk_database<Dimension>::const_record_iterator t = this->database.record_begin(); t != this->database.record_end(); ++t)
      {
        this->table->addSample(t->get_wavenumber() * Mpc_units::Mpc, t->get_Pk() / Pk_database_impl::DimensionTraits<Dimension>().unit());
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
template <bool P, typename std::enable_if<P>::type*>
Dimension generic_Pk<Tag, Dimension, protect>::operator()(const Mpc_units::energy& k) const
  {
    if(k > SPLINE_PK_DEFAULT_TOP_CLEARANCE * this->database.get_k_max())
      {
        std::ostringstream msg;
        msg << ERROR_POWERSPECTRUM_SPLINE_TOO_BIG << " (k = " << k * Mpc_units::Mpc << " h/Mpc, k_max = "
            << this->database.get_k_max() * Mpc_units::Mpc << " h/Mpc)";
        throw std::overflow_error(msg.str());
      }

    if(k < SPLINE_PK_DEFAULT_BOTTOM_CLEARANCE * this->database.get_k_min())
      {
        std::ostringstream msg;
        msg << ERROR_POWERSPECTRUM_SPLINE_TOO_SMALL << " (k = " << k * Mpc_units::Mpc << " h/Mpc, k_min = "
            << this->database.get_k_min() * Mpc_units::Mpc << " h/Mpc)";
        throw std::overflow_error(msg.str());
      }
    
    return this->evaluate(k);
  }


template <typename Tag, typename Dimension, bool protect>
template <bool P, typename std::enable_if<!P>::type*>
Dimension generic_Pk<Tag, Dimension, protect>::operator()(const Mpc_units::energy& k) const
  {
    return this->evaluate(k);
  }


template <typename Tag, typename Dimension, bool protect>
Dimension generic_Pk<Tag, Dimension, protect>::evaluate(const Mpc_units::energy& k) const
  {
    SPLINTER::DenseVector x(1);
    x(0) = k * Mpc_units::Mpc;
    
    return(this->spline->eval(x) * Pk_database_impl::DimensionTraits<Dimension>().unit());
  }


template <typename Tag, typename Dimension, bool protect>
bool generic_Pk<Tag, Dimension, protect>::is_valid(const Mpc_units::energy& k, double bottom_clearance, double top_clearance) const
  {
    if(this->database.size() == 0) return false;
    
    if(bottom_clearance < SPLINE_PK_DEFAULT_BOTTOM_CLEARANCE) bottom_clearance = SPLINE_PK_DEFAULT_BOTTOM_CLEARANCE;
    if(top_clearance > SPLINE_PK_DEFAULT_TOP_CLEARANCE)       top_clearance = SPLINE_PK_DEFAULT_TOP_CLEARANCE;
    
    // TODO: detect max-k < min-k after accounting for clearances
    
    if(k > top_clearance * this->database.get_k_max()) return false;
    if(k < bottom_clearance * this->database.get_k_min()) return false;
    
    return true;
  }


namespace boost
  {

    namespace serialization
      {

        template <typename Archive, typename Tag, typename Dimension, bool protect>
        inline void save_construct_data(Archive& ar, const generic_Pk<Tag, Dimension, protect>* t, const unsigned int file_version)
          {
            ar << t->get_db();
          }


        template <typename Archive, typename Tag, typename Dimension, bool protect>
        inline void load_construct_data(Archive& ar, generic_Pk<Tag, Dimension, protect>* t, const unsigned int file_version)
          {
            Pk_database<Dimension> db;
            ar >> db;

            ::new(t) generic_Pk<Tag, Dimension, protect>(db);
          }

      }   // namespace serialization

  }   // namespace boost


#endif //LSSEFT_GENERIC_POWER_SPECTRUM_H