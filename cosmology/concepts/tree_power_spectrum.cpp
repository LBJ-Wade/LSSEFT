//
// Created by David Seery on 29/10/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>

#include "exceptions.h"
#include "localizations/messages.h"

#include "tree_power_spectrum.h"

#include "datatable.h"
#include "bsplineapproximant.h"


tree_power_spectrum::tree_power_spectrum(const boost::filesystem::path& p)
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


tree_power_spectrum::tree_power_spectrum(const powerspectrum_database& db)
  : database(db)
  {
    this->recalculate_spline();
  }


void tree_power_spectrum::ingest_CAMB(const boost::filesystem::path& p)
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

        double k, Pk;
        line_stream >> k >> Pk;

        eV_units::energy k_in_eV = k / eV_units::Mpc;
        eV_units::inverse_energy3 Pk_in_inv_eV3 = Pk * eV_units::Mpc3;
        this->database.add_record(k_in_eV, Pk_in_inv_eV3);
      }

    in.close();
 }


void tree_power_spectrum::recalculate_spline()
  {
    this->table.release();
    this->spline.release();

    this->table = std::make_unique<SPLINTER::DataTable>();

    for(powerspectrum_database::const_record_iterator t = this->database.record_begin(); t != this->database.record_end(); ++t)
      {
        this->table->addSample(static_cast<double>(t->get_wavenumber()), static_cast<double>(t->get_Pk()));
      }

    this->spline = std::make_unique<SPLINTER::BSplineApproximant>(*this->table, SPLINTER::BSplineType::CUBIC);
  }


eV_units::inverse_energy3 tree_power_spectrum::operator()(const eV_units::energy& k) const
  {
    if(k > 0.9*this->database.get_k_max())
      {
        std::ostringstream msg;
        msg << ERROR_POWERSPECTRUM_SPLINE_TOO_BIG << " (k = " << k * eV_units::Mpc << " h/Mpc, k_max = " << this->database.get_k_max() * eV_units::Mpc << " h/Mpc)";
        throw std::overflow_error(msg.str());
      }
    if(k < 1.1*this->database.get_k_min())
      {
        std::ostringstream msg;
        msg << ERROR_POWERSPECTRUM_SPLINE_TOO_SMALL << " (k = " << k * eV_units::Mpc << " h/Mpc, k_min = " << this->database.get_k_min() * eV_units::Mpc << " h/Mpc)";
        throw std::overflow_error(msg.str());
      }

    SPLINTER::DenseVector x(1);
    x(0) = static_cast<double>(k);

    return eV_units::inverse_energy3(this->spline->eval(x));
  }
