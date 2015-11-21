//
// Created by David Seery on 29/10/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

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
        db.add_record(k_in_eV, Pk_in_inv_eV3);
      }

    in.close();
  }
