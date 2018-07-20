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

#include <iostream>
#include <fstream>

#include "boost/algorithm/string.hpp"

#include "FRW_model.h"

#include "localizations/messages.h"
#include "exceptions.h"


FRW_model::FRW_model(std::string nm, double om, double occ, double h_, Mpc_units::energy tc, double ne, double fb, double zs,
                     double zd, double ze, double Ac, double n, Mpc_units::energy kp)
  : name(nm),
    omega_m(om),
    omega_cc(occ),
    h(h_),
    T_CMB(tc),
    Neff(ne),
    f_baryon(fb),
    z_star(zs),
    z_drag(zd),
    z_eq(ze),
    A_curv(Ac),
    ns(n),
    k_piv(kp)
  {
  }


FRW_model::FRW_model(boost::filesystem::path p, double om, double occ, double h_, Mpc_units::energy tc, double ne,
                     double fb, double zs, double zd, double ze, double Ac, double n, Mpc_units::energy kp)
  : name(p.string()),
    omega_m(om),
    omega_cc(occ),
    h(h_),
    T_CMB(tc),
    Neff(ne),
    f_baryon(fb),
    z_star(zs),
    z_drag(zd),
    z_eq(ze),
    A_curv(Ac),
    ns(n),
    k_piv(kp)
  {
    std::ifstream in;
    in.open(p.string());

    if(!in.good())
      {
        std::ostringstream msg;
        msg << ERROR_PARAMETERS_FILE_NOT_READABLE_A << " " << p << " " << ERROR_PARAMETERS_FILE_NOT_READABLE_B;
        throw runtime_exception(exception_type::runtime_error, msg.str());
      }

      for(std::string line; std::getline(in, line); )
        {
          std::stringstream line_stream(line);

          if(line.front() != '#')   // hash is a comment
            {
              std::string label;
              double value;
              line_stream >> label >> value;

              // convert label to lower case
              boost::algorithm::to_lower(label);

              if(label == "h0" || label == "h")
                {
                  this->h = value;
                  std::cout << "Set h = " << this->h << '\n';
                }
              else if(label == "omegam" || label == "omega_m" || label == "omegam0" || label == "omega_m0")
                {
                  this->omega_m = value;
                  this->omega_cc = 1.0 - value;
                  std::cout << "Set Omega_m = " << this->omega_m << ", Omega_CC = " << this->omega_cc << '\n';
                }
              else if(label == "omegacc" || label == "omega_cc" || label == "omegacc0" || label == "omega_cc0")
                {
                  this->omega_cc = value;
                  this->omega_m = 1.0 - value;
                  std::cout << "Set Omega_m = " << this->omega_m << ", Omega_CC = " << this->omega_cc << '\n';
                }
              else if(label == "f" || label == "fbaryon" || label == "f_baryon")
                {
                  this->f_baryon = value;
                  std::cout << "Set f_baryon = " << this->f_baryon << '\n';
                }
              else
                {
                  std::ostringstream msg;
                  msg << ERROR_PARAMETERS_FILE_UNKNOWN_LABEL << " '" << label << "'";
                  throw runtime_exception(exception_type::runtime_error, msg.str());
                }
            }
        }
  }
