//
// Created by David Seery on 12/08/2015.
// --@@
// Copyright (c) 2017 University of Sussex. All rights reserved.
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

#ifndef LSSEFT_RANGE_H
#define LSSEFT_RANGE_H

#include <iostream>
#include <fstream>
#include <sstream>


#include "range_detail/range.h"
#include "range_detail/stepping.h"
#include "range_detail/aggregation.h"

#include "units/Mpc_units.h"

#include "exceptions.h"
#include "localizations/messages.h"

#include "boost/filesystem/path.hpp"


//! utility function to build a range from an input file
template <typename Value>
aggregation_range<Value> load_range_from_file(const boost::filesystem::path& p, Value units)
  {
    aggregation_range<Value> r;  // construct empty aggregation range to hold result

    std::ifstream in(p.string(), std::ios::in);

    if(!in.good())
      {
        std::ostringstream msg;
        msg << ERROR_KMODES_FILE_NOT_READABLE_A << " " << p << " " << ERROR_KMODES_FILE_NOT_READABLE_B;
        throw runtime_exception(exception_type::runtime_error, msg.str());
      }

    for(std::string line; std::getline(in, line); )
      {
        std::stringstream line_stream(line);

        if(line.front() != '#')   // hash is a comment
          {
            double val;
            line_stream >> val;

            r += stepping_range<Value>(val, val, 0, units, spacing_type::linear);
          }
      }

    return r;
  }


#endif //LSSEFT_RANGE_H
