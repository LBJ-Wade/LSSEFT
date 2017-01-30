//
// Created by David Seery on 17/08/2015.
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

#include "transfer_function.h"


transfer_function::transfer_function(const Mpc_units::energy& _k, const k_token& t, std::shared_ptr<z_database> z)
  : k(_k),
    token(t),
    z_db(std::move(z))
  {
    // if we were passed a non-null redshift database, set up containers
    // (perhaps should disallow construction with a null database?)
    if(z_db)
      {
        delta_m.reset(new std::vector<double>);
        delta_r.reset(new std::vector<double>);
        theta_m.reset(new std::vector<double>);
        theta_r.reset(new std::vector<double>);
        Phi.reset(new std::vector<double>);

        delta_m->reserve(z_db->size());
        delta_r->reserve(z_db->size());
        theta_m->reserve(z_db->size());
        theta_r->reserve(z_db->size());
        Phi->reserve(z_db->size());
      }
  }


void transfer_function::set_integration_metadata(boost::timer::nanosecond_type t, size_t s)
  {
    this->integration_time = t;
    this->steps = s;
  }


void transfer_function::push_back(double dm, double dr, double tm, double tr, double P)
  {
    this->delta_m->push_back(dm);
    this->delta_r->push_back(dr);
    this->theta_m->push_back(tm);
    this->theta_r->push_back(tr);
    this->Phi->push_back(P);
  }
