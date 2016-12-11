//
// Created by David Seery on 06/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#include "spline.h"


wiggle_Pk_raw_adapter::wiggle_Pk_raw_adapter(const wiggle_Pk& w)
  : container(w)
  {
  }


Mpc_units::inverse_energy3 wiggle_Pk_raw_adapter::operator()(const Mpc_units::energy& k) const
  {
    return this->container.Pk_raw(k);
  }


wiggle_Pk_wiggle_adapter::wiggle_Pk_wiggle_adapter(const wiggle_Pk& w)
  : container(w)
  {
  }


Mpc_units::inverse_energy3 wiggle_Pk_wiggle_adapter::operator()(const Mpc_units::energy& k) const
  {
    return this->container.Pk_wiggle(k);
  }


wiggle_Pk_nowiggle_adapter::wiggle_Pk_nowiggle_adapter(const wiggle_Pk& w)
  : container(w)
  {
  }


Mpc_units::inverse_energy3 wiggle_Pk_nowiggle_adapter::operator()(const Mpc_units::energy& k) const
  {
    return this->container.Pk_nowiggle(k);
  }
