//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include "powerspectrum_record.h"


Pk_record::Pk_record(const eV_units::energy& _k, const eV_units::inverse_energy3& _Pk)
  : k(_k),
    Pk(_Pk)
  {
  }
