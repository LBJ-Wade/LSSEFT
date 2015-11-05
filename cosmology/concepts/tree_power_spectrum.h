//
// Created by David Seery on 29/10/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_TREE_POWER_SPECTRUM_H
#define LSSEFT_TREE_POWER_SPECTRUM_H

#include "database/powerspectrum_database.h"

#include "boost/filesystem/operations.hpp"


class tree_power_spectrum
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    tree_power_spectrum(const boost::filesystem::path& p);

    //! destructor is default
    ~tree_power_spectrum() = default;


    // INTERNAL DATA

  private:

    //! power spectrum
    powerspectrum_database db;

  };


#endif //LSSEFT_TREE_POWER_SPECTRUM_H
