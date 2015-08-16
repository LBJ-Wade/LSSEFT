//
// Created by David Seery on 15/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_TRANSFER_INTEGRATOR_H
#define LSSEFT_TRANSFER_INTEGRATOR_H


#include <memory>

#include "cosmology/FRW_model.h"
#include "units/eV_units.h"
#include "database/tokens.h"
#include "database/redshift_database.h"

#include "defaults.h"


class transfer_function
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    transfer_function(const FRW_model& m, const eV_units::energy& _k, const wavenumber_token& t, std::shared_ptr<redshift_database>& z);

    //! destructor is default
    ~transfer_function() = default;


    // INTERNAL DATA

  private:

    // CONFIGURATION DATA

    //! FRW model
    FRW_model model;

    //! wavenumber
    eV_units::energy k;

    //! wavenumber token
    wavenumber_token token;

    //! redshift database; managed using a std::shared_ptr<>
    //! to avoid unnecessary duplication expense
    std::shared_ptr<redshift_database> z_db;


    // TRANSFER FUNCTIONS

    // these are managed using std::shared_ptr<>s to avoid
    // expensive duplication

    //! delta_m transfer function
    std::shared_ptr< std::vector<double> > delta_m;

    //! theta_m transfer function
    std::shared_ptr< std::vector<double> > theta_m;

    //! delta_r transfer function
    std::shared_ptr< std::vector<double> > delta_r;

    //! theta_r transfer function
    std::shared_ptr< std::vector<double> > theta_r;

    //! Phi transfer function
    std::shared_ptr< std::vector<double> > Phi;

  };


class transfer_integrator
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor is default
    transfer_integrator(double a=LSSEFT_DEFAULT_ABS_ERR, double r=LSSEFT_DEFAULT_REL_ERR);

    //! destructor is default
    ~transfer_integrator() = default;


    // TRANSFER FUNCTION

  public:

    //! integrate transfer function for a given k-mode and set of redshift samples
    transfer_function integrate(const FRW_model& model, const eV_units::energy& k, const wavenumber_token& tok,
                                std::shared_ptr<redshift_database>& z_db);


    // INTERNAL DATA

  private:

    // ERRORS

    //! required absolute error
    double abs_err;

    //! required relative error
    double rel_err;

  };


#endif //LSSEFT_TRANSFER_INTEGRATOR_H
