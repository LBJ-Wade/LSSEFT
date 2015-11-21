//
// Created by David Seery on 13/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_TYPES_H
#define LSSEFT_TYPES_H


#include <map>

#include "database/tokens.h"
#include "database/z_database.h"
#include "cosmology/concepts/tree_power_spectrum.h"

#include "units/eV_units.h"


//! work record for a transfer function calculation
class transfer_work_record
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor - this takes a shared pointer to the redshift database, because eventually
    //! we are going to share ownership with an object representing an MPI
    //! message. Ultimately we may want to look at this again (it doesn't seem to get the
    //! ownership concept right) but it avoids costly copies of the redshift database
    transfer_work_record(const eV_units::energy& _k, const k_token& kt, std::shared_ptr<z_database>& z)
      : k(_k),
        k_tok(kt),
        z_db(z)
      {
      }

    //! destructor is default
    ~transfer_work_record() = default;


    // INTERFACE

  public:

    //! dereference to get wavenumber
    const eV_units::energy& operator*() const { return(this->k); }

    //! get token
    const k_token& get_token() const { return(this->k_tok); }

    //! get redshift database
    const std::shared_ptr<z_database>& get_z_db() const { return(this->z_db); }


    // INTERNAL DATA

  private:

    //! wavenumber
    eV_units::energy k;

    //! wavenumber token
    k_token k_tok;

    //! redshift database
    std::shared_ptr<z_database> z_db;

  };

//! list of work for transfer function calculation
typedef std::list<transfer_work_record> transfer_work_list;


//! work record for a momentum integral calculation
class loop_momentum_work_record
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    loop_momentum_work_record(const eV_units::energy& _k, const k_token& kt,
                              const eV_units::energy& _UV, const UV_token& UVt,
                              const eV_units::energy& _IR, const IR_token& IRt,
                              std::shared_ptr<tree_power_spectrum>& _Pk)
      : k(_k),
        UV_cutoff(_UV),
        IR_cutoff(_IR),
        k_tok(kt),
        UV_tok(UVt),
        IR_tok(IRt),
        Pk(_Pk)
      {
      }

    //! destructor is default
    ~loop_momentum_work_record() = default;


    // INTERFACE

  public:

    //! deference to get wavenumber
    const eV_units::energy& operator*() const { return(this->k); }

    //! get token
    const k_token& kt_get_token() const { return(this->k_tok); }

    //! get tree-level power spectrum
    const std::shared_ptr<tree_power_spectrum>& get_Pk_db() const { return(this->Pk); }


    // INTERNAL DATA

  private:

    //! wavenumber
    eV_units::energy k;

    //! UV cutoff
    eV_units::energy UV_cutoff;

    //! IR cutoff
    eV_units::energy IR_cutoff;

    //! wavenumber token
    k_token k_tok;

    //! UV cutoff token
    UV_token UV_tok;

    //! IR cutoff token
    IR_token IR_tok;

    //! tree-level power spectrum
    std::shared_ptr<tree_power_spectrum> Pk;

  };

//! list of work for loop integral calculation
typedef std::list<loop_momentum_work_record> loop_momentum_work_list;


#endif //LSSEFT_TYPES_H
