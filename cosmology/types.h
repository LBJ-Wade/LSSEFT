//
// Created by David Seery on 13/08/2015.
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

#ifndef LSSEFT_TYPES_H
#define LSSEFT_TYPES_H


#include <map>

#include "database/tokens.h"
#include "database/z_database.h"
#include "cosmology/concepts/power_spectrum.h"
#include "cosmology/concepts/oneloop_growth.h"
#include "cosmology/concepts/loop_integral.h"
#include "cosmology/concepts/oneloop_Pk.h"
#include "cosmology/concepts/Matsubara_XY.h"

#include "units/Mpc_units.h"


//! work record for a transfer function calculation
class transfer_work_record
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor - this takes a shared pointer to the redshift database, because eventually
    //! we are going to share ownership with an object representing an MPI
    //! message. Ultimately we may want to look at this again (it doesn't seem to get the
    //! ownership concept right) but it avoids costly copies of the redshift database
    transfer_work_record(const Mpc_units::energy& _k, const k_token& kt, std::shared_ptr<z_database>& z)
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
    const Mpc_units::energy& operator*() const { return(this->k); }

    //! get token
    const k_token& get_token() const { return(this->k_tok); }

    //! get redshift database
    const std::shared_ptr<z_database>& get_z_db() const { return(this->z_db); }


    // INTERNAL DATA

  private:

    //! wavenumber
    Mpc_units::energy k;

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
    loop_momentum_work_record(const Mpc_units::energy& _k, const k_token& kt,
                              const Mpc_units::energy& _UV, const UV_cutoff_token& UVt,
                              const Mpc_units::energy& _IR, const IR_cutoff_token& IRt,
                              const std::shared_ptr<initial_filtered_Pk>& _Pk)
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
    const Mpc_units::energy& operator*() const { return(this->k); }

    //! get wavenumber token
    const k_token& get_k_token() const { return(this->k_tok); }

    //! get UV cutoff
    const Mpc_units::energy& get_UV_cutoff() const { return(this->UV_cutoff); }

    //! get UV cutoff token
    const UV_cutoff_token& get_UV_token() const { return(this->UV_tok); }

    //! get IR cutoff
    const Mpc_units::energy& get_IR_cutoff() const { return(this->IR_cutoff); }

    //! get IR cutoff token
    const IR_cutoff_token& get_IR_token() const { return(this->IR_tok); }

    //! get tree-level power spectrum
    const std::shared_ptr<initial_filtered_Pk>& get_tree_Pk_db() const { return(this->Pk); }


    // INTERNAL DATA

  private:

    //! wavenumber
    Mpc_units::energy k;

    //! UV cutoff
    Mpc_units::energy UV_cutoff;

    //! IR cutoff
    Mpc_units::energy IR_cutoff;

    //! wavenumber token
    k_token k_tok;

    //! UV cutoff token
    UV_cutoff_token UV_tok;

    //! IR cutoff token
    IR_cutoff_token IR_tok;

    //! tree-level power spectrum
    std::shared_ptr<initial_filtered_Pk> Pk;

  };

//! list of work for loop integral calculation
typedef std::list<loop_momentum_work_record> loop_momentum_work_list;


//! work record for a one-loop power spectrum calculation
class one_loop_Pk_work_record
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    one_loop_Pk_work_record(const Mpc_units::energy _k,
                            const std::shared_ptr<oneloop_growth>& gf, const std::shared_ptr<loop_integral>& k,
                            const std::shared_ptr<initial_filtered_Pk>& _Pk_init,
                            const std::shared_ptr<final_filtered_Pk>& _Pk_final)
      : k(_k),
        gf_factors(gf),
        loop_data(k),
        Pk_init(_Pk_init),
        Pk_final(_Pk_final)
      {
      }
    
    
    // INTERFACE
    
  public:
    
    //! get wavenumber
    const Mpc_units::energy& operator*() const { return this->k; }
    
    //! get growth factor database
    const std::shared_ptr<oneloop_growth>& get_gf_factors() const { return this->gf_factors; }
    
    //! get loop kernels
    const std::shared_ptr<loop_integral>& get_loop_data() const { return this->loop_data; }
    
    //! get initial linear power spectrum
    const std::shared_ptr<initial_filtered_Pk>& get_init_linear_Pk() const { return this->Pk_init; }
    
    //! get final linear power spectrum
    const std::shared_ptr<final_filtered_Pk>& get_final_linear_Pk() const { return this->Pk_final; }
    
    
    // INTERNAL DATA
  
  private:
    
    // Payload data
    
    //! physical scale k
    Mpc_units::energy k;
    
    //! growth factors
    std::shared_ptr<oneloop_growth> gf_factors;
    
    //! loop momentum kernels for this (k, IR, UV) combination
    std::shared_ptr<loop_integral> loop_data;
    
    //! initial linear power spectrum
    std::shared_ptr<initial_filtered_Pk> Pk_init;
    
    //! final linear power spectrum, if provided
    std::shared_ptr<final_filtered_Pk> Pk_final;

  };

//! list of work for power spectrum calculation
typedef std::list<one_loop_Pk_work_record> one_loop_Pk_work_list;


//! work record for Matsubara X & Y coefficients
class Matsubara_XY_work_record
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! constructor
    Matsubara_XY_work_record(const Mpc_units::energy& _IR, const IR_resum_token& _IRt,
                            const std::shared_ptr<initial_filtered_Pk>& _Pk)
      : IR_resum(_IR),
        IR_resum_tok(_IRt),
        Pk(_Pk)
      {
      }
    
    
    // INTERFACE
  
  public:
    
    //! get IR resummation scale
    const Mpc_units::energy& get_IR_resum() const { return this->IR_resum; }
    
    //! get IR resummation token
    const IR_resum_token& get_IR_resum_token() const { return this->IR_resum_tok; }
    
    //! get tree-level power spectrum
    const std::shared_ptr<initial_filtered_Pk>& get_linear_Pk() const { return this->Pk; }
    
    
    // INTERNAL DATA
  
  private:
    
    // Payload data
    
    //! IR resummation scale
    Mpc_units::energy IR_resum;
    
    //! IR resummation token
    IR_resum_token IR_resum_tok;
    
    //! tree-level power spectrum
    std::shared_ptr<initial_filtered_Pk> Pk;
    
  };

//! list of work
typedef std::list<Matsubara_XY_work_record> Matsubara_XY_work_list;


// work record for a one-loop resummed P(k) calculation
class one_loop_resum_Pk_work_record
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor
    one_loop_resum_Pk_work_record(const Mpc_units::energy& _k, const Matsubara_XY& _XY,
                                  const std::shared_ptr<oneloop_Pk>& _data,
                                  const oneloop_growth_record& _gf_data,
                                  const std::shared_ptr<initial_filtered_Pk>& _Pk_init,
                                  const std::shared_ptr<final_filtered_Pk>& _Pk_final)
      : k(_k),
        XY(_XY),
        data(_data),
        gf_data(_gf_data),
        Pk_init(_Pk_init),
        Pk_final(_Pk_final)
      {
      }
    
    //! destructor is default
    ~one_loop_resum_Pk_work_record() = default;
    
    
    // INTERFACE
  
  public:
    
    //! get wavenumber
    const Mpc_units::energy& operator*() const { return this->k; }
    
    //! get Matsubara X & Y coefficients
    const Matsubara_XY& get_Matsubara_XY() const { return this->XY; }
    
    //! get one-loop P(k) data
    const std::shared_ptr<oneloop_Pk>& get_Pk_data() const { return this->data; }
    
    //! get gf growth factors
    const oneloop_growth_record& get_gf_data() const { return this->gf_data; }
    
    //! get initial linear power spectrum
    const std::shared_ptr<initial_filtered_Pk>& get_init_linear_Pk() const { return this->Pk_init; }
    
    //! get final linear power spectrum
    const std::shared_ptr<final_filtered_Pk>& get_final_linear_Pk() const { return this->Pk_final; }
    
    
    // INTERNAL DATA
  
  private:
    
    // Payload data
    
    //! physical scale k
    Mpc_units::energy k;
    
    //! Matsubara X & Y coefficients
    Matsubara_XY XY;
    
    //! one-loop power spectrum data
    std::shared_ptr<oneloop_Pk> data;
    
    //! gf growth factors
    oneloop_growth_record gf_data;
    
    //! initial linear power spectrum
    std::shared_ptr<initial_filtered_Pk> Pk_init;
    
    //! final linear power spectrum
    std::shared_ptr<final_filtered_Pk> Pk_final;
    
  };

//! list of work
typedef std::list<one_loop_resum_Pk_work_record> one_loop_resum_Pk_work_list;


//! work record for a (one-loop) multipole power spectrum calculation
class multipole_Pk_work_record
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor
    multipole_Pk_work_record(const Mpc_units::energy& _k, const Matsubara_XY& _XY,
                             const std::shared_ptr<oneloop_Pk>& _data, const oneloop_growth_record& _gf_data,
                             const std::shared_ptr<initial_filtered_Pk>& _Pk_init,
                             const std::shared_ptr<final_filtered_Pk>& _Pk_final)
      : k(_k),
        XY(_XY),
        data(_data),
        gf_data(_gf_data),
        Pk_init(_Pk_init),
        Pk_final(_Pk_final)
      {
      }
    
    
    // INTERFACE
    
  public:
    
    //! get wavenumber
    const Mpc_units::energy& operator*() const { return this->k; }
    
    //! get Matsubara X & Y coefficients
    const Matsubara_XY& get_Matsubara_XY() const { return this->XY; }
    
    //! get one-loop P(k) data
    const std::shared_ptr<oneloop_Pk>& get_Pk_data() const { return this->data; }
    
    //! get gf growth factors
    const oneloop_growth_record& get_gf_data() const { return this->gf_data; }
    
    //! get initial linear power spectrum
    const std::shared_ptr<initial_filtered_Pk>& get_init_linear_Pk() const { return this->Pk_init; }
    
    //! get final linear power spectrum
    const std::shared_ptr<final_filtered_Pk>& get_final_linear_Pk() const { return this->Pk_final; }
    
    
    // INTERNAL DATA
    
  private:
    
    // Payload data
    
    //! physical scale k
    Mpc_units::energy k;
    
    //! Matsubara X & Y coefficients
    Matsubara_XY XY;
    
    //! one-loop power spectrum data
    std::shared_ptr<oneloop_Pk> data;
    
    //! gf growth factors
    oneloop_growth_record gf_data;
    
    //! initial linear power spectrum
    std::shared_ptr<initial_filtered_Pk> Pk_init;
    
    //! final linear power spectrum
    std::shared_ptr<final_filtered_Pk> Pk_final;
    
  };

//! list of work
typedef std::list<multipole_Pk_work_record> multipole_Pk_work_list;


//! work record for filtering the linear power spectrum into wiggle/no-wiggle components
class filter_Pk_work_record
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor
    filter_Pk_work_record(const Mpc_units::energy& _k, const k_token& kt,
                          std::shared_ptr<filterable_Pk>& Pk, const linear_Pk_token& Pt)
      : k(_k),
        k_tok(kt),
        Pk_lin(Pk),
        Pk_tok(Pt)
      {
      }
    
    
    // INTERFACE
    
  public:
    
    //! get wavenumber
    const Mpc_units::energy& operator*() const { return this->k; }
    
    //! get wavenumber token
    const k_token& get_k_token() const { return this->k_tok; }
    
    //! get unfiltered linear power spectrum
    const std::shared_ptr<filterable_Pk>& get_linear_Pk() const { return this->Pk_lin; }
    
    //! get power spectrum token
    const linear_Pk_token& get_Pk_token() const { return this->Pk_tok; }
    
    
    // INTERNAL DATA
    
  private:
    
    // Payload data
    
    //! physical scale k
    Mpc_units::energy k;
    
    //! wavenumber token
    k_token k_tok;
    
    //! unfiltered linear power spectrum
    std::shared_ptr<filterable_Pk> Pk_lin;
    
    //! token for unfiltered linear power spectrum
    linear_Pk_token Pk_tok;
    
  };

//! list of work
typedef std::list<filter_Pk_work_record> filter_Pk_work_list;


#endif //LSSEFT_TYPES_H
