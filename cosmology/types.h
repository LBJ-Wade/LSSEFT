//
// Created by David Seery on 13/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
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
#include "cosmology/concepts/Matsubara_A.h"

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
                              const std::shared_ptr<tree_power_spectrum>& _Pk)
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
    const std::shared_ptr<tree_power_spectrum>& get_tree_Pk_db() const { return(this->Pk); }


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
    std::shared_ptr<tree_power_spectrum> Pk;

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
                            const std::shared_ptr<tree_power_spectrum>& _Pk)
      : k(_k),
        gf_factors(gf),
        loop_data(k),
        Pk(_Pk)
      {
      }
    
    
    // INTERFACE
    
  public:
    
    //! get wavenubmer
    const Mpc_units::energy& operator*() const { return this->k; }
    
    //! get growth factor database
    const std::shared_ptr<oneloop_growth>& get_gf_factors() const { return this->gf_factors; }
    
    //! get loop kernels
    const std::shared_ptr<loop_integral>& get_loop_data() const { return this->loop_data; }
    
    //! get tree-level power spectrum
    const std::shared_ptr<tree_power_spectrum>& get_tree_Pk_db() const { return this->Pk; }
    
    
    // INTERNAL DATA
  
  private:
    
    // Payload data
    
    //! physical scale k
    Mpc_units::energy k;
    
    //! growth factors
    std::shared_ptr<oneloop_growth> gf_factors;
    
    //! loop momentum kernels for this (k, IR, UV) combination
    std::shared_ptr<loop_integral> loop_data;
    
    //! tree-level power spectrum
    std::shared_ptr<tree_power_spectrum> Pk;

  };

//! list of work for power spectrum calculation
typedef std::list<one_loop_Pk_work_record> one_loop_Pk_work_list;


//! work record for a Matsubara A-coefficient
class Matsubara_A_work_record
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! constructor
    Matsubara_A_work_record(const Mpc_units::energy& _IR, const IR_resum_token& _IRt,
                            const std::shared_ptr<tree_power_spectrum>& _Pk)
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
    const std::shared_ptr<tree_power_spectrum>& get_tree_Pk_db() const { return this->Pk; }
    
    
    // INTERNAL DATA
  
  private:
    
    // Payload data
    
    //! IR resummation scale
    Mpc_units::energy IR_resum;
    
    //! IR resummation token
    IR_resum_token IR_resum_tok;
    
    //! tree-level power spectrum
    std::shared_ptr<tree_power_spectrum> Pk;
    
  };

//! list of work
typedef std::list<Matsubara_A_work_record> Matsubara_A_work_list;


//! work record for a (one-loop) multipole power spectrum calculation
class multipole_Pk_work_record
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! constructor
    multipole_Pk_work_record(const Mpc_units::energy& _k, const Matsubara_A& _A,
                             const std::shared_ptr<oneloop_Pk>& _data,
                             const oneloop_growth_record& _gf_data,
                             const std::shared_ptr<tree_power_spectrum>& _Pk)
      : k(_k),
        A(_A),
        data(_data),
        gf_data(_gf_data),
        Pk(_Pk)
      {
      }
    
    
    // INTERFACE
    
  public:
    
    //! get wavenumber
    const Mpc_units::energy& operator*() const { return this->k; }
    
    //! get Matsubara A coefficient
    const Matsubara_A& get_Matsubara_A() const { return this->A; }
    
    //! get one-loop P(k) data
    const std::shared_ptr<oneloop_Pk>& get_Pk_data() const { return this->data; }
    
    //! get gf growth factors
    const oneloop_growth_record& get_gf_data() const { return this->gf_data; }
    
    //! get tree-level power spectrum
    const std::shared_ptr<tree_power_spectrum>& get_tree_Pk_db() const { return this->Pk; }
    
    
    // INTERNAL DATA
    
  private:
    
    // Payload data
    
    //! physical scale k
    Mpc_units::energy k;
    
    //! Matsubara-A coefficient
    Matsubara_A A;
    
    //! one-loop power spectrum data
    std::shared_ptr<oneloop_Pk> data;
    
    //! gf growth factors
    oneloop_growth_record gf_data;
    
    //! tree-level power spectrum
    std::shared_ptr<tree_power_spectrum> Pk;
    
  };

//! list of work
typedef std::list<multipole_Pk_work_record> multipole_Pk_work_list;


#endif //LSSEFT_TYPES_H
