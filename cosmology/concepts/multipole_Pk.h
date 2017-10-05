//
// Created by David Seery on 18/11/2016.
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

#ifndef LSSEFT_MULTIPOLE_PK_H
#define LSSEFT_MULTIPOLE_PK_H


#include <map>

#include "Pk_resum_value.h"

#include "database/tokens.h"
#include "units/Mpc_units.h"

#include "boost/timer/timer.hpp"
#include "boost/optional.hpp"
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/map.hpp"


template <typename ValueType>
class Legendre_multiplet
  {
    
    // TYPES
  
  public:
    
    typedef ValueType value_type;
    
    
    // CONSTRUCTOR, DESTSRUCTOR
    
  public:
    
    //! constructor
    Legendre_multiplet(value_type ell0_, value_type ell2_, value_type ell4_)
      : ell0(std::move(ell0_)),
        ell2(std::move(ell2_)),
        ell4(std::move(ell4_))
      {
      }
    
    //! empty constructor
    Legendre_multiplet()
      : ell0(),
        ell2(),
        ell4()
      {
      }
    
    //! destructor is default
    ~Legendre_multiplet() = default;
    
    
    // ACCESSORS
    
  public:
    
    //! get ell = 0 mode
    const value_type& get_ell0() const { return this->ell0; }
    
    //! get ell = 2 mode
    const value_type& get_ell2() const { return this->ell2; }
    
    //! get ell = 4 mode
    const value_type& get_ell4() const { return this->ell4; }
    
    
    // INTERNAL DATA
    
  private:
    
    //! ell = 0 mode
    value_type ell0;
    
    //! ell = 2 mode
    value_type ell2;
    
    //! ell = 4 mode
    value_type ell4;
    
    
    // enable boost::serialization support
    friend class boost::serialization::access;
    
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & ell0;
        ar & ell2;
        ar & ell4;
      }
    
  };


typedef Legendre_multiplet<k2_Pk_resum> k2_Pk_resum_multiplet;
typedef Legendre_multiplet<Pk_resum> Pk_resum_multiplet;


class Pk_ell
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    Pk_ell(const Pk_resum& _tree, const Pk_resum& _P13, const Pk_resum& _P22, const Pk_resum& _PSPT);
    
    //! destructor is default
    ~Pk_ell() = default;
    
    
    // INTERFACE

  public:

    //! get tree value
    const Pk_resum& get_tree() const { return Ptree; }
    
    //! get 13 term
    const Pk_resum& get_13() const { return P13; }
    
    //! get 22 term
    const Pk_resum& get_22() const { return P22; }

    //! get full 1-loop SPT
    const Pk_resum& get_1loop_SPT() const { return PSPT; }

    
    // INTERNAL DATA
    
  private:
    
    //! tree term
    Pk_resum Ptree;
    
    //! 13 term
    Pk_resum P13;
    
    //! 22 term
    Pk_resum P22;
    
    //! full 1-loop SPT power spectrum
    Pk_resum PSPT;

    
    // enable boost::serialization support
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & Ptree;
        ar & P13;
        ar & P22;
        ar & PSPT;
      }
    
  };


namespace boost
  {
    
    namespace serialization
      {
        
        template <typename Archive>
        inline void save_construct_data(Archive& ar, const Pk_ell* t, const unsigned int file_version)
          {
          }
        
        template <typename Archive>
        inline void load_construct_data(Archive& ar, Pk_ell* t, const unsigned int file_version)
          {
            // invoke in-place constructor with zero elements; will be overwritten during deserialization
            ::new(t) Pk_ell(Pk_resum(), Pk_resum(), Pk_resum(), Pk_resum());
          }
        
      }   // namespace serialization
    
  }   // namespace boost


class multipole_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    multipole_Pk(const k_token kt, const growth_params_token& gp, const loop_integral_params_token& lp,
                 const MatsubaraXY_params_token& XYp, const linear_Pk_token Pkt_i,
                 const boost::optional<linear_Pk_token> Pkt_f, const IR_cutoff_token IRt, const UV_cutoff_token UVt,
                 const z_token zt, const IR_resum_token IRrt, const Pk_ell _P0, const Pk_ell _P2, const Pk_ell _P4);
    
    //! empty constructor for us in MPI payloads
    multipole_Pk();
    
    //! destructor is default
    ~multipole_Pk() = default;
    
    
    // INTERFACE
  
  public:
    
    //! get wavenumber token
    const k_token& get_k_token() const { return this->k; }
    
    //! get growth parameters token
    const growth_params_token& get_growth_params_token() const { return this->growth_params; }
    
    //! get loop integral parameters token
    const loop_integral_params_token& get_loop_params_token() const { return this->loop_params; }
    
    //! get XY parameters token
    const MatsubaraXY_params_token& get_XY_params_token() const { return this->XY_params; }
    
    //! get UV cutoff token
    const UV_cutoff_token& get_UV_cutoff_token() const { return this->UV_cutoff; }
    
    //! get IR cutoff token
    const IR_cutoff_token& get_IR_cutoff_token() const { return this->IR_cutoff; }
    
    //! get z token
    const z_token& get_z_token() const { return this->z; }
    
    //! get IR resummation scale token
    const IR_resum_token& get_IR_resum_token() const { return this->IR_resum; }
    
    //! get initial linear power spectrum token
    const linear_Pk_token& get_init_Pk_token() const { return this->init_Pk; }
    
    //! get final power spectrum token, if provided
    const boost::optional<linear_Pk_token>& get_final_Pk_token() const { return this->final_Pk; }
    
    
    // ACCESSORS
    
  public:
    
    //! get monopole
    const Pk_ell& get_P0() const { return this->P0; }
    
    //! get quadrupole
    const Pk_ell& get_P2() const { return this->P2; }
    
    //! get hexadecapole
    const Pk_ell& get_P4() const { return this->P4; }
    
    
    // INTERNAL DATA
  
  private:
    
    // CONFIGURATION DATA
    
    //! wavenumber token
    k_token k;
    
    //! growth parameters token
    growth_params_token growth_params;
    
    //! loop parameters token
    loop_integral_params_token loop_params;
    
    //! XY parameters token
    MatsubaraXY_params_token XY_params;
    
    //! initial linear power spectrum token
    linear_Pk_token init_Pk;
    
    //! final linear power spectrum token, if provided
    boost::optional<linear_Pk_token> final_Pk;
    
    //! UV cutoff token
    UV_cutoff_token UV_cutoff;
    
    //! IR cutoff token
    IR_cutoff_token IR_cutoff;
    
    //! redshift token
    z_token z;
    
    //! IR resummation scale token
    IR_resum_token IR_resum;
    
    
    // VALUES
    
    //! P0 multipole
    Pk_ell P0;
    
    //! P2 multipole
    Pk_ell P2;
    
    //! P4 multipole
    Pk_ell P4;
    
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & k;
        ar & growth_params;
        ar & loop_params;
        ar & XY_params;
        ar & init_Pk;
        ar & final_Pk;
        ar & UV_cutoff;
        ar & IR_cutoff;
        ar & z;
        ar & IR_resum;
        ar & P0;
        ar & P2;
        ar & P4;
      }

  };


using multipole_Pk_set = std::map< std::string, multipole_Pk >;


namespace boost
  {
    
    namespace serialization
      {
        
        template <typename Archive>
        inline void save_construct_data(Archive& ar, const multipole_Pk* t, const unsigned int file_version)
          {
          }
        
        template <typename Archive>
        inline void load_construct_data(Archive& ar, multipole_Pk* t, const unsigned int file_version)
          {
            // invoke in-place constructor with zero elements; will be overwritten during deserialization
            Pk_resum empty_Pk;

            Pk_ell empty(empty_Pk, empty_Pk, empty_Pk, empty_Pk);
    
            ::new(t) multipole_Pk(k_token(0), growth_params_token(0), loop_integral_params_token(0),
                                  MatsubaraXY_params_token(0), linear_Pk_token(0),
                                  boost::none, IR_cutoff_token(0), UV_cutoff_token(0), z_token(0), IR_resum_token(0),
                                  empty, empty, empty);
          }
        
      }   // namespace serialization
    
  }   // namespace boost


#endif //LSSEFT_MULTIPOLE_PK_H
