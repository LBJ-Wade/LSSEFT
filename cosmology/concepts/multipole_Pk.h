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


#include "database/tokens.h"
#include "units/Mpc_units.h"

#include "boost/timer/timer.hpp"
#include "boost/serialization/serialization.hpp"


class Pk_ell
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    Pk_ell(const Mpc_units::inverse_energy3& _Pt, const Mpc_units::inverse_energy3& _Pt_resum,
           const Mpc_units::inverse_energy3& _P13, const Mpc_units::inverse_energy3& _P13_resum,
           const Mpc_units::inverse_energy3& _P22, const Mpc_units::inverse_energy3& _P22_resum,
           const Mpc_units::inverse_energy3& _PSPT, const Mpc_units::inverse_energy3& _PSPT_resum,
           const Mpc_units::inverse_energy& _Z2d, const Mpc_units::inverse_energy3& _Z0v,
           const Mpc_units::inverse_energy& _Z2v, const Mpc_units::inverse_energy3& _Z0vd,
           const Mpc_units::inverse_energy& _Z2vd, const Mpc_units::inverse_energy& _Z2vv,
           const Mpc_units::inverse_energy& _Z2vvd, const Mpc_units::inverse_energy& _Z2vvv);
    
    //! empty constructor for use when overwriting with MPI payloads
    Pk_ell();
    
    //! destructor is default
    ~Pk_ell() = default;
    
    
    // INTERFACE

  public:

    //! get tree value
    Mpc_units::inverse_energy3& get_tree() { return Ptree; }
    const Mpc_units::inverse_energy3& get_tree() const { return Ptree; }
    
    //! get resummed tree value
    Mpc_units::inverse_energy3& get_tree_resum() { return Ptree_resum; }
    const Mpc_units::inverse_energy3& get_tree_resum() const { return Ptree_resum; }
    
    //! get 13 term
    Mpc_units::inverse_energy3& get_13() { return P13; }
    const Mpc_units::inverse_energy3& get_13() const { return P13; }
    
    //! get resummed 13 term
    Mpc_units::inverse_energy3& get_13_resum() { return P13_resum; }
    const Mpc_units::inverse_energy3& get_13_resum() const { return P13_resum; }
    
    //! get 22 term
    Mpc_units::inverse_energy3& get_22() { return P22; }
    const Mpc_units::inverse_energy3& get_22() const { return P22; }
    
    //! get resummed 22 term
    Mpc_units::inverse_energy3& get_22_resum() { return P22_resum; }
    const Mpc_units::inverse_energy3& get_22_resum() const { return P22_resum; }

    //! get full 1-loop SPT
    Mpc_units::inverse_energy3& get_1loop_SPT() { return PSPT; }
    const Mpc_units::inverse_energy3& get_1loop_SPT() const { return PSPT; }
    
    //! get full 1-loop SPT resummed
    Mpc_units::inverse_energy3& get_1loop_SPT_resum() { return PSPT_resum; }
    const Mpc_units::inverse_energy3& get_1loop_SPT_resum() const { return PSPT_resum; }
    
    //! get Z2_delta counterterm
    Mpc_units::inverse_energy& get_Z2_delta() { return Z2_delta; }
    const Mpc_units::inverse_energy& get_Z2_delta() const { return Z2_delta; }
    
    //! get Z0_v counterterm
    Mpc_units::inverse_energy3& get_Z0_v() { return Z0_v; }
    const Mpc_units::inverse_energy3& get_Z0_v() const { return Z0_v; }
    
    //! get Z2_v counterterm
    Mpc_units::inverse_energy& get_Z2_v() { return Z2_v; }
    const Mpc_units::inverse_energy& get_Z2_v() const { return Z2_v; }
    
    //! get Z0_vdelta counterterm
    Mpc_units::inverse_energy3& get_Z0_vdelta() { return Z0_vdelta; }
    const Mpc_units::inverse_energy3& get_Z0_vdelta() const { return Z0_vdelta; }
    
    //! get Z2_vdelta counterterm
    Mpc_units::inverse_energy& get_Z2_vdelta() { return Z2_vdelta; }
    const Mpc_units::inverse_energy& get_Z2_vdelta() const { return Z2_vdelta; }
    
    //! get Z2_v counterterm
    Mpc_units::inverse_energy& get_Z2_vv() { return Z2_vv; }
    const Mpc_units::inverse_energy& get_Z2_vv() const { return Z2_vv; }
    
    //! get Z2_vvdelta counterterm
    Mpc_units::inverse_energy& get_Z2_vvdelta() { return Z2_vvdelta; }
    const Mpc_units::inverse_energy& get_Z2_vvdelta() const { return Z2_vvdelta; }
    
    //! get Z2_vvv counterterm
    Mpc_units::inverse_energy& get_Z2_vvv() { return Z2_vvv; }
    const Mpc_units::inverse_energy& get_Z2_vvv() const { return Z2_vvv; }
    
    
    // INTERNAL DATA
    
  private:
    
    //! tree term
    Mpc_units::inverse_energy3 Ptree;
  
    //! tree term - resummed
    Mpc_units::inverse_energy3 Ptree_resum;
    
    //! 13 term
    Mpc_units::inverse_energy3 P13;
    
    //! 13 term - resummed
    Mpc_units::inverse_energy3 P13_resum;
    
    //! 22 term
    Mpc_units::inverse_energy3 P22;
    
    //! 22 term - resummed
    Mpc_units::inverse_energy3 P22_resum;
    
    //! full 1-loop SPT power spectrum
    Mpc_units::inverse_energy3 PSPT;
    
    //! full 1-loop SPT power spectrum - resummed
    Mpc_units::inverse_energy3 PSPT_resum;
    
    
    // COUNTERTERMS
    
    //! Z2_delta
    Mpc_units::inverse_energy Z2_delta;
    
    //! Z0_v
    Mpc_units::inverse_energy3 Z0_v;
    
    //! Z2_v
    Mpc_units::inverse_energy Z2_v;
    
    //! Z0_vdelta
    Mpc_units::inverse_energy3 Z0_vdelta;
    
    //! Z2_vdelta
    Mpc_units::inverse_energy Z2_vdelta;
    
    //! Z2_vv
    Mpc_units::inverse_energy Z2_vv;
    
    //! Z2_vvdelta
    Mpc_units::inverse_energy Z2_vvdelta;
    
    //! Z2_vvv
    Mpc_units::inverse_energy Z2_vvv;
    
    
    // enable boost::serialization support
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & Ptree;
        ar & Ptree_resum;
        ar & P13;
        ar & P13_resum;
        ar & P22;
        ar & P22_resum;
        ar & PSPT;
        ar & PSPT_resum;
        ar & Z2_delta;
        ar & Z0_v;
        ar & Z2_v;
        ar & Z0_vdelta;
        ar & Z2_vdelta;
        ar & Z2_vv;
        ar & Z2_vvdelta;
        ar & Z2_vvv;
      }
    
  };


class multipole_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    multipole_Pk(const k_token& kt, const linear_Pk_token& Pkt_i, const boost::optional<linear_Pk_token>& Pkt_f,
                 const IR_cutoff_token& IRt, const UV_cutoff_token& UVt, const z_token& zt, const IR_resum_token& IRrt,
                 const Pk_ell& _P0, const Pk_ell& _P2, const Pk_ell& _P4);
    
    //! empty constructor used for receiving MPI payloads
    multipole_Pk();
    
    //! destructor is default
    ~multipole_Pk() = default;
    
    
    // INTERFACE
  
  public:
    
    //! get wavenumber token
    const k_token& get_k_token() const { return this->k; }
    
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
        ar & init_Pk;
        ar & UV_cutoff;
        ar & IR_cutoff;
        ar & z;
        ar & IR_resum;
        ar & P0;
        ar & P2;
        ar & P4;
      }

  };


#endif //LSSEFT_MULTIPOLE_PK_H
