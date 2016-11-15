//
// Created by David Seery on 14/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_ONE_LOOP_PK_H
#define LSSEFT_ONE_LOOP_PK_H


#include "database/tokens.h"
#include "units/Mpc_units.h"

#include "boost/timer/timer.hpp"
#include "boost/serialization/serialization.hpp"


class Pk_value
  {
    
  public:
    
    typedef Mpc_units::inverse_energy3 value_type;
    
    //! value constructor
    Pk_value(value_type v)
      : value(std::move(v))
      {
      }
    
    //! empty constructor
    Pk_value()
      : value(value_type(0.0))
      {
      }
    
    
    operator value_type() const { return this->value; }
    
    value_type value;
    
  private:
    
    // enable boost::serialization support
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & value;
      }
    
  };


class k2_Pk_value
  {
  
  public:
    
    typedef Mpc_units::inverse_energy value_type;
    
    //! value constructor
    k2_Pk_value(value_type v)
      : value(std::move(v))
      {
      }
    
    //! empty constructor
    k2_Pk_value()
      : value(value_type(0.0))
      {
      }
    
    
    operator value_type() const { return this->value; }
    
    value_type value;
  
  private:
    
    // enable boost::serialization support
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & value;
      }
    
  };


//! overload + so that power spectrum and counterterm values can be added
inline Pk_value operator+(const Pk_value& A, const Pk_value& B)
  {
    Pk_value res;
    
    res.value = A.value + B.value;
    
    return std::move(res);
  }


//! overload + so that power spectrum and counterterm values can be added
inline k2_Pk_value operator+(const k2_Pk_value& A, const k2_Pk_value& B)
  {
    k2_Pk_value res;
    
    res.value = A.value + B.value;
    
    return std::move(res);
  }


class dd_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    dd_Pk(const Pk_value& _Pt, const Pk_value& _P13, const Pk_value& _P22, const k2_Pk_value& _Z2d);
    
    //! empty constructor for use when overwriting with MPI payloads
    dd_Pk();
    
    //! destructor is default
    ~dd_Pk() = default;
    
    
    // INTERFACE
    
  public:
    
    //! get tree value
    Pk_value& get_tree() { return this->Ptree; }
    const Pk_value& get_tree() const { return this->Ptree; }
    
    //! get 13 value
    Pk_value& get_13() { return this->P13; }
    const Pk_value& get_13() const { return this->P13; }
    
    //! get 22 value
    Pk_value& get_22() { return this->P22; }
    const Pk_value& get_22() const { return this->P22; }
    
    //! get total SPT power spectrum = 13 + 22
    Pk_value& get_1loop_SPT() { return this->P1loopSPT; }
    const Pk_value& get_1loop_SPT() const { return this->P1loopSPT; }
    
    
    // COUNTERTERMS
    
    //! get EFT counterterm
    k2_Pk_value& get_Z2_delta() { return this->Z2_delta; }
    const k2_Pk_value& get_Z2_delta() const { return this->Z2_delta; }
    
    
    // INTERNAL DATA
    
  private:
    
    //! tree power spectrum
    Pk_value Ptree;
    
    //! 13 terms
    Pk_value P13;
    
    //! 22 terms
    Pk_value P22;
    
    //! total 1-loop SPT value
    Pk_value P1loopSPT;
    
    //! coefficient of the counterterm Z2_delta
    k2_Pk_value Z2_delta;
    
    
    // enable boost::serialization support
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & Ptree;
        ar & P13;
        ar & P22;
        ar & P1loopSPT;
        ar & Z2_delta;
      }
    
  };


class rsd_dd_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
  
  public:
    
    //! value constructor
    rsd_dd_Pk(const Pk_value& _Pt, const Pk_value& _P13, const Pk_value& _P22,
              const k2_Pk_value& _Z2d, const Pk_value& _Z0v, const k2_Pk_value& _Z2v,
              const Pk_value& _Z0vd, const k2_Pk_value& _Z2vd,
              const k2_Pk_value& _Z2vv, const k2_Pk_value& _Z2vvd, const k2_Pk_value& _Z2vvv);
    
    //! empty constructor for use when overwriting with MPI payloads
    rsd_dd_Pk();
    
    //! destructor is default
    ~rsd_dd_Pk() = default;
    
    
    // INTERFACE
  
  public:
    
    //! get tree value
    Pk_value& get_tree() { return this->Ptree; }
    const Pk_value& get_tree() const { return this->Ptree; }
    
    //! get 13 value
    Pk_value& get_13() { return this->P13; }
    const Pk_value& get_13() const { return this->P13; }
    
    //! get 22 value
    Pk_value& get_22() { return this->P22; }
    const Pk_value& get_22() const { return this->P22; }
    
    //! get total 13 + 22
    Pk_value& get_1loop_SPT() { return this->P1loopSPT; }
    const Pk_value& get_1loop_SPT() const { return this->P1loopSPT; }
    
    
    // COUNTERTERMS
    
    //! get Z2_delta counterterm
    k2_Pk_value& get_Z2_delta() { return this->Z2_delta; }
    const k2_Pk_value& get_Z2_delta() const { return this->Z2_delta; }
    
    //! get Z0_v counterterm
    Pk_value& get_Z0_v() { return this->Z0_v; }
    const Pk_value& get_Z0_v() const { return this->Z0_v; }
    
    //! get Z2_v counterterm
    k2_Pk_value& get_Z2_v() { return this->Z2_v; }
    const k2_Pk_value& get_Z2_v() const { return this->Z2_v; }
    
    //! get Z0_vdelta counterterm
    Pk_value& get_Z0_vdelta() { return this->Z0_vdelta; }
    const Pk_value& get_Z0_vdelta() const { return this->Z0_vdelta; }
    
    //! get Z2_vdelta counterterm
    k2_Pk_value& get_Z2_vdelta() { return this->Z2_vdelta; }
    const k2_Pk_value& get_Z2_vdelta() const { return this->Z2_vdelta; }
    
    //! get Z2_vv counterterm
    k2_Pk_value& get_Z2_vv() { return this->Z2_vv; }
    const k2_Pk_value& get_Z2_vv() const { return this->Z2_vv; }
    
    //! get Z2_vvdelta counterterm
    k2_Pk_value& get_Z2_vvdelta() { return this->Z2_vvdelta; }
    const k2_Pk_value& get_Z2_vvdelta() const { return this->Z2_vvdelta; }
    
    //! get Z2_vvv counterterm
    k2_Pk_value& get_Z2_vvv() { return this->Z2_vvv; }
    const k2_Pk_value& get_Z2_vvv() const { return this->Z2_vvv; }
    
    
    // INTERNAL DATA
  
  private:
    
    //! tree power spectrum
    Pk_value Ptree;
    
    //! 13 terms
    Pk_value P13;
    
    //! 22 terms
    Pk_value P22;
    
    //! total 1-loop SPT value
    Pk_value P1loopSPT;
    
    //! coefficient of the counterterm Z2_delta
    k2_Pk_value Z2_delta;
    
    //! coefficient of the counterterm Z0_v
    Pk_value Z0_v;
    
    //! coefficient of the counterterm Z2_v
    k2_Pk_value Z2_v;
    
    //! coefficient of the counterterm Z0_vd
    Pk_value Z0_vdelta;
    
    //! coefficient of the counterterm Z2_vd
    k2_Pk_value Z2_vdelta;
    
    //! coefficient of the counterterm Z2_vv
    k2_Pk_value Z2_vv;
    
    //! coefficient of the counterterm Z2_vvdelta
    k2_Pk_value Z2_vvdelta;
    
    //! coefficient of the counterterm Z2_vvv
    k2_Pk_value Z2_vvv;
    
    
    // enable boost::serialization support
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & Ptree;
        ar & P13;
        ar & P22;
        ar & P1loopSPT;
      }
    
  };


class one_loop_Pk
  {
    
    // CONSTRUCTOR, DESTRUCTOR
    
  public:
    
    //! value constructor
    one_loop_Pk(const k_token& kt, const UV_token& UVt, const IR_token& IRt, const z_token& zt,
                const dd_Pk& _dd, const rsd_dd_Pk& _rsd_mu0, const rsd_dd_Pk& _rsd_mu2,
                const rsd_dd_Pk& _rsd_mu4, const rsd_dd_Pk& _rsd_mu6, const rsd_dd_Pk& _rsd_mu8);
    
    //! empty constructor, used when receiving an MPI payload
    one_loop_Pk();
    
    //! destructor is default
    ~one_loop_Pk() = default;
    
    
    // INTERFACE
    
  public:
    
    //! get wavenumber token
    const k_token& get_k_token() const { return this->k; }
    
    //! get UV cutoff token
    const UV_token& get_UV_token() const { return this->UV_cutoff; }
    
    //! get IR cutoff token
    const IR_token& get_IR_token() const { return this->IR_cutoff; }
    
    //! get z token
    const z_token& get_z_token() const { return this->z; }
    
    
    //! get delta-delta power spectrum
    const dd_Pk& get_dd() const { return this->dd; }
    
    //! get delta-delta RSD power spectrum mu^0 coefficient
    const rsd_dd_Pk& get_dd_rsd_mu0() const { return this->rsd_dd_mu0; }
    
    //! get delta-delta RSD power spectrum mu^2 coefficient
    const rsd_dd_Pk& get_dd_rsd_mu2() const { return this->rsd_dd_mu2; }
    
    //! get delta-delta RSD power spectrum mu^4 coefficient
    const rsd_dd_Pk& get_dd_rsd_mu4() const { return this->rsd_dd_mu4; }
    
    //! get delta-delta RSD power spectrum mu^6 coefficient
    const rsd_dd_Pk& get_dd_rsd_mu6() const { return this->rsd_dd_mu6; }
    
    //! get delta-delta RSD power spectrum mu^8 coefficient
    const rsd_dd_Pk& get_dd_rsd_mu8() const { return this->rsd_dd_mu8; }
    
    
    // INTERNAL DATA
  
  private:
    
    // CONFIGURATION DATA
    
    //! wavenumber token
    k_token k;
    
    //! UV cutoff token
    UV_token UV_cutoff;
    
    //! IR cutoff token
    IR_token IR_cutoff;
    
    //! redshift token
    z_token z;
    
    
    // VALUES
    
    //! delta-delta power spectrum
    dd_Pk dd;
    
    //! mu^0 term in delta_s-delta_s power spectrum
    rsd_dd_Pk rsd_dd_mu0;
    
    //! mu^2 term in delta_s-delta_s power spectrum
    rsd_dd_Pk rsd_dd_mu2;
    
    //! mu^4 term in delta_s-delta_s power spectrum
    rsd_dd_Pk rsd_dd_mu4;
    
    //! mu^6 term in delta_s-delta_s power spectrum
    rsd_dd_Pk rsd_dd_mu6;
    
    //! mu^8 term in delta_s-delta_s power spectrum
    rsd_dd_Pk rsd_dd_mu8;
    
    
    // enable boost::serialization support, and hence automated packing for transmission over MPI
    friend class boost::serialization::access;
    
    template <typename Archive>
    void serialize(Archive& ar, unsigned int version)
      {
        ar & k;
        ar & UV_cutoff;
        ar & IR_cutoff;
        ar & z;
        ar & dd;
        ar & rsd_dd_mu0;
        ar & rsd_dd_mu2;
        ar & rsd_dd_mu4;
        ar & rsd_dd_mu6;
        ar & rsd_dd_mu8;
      }
    
    
  };


#endif //LSSEFT_ONE_LOOP_PK_H
