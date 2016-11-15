//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_MPI_OPERATIONS_H
#define LSSEFT_MPI_OPERATIONS_H


#include <memory>

#include "cosmology/FRW_model.h"
#include "cosmology/concepts/transfer_function.h"
#include "cosmology/concepts/loop_integral.h"
#include "cosmology/concepts/oneloop_Pk.h"
#include "units/Mpc_units.h"
#include "database/z_database.h"

#include "boost/serialization/serialization.hpp"
#include "boost/serialization/shared_ptr.hpp"


namespace MPI_detail
  {

    constexpr unsigned int RANK_MASTER = 0;


    // MESSAGE TAGS

    constexpr unsigned int MESSAGE_NEW_TRANSFER_TASK          = 0;
    constexpr unsigned int MESSAGE_NEW_TRANSFER_INTEGRATION   = 1;

    constexpr unsigned int MESSAGE_NEW_LOOP_INTEGRAL_TASK     = 10;
    constexpr unsigned int MESSAGE_NEW_LOOP_INTEGRATION       = 11;
    
    constexpr unsigned int MESSAGE_NEW_ONE_LOOP_PK_TASK       = 20;
    constexpr unsigned int MESSAGE_NEW_ONE_LOOP_PK            = 21;

    constexpr unsigned int MESSAGE_WORKER_READY               = 90;
    constexpr unsigned int MESSAGE_WORK_PRODUCT_READY         = 91;

    constexpr unsigned int MESSAGE_END_OF_WORK                = 98;
    constexpr unsigned int MESSAGE_END_OF_WORK_ACK            = 99;

    constexpr unsigned int MESSAGE_TERMINATE                  = 999;


    // TRANSFER INTEGRATION PAYLOADS


    class new_transfer_integration
      {

        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! empty constructor: used to receive a payload
        new_transfer_integration()
          : model(),
            k(0),       // note k has no default constructor
            token(0),   // note token has no default constructor
            z_db()
          {
          }

        //! value constructor: used to construct and send a payload
        new_transfer_integration(const FRW_model& m, const Mpc_units::energy& _k, const k_token& t, const std::shared_ptr<z_database>& z)
          : model(m),
            k(_k),
            token(t),
            z_db(z)
          {
          }

        //! destructor is default
        ~new_transfer_integration() = default;


        // ACCESS PAYLOAD

      public:

        //! get model
        const FRW_model& get_model() const { return(this->model); }

        //! get wavenumber
        const Mpc_units::energy& get_k() const { return(this->k); }

        //! get wavenumber token
        const k_token& get_token() const { return(this->token); }

        //! get redshift database
        const z_database& get_z_db() const { return *this->z_db; }


        // INTERNAL DATA

      private:

        //! FRW model to use for the integration
        FRW_model model;

        //! wavenumber to integrate
        Mpc_units::energy k;

        //! wavenumber token
        k_token token;

        //! redshifts to sample; use shared_ptr to avoid costly copies in case
        //! z_db is large
        std::shared_ptr<z_database> z_db;


        // enable boost::serialization support, and hence automated packing for transmission over MPI
        friend class boost::serialization::access;

        template <typename Archive>
        void serialize(Archive& ar, unsigned int version)
          {
            ar & model;
            ar & k;
            ar & token;
            ar & z_db;
          }

      };


    class transfer_integration_ready
      {

        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! empty constructor: used to receive a payload
        //! transfer_function, Mpc_units::energy and k_token have no default constructors
        transfer_integration_ready()
          : data(Mpc_units::energy(0), k_token(0), std::make_shared<z_database>())
          {
          }

        //! value constructor: used to construct and send a payload
        transfer_integration_ready(const transfer_function& f)
          : data(f)
          {
          }

        //! destructor is default
        ~transfer_integration_ready() = default;


        // INTERFACE

      public:

        const transfer_function& get_data() const { return(this->data); }


        // INTERNAL DATA

      private:

        //! transfer function container
        transfer_function data;


        // enable boost::serialization support, and hence automated packing for transmission over MPI
        friend class boost::serialization::access;

        template <typename Archive>
        void serialize(Archive& ar, unsigned int version)
          {
            ar & data;
          }

      };


    // LOOP INTEGRAL PAYLOADS


    class new_loop_momentum_integration
      {

        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! empty constructor: used to receive a payload
        new_loop_momentum_integration()
          : model(),
            k(0),
            UV_cutoff(0),
            IR_cutoff(0),
            k_tok(0),
            UV_tok(0),
            IR_tok(0),
            Pk()
          {
          }

        //! value constructor: used to construct and send a payload
        new_loop_momentum_integration(const FRW_model& m, const Mpc_units::energy& _k, const k_token& kt,
                                      const Mpc_units::energy& UV, const UV_token& UVt,
                                      const Mpc_units::energy& IR, const IR_token& IRt,
                                      const std::shared_ptr<tree_power_spectrum>& _Pk)
          : model(m),
            k(_k),
            UV_cutoff(UV),
            IR_cutoff(IR),
            k_tok(kt),
            UV_tok(UVt),
            IR_tok(IRt),
            Pk(_Pk)
          {
          }

        //! destructor is default
        ~new_loop_momentum_integration() = default;


        // ACCESS PAYLOAD

      public:

        //! get model
        const FRW_model& get_model() const { return(this->model); }

        //! get wavenumber
        const Mpc_units::energy& get_k() const { return(this->k); }

        //! get wavenumber token
        const k_token& get_k_token() const { return(this->k_tok); }

        //! get UV cutoff
        const Mpc_units::energy& get_UV_cutoff() const { return(this->UV_cutoff); }

        //! get UV cutoff token
        const UV_token& get_UV_token() const { return(this->UV_tok); }

        //! get IR cutoff
        const Mpc_units::energy& get_IR_cutoff() const { return(this->IR_cutoff); }

        //! get IR cutoff token
        const IR_token& get_IR_token() const { return(this->IR_tok); }

        //! get tree-level power spectrum
        const tree_power_spectrum& get_tree_power_spectrum() const { return *this->Pk; }


        // INTERNAL DATA

      private:

        //! FRW model to use for the integration
        FRW_model model;

        //! wavenumber to integrate
        Mpc_units::energy k;

        //! wavenumber token
        k_token k_tok;

        //! UV cutoff to use
        Mpc_units::energy UV_cutoff;

        //! UV cutoff token
        UV_token UV_tok;

        //! IR cutoff to USE
        Mpc_units::energy IR_cutoff;

        //! IR cutoff token
        IR_token IR_tok;

        //! tree-level power spectrum
        std::shared_ptr<tree_power_spectrum> Pk;


        // enable boost::serialization support, and hence automated packing for transmission over MPI
        friend class boost::serialization::access;

        template <typename Archive>
        void serialize(Archive& ar, unsigned int version)
          {
            ar & model;
            ar & k;
            ar & k_tok;
            ar & UV_cutoff;
            ar & UV_tok;
            ar & IR_cutoff;
            ar & IR_tok;
            ar & Pk;
          }

      };


    class loop_momentum_integration_ready
      {

        // CONSTRUCTOR, DESTRUCTOR

      public:

        //! empty constructor: used to receive a payload;
        //! note Mpc_units::energy, k_token, IR_token and UV_token have no default constructor
        loop_momentum_integration_ready()
          : data()
          {
          }

        //! value constructor: used to construct and send a payload
        loop_momentum_integration_ready(const loop_integral& l)
          : data(l)
          {
          }

        //! destructor is default
        ~loop_momentum_integration_ready() = default;


        // INTERFACE

      public:

        const loop_integral& get_data() const { return(this->data); }


        // INTERNAL DATA

      private:

        //! loop integral container
        loop_integral data;


        // enable boost::serialization support, and hence automated packing for transmission over MPI
        friend class boost::serialization::access;

        template <typename Archive>
        void serialize(Archive& ar, unsigned int version)
          {
            ar & data;
          }

      };
    
    
    // POWER SPECTRUM PAYLOADS
    
    
    class new_one_loop_Pk
      {
        
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! empty constructor: used to receive a payload
        new_one_loop_Pk()
          : k(0.0),
            gf_factors(),
            loop_data(),
            Pk()
          {
          }
        
        //! value constructor: used to construct and send a payload
        new_one_loop_Pk(const Mpc_units::energy& _k,
                        const std::shared_ptr<oneloop_growth>& gf, const std::shared_ptr<loop_integral>& k,
                        const std::shared_ptr<tree_power_spectrum>& _Pk)
          : k(_k),
            gf_factors(gf),
            loop_data(k),
            Pk(_Pk)
          {
          }
        
        //! destructor is default
        ~new_one_loop_Pk() = default;
        
        
        // ACCESS PAYLOAD
        
      public:
    
        //! get wavenumber
        const Mpc_units::energy& get_k() const { return(this->k); }
        
        //! get growth data
        const oneloop_growth& get_gf_factors() const { return *this->gf_factors; }
        
        //! get one-loop kernel data
        const loop_integral& get_loop_data() const { return *this->loop_data; }
    
        //! get tree-level power spectrum
        const tree_power_spectrum& get_tree_power_spectrum() const { return *this->Pk; }

        
        // INTERNAL DATA
        
        private:
    
        //! wavenumber being computed
        Mpc_units::energy k;
    
        //! growth data
        std::shared_ptr<oneloop_growth> gf_factors;
        
        //! loop kernel data
        std::shared_ptr<loop_integral> loop_data;
    
        //! tree-level power spectrum
        std::shared_ptr<tree_power_spectrum> Pk;
    
    
        // enable boost::serialization support, and hence automated packing for transmission over MPI
        friend class boost::serialization::access;
    
        template <typename Archive>
        void serialize(Archive& ar, unsigned int version)
          {
            ar & k;
            ar & gf_factors;
            ar & loop_data;
            ar & Pk;
          }
    
      };
    
    
    class one_loop_Pk_ready
      {
        
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! empty constructor: used to receive a payload
        one_loop_Pk_ready()
          : data()
          {
          }
        
        //! value constructor: used to construct and send a payload
        one_loop_Pk_ready(const oneloop_Pk& Pk)
          : data(Pk)
          {
          }
        
        //! destructor is default
        ~one_loop_Pk_ready() = default;
        
        
        // INTERFACE
        
      public:
        
        const oneloop_Pk& get_data() const { return this->data; }
        
        
        // INTERNAL DATA
        
      private:
        
        //! one-loop Pk container
        oneloop_Pk data;
    
    
        // enable boost::serialization support, and hence automated packing for transmission over MPI
        friend class boost::serialization::access;
    
        template <typename Archive>
        void serialize(Archive& ar, unsigned int version)
          {
            ar & data;
          }
        
      };

  }   // namespace MPI


#endif //LSSEFT_MPI_OPERATIONS_H
