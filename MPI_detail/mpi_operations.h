//
// Created by David Seery on 10/08/2015.
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

#ifndef LSSEFT_MPI_OPERATIONS_H
#define LSSEFT_MPI_OPERATIONS_H


#include <memory>
#include <cosmology/concepts/oneloop_growth.h>

#include "cosmology/FRW_model.h"
#include "cosmology/concepts/transfer_function.h"
#include "cosmology/concepts/loop_integral.h"
#include "cosmology/concepts/oneloop_Pk.h"
#include "cosmology/concepts/multipole_Pk.h"
#include "cosmology/concepts/Matsubara_XY.h"
#include "cosmology/concepts/filtered_Pk_value.h"

#include "units/Mpc_units.h"
#include "database/z_database.h"

#include "boost/serialization/serialization.hpp"
#include "boost/serialization/shared_ptr.hpp"
#include "boost/serialization/map.hpp"


namespace MPI_detail
  {

    constexpr unsigned int RANK_MASTER = 0;


    // MESSAGE TAGS

    constexpr unsigned int MESSAGE_NEW_TRANSFER_TASK          = 0;
    constexpr unsigned int MESSAGE_NEW_TRANSFER_INTEGRATION   = 1;
    
    constexpr unsigned int MESSAGE_NEW_FILTER_PK_TASK         = 10;
    constexpr unsigned int MESSAGE_NEW_FILTER_PK              = 21;

    constexpr unsigned int MESSAGE_NEW_LOOP_INTEGRAL_TASK     = 20;
    constexpr unsigned int MESSAGE_NEW_LOOP_INTEGRATION       = 21;
    
    constexpr unsigned int MESSAGE_NEW_ONE_LOOP_PK_TASK       = 30;
    constexpr unsigned int MESSAGE_NEW_ONE_LOOP_PK            = 31;
    
    constexpr unsigned int MESSAGE_NEW_MATSUBARA_XY_TASK      = 40;
    constexpr unsigned int MESSAGE_NEW_MATSUBARA_XY           = 51;
    
    constexpr unsigned int MESSAGE_NEW_MULTIPOLE_PK_TASK      = 50;
    constexpr unsigned int MESSAGE_NEW_MULTIPOLE_PK           = 51;

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
          : model("", 0, 0, 0, Mpc_units::energy(0), 0, 0, 0, 0, 0, 0, 0, Mpc_units::energy(0)),
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
    
    
    // FILTERING OF LINEAR POWER SPECTRUM


    class new_filter_Pk
      {
        
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! empty constructor, used to receive a payload
        new_filter_Pk()
          : model("", 0, 0, 0, Mpc_units::energy(0), 0, 0, 0, 0, 0, 0, 0, Mpc_units::energy(0)),
            k(0),
            k_tok(0),
            Pk_tok(0),
            Pk_lin(),
            params_tok(0),
            params()
          {
          }
           
        //! value constructor, used to construct and send a payload
        new_filter_Pk(const FRW_model& m, const Mpc_units::energy& _k, const k_token& kt, const linear_Pk_token& Pt,
                      const std::shared_ptr<filterable_Pk>& _Pk, const filter_params_token& pt, const Pk_filter_params& p)
          : model(m),
            k(_k),
            k_tok(kt),
            Pk_tok(Pt),
            Pk_lin(_Pk),
            params_tok(pt),
            params(p)
          {
          }
        
        //! destructor is default
        ~new_filter_Pk() = default;
        
        
        // ACCESS PAYLOAD
        
      public:
    
        //! get model
        const FRW_model& get_model() const { return(this->model); }
    
        //! get wavenumber
        const Mpc_units::energy& get_k() const { return(this->k); }
    
        //! get wavenumber token
        const k_token& get_k_token() const { return(this->k_tok); }
    
        //! get power specturm token
        const linear_Pk_token& get_Pk_token() const { return this->Pk_tok; }
        
        //! get linear power specturm container
        const filterable_Pk& get_Pk_linear() const { return *this->Pk_lin; }
        
        //! get filtering parameters token
        const filter_params_token& get_params_token() const { return this->params_tok; }
        
        //! get filtering parameters
        const Pk_filter_params& get_params() const { return this->params; }
    
    
        // INTERNAL DATA
        
      private:
    
        //! FRW model to use for this calculation
        FRW_model model;
    
        //! wavenumber to integrate
        Mpc_units::energy k;
    
        //! wavenumber token
        k_token k_tok;
        
        //! power spectrum token
        linear_Pk_token Pk_tok;
        
        //! linear power spectrum container
        std::shared_ptr<filterable_Pk> Pk_lin;
        
        //! token for filtering parameters
        filter_params_token params_tok;
        
        //! filtering parameters
        Pk_filter_params params;
    
    
        // enable boost::serialization support, and hence automated packing for transmission over MPI
        friend class boost::serialization::access;
    
        template <typename Archive>
        void serialize(Archive& ar, unsigned int version)
          {
            ar & model;
            ar & k;
            ar & k_tok;
            ar & Pk_tok;
            ar & Pk_lin;
            ar & params_tok;
            ar & params;
          }
    
      };
    
    
    class filter_Pk_ready
      {
        
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! empty constructor: used to receive a payload
        filter_Pk_ready()
          : data()
          {
          }
        
        //! value constructor: used to construct and send a payload
        filter_Pk_ready(const filtered_Pk_value& f)
          : data(f)
          {
          }
    
        //! destructor is default
        ~filter_Pk_ready() = default;
    
    
        // INTERFACE
  
      public:
    
        const filtered_Pk_value& get_data() const { return this->data; }
    
    
        // INTERNAL DATA
  
      private:
    
        //! loop integral container
        filtered_Pk_value data;
    
    
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
          : model("", 0, 0, 0, Mpc_units::energy(0), 0, 0, 0, 0, 0, 0, 0, Mpc_units::energy(0)),
            k(0),
            UV_cutoff(0),
            IR_cutoff(0),
            k_tok(0),
            UV_tok(0),
            IR_tok(0),
            Pk(),
            params_tok(0),
            params()
          {
          }

        //! value constructor: used to construct and send a payload
        new_loop_momentum_integration(const FRW_model& m, const Mpc_units::energy& _k, const k_token& kt,
                                      const Mpc_units::energy& UV, const UV_cutoff_token& UVt, const Mpc_units::energy& IR,
                                      const IR_cutoff_token& IRt, const std::shared_ptr<initial_filtered_Pk>& _Pk,
                                      const loop_integral_params_token& pt, const loop_integral_params& p)
          : model(m),
            k(_k),
            UV_cutoff(UV),
            IR_cutoff(IR),
            k_tok(kt),
            UV_tok(UVt),
            IR_tok(IRt),
            Pk(_Pk),
            params_tok(pt),
            params(p)
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
        const UV_cutoff_token& get_UV_token() const { return(this->UV_tok); }

        //! get IR cutoff
        const Mpc_units::energy& get_IR_cutoff() const { return(this->IR_cutoff); }

        //! get IR cutoff token
        const IR_cutoff_token& get_IR_token() const { return(this->IR_tok); }

        //! get tree-level power spectrum
        const initial_filtered_Pk& get_tree_power_spectrum() const { return *this->Pk; }
        
        //! get parameters token
        const loop_integral_params_token& get_params_token() const { return this->params_tok; }
        
        //! get parameters block
        const loop_integral_params& get_params() const { return this->params; }


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
        UV_cutoff_token UV_tok;

        //! IR cutoff to USE
        Mpc_units::energy IR_cutoff;

        //! IR cutoff token
        IR_cutoff_token IR_tok;

        //! tree-level power spectrum
        std::shared_ptr<initial_filtered_Pk> Pk;
        
        //! parameters token
        loop_integral_params_token params_tok;
        
        //! parameters block
        loop_integral_params params;


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

        //! empty constructor: used to receive a payload
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
    
    
    // MATSUBARA XY-COEFFICIENT PAYLOADS
    
    class new_Matsubara_XY
      {
        
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! empty constructor: used to receive a payload
        new_Matsubara_XY()
          : IR_resum(0.0),
            IR_resum_tok(0),
            Pk(),
            params_tok(0),
            params()
          {
          }
        
        //! value constructor: used to construct and send a payload
        new_Matsubara_XY(const Mpc_units::energy& _IR, const IR_resum_token& _IRt,
                         const std::shared_ptr<initial_filtered_Pk>& _Pk, const MatsubaraXY_params_token& _pt,
                         const MatsubaraXY_params& _pm)
          : IR_resum(_IR),
            IR_resum_tok(_IRt),
            Pk(_Pk),
            params_tok(_pt),
            params(_pm)
          {
          }
        
        //! destructor is default
        ~new_Matsubara_XY() = default;
        
        
        // ACCESS PAYLOAD
        
      public:
    
        //! get IR resummation scale
        const Mpc_units::energy& get_IR_resum() const { return this->IR_resum; }
    
        //! get IR resummation token
        const IR_resum_token& get_IR_resum_token() const { return this->IR_resum_tok; }
    
        //! get tree-level power spectrum
        const initial_filtered_Pk& get_tree_power_spectrum() const { return *this->Pk; }
        
        //! get parameters token
        const MatsubaraXY_params_token& get_params_token() const { return this->params_tok; }
        
        //! get parameters block
        const MatsubaraXY_params& get_params() const { return this->params; }
    
    
        // INTERNAL DATA
  
      private:
    
        //! IR resummation scale
        Mpc_units::energy IR_resum;
    
        //! IR resummation token
        IR_resum_token IR_resum_tok;
    
        //! tree-level power spectrum
        std::shared_ptr<initial_filtered_Pk> Pk;
        
        //! parameters token
        MatsubaraXY_params_token params_tok;
        
        //! parameters
        MatsubaraXY_params params;
    
    
        // enable boost::serialization support, and hence automated packing for transmission over MPI
        friend class boost::serialization::access;
    
        template <typename Archive>
        void serialize(Archive& ar, unsigned int version)
          {
            ar & IR_resum;
            ar & IR_resum_tok;
            ar & Pk;
            ar & params_tok;
            ar & params;
          }
    
      };
    
    
    class Matsubara_XY_ready
      {
        
        // CONSTRUCTOR, DESTRUCTOR
      
      public:
        
        //! empty constructor: used to receive a payload
        Matsubara_XY_ready()
          : data()
          {
          }
        
        //! value constructor: used to send a payload
        Matsubara_XY_ready(const Matsubara_XY _d)
          : data(_d)
          {
          }
        
        //! destructor is default
        ~Matsubara_XY_ready() = default;
        
        
        // INTERFACE
      
      public:
        
        //! accessor for payload
        const Matsubara_XY& get_data() const { return this->data; }
        
        
        // INTERNAL DATA
      
      private:
        
        //! payload
        Matsubara_XY data;
        
        
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
            Pk_init(),
            Pk_final()
          {
          }
        
        //! value constructor: used to construct and send a payload
        new_one_loop_Pk(const Mpc_units::energy& _k,
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
    
        //! get initial linear power spectrum
        const initial_filtered_Pk& get_init_linear_Pk() const { return *this->Pk_init; }
        
        //! get final linear power spectrum, if provided
        boost::optional<const final_filtered_Pk&> get_final_linear_Pk() const
          {
            if(this->Pk_final) return *this->Pk_final;
            return boost::none;
          }

        
        // INTERNAL DATA
        
        private:
    
        //! wavenumber being computed
        Mpc_units::energy k;
    
        //! growth data
        std::shared_ptr<oneloop_growth> gf_factors;
        
        //! loop kernel data
        std::shared_ptr<loop_integral> loop_data;
    
        //! initial linear power spectrum
        std::shared_ptr<initial_filtered_Pk> Pk_init;
        
        //! final linear power spectrum, if provided
        std::shared_ptr<final_filtered_Pk> Pk_final;
    
    
        // enable boost::serialization support, and hence automated packing for transmission over MPI
        friend class boost::serialization::access;
    
        template <typename Archive>
        void serialize(Archive& ar, unsigned int version)
          {
            ar & k;
            ar & gf_factors;
            ar & loop_data;
            ar & Pk_init;
            ar & Pk_final;
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
        one_loop_Pk_ready(const std::list<oneloop_Pk_set>& Pks)
          : data(Pks)
          {
          }
        
        //! destructor is default
        ~one_loop_Pk_ready() = default;
        
        
        // INTERFACE
        
      public:
        
        const std::list<oneloop_Pk_set>& get_data() const { return this->data; }
        
        
        // INTERNAL DATA
        
      private:
        
        //! one-loop Pk container
        std::list<oneloop_Pk_set> data;
    
    
        // enable boost::serialization support, and hence automated packing for transmission over MPI
        friend class boost::serialization::access;
    
        template <typename Archive>
        void serialize(Archive& ar, unsigned int version)
          {
            ar & data;
          }
        
      };

    
    class new_multipole_Pk
      {
        
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! empty constructor: used to receive a payload
        new_multipole_Pk()
          : k(0.0),
            XY(),
            data(),
            Df_data(),
            Pk_init(),
            Pk_final()
          {
          }
        
        //! value constructor: used to construct and send a payload
        new_multipole_Pk(const Mpc_units::energy& _k, const Matsubara_XY& _XY, const std::shared_ptr<oneloop_Pk_set>& _data,
                         const oneloop_growth_record& _Df_data, const std::shared_ptr<initial_filtered_Pk>& _Pk_init,
                         const std::shared_ptr<final_filtered_Pk>& _Pk_final)
          : k(_k),
            XY(_XY),
            data(_data),
            Df_data(_Df_data),
            Pk_init(_Pk_init),
            Pk_final(_Pk_final)
          {
          }
        
        //! destructor is default
        ~new_multipole_Pk() = default;
        
        
        // ACCESS PAYLOAD
        
      public:
        
        //! get wavenumber
        const Mpc_units::energy& get_k() const { return this->k; }
        
        //! get Matsubara X & Y coefficients
        const Matsubara_XY& get_Matsubara_XY() const { return this->XY; }
        
        //! get one-loop data
        const oneloop_Pk_set& get_oneloop_Pk_data() const { return *this->data; }
        
        //! get gf growth factors
        const oneloop_growth_record& get_Df_data() const { return this->Df_data; }
        
        //! get initial linear power spectrum
        const initial_filtered_Pk& get_init_linear_Pk() const { return *this->Pk_init; }
        
        //! get final linear power spectrum, if provided
        boost::optional<const final_filtered_Pk&> get_final_linear_Pk() const
          {
            if(this->Pk_final) return *this->Pk_final;
            return boost::none;
          }
    
    
        // INTERNAL DATA
  
      private:
    
        // Payload data
    
        //! physical scale k
        Mpc_units::energy k;
        
        //! Matsubara X & Y coefficients
        Matsubara_XY XY;
    
        //! one-loop power spectrum data
        std::shared_ptr<oneloop_Pk_set> data;
        
        //! gf growth factors
        oneloop_growth_record Df_data;
    
        //! initial linear power spectrum
        std::shared_ptr<initial_filtered_Pk> Pk_init;
        
        //! final linear power spectrum, if provided
        std::shared_ptr<final_filtered_Pk> Pk_final;
    
    
        // enable boost::serialization support, and hence automated packing for transmission over MPI
        friend class boost::serialization::access;
    
        template <typename Archive>
        void serialize(Archive& ar, unsigned int version)
          {
            ar & k;
            ar & XY;
            ar & data;
            ar & Df_data;
            ar & Pk_init;
            ar & Pk_final;
          }
        
      };
    
    
    class multipole_Pk_ready
      {
        
        // CONSTRUCTOR, DESTRUCTOR
        
      public:
        
        //! empty constructor: used to receive a payload
        multipole_Pk_ready()
          {
          }
        
        //! value constructor: used to construct and send a payload
        multipole_Pk_ready(const multipole_Pk_set& Pk)
          : data(Pk)
          {
          }
    
        //! destructor is default
        ~multipole_Pk_ready() = default;
        
        
        // INTERFACE
        
      public:
        
        const multipole_Pk_set& get_data() const { return this->data; }
        
        
        // INTERNAL DATA
        
      private:
        
        //! multipole Pk container
        multipole_Pk_set data;
        

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
