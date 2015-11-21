//
// Created by David Seery on 10/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_MPI_OPERATIONS_H
#define LSSEFT_MPI_OPERATIONS_H


#include <memory>

#include "cosmology/FRW_model.h"
#include "cosmology/transfer_integrator.h"
#include "units/eV_units.h"
#include "database/z_database.h"

#include "boost/serialization/serialization.hpp"
#include "boost/serialization/shared_ptr.hpp"


namespace MPI_detail
  {

    constexpr unsigned int RANK_MASTER = 0;


    // MESSAGE TAGS

    constexpr unsigned int MESSAGE_NEW_TRANSFER_TASK          = 0;
    constexpr unsigned int MESSAGE_NEW_TRANSFER_INTEGRATION   = 1;
    constexpr unsigned int MESSAGE_TRANSFER_INTEGRATION_READY = 2;

    constexpr unsigned int MESSAGE_WORKER_READY               = 90;
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
        new_transfer_integration(const FRW_model& m, const eV_units::energy& _k, const k_token& t, const std::shared_ptr<z_database>& z)
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
        const eV_units::energy& get_k() const { return(this->k); }

        //! get wavenumber token
        const k_token& get_token() const { return(this->token); }

        //! get redshift database
        std::shared_ptr<z_database> get_z_db() const { return(this->z_db); }


        // INTERNAL DATA

      private:

        //! FRW model to use for the integration
        FRW_model model;

        //! wavenumber to integrate
        eV_units::energy k;

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
        //! transfer_function, eV_units::energy and k_token have no default constructors
        transfer_integration_ready()
          : data(eV_units::energy(0), k_token(0), std::shared_ptr<z_database>())
          {
          }

        //! value constructor: used to construct and send a payload
        transfer_integration_ready(const transfer_function& f)
          : data(f)
          {
          }


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

  }   // namespace MPI


#endif //LSSEFT_MPI_OPERATIONS_H
