//
// Created by David Seery on 22/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_DATA_MANAGER_IMPL_H
#define LSSEFT_DATA_MANAGER_IMPL_H


#include "k_database.h"
#include "IR_database.h"
#include "UV_database.h"


namespace data_manager_impl
  {

    class loop_momentum_configuration
      {

      public:

        loop_momentum_configuration(k_database::const_record_iterator& _k, UV_database::const_record_iterator& _UV, IR_database::const_record_iterator& _IR)
          : k(_k),
            UV(_UV),
            IR(_IR),
            include(true)
          {
          }

        ~loop_momentum_configuration() = default;

        k_database::const_record_iterator  k;
        UV_database::const_record_iterator UV;
        IR_database::const_record_iterator IR;

        bool include;
      };

  }


#endif //LSSEFT_DATA_MANAGER_IMPL_H
