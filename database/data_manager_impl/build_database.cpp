//
// Created by David Seery on 09/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#include <iostream>
#include <sstream>
#include <assert.h>

#include <set>
#include <unordered_set>

#include "database/data_manager.h"
#include "database/data_manager_impl/types.h"

#include "sqlite3_detail/utilities.h"
#include "sqlite3_detail/operations.h"

#include "utilities/formatter.h"

#include "defaults.h"

#include "boost/timer/timer.hpp"


// GENERATE REDSHIFT AND WAVENUMBER DATABASES


std::unique_ptr<z_database> data_manager::build_redshift_db(range<double>& sample)
  {
    // construct an empty redshift database
    std::unique_ptr<z_database> z_db(new z_database);
    
    // grab the grid of redshift samples
    const std::vector<double>& z_samples = sample.grid();
    
    for(std::vector<double>::const_iterator t = z_samples.begin(); t != z_samples.end(); ++t)
      {
        std::unique_ptr<z_token> tok = this->tokenize(*t);
        z_db->add_record(*t, *tok);
      }
    
    return(z_db);
  }


std::unique_ptr<k_database> data_manager::build_k_db(range<Mpc_units::energy>& sample)
  {
    return this->build_wavenumber_db<k_token>(sample);
  }


std::unique_ptr<IR_cutoff_database> data_manager::build_IR_cutoff_db(range<Mpc_units::energy>& sample)
  {
    return this->build_wavenumber_db<IR_cutoff_token>(sample);
  }


std::unique_ptr<UV_cutoff_database> data_manager::build_UV_cutoff_db(range<Mpc_units::energy>& sample)
  {
    return this->build_wavenumber_db<UV_cutoff_token>(sample);
  }


std::unique_ptr<IR_resum_database> data_manager::build_IR_resum_db(range<Mpc_units::energy>& sample)
  {
    return this->build_wavenumber_db<IR_resum_token>(sample);
  }


std::unique_ptr<k_database> data_manager::build_k_db(transaction_manager& mgr, const linear_Pk& Pk_lin,
                                                     double bottom_clearance, double top_clearance)
  {
    // construct an empty wavenumber database
    std::unique_ptr<k_database> k_db = std::make_unique<k_database>();
    
    // get power spectrum database underlying this container
    const tree_Pk::database_type& Pk_db = Pk_lin.get_db();
    
    for(tree_Pk::database_type::const_record_iterator t = Pk_db.record_cbegin(); t != Pk_db.record_cend(); ++t)
      {
        // ask linear_Pk container whether this P(k) value is acceptable
        const Mpc_units::energy& k = t->get_wavenumber();
        if(Pk_lin.is_valid(k, bottom_clearance, top_clearance))
          {
            std::unique_ptr<k_token> tok = this->tokenize<k_token>(mgr, k);
            k_db->add_record(k, *tok);
          }
      }
    
    return k_db;
  }
