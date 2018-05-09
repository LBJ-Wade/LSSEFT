//
// Created by David Seery on 10/11/2016.
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

#ifndef LSSEFT_SQLITE3_FIND_H
#define LSSEFT_SQLITE3_FIND_H


#include "sqlite3_policy.h"

#include "utilities.h"
#include "temporary_tables.h"

#include "database/transaction_manager.h"
#include "database/tokens.h"
#include "database/z_database.h"
#include "database/k_database.h"

#include "cosmology/concepts/oneloop_growth.h"
#include "cosmology/concepts/loop_integral.h"
#include "cosmology/concepts/oneloop_Pk.h"
#include "cosmology/concepts/Matsubara_XY.h"
#include "cosmology/concepts/power_spectrum.h"

#include "sqlite3.h"


namespace sqlite3_operations
  {
    
    //! extract one-loop growth g- and f-functions for a given set of redshifts
    std::unique_ptr<oneloop_growth>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& token,
             const growth_params_token& params, const z_database& z_db);
    
    //! extract loop integrals for a given wavenumber, linear power spectrum, UV-cutoff and IR-cutoff combination
    std::unique_ptr<loop_integral>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
         const loop_integral_params_token& params, const k_token& k, const linear_Pk_token& Pk,
         const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff);
    
    //! extract P(k) data for a given wavenumber, z-value, linear power spectrum, UV-cutoff and IR-cutoff combination
    std::unique_ptr<oneloop_Pk_set>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
             const growth_params_token& growth_params, const loop_integral_params_token& loop_params, const k_token& k,
             const z_token& z, const linear_Pk_token& init_Pk_lin, const boost::optional<linear_Pk_token>& final_Pk_lin,
             const IR_cutoff_token& IR_cutoff, const UV_cutoff_token& UV_cutoff);
    
    //! extract Matsubara X & Y coefficient sfor a given linear power spectrum and IR resummation scale
    std::unique_ptr<Matsubara_XY>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
         const MatsubaraXY_params_token& params, const linear_Pk_token& Pk, const IR_resum_token& IR_resum);
    
    
    //! extract filtered linear power spectrum (of given type Payload) for a given set of k-modes
    template <typename Payload>
    std::unique_ptr<Payload>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const linear_Pk_token& token,
         const k_database& k_db)
      {
        // set up temporary table of desired k identifiers
        std::string ktab = k_table(db, mgr, policy, k_db);
        
        std::ostringstream read_stmt;
        read_stmt << "SELECT sample.Pk_raw, sample.Pk_nw "
                  << "FROM " << ktab << " "
                  << "INNER JOIN (SELECT * FROM " << policy.Pk_linear_table() << " WHERE Pk_id=@Pk_id) sample "
                  << "ON " << ktab << ".id = sample.kid "
                  << "ORDER BY " << ktab << ".ROWID ASC;";
        
        // prepare statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, read_stmt.str().c_str(), read_stmt.str().length()+1, &stmt, nullptr));
        
        // bind parameter values
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@Pk_id"), token.get_id()));
        
        // set up databases to hold the result
        tree_Pk::database_type raw_db;
        tree_Pk_w::database_type nw_db;
        
        // perform read
        int result = 0;
        k_database::const_record_iterator t = k_db.record_cbegin();
        while((result = sqlite3_step(stmt)) != SQLITE_DONE && t != k_db.record_cend())
          {
            if(result == SQLITE_ROW)
              {
                Mpc_units::inverse_energy3 raw = sqlite3_column_double(stmt, 0) * dimensionful_unit<Mpc_units::inverse_energy3>();
                Mpc_units::inverse_energy3 nw  = sqlite3_column_double(stmt, 1) * dimensionful_unit<Mpc_units::inverse_energy3>();
                
                raw_db.add_record(*(*t), raw);
                nw_db.add_record(*(*t), nw);
                
                ++t;
              }
            else
              {
                check_stmt(db, sqlite3_clear_bindings(stmt));
                check_stmt(db, sqlite3_finalize(stmt));
                
                throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_READ_FILTERED_PK_FAIL);
              }
          }
        
        // clear bindings and release
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));
        
        // drop temporary table
        drop_temp(db, mgr, ktab);
        
        if(t != k_db.record_cend()) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_FILTERED_PK_MISREAD);
        
        auto payload = std::make_unique<Payload>(token, nw_db, raw_db);
        
        return std::move(payload);
      }
    
  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_FIND_H
