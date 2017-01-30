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
         const z_database& z_db);
    
    //! extract filtered linear power spectrum for a given set of k-moes
    std::unique_ptr<wiggle_Pk>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const linear_Pk_token& token,
         const k_database& k_db);
    
    //! extract loop integrals for a given wavenumber, linear power spectrum, UV-cutoff and IR-cutoff combination
    std::unique_ptr<loop_integral>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
         const k_token& k, const linear_Pk_token& Pk, const IR_cutoff_token& IR_cutoff,
         const UV_cutoff_token& UV_cutoff);
    
    //! extract P(k) data for a given wavenumber, z-value, linear power spectrum, UV-cutoff and IR-cutoff combination
    std::unique_ptr<oneloop_Pk>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
             const k_token& k, const z_token& z, const linear_Pk_token& Pk_lin, const IR_cutoff_token& IR_cutoff,
             const UV_cutoff_token& UV_cutoff);
    
    //! extract Matsubara X & Y coefficient sfor a given linear power spectrum and IR resummation scale
    std::unique_ptr<Matsubara_XY>
    find(sqlite3* db, transaction_manager& mgr, const sqlite3_policy& policy, const FRW_model_token& model,
         const linear_Pk_token& Pk, const IR_resum_token& IR_resum);
    
  }   // namespace sqlite3_operations


#endif //LSSEFT_SQLITE3_FIND_H
