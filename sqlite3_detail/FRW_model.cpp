//
// Created by David Seery on 12/08/2015.
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


#include <sstream>
#include <assert.h>

#include "utilities.h"
#include "FRW_model.h"

#include "exceptions.h"
#include "localizations/messages.h"


namespace sqlite3_operations
  {

    boost::optional<unsigned int> lookup_FRW_model(sqlite3* db, transaction_manager& mgr, const FRW_model& obj, const sqlite3_policy& policy, double tol)
      {
        assert(db != nullptr);

        double            omega_m  = obj.get_omega_m();
        double            omega_cc = obj.get_omega_cc();
        double            h        = obj.get_h();
        Mpc_units::energy T_CMB    = obj.get_T_CMB();
        double            Neff     = obj.get_Neff();
        double            f_baryon = obj.get_f_baryon();
        double            z_star   = obj.get_z_star();
        double            z_drag   = obj.get_z_drag();
        double            z_eq     = obj.get_z_eq();
        double            A_curv   = obj.get_A_curv();
        double            ns       = obj.get_ns();
        Mpc_units::energy k_piv    = obj.get_k_piv();
    
        double T_CMB_in_Kelvin = T_CMB / Mpc_units::Kelvin;
        double k_piv_in_invMpc = k_piv * Mpc_units::Mpc;

        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id FROM " << policy.FRW_model_table() << " WHERE "
          << "ABS((omega_m-@om)/omega_m)<@tol "
          << "AND ABS((omega_cc-@occ)/omega_cc)<@tol "
          << "AND ABS((h-@h)/h)<@tol "
          << "AND ABS((T_CMB-@Tcmb)/T_CMB)<@tol "
          << "AND ABS((Neff-@Neff)/Neff)<@tol "
          << "AND ABS((f_baryon-@f_baryon)/f_baryon)<@tol "
          << "AND ABS((z_star-@z_star)/z_star)<@tol "
          << "AND ABS((z_drag-@z_drag)/z_drag)<@tol "
          << "AND ABS((z_eq-@z_eq)/z_eq)<@tol "
          << "AND ABS((A_curv-@A_curv)/A_curv)<@tol "
          << "AND ABS((ns-@ns)/ns)<@tol "
          << "AND ABS((k_piv-@k_piv)/k_piv)<@tol;";

        // prepare SQL statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));

        // bind values to the parameters
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@tol"), tol));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@om"), omega_m));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@occ"), omega_cc));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@h"), h));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Tcmb"), T_CMB_in_Kelvin));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Neff"), Neff));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@f_baryon"), f_baryon));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@z_star"), z_star));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@z_drag"), z_drag));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@z_eq"), z_eq));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@A_curv"), A_curv));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@ns"), ns));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@k_piv"), k_piv_in_invMpc));

        // execute statement and step through results
        int status = 0;
        boost::optional<unsigned int> id = boost::none;
        while((status = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(status == SQLITE_ROW)
              {
                if(id) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_MULTIPLE_FRW_MODELS);
                id = static_cast<unsigned int>(sqlite3_column_int(stmt, 0));
              }
          }

        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));

        return(id);
      }


    unsigned int insert_FRW_model(sqlite3* db, transaction_manager& mgr, const FRW_model& obj, const sqlite3_policy& policy)
      {
        assert(db != nullptr);

        // get number of rows in table; this will be the identifier for the new model
        unsigned int new_id = count(db, policy.FRW_model_table());

        double            omega_m  = obj.get_omega_m();
        double            omega_cc = obj.get_omega_cc();
        double            h        = obj.get_h();
        Mpc_units::energy T_CMB    = obj.get_T_CMB();
        double            Neff     = obj.get_Neff();
        double            f_baryon = obj.get_f_baryon();
        double            z_star   = obj.get_z_star();
        double            z_drag   = obj.get_z_drag();
        double            z_eq     = obj.get_z_eq();
        double            A_curv   = obj.get_A_curv();
        double            ns       = obj.get_ns();
        Mpc_units::energy k_piv    = obj.get_k_piv();

        double T_CMB_in_Kelvin = T_CMB / Mpc_units::Kelvin;
        double k_piv_in_invMpc = k_piv * Mpc_units::Mpc;

        std::ostringstream insert_stmt;
        insert_stmt
          << "INSERT INTO " << policy.FRW_model_table() << " VALUES (@id, @omega_m, @omega_cc, @h, @T_CMB, @Neff, @f_baryon, @z_star, @z_drag, @z_eq, @A_curv, @ns, @k_piv);";

        // prepare SQL statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));

        // bind values to the parameters
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@id"), new_id));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@omega_m"), omega_m));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@omega_cc"), omega_cc));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@h"), h));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@T_CMB"), T_CMB_in_Kelvin));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Neff"), Neff));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@f_baryon"), f_baryon));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@z_star"), z_star));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@z_drag"), z_drag));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@z_eq"), z_eq));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@A_curv"), A_curv));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@ns"), ns));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@k_piv"), k_piv_in_invMpc));

        // perform insertion
        check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_FRW_MODEL_FAIL, SQLITE_DONE);

        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));

        return(new_id);
      }

  }   // namespace sqlite3_operations
