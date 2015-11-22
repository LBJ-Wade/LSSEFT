//
// Created by David Seery on 12/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
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

        double           omega_m  = obj.get_omega_m();
        double           omega_cc = obj.get_omega_cc();
        double           h        = obj.get_h();
        Mpc_units::energy T_CMB    = obj.get_T_CMB();

        double T_CMB_in_Kelvin = T_CMB / Mpc_units::Kelvin;

        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id FROM " << policy.FRW_model_table() << " WHERE "
          << "ABS((omega_m-@om)/omega_m)<@tol "
          << "AND ABS((omega_cc-@occ)/omega_cc)<@tol "
          << "AND ABS((h-@h)/h)<@tol "
          << "AND ABS((T_CMB-@Tcmb)/T_CMB)<@tol;";

        // prepare SQL statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));

        // bind values to the parameters
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@tol"), tol));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@om"), omega_m));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@occ"), omega_cc));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@h"), h));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@Tcmb"), T_CMB_in_Kelvin));

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

        double           omega_m  = obj.get_omega_m();
        double           omega_cc = obj.get_omega_cc();
        double           h        = obj.get_h();
        Mpc_units::energy T_CMB    = obj.get_T_CMB();

        double T_CMB_in_Kelvin = T_CMB / Mpc_units::Kelvin;

        std::ostringstream insert_stmt;
        insert_stmt
          << "INSERT INTO " << policy.FRW_model_table() << " VALUES (@id, @omega_m, @omega_cc, @h, @T_CMB);";

        // prepare SQL statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));

        // bind values to the parameters
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@id"), new_id));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@omega_m"), omega_m));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@omega_cc"), omega_cc));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@h"), h));
        check_stmt(db, sqlite3_bind_double(stmt, sqlite3_bind_parameter_index(stmt, "@T_CMB"), T_CMB_in_Kelvin));

        // perform insertion
        check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_FRW_MODEL_FAIL, SQLITE_DONE);

        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));

        return(new_id);
      }

  }   // namespace sqlite3_operations
