//
// Created by David Seery on 05/12/2016.
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

#ifndef LSSEFT_SQLITE3_POWER_SPECTRUM_H
#define LSSEFT_SQLITE3_POWER_SPECTRUM_H


#include "database/transaction_manager.h"
#include "database/tokens.h"

#include "cosmology/concepts/power_spectrum.h"

#include "units/Mpc_units.h"

#include "sqlite3_policy.h"

#include "utilities.h"
#include "exceptions.h"

#include "localizations/messages.h"

#include "boost/optional.hpp"

#include "sqlite3.h"


namespace sqlite3_operations
  {
    
    //! lookup ID for linear power spectrum
    template <typename PkContainer>
    boost::optional<unsigned int> lookup_Pk_linear(sqlite3* db, transaction_manager& mgr, const FRW_model_token& model,
                                                   const PkContainer& Pk_lin, const sqlite3_policy& policy)
      {
        assert(db != nullptr);
        
        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id, mid, md5_hash FROM " << tokenization_table<linear_Pk_token>(policy) << " WHERE "
          << "path = @p;";
        
        // prepare SQL statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, select_stmt.str().c_str(), select_stmt.str().length()+1, &stmt, nullptr));
        
        // bind values to the parameters
        const std::string path = Pk_lin.get_path().string();
        check_stmt(db, sqlite3_bind_text(stmt, sqlite3_bind_parameter_index(stmt, "@p"), path.c_str(), path.length(), SQLITE_STATIC));
        
        // execute statement and step through results
        int status;
        boost::optional<unsigned int> id = boost::none;
        while((status = sqlite3_step(stmt)) != SQLITE_DONE)
          {
            if(status == SQLITE_ROW)
              {
                if(id) throw runtime_exception(exception_type::database_error, ERROR_SQLITE3_MULTIPLE_PK_LINEAR);
                
                unsigned int model_id = static_cast<unsigned int>(sqlite3_column_int(stmt, 1));
                std::string md5(reinterpret_cast<const char*>(sqlite3_column_text(stmt, 2)));
                
                if(model_id != model.get_id())
                  {
                    std::ostringstream msg;
                    msg << ERROR_SQLITE3_PK_LINEAR << " '" << path << "' "
                        << ERROR_SQLITE3_PK_LINEAR_WRONG_MODEL << " " << model.get_id();
                    throw runtime_exception(exception_type::runtime_error, msg.str());
                  }
                
                if(md5 != Pk_lin.get_MD5_hash())
                  {
                    std::ostringstream msg;
                    msg << ERROR_SQLITE3_PK_LINEAR << " '" << path << "' "
                        << ERROR_SQLITE3_PK_LINEAR_WRONG_MD5;
                    throw runtime_exception(exception_type::runtime_error, msg.str());
                  }
                
                id = static_cast<unsigned int>(sqlite3_column_int(stmt, 0));
              }
          }
        
        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));
        
        return id;
      }
    
    
    //! insert a linear power spectrum configuration
    template <typename PkContainer>
    unsigned int insert_Pk_linear(sqlite3* db, transaction_manager& mgr, const FRW_model_token& model,
                                  const PkContainer& Pk_lin, const sqlite3_policy& policy)
      {
        assert(db != nullptr);
        
        // get number of rows in table; this will be the identifier for the new power spectrum record
        unsigned int new_id = count(db, tokenization_table<linear_Pk_token>(policy));
        
        std::ostringstream insert_stmt;
        insert_stmt
          << "INSERT INTO " << tokenization_table<linear_Pk_token>(policy) << " VALUES (@id, @mid, @path, @md5_hash);";
        
        // prepare SQL statement
        sqlite3_stmt* stmt;
        check_stmt(db, sqlite3_prepare_v2(db, insert_stmt.str().c_str(), insert_stmt.str().length()+1, &stmt, nullptr));
        
        // bind values to the parameters
        const std::string path = Pk_lin.get_path().string();
        const std::string md5 = Pk_lin.get_MD5_hash();
        
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@id"), new_id));
        check_stmt(db, sqlite3_bind_int(stmt, sqlite3_bind_parameter_index(stmt, "@mid"), model.get_id()));
        check_stmt(db, sqlite3_bind_text(stmt, sqlite3_bind_parameter_index(stmt, "@path"), path.c_str(), path.length(), SQLITE_STATIC));
        check_stmt(db, sqlite3_bind_text(stmt, sqlite3_bind_parameter_index(stmt, "@md5_hash"), md5.c_str(), md5.length(), SQLITE_STATIC));
        
        // perform insertion
        check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_PK_LINEAR_CONFIG_FAIL, SQLITE_DONE);
        
        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));
        
        return(new_id);
      }
    
    
  }   // namespace sqlite3_operations

#endif //LSSEFT_SQLITE3_POWER_SPECTRUM_H
