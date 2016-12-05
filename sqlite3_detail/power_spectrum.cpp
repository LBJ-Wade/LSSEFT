//
// Created by David Seery on 05/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#include "power_spectrum.h"


namespace sqlite3_operations
  {
    
    boost::optional<unsigned int> lookup_Pk_linear(sqlite3* db, transaction_manager& mgr, const FRW_model_token& model,
                                                   const linear_power_spectrum& Pk_lin, const sqlite3_policy& policy)
      {
        assert(db != nullptr);
        
        std::ostringstream select_stmt;
        select_stmt
          << "SELECT id, mid, md5_hash FROM " << tokenization_table<Pk_linear_token>(policy) << " WHERE "
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
        
        // finalize statement and release resoureces
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));

        return id;
      }
    
    
    unsigned int insert_Pk_linear(sqlite3* db, transaction_manager& mgr, const FRW_model_token& model,
                                  const linear_power_spectrum& Pk_lin, const sqlite3_policy& policy)
      {
        assert(db != nullptr);
        
        // get number of rows in table; this will be the identifier for the new power spectrum record
        unsigned int new_id = count(db, tokenization_table<Pk_linear_token>(policy));
        
        std::ostringstream insert_stmt;
        insert_stmt
          << "INSERT INTO " << tokenization_table<Pk_linear_token>(policy) << " VALUES (@id, @mid, @path, @md5_hash);";
    
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
        check_stmt(db, sqlite3_step(stmt), ERROR_SQLITE3_INSERT_PK_LINEAR_FAIL, SQLITE_DONE);
    
        // finalize statement and release resources
        check_stmt(db, sqlite3_clear_bindings(stmt));
        check_stmt(db, sqlite3_finalize(stmt));
    
        return(new_id);
      }

  }   // namespace sqlite3_operations
