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

#include <fstream>

#include "linear.h"

#include "openssl/md5.h"


linear_Pk::linear_Pk(const boost::filesystem::path& p)
  : path(p.is_absolute() ? p : boost::filesystem::absolute(p)),
    container(path)    // need to be sure path is initialized before container
  {
    this->md5_hash = this->hash(p);
  }


linear_Pk::linear_Pk(const std::string& p, const tree_Pk::database_type& d, const std::string& h)
  : path(p),
    container(d),
    md5_hash(h)
  {
  }


std::string linear_Pk::hash(const boost::filesystem::path& p)
  {
    unsigned char result[MD5_DIGEST_LENGTH];

    std::ifstream in(p.string().c_str(), std::ios_base::in);

    if(!in)
      {
        in.close();
        return std::string();   // return empty hash
      }
    
    // read in file in 1k chunks
    constexpr unsigned int BUFFER_SIZE = 64*1024;
    std::unique_ptr<char> buffer(new char[BUFFER_SIZE]);
    
    MD5_CTX md5_context;
    MD5_Init(&md5_context);
    
    while(!in.eof())
      {
        in.read(buffer.get(), BUFFER_SIZE);
        MD5_Update(&md5_context, buffer.get(), in.gcount());
      }
    
    MD5_Final(result, &md5_context);
    
    in.close();
    
    std::ostringstream hash;
    for(unsigned int i = 0; i < MD5_DIGEST_LENGTH; ++i)
      {
        hash << std::setw(2) << std::hex << static_cast<int>(result[i]);
      }
    
    return hash.str();
  }


bool linear_Pk::is_valid(const Mpc_units::energy& k, double bottom_clearance, double top_clearance) const
  {
    return this->container.is_valid(k, bottom_clearance, top_clearance);
  }
