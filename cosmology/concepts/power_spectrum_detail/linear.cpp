//
// Created by David Seery on 05/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#include <fstream>

#include "linear.h"

#include "openssl/md5.h"


linear_power_spectrum::linear_power_spectrum(const boost::filesystem::path& p)
  : path(p.is_absolute() ? p : boost::filesystem::absolute(p)),
    container(path)    // need to be sure path is initialized before container
  {
    this->md5_hash = this->hash(p);
  }


std::string linear_power_spectrum::hash(const boost::filesystem::path& p)
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
