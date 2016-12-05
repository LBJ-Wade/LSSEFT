//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "tokens.h"


generic_token::generic_token(unsigned int i)
  : id(i)
  {
  }


bool operator<(const generic_token& a, const generic_token& b)
  {
    return(a.id < b.id);
  }


bool operator==(const generic_token& a, const generic_token& b)
  {
    return(a.id == b.id);
  }


FRW_model_token::FRW_model_token(unsigned int i)
  : generic_token(i)
  {
  }


z_token::z_token(unsigned int i)
  : generic_token(i)
  {
  }


k_token::k_token(unsigned int i)
  : generic_token(i)
  {
  }


IR_cutoff_token::IR_cutoff_token(unsigned int i)
  : generic_token(i)
  {
  }


UV_cutoff_token::UV_cutoff_token(unsigned int i)
  : generic_token(i)
  {
  }


IR_resum_token::IR_resum_token(unsigned int i)
  : generic_token(i)
  {
  }


linear_Pk_token::linear_Pk_token(unsigned int i)
  : generic_token(i)
  {
  }
