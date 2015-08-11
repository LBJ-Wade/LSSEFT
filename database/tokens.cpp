//
// Created by David Seery on 11/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "tokens.h"


generic_token::generic_token(unsigned int i)
  : id(i)
  {
  }


FRW_model_token::FRW_model_token(unsigned int i)
  : generic_token(id)
  {
  }
