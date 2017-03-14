//
// Created by David Seery on 11/08/2015.
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


filter_data_token::filter_data_token(unsigned int i)
  : generic_token(i)
  {
  }