//
// Created by David Seery on 21/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#include "Matsubara_A.h"


Matsubara_A::Matsubara_A(const IR_resum_token& _IRt, const Mpc_units::inverse_energy2& _A)
  : IR_resum_tok(_IRt),
    A(_A)
  {
  }


Matsubara_A::Matsubara_A()
  : IR_resum_tok(0),
    A(0.0)
  {
  }
