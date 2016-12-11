//
// Created by David Seery on 21/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#include "Matsubara_XY.h"


Matsubara_XY::Matsubara_XY(const linear_Pk_token& Pkt, const IR_resum_token& IRt, const Mpc_units::inverse_energy2& _X,
                           const Mpc_units::inverse_energy2& _Y)
  : IR_resum_tok(IRt),
    Pk_lin(Pkt),
    X(_X),
    Y(_Y)
  {
  }


Matsubara_XY::Matsubara_XY()
  : IR_resum_tok(0),
    Pk_lin(0),
    X(0.0),
    Y(0.0)
  {
  }
