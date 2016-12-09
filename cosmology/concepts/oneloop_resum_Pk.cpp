//
// Created by David Seery on 09/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#include "oneloop_resum_Pk.h"


oneloop_resum_Pk::oneloop_resum_Pk(const k_token& kt, const linear_Pk_token& Pkt, const IR_cutoff_token& IRt,
                                   const UV_cutoff_token& UVt, const z_token& zt, const IR_resum_token& IRrt,
                                   const resum_dd_Pk& Pkr)
  : k(kt),
    Pk_lin(Pkt),
    IR_cutoff(IRt),
    UV_cutoff(UVt),
    z(zt),
    IR_resum(IRrt),
    Pk_resum(Pkr)
  {
  }


oneloop_resum_Pk::oneloop_resum_Pk()
  : k(0.0),
    Pk_lin(0),
    IR_cutoff(0),
    UV_cutoff(0),
    z(0),
    IR_resum(0),
    Pk_resum()
  {
  }
