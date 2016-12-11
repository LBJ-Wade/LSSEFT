//
// Created by David Seery on 05/12/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#include "filtered_Pk.h"


filtered_Pk::filtered_Pk(const k_token& kt, const linear_Pk_token& Pt, Mpc_units::inverse_energy3 _Pk_nw,
                         Mpc_units::inverse_energy3 _Pk_raw, Mpc_units::inverse_energy3 _Pk_ref)
  : k_tok(kt),
    Pk_tok(Pt),
    Pk_nw(std::move(_Pk_nw)),
    Pk_raw(std::move(_Pk_raw)),
    Pk_ref(std::move(_Pk_ref))
  {
  }


filtered_Pk::filtered_Pk()
  : k_tok(0),
    Pk_tok(0),
    Pk_nw(0.0),
    Pk_raw(0.0),
    Pk_ref(0.0)
  {
  }
