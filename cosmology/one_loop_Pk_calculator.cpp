//
// Created by David Seery on 14/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#include "one_loop_Pk_calculator.h"


std::list<one_loop_Pk> one_loop_Pk_calculator::calculate(const Mpc_units::energy& k, const k_token& k_tok, const IR_token& IR_tok,
                                                         const UV_token& UV_tok, const oneloop_growth& gf_factors,
                                                         const loop_integral& loop_data, const tree_power_spectrum& Ptree)
  {
    std::list<one_loop_Pk> container;
    
    for(const oneloop_value& val : gf_factors)
      {
        dd_Pk dd = this->compute_dd(k, val, loop_data, Ptree);

        rsd_dd_Pk rsd_mu0 = this->compute_rsd_dd_mu0(k, val, loop_data, Ptree);
        rsd_dd_Pk rsd_mu2 = this->compute_rsd_dd_mu2(k, val, loop_data, Ptree);
        rsd_dd_Pk rsd_mu4 = this->compute_rsd_dd_mu4(k, val, loop_data, Ptree);
        rsd_dd_Pk rsd_mu6 = this->compute_rsd_dd_mu6(k, val, loop_data, Ptree);
        rsd_dd_Pk rsd_mu8 = this->compute_rsd_dd_mu8(k, val, loop_data, Ptree);
    
        container.emplace_back(k_tok, UV_tok, IR_tok, val.first.get_id(), dd,
                               rsd_mu0, rsd_mu2, rsd_mu4, rsd_mu6, rsd_mu8);
      }
    
    return container;
  }


dd_Pk
one_loop_Pk_calculator::compute_dd(const Mpc_units::energy& k, const oneloop_value& val, const loop_integral& loop_data,
                                   const tree_power_spectrum& Ptree)
  {
    return dd_Pk();
  }


rsd_dd_Pk one_loop_Pk_calculator::compute_rsd_dd_mu0(const Mpc_units::energy& k, const oneloop_value& val,
                                                     const loop_integral& loop_data, const tree_power_spectrum& Ptree)
  {
    return rsd_dd_Pk();
  }


rsd_dd_Pk one_loop_Pk_calculator::compute_rsd_dd_mu2(const Mpc_units::energy& k, const oneloop_value& val,
                                                     const loop_integral& loop_data, const tree_power_spectrum& Ptree)
  {
    return rsd_dd_Pk();
  }


rsd_dd_Pk one_loop_Pk_calculator::compute_rsd_dd_mu4(const Mpc_units::energy& k, const oneloop_value& val,
                                                     const loop_integral& loop_data, const tree_power_spectrum& Ptree)
  {
    return rsd_dd_Pk();
  }


rsd_dd_Pk one_loop_Pk_calculator::compute_rsd_dd_mu6(const Mpc_units::energy& k, const oneloop_value& val,
                                                     const loop_integral& loop_data, const tree_power_spectrum& Ptree)
  {
    return rsd_dd_Pk();
  }


rsd_dd_Pk one_loop_Pk_calculator::compute_rsd_dd_mu8(const Mpc_units::energy& k, const oneloop_value& val,
                                                     const loop_integral& loop_data, const tree_power_spectrum& Ptree)
  {
    return rsd_dd_Pk();
  }

