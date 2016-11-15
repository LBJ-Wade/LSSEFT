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
        dd_Pk dd = this->compute_dd(k, val.second, loop_data, Ptree);

        rsd_dd_Pk rsd_mu0 = this->compute_rsd_dd_mu0(k, val.second, loop_data, Ptree);
        rsd_dd_Pk rsd_mu2 = this->compute_rsd_dd_mu2(k, val.second, loop_data, Ptree);
        rsd_dd_Pk rsd_mu4 = this->compute_rsd_dd_mu4(k, val.second, loop_data, Ptree);
        rsd_dd_Pk rsd_mu6 = this->compute_rsd_dd_mu6(k, val.second, loop_data, Ptree);
        rsd_dd_Pk rsd_mu8 = this->compute_rsd_dd_mu8(k, val.second, loop_data, Ptree);
    
        container.emplace_back(k_tok, UV_tok, IR_tok, val.first.get_id(), dd,
                               rsd_mu0, rsd_mu2, rsd_mu4, rsd_mu6, rsd_mu8);
      }
    
    return container;
  }


dd_Pk
one_loop_Pk_calculator::compute_dd(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                   const loop_integral& loop_data, const tree_power_spectrum& Ptree)
  {
    delta_22_integrals d22 = loop_data.get_delta22();
    delta_13_integrals d13 = loop_data.get_delta13();
    
    Pk_value tree = val.g * val.g * Ptree(k);
    
    Pk_value P13 = Ptree(k) * val.g * ( (val.D - val.J)   * d13.get_D()
                                        + val.E           * d13.get_E()
                                        + (val.F + val.J) * d13.get_F()
                                        + val.G           * d13.get_G()
                                        + (val.J / 2.0)   * (d13.get_J2() - 2.0 * d13.get_J1()) );
    
    Pk_value P22 = val.A * val.A   * d22.get_AA()
                   + val.A * val.B * d22.get_AB()
                   + val.B * val.B * d22.get_BB();
    
    k2_Pk_value Z2_delta = -2.0
                           * (-18.0*val.D - 28.0*val.E + 7.0*val.F + 2.0*val.G + 13.0*val.J)
                           * val.g * k*k * Ptree(k);
    
    return std::move(dd_Pk(tree, P13, P22, Z2_delta));
  }


rsd_dd_Pk one_loop_Pk_calculator::compute_rsd_dd_mu0(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                                     const loop_integral& loop_data, const tree_power_spectrum& Ptree)
  {
    return rsd_dd_Pk();
  }


rsd_dd_Pk one_loop_Pk_calculator::compute_rsd_dd_mu2(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                                     const loop_integral& loop_data, const tree_power_spectrum& Ptree)
  {
    return rsd_dd_Pk();
  }


rsd_dd_Pk one_loop_Pk_calculator::compute_rsd_dd_mu4(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                                     const loop_integral& loop_data, const tree_power_spectrum& Ptree)
  {
    return rsd_dd_Pk();
  }


rsd_dd_Pk one_loop_Pk_calculator::compute_rsd_dd_mu6(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                                     const loop_integral& loop_data, const tree_power_spectrum& Ptree)
  {
    return rsd_dd_Pk();
  }


rsd_dd_Pk one_loop_Pk_calculator::compute_rsd_dd_mu8(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                                     const loop_integral& loop_data, const tree_power_spectrum& Ptree)
  {
    return rsd_dd_Pk();
  }

