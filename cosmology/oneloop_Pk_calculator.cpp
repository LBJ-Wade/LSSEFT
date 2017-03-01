//
// Created by David Seery on 14/11/2016.
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

#include "oneloop_Pk_calculator.h"


std::list<oneloop_Pk>
oneloop_Pk_calculator::calculate_dd(const Mpc_units::energy& k, const k_token& k_tok, const IR_cutoff_token& IR_tok,
                                    const UV_cutoff_token& UV_tok, const oneloop_growth& gf_factors,
                                    const loop_integral& loop_data, const initial_filtered_Pk& Pk_init,
                                    const boost::optional<const final_filtered_Pk&>& Pk_final)
  {
    std::list<oneloop_Pk> container;
    
    // use spline to evaluate initial and final linear power spectra at scale k
    Pk_value Ptr_init = build_Pk_value(k, Pk_init);
    
    boost::optional<linear_Pk_token> final_tok;
    if(Pk_final) final_tok = Pk_final->get_token();
    
    boost::optional<Pk_value> Ptr_final;
    if(Pk_final) Ptr_final = build_Pk_value(k, *Pk_final);
    
    for(const oneloop_value& val : gf_factors)
      {
        dd_Pk dd = this->compute_dd(k, val.second, loop_data, Ptr_init, Ptr_final);

        rsd_dd_Pk rsd_mu0 = this->compute_rsd_dd_mu0(k, val.second, loop_data, Ptr_init, Ptr_final);
        rsd_dd_Pk rsd_mu2 = this->compute_rsd_dd_mu2(k, val.second, loop_data, Ptr_init, Ptr_final);
        rsd_dd_Pk rsd_mu4 = this->compute_rsd_dd_mu4(k, val.second, loop_data, Ptr_init, Ptr_final);
        rsd_dd_Pk rsd_mu6 = this->compute_rsd_dd_mu6(k, val.second, loop_data, Ptr_init, Ptr_final);
        rsd_dd_Pk rsd_mu8 = this->compute_rsd_dd_mu8(k, val.second, loop_data, Ptr_init, Ptr_final);

        container.emplace_back(k_tok, Pk_init.get_token(), final_tok, IR_tok, UV_tok, val.first.get_id(), dd,
                               rsd_mu0, rsd_mu2, rsd_mu4, rsd_mu6, rsd_mu8);
      }
    
    return container;
  }


dd_Pk
oneloop_Pk_calculator::compute_dd(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                  const loop_integral& loop_data, const Pk_value& Ptr_init,
                                  const boost::optional<Pk_value>& Ptr_final)
  {
    delta_22_integrals d22 = loop_data.get_delta22();
    delta_13_integrals d13 = loop_data.get_delta13();
    
    // TREE TERMS
    
    Pk_value tree = val.g * val.g * (Ptr_final ? *Ptr_final : Ptr_init);
    
    // 13 TERMS
    
    Pk_value P13 = Ptr_init * val.g * ( (val.D - val.J) * d13.get_D()
                                      + val.E           * d13.get_E()
                                      + (val.F + val.J) * d13.get_F()
                                      + val.G           * d13.get_G()
                                      - val.J           * d13.get_J1()
                                      + val.J           * d13.get_J2() / 2.0 );
    
    // 22 TERMS
    
    Pk_value P22 = val.A*val.A   * d22.get_AA()
                   + val.A*val.B * d22.get_AB()
                   + val.B*val.B * d22.get_BB();
    
    // COUNTERTERMS
    
    k2_Pk_value Z2_delta = -2.0
                           * (-18.0*val.D - 28.0*val.E + 7.0*val.F + 2.0*val.G + 13.0*val.J)
                           * val.g * k*k * Ptr_init;
    
    return std::move(dd_Pk(tree, P13, P22, Z2_delta));
  }


rsd_dd_Pk
oneloop_Pk_calculator::compute_rsd_dd_mu0(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                          const loop_integral& loop_data, const Pk_value& Ptr_init,
                                          const boost::optional<Pk_value>& Ptr_final)
  {
    delta_22_integrals d22 = loop_data.get_delta22();
    delta_13_integrals d13 = loop_data.get_delta13();
    
    rsd_22_integrals rsd22 = loop_data.get_rsd22();
    rsd_13_integrals rsd13 = loop_data.get_rsd13();
    
    // TREE TERMS
    
    Pk_value tree = val.g * val.g * (Ptr_final ? *Ptr_final : Ptr_init);
    
    // 13 TERMS
    
    Pk_value P13_simple = Ptr_init * val.g * ( (val.D - val.J) * d13.get_D()
                                             + val.E           * d13.get_E()
                                             + (val.F + val.J) * d13.get_F()
                                             + val.G           * d13.get_G()
                                             - val.J           * d13.get_J1()
                                             + val.J           * d13.get_J2() / 2.0 );
    
    Pk_value P13_tensor;    // no tensor term at mu^0
    
    Pk_value P13 = P13_simple + P13_tensor;
    
    // 22 TERMS
    
    Pk_value P22_simple = val.A*val.A   * d22.get_AA()
                          + val.A*val.B * d22.get_AB()
                          + val.B*val.B * d22.get_BB();
    
    Pk_value P22_tensor;    // no tensor term at mu^0
    
    Pk_value P22 = P22_simple + P22_tensor;
    
    // COUNTERTERMS
    
    k2_Pk_value Z2_delta = 2.0 * val.g * (18.0*val.D + 28.0*val.E - 7.0*val.F - 2.0*val.G - 13.0*val.J) * k*k * Ptr_init;

    Pk_value Z0_v;    // no contribution at mu^0
    
    k2_Pk_value Z2_v;   // no contribution at mu^0
    
    Pk_value Z0_vdelta;   // no contribution at mu^0
    
    k2_Pk_value Z2_vdelta;    // no contribution at mu^0
    
    k2_Pk_value Z2_vv;    // no contribution at mu^0
    
    k2_Pk_value Z2_vvdelta;   // no contribution at mu^0
    
    k2_Pk_value Z2_vvv;   // no contribution at mu^0
    
    return rsd_dd_Pk(tree, P13, P22, Z2_delta, Z0_v, Z2_v, Z0_vdelta, Z2_vdelta, Z2_vv, Z2_vvdelta, Z2_vvv);
  }


rsd_dd_Pk
oneloop_Pk_calculator::compute_rsd_dd_mu2(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                          const loop_integral& loop_data, const Pk_value& Ptr_init,
                                          const boost::optional<Pk_value>& Ptr_final)
  {
    delta_22_integrals d22 = loop_data.get_delta22();
    delta_13_integrals d13 = loop_data.get_delta13();
    
    rsd_22_integrals rsd22 = loop_data.get_rsd22();
    rsd_13_integrals rsd13 = loop_data.get_rsd13();
    
    // TREE TERMS
    
    Pk_value tree = 2.0 * val.f * (val.g * val.g * (Ptr_final ? *Ptr_final : Ptr_init));
    
    // 13 TERMS
    
    Pk_value P13_simple = Ptr_init * val.g * ( (val.D*(val.fD + val.f) - (val.fJ + val.f)*val.J)                                             * d13.get_D()
                                           + val.E*(val.fE + val.f)                                                                     * d13.get_E()
                                           + (val.fF*val.F + val.fJ*val.J + val.f*(val.F - val.A*val.g + val.J))                        * d13.get_F()
                                           + (-val.B*val.f*val.g + (val.fG + val.f)*val.G)                                              * d13.get_G()
                                           + (val.A*(-val.fA + val.f)*val.g + val.f*(val.g*val.g*val.g - 2.0*val.J) - 2.0*val.fJ*val.J) * d13.get_J1() / 2.0
                                           + (val.B*(-val.fB + val.f)*val.g + (val.fJ + val.f)*val.J)                                   * d13.get_J2() / 2.0 );
    
    Pk_value P13_tensor_a = Ptr_init * val.f * val.g*val.g * (val.fA*val.A + val.fB*val.B - val.f*val.g*val.g) * rsd13.get_a() / 2.0;
    Pk_value P13_tensor_b = Ptr_init * 2.0 * (-1.0+val.f) * val.g*val.g * (val.fA*val.A + val.fB*val.B - val.f*val.g*val.g) * rsd13.get_b() / 6.0;
    Pk_value P13_tensor_c = Ptr_init * val.g*val.g * (-12.0*val.fB*val.B + val.B*(32.0+11.0*val.fB)*val.f + val.A*(16.0*val.f + val.fA*(-12.0+11.0*val.f))
                                                 - 3.0*(-4.0+val.f)*val.f*val.g*val.g) * rsd13.get_c() / (-6.0);
    Pk_value P13_tensor_d = Ptr_init * 3.0 * (-2.0+val.f) * val.g*val.g * (-val.A*val.fA - val.B*val.fB + val.f*val.g*val.g) * rsd13.get_d();
    Pk_value P13_tensor_e = Ptr_init * val.g*val.g * (64.0*val.B*val.fB + val.A*(val.fA*(48.0-11.0*val.f) - 16.0*val.f) - val.B*(32.0+11.0*val.fB)*val.f
                                                 + val.f*(-48.0+11.0*val.f)*val.g*val.g) * rsd13.get_e() / 6.0;
    Pk_value P13_tensor_f = Ptr_init * 2.0 * (-3.0*val.f) * val.g*val.g * (val.A*val.fA + val.B*val.fB - val.f*val.g*val.g) * rsd13.get_f();
    Pk_value P13_tensor_g = Ptr_init * (-4.0+val.f) * val.g*val.g * (val.A*val.fA + val.B*val.fB - val.f*val.g*val.g) * rsd13.get_g() / 2.0;
    
    Pk_value P13 = P13_simple + P13_tensor_a + P13_tensor_b + P13_tensor_c + P13_tensor_d + P13_tensor_e + P13_tensor_f + P13_tensor_g;
    
    // 22 TERMS
    
    Pk_value P22_simple = 2.0*val.A*(val.fA*val.A - val.f*val.g*val.g)                * d22.get_AA()
                          + (val.A*val.B*(val.fA + val.fB) - val.f*val.g*val.g*val.B) * d22.get_AB()
                          + 2.0*val.B*val.B*val.fB                                    * d22.get_BB();
    
    Pk_value P22_tensor_1 = (-1.0/2.0) * val.f*val.f * val.g*val.g*val.g*val.g * rsd22.get_A1();
    Pk_value P22_tensor_2 = 1.0 * val.A * val.f*val.f * val.g*val.g * rsd22.get_A2();
    Pk_value P22_tensor_3 = 2.0 * val.A * val.f * val.g*val.g * rsd22.get_A3();
    Pk_value P22_tensor_4 = 2.0 * val.B * val.f*val.f * val.g*val.g * rsd22.get_A4();
    Pk_value P22_tensor_5 = 4.0 * val.B * val.f * val.g*val.g * rsd22.get_A5();
    
    Pk_value P22 = P22_simple + P22_tensor_1 + P22_tensor_2 + P22_tensor_3 + P22_tensor_4 + P22_tensor_5;
    
    // COUNTERTERMS
    
    k2_Pk_value Z2_delta = 2.0 * val.g * (18.0*val.D*(val.fD + val.f) + 28.0*val.E*(val.fE + val.f) - 7.0*(val.fF +val.f)*val.F - 2.0*(val.fG + val.f)*val.G - 13.0*(val.fJ + val.f)*val.J) * k*k * Ptr_init;
    
    Pk_value Z0_v = -2.0 * val.g*val.g * (-2.0*val.B*val.fB + 2.0*val.B*val.f + val.A*(-val.fA + val.f) + val.f*val.g*val.g) * Ptr_init;
    
    k2_Pk_value Z2_v = 2.0 * val.g*val.g * (-12.0*val.A*val.fA - 12.0*val.B*val.fB + 5.0*val.A*val.f + 10.0*val.B*val.f + 12.0*val.f*val.g*val.g) * k*k * Ptr_init;
    
    Pk_value Z0_vdelta = -2.0 * val.g*val.g * (-2.0*val.B*val.fB + 2.0*val.B*val.f + val.A*(-val.fA + val.f) + val.f*val.g*val.g) * Ptr_init;
    
    k2_Pk_value Z2_vdelta = 2.0 * val.g*val.g * (-12.0*val.A*val.fA - 12.0*val.B*val.fB + 5.0*val.A*val.f + 10.0*val.B*val.f + 12.0*val.f*val.g*val.g) * k*k * Ptr_init;
    
    k2_Pk_value Z2_vv = -16.0 * val.f*val.g*val.g * (-val.A*val.fA - val.B*val.fB + val.f*val.g*val.g) * k*k * Ptr_init;
    
    k2_Pk_value Z2_vvdelta = 5.0 * val.f*val.f * val.g*val.g*val.g*val.g * k*k * Ptr_init;
    
    k2_Pk_value Z2_vvv;   // no contribution at mu^2
    
    return rsd_dd_Pk(tree, P13, P22, Z2_delta, Z0_v, Z2_v, Z0_vdelta, Z2_vdelta, Z2_vv, Z2_vvdelta, Z2_vvv);
  }


rsd_dd_Pk
oneloop_Pk_calculator::compute_rsd_dd_mu4(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                          const loop_integral& loop_data, const Pk_value& Ptr_init,
                                          const boost::optional<Pk_value>& Ptr_final)
  {
    delta_22_integrals d22 = loop_data.get_delta22();
    delta_13_integrals d13 = loop_data.get_delta13();
    
    rsd_22_integrals rsd22 = loop_data.get_rsd22();
    rsd_13_integrals rsd13 = loop_data.get_rsd13();
    
    // TREE TERMS
    
    Pk_value tree = val.f * val.f * (val.g * val.g * (Ptr_final ? *Ptr_final : Ptr_init));
    
    // 13 TERMS
    
    Pk_value P13_simple = Ptr_init * val.g * ( val.f*(val.D*val.fD - val.fJ*val.J)                                              * d13.get_D()
                                             + val.E*val.fE*val.f                                                               * d13.get_E()
                                             + val.f*(val.fF*val.F - val.A*val.f*val.g + val.fJ*val.J)                          * d13.get_F()
                                             + val.f*(-val.B*val.f*val.g + val.fG*val.G)                                        * d13.get_G()
                                             + val.f*(val.A*(-val.fA+val.f)*val.g + val.f*val.g*val.g*val.g - 2.0*val.fJ*val.J) * d13.get_J1() / 2.0
                                             + val.f*(val.B*(-val.fB+val.f)*val.g + val.fJ*val.J)                               * d13.get_J2() / 2.0 );
    
    Pk_value P13_tensor_a = Ptr_init * val.f * (1.0+val.f) * val.g*val.g * (val.A*val.fA + val.B*val.fB - val.f*val.g*val.g) * rsd13.get_a() / 2.0;
    Pk_value P13_tensor_b = Ptr_init * 2.0 * (-1.0+val.f) * val.f * val.g*val.g * (val.A*val.fA + val.B*val.fB - val.f*val.g*val.g) * rsd13.get_b();
    Pk_value P13_tensor_c = Ptr_init * val.f * val.g*val.g * (19.0*val.B*val.fB + val.B*(32.0+11.0*val.fB)*val.f
                                                         + val.A*(16.0*val.f+val.fA*(3.0+11.0*val.f)) + val.f*(-3.0+5.0*val.f)*val.g*val.g) * rsd13.get_c() / (-6.0);
    Pk_value P13_tensor_d = Ptr_init * 3.0 * (-3.0+val.f) * val.f * val.g*val.g * (-val.A*val.fA - val.B*val.fB + val.f*val.g*val.g) * rsd13.get_d();
    Pk_value P13_tensor_e = Ptr_init * val.f * val.g*val.g * (val.B*(-85.0*val.fB + 32.0*val.f + 11.0*val.fB*val.f)
                                                         + val.A*(16.0*val.f + val.fA*(-69.0+11.0*val.f))
                                                         + (69.0-11.0*val.f)*val.f*val.g*val.g) * rsd13.get_e() / (-6.0);
    Pk_value P13_tensor_f = Ptr_init * 2.0 * (-5.0+val.f) * val.f * val.g*val.g * (val.A*val.fA + val.B*val.fB - val.f*val.g*val.g) * rsd13.get_f();
    Pk_value P13_tensor_g = Ptr_init * (-7.0+val.f) * val.f * val.g*val.g * (val.A*val.fA + val.B*val.fB - val.f*val.g*val.g) * rsd13.get_g() / 2.0;
    
    Pk_value P13 = P13_simple + P13_tensor_a + P13_tensor_b + P13_tensor_c + P13_tensor_d + P13_tensor_e + P13_tensor_f + P13_tensor_g;
    
    // 22 TERMS
    
    Pk_value P22_simple = (val.fA*val.A - val.f*val.g*val.g)*(val.fA*val.A - val.f*val.g*val.g) * d22.get_AA()
                          + val.fB*val.B*(val.fA*val.A - val.f*val.g*val.g)                     * d22.get_AB()
                          + val.fB*val.fB*val.B*val.B                                           * d22.get_BB();
    
    Pk_value P22_tensor_1 = -1.0 * val.g*val.g*val.g*val.g * val.f*val.f*val.f * rsd22.get_A1();
    Pk_value P22_tensor_2 = (3.0/8.0) * val.g*val.g*val.g*val.g * val.f*val.f*val.f*val.f * rsd22.get_B2();
    Pk_value P22_tensor_3 = (-1.0/2.0) * val.g*val.g*val.g*val.g * val.f*val.f * rsd22.get_B3();
    Pk_value P22_tensor_4 = 2.0 * val.g*val.g * val.A * val.fA * val.f * rsd22.get_A3();
    Pk_value P22_tensor_5 = 1.0 * val.g*val.g * val.A * val.fA * val.f*val.f * rsd22.get_A2();
    Pk_value P22_tensor_6 = 4.0 * val.g*val.g * val.B * val.fB * val.f * rsd22.get_B6();
    Pk_value P22_tensor_7 = 2.0 * val.g*val.g * val.B * val.fB * val.f*val.f * rsd22.get_A4();
    Pk_value P22_tensor_8 = 1.0 * val.g*val.g * val.A * val.f*val.f * rsd22.get_B8();
    Pk_value P22_tensor_9 = 2.0 * val.g*val.g * val.B * val.f*val.f * rsd22.get_B9();
    
    Pk_value P22 = P22_simple + P22_tensor_1 + P22_tensor_2 + P22_tensor_3 + P22_tensor_4 + P22_tensor_5 + P22_tensor_6 + P22_tensor_7 + P22_tensor_8 + P22_tensor_9;
    
    // COUNTERTERMS
    
    k2_Pk_value Z2_delta = 2.0 * val.f * val.g * (18.0*val.D*val.fD + 28.0*val.E*val.fE - 7.0*val.fF*val.F - 2.0*val.fG*val.G - 13.0*val.fJ*val.J) * k*k * Ptr_init;
    
    Pk_value Z0_v = -2.0 * val.f * val.g*val.g * (-2.0*val.B*val.fB + 2.0*val.B*val.f + val.A*(-val.fA + val.f) + val.f*val.g*val.g) * Ptr_init;
    
    k2_Pk_value Z2_v = 2.0 * val.f * val.g*val.g * (-12.0*val.A*val.fA - 12.0*val.B*val.fB + 5.0*val.A*val.f + 10.0*val.B*val.f + 12.0*val.f*val.g*val.g) * k*k * Ptr_init;
    
    Pk_value Z0_vdelta = -2.0 * val.f * val.g*val.g * (-2.0*val.B*val.fB + 2.0*val.B*val.f + val.A*(-val.fA + val.f) + val.f*val.g*val.g) * Ptr_init;
    
    k2_Pk_value Z2_vdelta = 2.0 * val.f * val.g*val.g * (-12.0*val.A*val.fA - 12.0*val.B*val.fB + 5.0*val.A*val.f + 10.0*val.B*val.f + 12.0*val.f*val.g*val.g) * k*k * Ptr_init;
    
    k2_Pk_value Z2_vv = 2.0 * val.f * val.g*val.g * (2.0*val.B*val.fB*(3.0 + 4.0*val.f) + val.A*(val.fA + 8.0*val.fA*val.f)
                                                     - val.f*(1.0 + 8.0*val.f)*val.g*val.g) * k*k * Ptr_init;
    
    k2_Pk_value Z2_vvdelta = 5.0 * val.f*val.f*val.f * val.g*val.g*val.g*val.g * k*k * Ptr_init;
    
    k2_Pk_value Z2_vvv = 5.0 * val.f*val.f*val.f * val.g*val.g*val.g*val.g * k*k * Ptr_init;
    
    return rsd_dd_Pk(tree, P13, P22, Z2_delta, Z0_v, Z2_v, Z0_vdelta, Z2_vdelta, Z2_vv, Z2_vvdelta, Z2_vvv);
  }


rsd_dd_Pk
oneloop_Pk_calculator::compute_rsd_dd_mu6(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                          const loop_integral& loop_data, const Pk_value& Ptr_init,
                                          const boost::optional<Pk_value>& Ptr_final)
  {
    delta_22_integrals d22 = loop_data.get_delta22();
    delta_13_integrals d13 = loop_data.get_delta13();
    
    rsd_22_integrals rsd22 = loop_data.get_rsd22();
    rsd_13_integrals rsd13 = loop_data.get_rsd13();
    
    // TREE TERMS
    
    Pk_value tree;  // no tree contribution at mu^6
    
    // 13 TERMS
    
    Pk_value P13_simple;  // no simple contributions at mu^6
    
    auto t_factor = val.f*val.f * val.g*val.g * (val.A*val.fA + val.B*val.fB - val.f*val.g*val.g);
    
    Pk_value P13_tensor_a = Ptr_init * t_factor * rsd13.get_a() / 2.0;
    Pk_value P13_tensor_d = Ptr_init * 3.0 * t_factor * rsd13.get_d();
    Pk_value P13_tensor_e = Ptr_init * (7.0/2.0) * t_factor * rsd13.get_e();
    Pk_value P13_tensor_f = Ptr_init * (-4.0) * t_factor * rsd13.get_f();
    Pk_value P13_tensor_g = Ptr_init * (-3.0/2.0) * t_factor * rsd13.get_g();
    
    Pk_value P13_tensor_c = Ptr_init * val.f*val.f * val.g*val.g * (15.0*val.A*val.fA + 31.0*val.B*val.fB + val.f*(-15.0+8.0*val.f)*val.g*val.g) * rsd13.get_c() / (-6.0);
    
    Pk_value P13 = P13_simple + P13_tensor_a + P13_tensor_c + P13_tensor_d + P13_tensor_e + P13_tensor_f + P13_tensor_g;
    
    // 22 TERMS
    
    Pk_value P22_simple;  // no simple terms at mu^6
    
    Pk_value P22_tensor1 = 1.0 * val.A * val.fA * val.f*val.f * val.g*val.g * rsd22.get_C1();
    Pk_value P22_tensor2 = 2.0 * val.B * val.fB * val.f*val.f * val.g*val.g * rsd22.get_C2();
    Pk_value P22_tensor3 = 1.0 * val.f*val.f*val.f * val.g*val.g*val.g*val.g * rsd22.get_A1();
    Pk_value P22_tensor4 = (-1.0/4.0) * val.f*val.f*val.f*val.f * val.g*val.g*val.g*val.g * rsd22.get_C4();
    
    Pk_value P22 = P22_simple + P22_tensor1 + P22_tensor2 + P22_tensor3 + P22_tensor4;
    
    // COUNTERTERMS
    
    k2_Pk_value Z2_delta;   // no contribution at mu^6
    
    Pk_value Z0_v;    // no contribution at mu^6
    
    k2_Pk_value Z2_v;   // no contribution at mu^6
    
    Pk_value Z0_vdelta;   // no contribution at mu^6
    
    k2_Pk_value Z2_vdelta;    // no contribution at mu^6
    
    k2_Pk_value Z2_vv = 2.0 * val.f*val.f * val.g*val.g * (val.A*val.fA + 6.0*val.B*val.fB - val.f*val.g*val.g) * k*k * Ptr_init;
    
    k2_Pk_value Z2_vvdelta;   // no contribution at mu^6
    
    k2_Pk_value Z2_vvv = 5.0 * val.f*val.f*val.f*val.f * val.g*val.g*val.g*val.g * k*k * Ptr_init;
    
    return rsd_dd_Pk(tree, P13, P22, Z2_delta, Z0_v, Z2_v, Z0_vdelta, Z2_vdelta, Z2_vv, Z2_vvdelta, Z2_vvv);
  }


rsd_dd_Pk
oneloop_Pk_calculator::compute_rsd_dd_mu8(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                          const loop_integral& loop_data, const Pk_value& Ptr_init,
                                          const boost::optional<Pk_value>& Ptr_final)
  {
    delta_22_integrals d22 = loop_data.get_delta22();
    delta_13_integrals d13 = loop_data.get_delta13();
    
    rsd_22_integrals rsd22 = loop_data.get_rsd22();
    rsd_13_integrals rsd13 = loop_data.get_rsd13();
    
    // TREE TERMS
    
    Pk_value tree;  // no tree contribution at mu^8
    
    // 13 TERMS
    
    Pk_value P13;   // no 13 contribution at mu^8
    
    // 22 TERMS
    
    Pk_value P22_simple;    // no simple terms at mu^8
    
    Pk_value P22_tensor = (1.0/8.0) * val.f*val.f*val.f*val.f * val.g*val.g*val.g*val.g * rsd22.get_D1();
    
    Pk_value P22 = P22_simple + P22_tensor;
    
    // COUNTERTERMS
    
    k2_Pk_value Z2_delta;   // no contribution at mu^8
    
    Pk_value Z0_v;    // no contribution at mu^8
    
    k2_Pk_value Z2_v;   // no contribution at mu^8
    
    Pk_value Z0_vdelta;   // no contribution at mu^8
    
    k2_Pk_value Z2_vdelta;    // no contribution at mu^8
    
    k2_Pk_value Z2_vv;    // no contribution at mu^8
    
    k2_Pk_value Z2_vvdelta;   // no contribution at mu^8
    
    k2_Pk_value Z2_vvv;   // no contribution at mu^8
    
    return rsd_dd_Pk(tree, P13, P22, Z2_delta, Z0_v, Z2_v, Z0_vdelta, Z2_vdelta, Z2_vv, Z2_vvdelta, Z2_vvv);
  }


oneloop_resum_Pk
oneloop_Pk_calculator::calculate_resum_dd(const Mpc_units::energy& k, const Matsubara_XY& XY, const oneloop_Pk& data,
                                          const oneloop_growth_record& gf_data, const initial_filtered_Pk& init_Pk,
                                          const boost::optional<const final_filtered_Pk&>& final_Pk)
  {
    const auto& input_dd = data.get_dd();
    const auto& input_tree = input_dd.get_tree();
    const auto& input_13 = input_dd.get_13();
    const auto& input_22 = input_dd.get_22();
    const auto& input_Z2_delta = input_dd.get_Z2_delta();
    
    // compute Matsubara suppression factor
    double MatsubaraA   = k*k * gf_data.g*gf_data.g * XY;
    double MatsubaraExp = std::exp(-MatsubaraA);
    
    resum_Pk_value Ptree       = input_tree.get_nowiggle()     + MatsubaraExp*input_tree.get_wiggle();
    resum_Pk_value P13         = input_13.get_nowiggle()       + MatsubaraExp*input_13.get_wiggle();
    resum_Pk_value P22         = input_22.get_nowiggle()       + MatsubaraExp*input_22.get_wiggle();
    k2_resum_Pk_value Z2_delta = input_Z2_delta.get_nowiggle() + MatsubaraExp*input_Z2_delta.get_wiggle();

    resum_Pk_value SPT_nowiggle = input_tree.get_nowiggle() + input_13.get_nowiggle() + input_22.get_nowiggle();
    resum_Pk_value SPT_wiggle   = input_tree.get_wiggle()   + input_13.get_wiggle()   + input_22.get_wiggle();
    resum_Pk_value subtraction  = MatsubaraA * input_tree.get_wiggle();
    
    resum_Pk_value P1loop_SPT   = SPT_nowiggle + MatsubaraExp * (SPT_wiggle + subtraction);
    
    resum_dd_Pk Pk_resum(Ptree, P13, P22, P1loop_SPT, Z2_delta);
    
    return oneloop_resum_Pk(data.get_k_token(), data.get_init_Pk_token(), data.get_final_Pk_token(),
                            data.get_IR_token(), data.get_UV_token(), data.get_z_token(), XY.get_IR_resum_token(),
                            Pk_resum);
  }
