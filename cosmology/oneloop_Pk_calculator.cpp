//
// Created by David Seery on 14/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#include "oneloop_Pk_calculator.h"


std::list<oneloop_Pk> oneloop_Pk_calculator::calculate(const Mpc_units::energy& k, const k_token& k_tok, const IR_cutoff_token& IR_tok,
                                                       const UV_cutoff_token& UV_tok, const oneloop_growth& gf_factors,
                                                       const loop_integral& loop_data, const tree_Pk& Ptree)
  {
    std::list<oneloop_Pk> container;
    
    for(const oneloop_value& val : gf_factors)
      {
        dd_Pk dd = this->compute_dd(k, val.second, loop_data, Ptree);

        rsd_dd_Pk rsd_mu0 = this->compute_rsd_dd_mu0(k, val.second, loop_data, Ptree);
        rsd_dd_Pk rsd_mu2 = this->compute_rsd_dd_mu2(k, val.second, loop_data, Ptree);
        rsd_dd_Pk rsd_mu4 = this->compute_rsd_dd_mu4(k, val.second, loop_data, Ptree);
        rsd_dd_Pk rsd_mu6 = this->compute_rsd_dd_mu6(k, val.second, loop_data, Ptree);
        rsd_dd_Pk rsd_mu8 = this->compute_rsd_dd_mu8(k, val.second, loop_data, Ptree);
    
        container.emplace_back(k_tok, IR_tok, UV_tok, val.first.get_id(), dd,
                               rsd_mu0, rsd_mu2, rsd_mu4, rsd_mu6, rsd_mu8);
      }
    
    return container;
  }


dd_Pk
oneloop_Pk_calculator::compute_dd(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                  const loop_integral& loop_data, const tree_Pk& Ptree)
  {
    delta_22_integrals d22 = loop_data.get_delta22();
    delta_13_integrals d13 = loop_data.get_delta13();
    
    // TREE TERMS
    
    auto Ptr = Ptree(k);
    
    Pk_value tree = val.g * val.g * Ptr;
    
    // 13 TERMS
    
    Pk_value P13 = Ptr * val.g * ( (val.D - val.J)    * d13.get_D()
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
                           * val.g * k*k * Ptr;
    
    return std::move(dd_Pk(tree, P13, P22, Z2_delta));
  }


rsd_dd_Pk oneloop_Pk_calculator::compute_rsd_dd_mu0(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                                    const loop_integral& loop_data, const tree_Pk& Ptree)
  {
    delta_22_integrals d22 = loop_data.get_delta22();
    delta_13_integrals d13 = loop_data.get_delta13();
    
    rsd_22_integrals rsd22 = loop_data.get_rsd22();
    rsd_13_integrals rsd13 = loop_data.get_rsd13();
    
    // TREE TERMS
    
    auto Ptr = Ptree(k);
    
    Pk_value tree = val.g * val.g * Ptr;
    
    // 13 TERMS
    
    Pk_value P13_simple = Ptr * val.g * ( (val.D - val.J)   * d13.get_D()
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
    
    k2_Pk_value Z2_delta = 2.0 * val.g * (18.0*val.D + 28.0*val.E - 7.0*val.F - 2.0*val.G - 13.0*val.J) * k*k * Ptr;

    Pk_value Z0_v;    // no contribution at mu^0
    
    k2_Pk_value Z2_v;   // no contribution at mu^0
    
    Pk_value Z0_vdelta;   // no contribution at mu^0
    
    k2_Pk_value Z2_vdelta;    // no contribution at mu^0
    
    k2_Pk_value Z2_vv;    // no contribution at mu^0
    
    k2_Pk_value Z2_vvdelta;   // no contribution at mu^0
    
    k2_Pk_value Z2_vvv;   // no contribution at mu^0
    
    return rsd_dd_Pk(tree, P13, P22, Z2_delta, Z0_v, Z2_v, Z0_vdelta, Z2_vdelta, Z2_vv, Z2_vvdelta, Z2_vvv);
  }


rsd_dd_Pk oneloop_Pk_calculator::compute_rsd_dd_mu2(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                                    const loop_integral& loop_data, const tree_Pk& Ptree)
  {
    delta_22_integrals d22 = loop_data.get_delta22();
    delta_13_integrals d13 = loop_data.get_delta13();
    
    rsd_22_integrals rsd22 = loop_data.get_rsd22();
    rsd_13_integrals rsd13 = loop_data.get_rsd13();
    
    // TREE TERMS
    
    auto Ptr = Ptree(k);

    Pk_value tree = 2.0 * val.g * val.g * val.f * Ptr;
    
    // 13 TERMS
    
    Pk_value P13_simple = Ptr * val.g * ( (val.D*(val.fD + val.f) - (val.fJ + val.f)*val.J)                                             * d13.get_D()
                                           + val.E*(val.fE + val.f)                                                                     * d13.get_E()
                                           + (val.fF*val.F + val.fJ*val.J + val.f*(val.F - val.A*val.g + val.J))                        * d13.get_F()
                                           + (-val.B*val.f*val.g + (val.fG + val.f)*val.G)                                              * d13.get_G()
                                           + (val.A*(-val.fA + val.f)*val.g + val.f*(val.g*val.g*val.g - 2.0*val.J) - 2.0*val.fJ*val.J) * d13.get_J1() / 2.0
                                           + (val.B*(-val.fB + val.f)*val.g + (val.fJ + val.f)*val.J)                                   * d13.get_J2() / 2.0 );
    
    Pk_value P13_tensor_a = Ptr * val.f * val.g*val.g * (val.fA*val.A + val.fB*val.B - val.f*val.g*val.g) * rsd13.get_a() / 2.0;
    Pk_value P13_tensor_b = Ptr * 2.0 * (-1.0+val.f) * val.g*val.g * (val.fA*val.A + val.fB*val.B - val.f*val.g*val.g) * rsd13.get_b() / 6.0;
    Pk_value P13_tensor_c = Ptr * val.g*val.g * (-12.0*val.fB*val.B + val.B*(32.0+11.0*val.fB)*val.f + val.A*(16.0*val.f + val.fA*(-12.0+11.0*val.f))
                                                 - 3.0*(-4.0+val.f)*val.f*val.g*val.g) * rsd13.get_c() / (-6.0);
    Pk_value P13_tensor_d = Ptr * 3.0 * (-2.0+val.f) * val.g*val.g * (-val.A*val.fA - val.B*val.fB + val.f*val.g*val.g) * rsd13.get_d();
    Pk_value P13_tensor_e = Ptr * val.g*val.g * (64.0*val.B*val.fB + val.A*(val.fA*(48.0-11.0*val.f) - 16.0*val.f) - val.B*(32.0+11.0*val.fB)*val.f
                                                 + val.f*(-48.0+11.0*val.f)*val.g*val.g) * rsd13.get_e() / 6.0;
    Pk_value P13_tensor_f = Ptr * 2.0 * (-3.0*val.f) * val.g*val.g * (val.A*val.fA + val.B*val.fB - val.f*val.g*val.g) * rsd13.get_f();
    Pk_value P13_tensor_g = Ptr * (-4.0+val.f) * val.g*val.g * (val.A*val.fA + val.B*val.fB - val.f*val.g*val.g) * rsd13.get_g() / 2.0;
    
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
    
    k2_Pk_value Z2_delta = 2.0 * val.g * (18.0*val.D*(val.fD + val.f) + 28.0*val.E*(val.fE + val.f) - 7.0*(val.fF +val.f)*val.F - 2.0*(val.fG + val.f)*val.G - 13.0*(val.fJ + val.f)*val.J) * k*k * Ptr;
    
    Pk_value Z0_v = -2.0 * val.g*val.g * (-2.0*val.B*val.fB + 2.0*val.B*val.f + val.A*(-val.fA + val.f) + val.f*val.g*val.g) * Ptr;
    
    k2_Pk_value Z2_v = 2.0 * val.g*val.g * (-12.0*val.A*val.fA - 12.0*val.B*val.fB + 5.0*val.A*val.f + 10.0*val.B*val.f + 12.0*val.f*val.g*val.g) * k*k * Ptr;
    
    Pk_value Z0_vdelta = -2.0 * val.g*val.g * (-2.0*val.B*val.fB + 2.0*val.B*val.f + val.A*(-val.fA + val.f) + val.f*val.g*val.g) * Ptr;
    
    k2_Pk_value Z2_vdelta = 2.0 * val.g*val.g * (-12.0*val.A*val.fA - 12.0*val.B*val.fB + 5.0*val.A*val.f + 10.0*val.B*val.f + 12.0*val.f*val.g*val.g) * k*k * Ptr;
    
    k2_Pk_value Z2_vv = -16.0 * val.f*val.g*val.g * (-val.A*val.fA - val.B*val.fB + val.f*val.g*val.g) * k*k * Ptr;
    
    k2_Pk_value Z2_vvdelta = 5.0 * val.f*val.f * val.g*val.g*val.g*val.g * k*k * Ptr;
    
    k2_Pk_value Z2_vvv;   // no contribution at mu^2
    
    return rsd_dd_Pk(tree, P13, P22, Z2_delta, Z0_v, Z2_v, Z0_vdelta, Z2_vdelta, Z2_vv, Z2_vvdelta, Z2_vvv);
  }


rsd_dd_Pk oneloop_Pk_calculator::compute_rsd_dd_mu4(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                                    const loop_integral& loop_data, const tree_Pk& Ptree)
  {
    delta_22_integrals d22 = loop_data.get_delta22();
    delta_13_integrals d13 = loop_data.get_delta13();
    
    rsd_22_integrals rsd22 = loop_data.get_rsd22();
    rsd_13_integrals rsd13 = loop_data.get_rsd13();
    
    // TREE TERMS
    auto Ptr = Ptree(k);
    
    
    Pk_value tree = val.g * val.g * val.f * val.f * Ptr;
    
    // 13 TERMS
    
    Pk_value P13_simple = Ptr * val.g * ( val.f*(val.D*val.fD - val.fJ*val.J)                                                * d13.get_D()
                                          + val.E*val.fE*val.f                                                               * d13.get_E()
                                          + val.f*(val.fF*val.F - val.A*val.f*val.g + val.fJ*val.J)                          * d13.get_F()
                                          + val.f*(-val.B*val.f*val.g + val.fG*val.G)                                        * d13.get_G()
                                          + val.f*(val.A*(-val.fA+val.f)*val.g + val.f*val.g*val.g*val.g - 2.0*val.fJ*val.J) * d13.get_J1() / 2.0
                                          + val.f*(val.B*(-val.fB+val.f)*val.g + val.fJ*val.J)                               * d13.get_J2() / 2.0 );
    
    Pk_value P13_tensor_a = Ptr * val.f * (1.0+val.f) * val.g*val.g * (val.A*val.fA + val.B*val.fB - val.f*val.g*val.g) * rsd13.get_a() / 2.0;
    Pk_value P13_tensor_b = Ptr * 2.0 * (-1.0+val.f) * val.f * val.g*val.g * (val.A*val.fA + val.B*val.fB - val.f*val.g*val.g) * rsd13.get_b();
    Pk_value P13_tensor_c = Ptr * val.f * val.g*val.g * (19.0*val.B*val.fB + val.B*(32.0+11.0*val.fB)*val.f
                                                         + val.A*(16.0*val.f+val.fA*(3.0+11.0*val.f)) + val.f*(-3.0+5.0*val.f)*val.g*val.g) * rsd13.get_c() / (-6.0);
    Pk_value P13_tensor_d = Ptr * 3.0 * (-3.0+val.f) * val.f * val.g*val.g * (-val.A*val.fA - val.B*val.fB + val.f*val.g*val.g) * rsd13.get_d();
    Pk_value P13_tensor_e = Ptr * val.f * val.g*val.g * (val.B*(-85.0*val.fB + 32.0*val.f + 11.0*val.fB*val.f)
                                                         + val.A*(16.0*val.f + val.fA*(-69.0+11.0*val.f))
                                                         + (69.0-11.0*val.f)*val.f*val.g*val.g) * rsd13.get_e() / (-6.0);
    Pk_value P13_tensor_f = Ptr * 2.0 * (-5.0+val.f) * val.f * val.g*val.g * (val.A*val.fA + val.B*val.fB - val.f*val.g*val.g) * rsd13.get_f();
    Pk_value P13_tensor_g = Ptr * (-7.0+val.f) * val.f * val.g*val.g * (val.A*val.fA + val.B*val.fB - val.f*val.g*val.g) * rsd13.get_g() / 2.0;
    
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
    
    k2_Pk_value Z2_delta = 2.0 * val.f * val.g * (18.0*val.D*val.fD + 28.0*val.E*val.fE - 7.0*val.fF*val.F - 2.0*val.fG*val.G - 13.0*val.fJ*val.J) * k*k * Ptr;
    
    Pk_value Z0_v = -2.0 * val.f * val.g*val.g * (-2.0*val.B*val.fB + 2.0*val.B*val.f + val.A*(-val.fA + val.f) + val.f*val.g*val.g) * Ptr;
    
    k2_Pk_value Z2_v = 2.0 * val.f * val.g*val.g * (-12.0*val.A*val.fA - 12.0*val.B*val.fB + 5.0*val.A*val.f + 10.0*val.B*val.f + 12.0*val.f*val.g*val.g) * k*k * Ptr;
    
    Pk_value Z0_vdelta = -2.0 * val.f * val.g*val.g * (-2.0*val.B*val.fB + 2.0*val.B*val.f + val.A*(-val.fA + val.f) + val.f*val.g*val.g) * Ptr;
    
    k2_Pk_value Z2_vdelta = 2.0 * val.f * val.g*val.g * (-12.0*val.A*val.fA - 12.0*val.B*val.fB + 5.0*val.A*val.f + 10.0*val.B*val.f + 12.0*val.f*val.g*val.g) * k*k * Ptr;
    
    k2_Pk_value Z2_vv = 2.0 * val.f * val.g*val.g * (2.0*val.B*val.fB*(3.0 + 4.0*val.f) + val.A*(val.fA + 8.0*val.fA*val.f)
                                                     - val.f*(1.0 + 8.0*val.f)*val.g*val.g) * k*k * Ptr;
    
    k2_Pk_value Z2_vvdelta = 5.0 * val.f*val.f*val.f * val.g*val.g*val.g*val.g * k*k * Ptr;
    
    k2_Pk_value Z2_vvv = 5.0 * val.f*val.f*val.f * val.g*val.g*val.g*val.g * k*k * Ptr;
    
    return rsd_dd_Pk(tree, P13, P22, Z2_delta, Z0_v, Z2_v, Z0_vdelta, Z2_vdelta, Z2_vv, Z2_vvdelta, Z2_vvv);
  }


rsd_dd_Pk oneloop_Pk_calculator::compute_rsd_dd_mu6(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                                    const loop_integral& loop_data, const tree_Pk& Ptree)
  {
    delta_22_integrals d22 = loop_data.get_delta22();
    delta_13_integrals d13 = loop_data.get_delta13();
    
    rsd_22_integrals rsd22 = loop_data.get_rsd22();
    rsd_13_integrals rsd13 = loop_data.get_rsd13();
    
    // TREE TERMS
    
    auto Ptr = Ptree(k);
    
    Pk_value tree;  // no tree contribution at mu^6
    
    // 13 TERMS
    
    Pk_value P13_simple;  // no simple contributions at mu^6
    
    auto t_factor = val.f*val.f * val.g*val.g * (val.A*val.fA + val.B*val.fB - val.f*val.g*val.g);
    
    Pk_value P13_tensor_a = Ptr * t_factor * rsd13.get_a() / 2.0;
    Pk_value P13_tensor_d = Ptr * 3.0 * t_factor * rsd13.get_d();
    Pk_value P13_tensor_e = Ptr * (7.0/2.0) * t_factor * rsd13.get_e();
    Pk_value P13_tensor_f = Ptr * (-4.0) * t_factor * rsd13.get_f();
    Pk_value P13_tensor_g = Ptr * (-3.0/2.0) * t_factor * rsd13.get_g();
    
    Pk_value P13_tensor_c = Ptr * val.f*val.f * val.g*val.g * (15.0*val.A*val.fA + 31.0*val.B*val.fB + val.f*(-15.0+8.0*val.f)*val.g*val.g) * rsd13.get_c() / (-6.0);
    
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
    
    k2_Pk_value Z2_vv = 2.0 * val.f*val.f * val.g*val.g * (val.A*val.fA + 6.0*val.B*val.fB - val.f*val.g*val.g) * k*k * Ptr;
    
    k2_Pk_value Z2_vvdelta;   // no contribution at mu^6
    
    k2_Pk_value Z2_vvv = 5.0 * val.f*val.f*val.f*val.f * val.g*val.g*val.g*val.g * k*k * Ptr;
    
    return rsd_dd_Pk(tree, P13, P22, Z2_delta, Z0_v, Z2_v, Z0_vdelta, Z2_vdelta, Z2_vv, Z2_vvdelta, Z2_vvv);
  }


rsd_dd_Pk oneloop_Pk_calculator::compute_rsd_dd_mu8(const Mpc_units::energy& k, const oneloop_growth_record& val,
                                                    const loop_integral& loop_data, const tree_Pk& Ptree)
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
