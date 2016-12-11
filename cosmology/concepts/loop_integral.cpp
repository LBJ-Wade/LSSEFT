//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "loop_integral.h"


delta_22_integrals::delta_22_integrals(const inverse_energy3_integral& _AA, const inverse_energy3_integral& _AB,
                                       const inverse_energy3_integral& _BB)
  : fail(false),
    AA(_AA),
    AB(_AB),
    BB(_BB)
  {
  }


delta_22_integrals::delta_22_integrals()
  : fail(false),
    AA(),
    AB(),
    BB()
  {
  }


delta_13_integrals::delta_13_integrals(const dimless_integral& _D, const dimless_integral& _E, const dimless_integral& _F, const dimless_integral& _G,
                                       const dimless_integral& _J1, const dimless_integral& _J2)
  : fail(false),
    D(_D),
    E(_E),
    F(_F),
    G(_G),
    J1(_J1),
    J2(_J2)
  {
  }


delta_13_integrals::delta_13_integrals()
  : fail(false),
    D(),
    E(),
    F(),
    G(),
    J1(),
    J2()
  {
  }


rsd_22_integrals::rsd_22_integrals()
  : fail(false),
    A1(),
    A2(),
    A3(),
    A4(),
    A5(),
    B2(),
    B3(),
    B6(),
    B8(),
    B9(),
    C1(),
    C2(),
    C4(),
    D1()

{
  }


rsd_22_integrals::rsd_22_integrals(const inverse_energy3_integral& _A1, const inverse_energy3_integral& _A2,
                                   const inverse_energy3_integral& _A3, const inverse_energy3_integral& _A4,
                                   const inverse_energy3_integral& _A5, const inverse_energy3_integral& _B2,
                                   const inverse_energy3_integral& _B3, const inverse_energy3_integral& _B6,
                                   const inverse_energy3_integral& _B8, const inverse_energy3_integral& _B9,
                                   const inverse_energy3_integral& _C1, const inverse_energy3_integral& _C2,
                                   const inverse_energy3_integral& _C4, const inverse_energy3_integral& _D1)
  : fail(false),
    A1(_A1),
    A2(_A2),
    A3(_A3),
    A4(_A4),
    A5(_A5),
    B2(_B2),
    B3(_B3),
    B6(_B6),
    B8(_B8),
    B9(_B9),
    C1(_C1),
    C2(_C2),
    C4(_C4),
    D1(_D1)
  {
  }


rsd_13_integrals::rsd_13_integrals()
  : fail(false),
    a(),
    b(),
    c(),
    d(),
    e(),
    f(),
    g()
  {
  }


rsd_13_integrals::rsd_13_integrals(const dimless_integral& _a, const dimless_integral& _b, const dimless_integral& _c,
                                   const dimless_integral& _d, const dimless_integral& _e, const dimless_integral& _f,
                                   const dimless_integral& _g)
  : fail(false),
    a(_a),
    b(_b),
    c(_c),
    d(_d),
    e(_e),
    f(_f),
    g(_g)
  {
  }


loop_integral::loop_integral(const k_token& kt, const linear_Pk_token& Pkt, const UV_cutoff_token& UVt, const IR_cutoff_token& IRt,
                             const delta_22_integrals& d22, const delta_13_integrals& d13,
                             const rsd_22_integrals& r22, const rsd_13_integrals& r13)
  : k(kt),
    Pk_lin(Pkt),
    UV_cutoff(UVt),
    IR_cutoff(IRt),
    delta22(d22),
    delta13(d13),
    rsd22(r22),
    rsd13(r13)
  {
  }


loop_integral::loop_integral()
  : k(0),
    Pk_lin(0),
    UV_cutoff(0),
    IR_cutoff(0),
    delta22(),
    delta13(),
    rsd22(),
    rsd13()
  {
  }

