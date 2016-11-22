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
  : fail(false)
  {
  }


rsd_13_integrals::rsd_13_integrals()
  : fail(false)
  {
  }

  
loop_integral::loop_integral(const k_token& kt, const UV_cutoff_token& UVt, const IR_cutoff_token& IRt,
                             const delta_22_integrals& d22, const delta_13_integrals& d13,
                             const rsd_22_integrals& r22, const rsd_13_integrals& r13)
  : k(kt),
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
    UV_cutoff(0),
    IR_cutoff(0),
    delta22(),
    delta13(),
    rsd22(),
    rsd13()
  {
  }