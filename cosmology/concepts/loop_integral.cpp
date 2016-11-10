//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "loop_integral.h"


delta_22_integrals::delta_22_integrals(const inverse_energy3_kernel& _AA, const inverse_energy3_kernel& _AB,
                                       const inverse_energy3_kernel& _BB)
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


delta_13_integrals::delta_13_integrals(const dimless_kernel& _D, const dimless_kernel& _E, const dimless_kernel& _F, const dimless_kernel& _G,
                                       const dimless_kernel& _J1, const dimless_kernel& _J2)
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

  
loop_integral::loop_integral(const Mpc_units::energy& _k, const k_token& kt,
                             const Mpc_units::energy& UV,  const UV_token& UVt,
                             const Mpc_units::energy& IR, const IR_token& IRt,
                             const delta_22_integrals& d22, const delta_13_integrals& d13,
                             const rsd_22_integrals& r22, const rsd_13_integrals& r13)
  : k(_k),
    k_tok(kt),
    UV_cutoff(UV),
    UV_tok(UVt),
    IR_cutoff(IR),
    IR_tok(IRt),
    delta22(d22),
    delta13(d13),
    rsd22(r22),
    rsd13(r13)
  {
  }


loop_integral::loop_integral()
  : k(Mpc_units::energy(0)),
    k_tok(0),
    UV_cutoff(Mpc_units::energy(0)),
    UV_tok(0),
    IR_cutoff(Mpc_units::energy(0)),
    IR_tok(0),
    delta22(),
    delta13(),
    rsd22(),
    rsd13()
  {
  }