//
// Created by David Seery on 21/11/2015.
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


#include <iostream>
#include <fstream>

#include "oneloop_momentum_integrator.h"
#include "oneloop_integrands/delta_13.h"
#include "oneloop_integrands/delta_22.h"
#include "oneloop_integrands/delta_rsd_13.h"
#include "oneloop_integrands/delta_rsd_22.h"

#include "cuba.h"

#include "boost/timer/timer.hpp"



oneloop_momentum_integrator::oneloop_momentum_integrator(double a_13, double r_13, double a_22, double r_22)
  : abs_err_13(std::abs(a_13)),
    rel_err_13(std::abs(r_13)),
    abs_err_22(std::abs(a_22)),
    rel_err_22(std::abs(r_22))
  {
    // seed random number generator
    mersenne_twister.seed(random_device());
  }


loop_integral oneloop_momentum_integrator::integrate(const FRW_model& model, const Mpc_units::energy& k,
                                                     const k_token& k_tok, const Mpc_units::energy& UV_cutoff,
                                                     const UV_cutoff_token& UV_tok, const Mpc_units::energy& IR_cutoff,
                                                     const IR_cutoff_token& IR_tok, const initial_filtered_Pk& Pk)
  {
    delta_13_integrals delta13;
    delta_22_integrals delta22;
    rsd_13_integrals rsd13;
    rsd_22_integrals rsd22;

    // delta 22 kernels
    bool failAA = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::AA_integrand, delta22.get_AA(), loop_integral_type::P22, "AA");
    bool failAB = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::AB_integrand, delta22.get_AB(), loop_integral_type::P22, "AB");
    bool failBB = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::BB_integrand, delta22.get_BB(), loop_integral_type::P22, "BB");
    
    if(failAA || failAB || failBB) delta22.mark_failed();

    // delta 13 kernels
    bool failD  = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::D_integrand, delta13.get_D(), loop_integral_type::P13, "D");
    bool failE  = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::E_integrand, delta13.get_E(), loop_integral_type::P13, "E");
    bool failF  = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::F_integrand, delta13.get_F(), loop_integral_type::P13, "F");
    bool failG  = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::G_integrand, delta13.get_G(), loop_integral_type::P13, "G");
    bool failJ1 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::J1_integrand, delta13.get_J1(), loop_integral_type::P13, "J1");
    bool failJ2 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::J2_integrand, delta13.get_J2(), loop_integral_type::P13, "J2");
    
    if(failD || failE || failF || failG || failJ1 || failJ2) delta13.mark_failed();
    
    // RSD 13 kernels
    bool fail_RSD13_a = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD13_a_integrand, rsd13.get_a(), loop_integral_type::P13, "a");
    bool fail_RSD13_b = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD13_b_integrand, rsd13.get_b(), loop_integral_type::P13, "b");
    bool fail_RSD13_c = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD13_c_integrand, rsd13.get_c(), loop_integral_type::P13, "c");
    bool fail_RSD13_d = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD13_d_integrand, rsd13.get_d(), loop_integral_type::P13, "d");
    bool fail_RSD13_e = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD13_e_integrand, rsd13.get_e(), loop_integral_type::P13, "e");
    bool fail_RSD13_f = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD13_f_integrand, rsd13.get_f(), loop_integral_type::P13, "f");
    bool fail_RSD13_g = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD13_g_integrand, rsd13.get_g(), loop_integral_type::P13, "g");
    
    if(fail_RSD13_a || fail_RSD13_b || fail_RSD13_c || fail_RSD13_d || fail_RSD13_e || fail_RSD13_f || fail_RSD13_g)
      rsd13.mark_failed();
    
    // RSD 22 kernels
    bool fail_RSD22_A1 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD22_A1_integrand, rsd22.get_A1(), loop_integral_type::P22, "A1");
    bool fail_RSD22_A2 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD22_A2_integrand, rsd22.get_A2(), loop_integral_type::P22, "A2");
    bool fail_RSD22_A3 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD22_A3_integrand, rsd22.get_A3(), loop_integral_type::P22, "A3");
    bool fail_RSD22_A4 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD22_A4_integrand, rsd22.get_A4(), loop_integral_type::P22, "A4");
    bool fail_RSD22_A5 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD22_A5_integrand, rsd22.get_A5(), loop_integral_type::P22, "A5");
    bool fail_RSD22_B2 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD22_B2_integrand, rsd22.get_B2(), loop_integral_type::P22, "B2");
    bool fail_RSD22_B3 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD22_B3_integrand, rsd22.get_B3(), loop_integral_type::P22, "B3");
    bool fail_RSD22_B6 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD22_B6_integrand, rsd22.get_B6(), loop_integral_type::P22, "B6");
    bool fail_RSD22_B8 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD22_B8_integrand, rsd22.get_B8(), loop_integral_type::P22, "B8");
    bool fail_RSD22_B9 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD22_B9_integrand, rsd22.get_B9(), loop_integral_type::P22, "B9");
    bool fail_RSD22_C1 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD22_C1_integrand, rsd22.get_C1(), loop_integral_type::P22, "C1");
    bool fail_RSD22_C2 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD22_C2_integrand, rsd22.get_C2(), loop_integral_type::P22, "C2");
    bool fail_RSD22_C4 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD22_C4_integrand, rsd22.get_C4(), loop_integral_type::P22, "C4");
    bool fail_RSD22_D1 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::RSD22_D1_integrand, rsd22.get_D1(), loop_integral_type::P22, "D1");
    
    if(fail_RSD22_A1 || fail_RSD22_A2 || fail_RSD22_A3 || fail_RSD22_A4 || fail_RSD22_A5 || fail_RSD22_B2 ||
       fail_RSD22_B3 || fail_RSD22_B6 || fail_RSD22_B8 || fail_RSD22_B9 || fail_RSD22_C1 || fail_RSD22_C2 ||
       fail_RSD22_C4 || fail_RSD22_D1)
      rsd22.mark_failed();
    
    loop_integral container(k_tok, Pk.get_token(), UV_tok, IR_tok, delta22, delta13, rsd22, rsd13);

    return container;
  }


template <typename KernelRecord>
bool oneloop_momentum_integrator::kernel_integral(const FRW_model& model, const Mpc_units::energy& k,
                                                  const Mpc_units::energy& UV_cutoff,
                                                  const Mpc_units::energy& IR_cutoff, const initial_filtered_Pk& Pk,
                                                  integrand_t integrand, KernelRecord& result, loop_integral_type type,
                                                  const std::string& name)
  {
    // disable CUBA's internal auto-parallelization
    // we're handling multiprocessor activity ourselves via the scheduler,
    // so it's preferable to keep each core fully active rather than have threads
    // trying to manage Cuba's subworkers
    cubacores(0, oneloop_momentum_impl::pcores);
    
    wiggle_Pk_raw_adapter raw(Pk);
    wiggle_Pk_nowiggle_adapter nw(Pk);
    
    bool fail_raw = this->evaluate_integral(model, k, UV_cutoff, IR_cutoff, raw, integrand, result.get_raw(), type, name, "raw");
    bool fail_nw  = this->evaluate_integral(model, k, UV_cutoff, IR_cutoff, nw, integrand, result.get_nowiggle(), type, name, "no-wiggle");
    
    return fail_raw || fail_nw;
  }


template <typename IntegralRecord>
bool oneloop_momentum_integrator::evaluate_integral(const FRW_model& model, const Mpc_units::energy& k,
                                                    const Mpc_units::energy& UV_cutoff,
                                                    const Mpc_units::energy& IR_cutoff, const spline_Pk& Pk,
                                                    integrand_t integrand, IntegralRecord& result,
                                                    loop_integral_type type, const std::string& name,
                                                    const std::string& component)
  {
    cubareal integral[oneloop_momentum_impl::dimensions];
    cubareal error[oneloop_momentum_impl::dimensions];
    cubareal prob[oneloop_momentum_impl::dimensions];
    
    int regions;
    int evaluations;
    int fail;

    boost::timer::cpu_timer raw_timer;
    
    std::unique_ptr<oneloop_momentum_impl::integrand_data> data =
      std::make_unique<oneloop_momentum_impl::integrand_data>(model, k, UV_cutoff, IR_cutoff, Pk);
    
    double re = (type == loop_integral_type::P13 ? this->rel_err_13 : this->rel_err_22);
    double ae = (type == loop_integral_type::P13 ? this->abs_err_13 : this->abs_err_22);

    constexpr unsigned int MAX_13_TRIES = 5;
    constexpr unsigned int MAX_22_TRIES = 3;
    unsigned int max_tries = (type == loop_integral_type::P13 ? MAX_13_TRIES : MAX_22_TRIES);
    unsigned int tries = 0;
    
    while(tries < max_tries)
      {
        if(tries > 0)
          {
            re = re * 4.0;
            std::cout << "lsseft: relaxing error tolerance for kernel = " << name << " (" << component << "), attempt "
                      << tries << ", now abstol = " << ae << ", reltol = " << re << '\n';
          }

        Cuhre(oneloop_momentum_impl::dimensions,
              oneloop_momentum_impl::components,
              integrand, data.get(),
              oneloop_momentum_impl::points_per_invocation,
              re, ae,
              oneloop_momentum_impl::verbosity_none | oneloop_momentum_impl::samples_last,
              oneloop_momentum_impl::min_eval, oneloop_momentum_impl::max_eval,
              oneloop_momentum_impl::cuhre_key,
              nullptr, nullptr,
              &regions, &evaluations, &fail,
              integral, error, prob);
        
        if(fail == 0)
          {
            raw_timer.stop();
    
            // an overall factor 1/8pi^2 is taken out of the integrand, so remember to put it back here
            result.value = typename IntegralRecord::value_type(integral[0] / (8.0 * M_PI * M_PI));
            result.regions = regions;
            result.evaluations = evaluations;
            result.error = error[0];
            result.time = raw_timer.elapsed().wall;
            
            break;
          }
        
        ++tries;
      }
    
    if(tries >= max_tries)
      std::cerr << "lsseft: integration failure: kernel = " << name << " (" << component << "), "
                << "regions = " << regions << ", evaluations = " << evaluations << ", fail = "
                << fail << ", value = " << integral[0] << ", error = " << error[0] << ", probability = " << prob[0]
                << '\n';
    
    
    return (tries >= max_tries);
  }

// Alternative Divonne integrator

//    Divonne(oneloop_momentum_impl::dimensions, oneloop_momentum_impl::components,
//            integrand, data.get(),
//            oneloop_momentum_impl::points_per_invocation,
//            this->rel_err, this->abs_err,
//            oneloop_momentum_impl::verbosity_none | oneloop_momentum_impl::samples_last,
//            this->mersenne_twister(),                                                          // seed for internal Cuba random number generator
//            oneloop_momentum_impl::min_eval, oneloop_momentum_impl::max_eval,
//            oneloop_momentum_impl::divonne_key1, oneloop_momentum_impl::divonne_key2, oneloop_momentum_impl::divonne_key3,
//            oneloop_momentum_impl::divonne_maxpass,
//            oneloop_momentum_impl::divonne_border, oneloop_momentum_impl::divonne_maxchisq, oneloop_momentum_impl::divonne_minchisq,
//            oneloop_momentum_impl::ngiven, oneloop_momentum_impl::ldxgiven, nullptr,
//            oneloop_momentum_impl::nextra, nullptr,
//            nullptr, nullptr,
//            &regions, &evaluations, &fail,
//            integral, error, prob);


void oneloop_momentum_integrator::write_integrands(const FRW_model& model, const Mpc_units::energy& k,
                                                   const Mpc_units::energy& UV_cutoff,
                                                   const Mpc_units::energy& IR_cutoff,
                                                   const spline_Pk& Pk, unsigned int Npoints)
  {
    std::ofstream AA;
    std::ofstream AB;
    std::ofstream BB;
    std::ofstream D;
    std::ofstream E;
    std::ofstream F;
    std::ofstream G;
    std::ofstream J1;
    std::ofstream J2;

    AA.open("AA.csv", std::ofstream::out | std::ofstream::trunc);
    AB.open("AB.csv", std::ofstream::out | std::ofstream::trunc);
    BB.open("BB.csv", std::ofstream::out | std::ofstream::trunc);
    D.open("D.csv", std::ofstream::out | std::ofstream::trunc);
    E.open("E.csv", std::ofstream::out | std::ofstream::trunc);
    F.open("F.csv", std::ofstream::out | std::ofstream::trunc);
    G.open("G.csv", std::ofstream::out | std::ofstream::trunc);
    J1.open("J1.csv", std::ofstream::out | std::ofstream::trunc);
    J2.open("J2.csv", std::ofstream::out | std::ofstream::trunc);

    AA.precision(12);
    AB.precision(12);
    BB.precision(12);
    D.precision(12);
    E.precision(12);
    F.precision(12);
    G.precision(12);
    J1.precision(12);
    J2.precision(12);

    std::shared_ptr<oneloop_momentum_impl::integrand_data> data = std::make_shared<oneloop_momentum_impl::integrand_data>(model, k, UV_cutoff, IR_cutoff, Pk);

    for(unsigned int l = 0; l <= Npoints; ++l)
      {
        for(unsigned int m = 0; m <= Npoints; ++m)
          {
            cubareal x[oneloop_momentum_impl::dimensions];
            cubareal f[oneloop_momentum_impl::components];

            x[0] = static_cast<cubareal>(l) / static_cast<cubareal>(Npoints);
            x[1] = static_cast<cubareal>(m) / static_cast<cubareal>(Npoints);

            f[0] = -1000.0;
            oneloop_momentum_impl::AA_integrand(nullptr, x, nullptr, f, data.get());
            AA << l << "," << m << "," << x[0] << "," << x[1] << "," << f[0] << '\n';

            f[0] = -1000.0;
            oneloop_momentum_impl::AB_integrand(nullptr, x, nullptr, f, data.get());
            AB << l << "," << m << "," << x[0] << "," << x[1] << "," << f[0] << '\n';

            f[0] = -1000.0;
            oneloop_momentum_impl::BB_integrand(nullptr, x, nullptr, f, data.get());
            BB << l << "," << m << "," << x[0] << "," << x[1] << "," << f[0] << '\n';

            f[0] = -1000.0;
            oneloop_momentum_impl::D_integrand(nullptr, x, nullptr, f, data.get());
            D << l << "," << m << "," << x[0] << "," << x[1] << "," << f[0] << '\n';

            f[0] = -1000.0;
            oneloop_momentum_impl::E_integrand(nullptr, x, nullptr, f, data.get());
            E << l << "," << m << "," << x[0] << "," << x[1] << "," << f[0] << '\n';

            f[0] = -1000.0;
            oneloop_momentum_impl::F_integrand(nullptr, x, nullptr, f, data.get());
            F << l << "," << m << "," << x[0] << "," << x[1] << "," << f[0] << '\n';

            f[0] = -1000.0;
            oneloop_momentum_impl::G_integrand(nullptr, x, nullptr, f, data.get());
            G << l << "," << m << "," << x[0] << "," << x[1] << "," << f[0] << '\n';

            f[0] = -1000.0;
            oneloop_momentum_impl::J1_integrand(nullptr, x, nullptr, f, data.get());
            J1 << l << "," << m << "," << x[0] << "," << x[1] << "," << f[0] << '\n';
            
            f[0] = -1000.0;
            oneloop_momentum_impl::J2_integrand(nullptr, x, nullptr, f, data.get());
            J2 << l << "," << m << "," << x[0] << "," << x[1] << "," << f[0] << '\n';
          }
      }

    AA.close();
    AB.close();
    BB.close();
    D.close();
    E.close();
    F.close();
    G.close();
    J1.close();
    J2.close();
  }
