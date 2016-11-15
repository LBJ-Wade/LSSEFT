//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
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
                                                     const UV_token& UV_tok, const Mpc_units::energy& IR_cutoff,
                                                     const IR_token& IR_tok, const tree_power_spectrum& Pk)
  {
    delta_13_integrals delta13;
    delta_22_integrals delta22;
    rsd_13_integrals rsd13;
    rsd_22_integrals rsd22;

    // delta 22 kernels
    bool failAA = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::AA_integrand, delta22.get_AA(), loop_integral_type::P22);
    bool failAB = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::AB_integrand, delta22.get_AB(), loop_integral_type::P22);
    bool failBB = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::BB_integrand, delta22.get_BB(), loop_integral_type::P22);
    
    if(failAA || failAB || failBB) delta22.mark_failed();

    // delta 13 kernels
    bool failD  = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::D_integrand, delta13.get_D(), loop_integral_type::P13);
    bool failE  = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::E_integrand, delta13.get_E(), loop_integral_type::P13);
    bool failF  = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::F_integrand, delta13.get_F(), loop_integral_type::P13);
    bool failG  = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::G_integrand, delta13.get_G(), loop_integral_type::P13);
    bool failJ1 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::J1_integrand, delta13.get_J1(), loop_integral_type::P13);
    bool failJ2 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::J2_integrand, delta13.get_J2(), loop_integral_type::P13);
    
    if(failD || failE || failF || failG || failJ1 || failJ2) delta13.mark_failed();
    
    // rsd 13 kernels
    
    loop_integral container(k_tok, UV_tok, IR_tok, delta22, delta13, rsd22, rsd13);

    return container;
  }


template <typename KernelRecord>
bool oneloop_momentum_integrator::kernel_integral(const FRW_model& model, const Mpc_units::energy& k, const Mpc_units::energy& UV_cutoff,
                                                  const Mpc_units::energy& IR_cutoff, const tree_power_spectrum& Pk, integrand_t integrand,
                                                  KernelRecord& result, loop_integral_type type)
  {
    boost::timer::cpu_timer timer;

    cubareal integral[oneloop_momentum_impl::dimensions];
    cubareal error[oneloop_momentum_impl::dimensions];
    cubareal prob[oneloop_momentum_impl::dimensions];

    int regions;
    int evaluations;
    int fail;

    std::unique_ptr<oneloop_momentum_impl::integrand_data> data = std::make_unique<oneloop_momentum_impl::integrand_data>(model, k, UV_cutoff, IR_cutoff, Pk);

    // disable CUBA's internal auto-parallelization
    // we're handling multiprocessor activity ourselves via the scheduler,
    // so it's preferable to keep each core fully active rather than have threads
    // trying to manage Cuba's subworkers
    cubacores(0, oneloop_momentum_impl::pcores);

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

    Cuhre(oneloop_momentum_impl::dimensions, oneloop_momentum_impl::components,
          integrand, data.get(),
          oneloop_momentum_impl::points_per_invocation,
          (type == loop_integral_type::P13 ? this->rel_err_13 : this->rel_err_22),
          (type == loop_integral_type::P13 ? this->abs_err_13 : this->abs_err_22),
          oneloop_momentum_impl::verbosity_none | oneloop_momentum_impl::samples_last,
          oneloop_momentum_impl::min_eval, oneloop_momentum_impl::max_eval,
          oneloop_momentum_impl::cuhre_key,
          nullptr, nullptr,
          &regions, &evaluations, &fail,
          integral, error, prob);

    timer.stop();

    if(fail != 0) std::cerr << "Integration failure: regions = " << regions << ", evaluations = " << evaluations << ", fail = " << fail << ", value = " << integral[0] << ", error = " << error[0] << ", probability = " << prob[0] << '\n';
//    else          std::cerr << "Integration success: regions = " << regions << ", evaluations = " << evaluations << ", fail = " << fail << ", value = " << integral[0] << ", error = " << error[0] << ", probability = " << prob[0] << '\n';

    // an overall factor 1 / (2pi)^3 is taken out of the integrand, so remember to put it back here
    result.value       = typename KernelRecord::value_type(integral[0] / (8.0 * M_PI * M_PI * M_PI));
    result.regions     = regions;
    result.evaluations = evaluations;
    result.error       = error[0];
    result.time        = timer.elapsed().wall;

    return(fail != 0);
  }


void oneloop_momentum_integrator::write_integrands(const FRW_model& model, const Mpc_units::energy& k,
                                                   const Mpc_units::energy& UV_cutoff, const Mpc_units::energy& IR_cutoff,
                                                   const tree_power_spectrum& Pk, unsigned int Npoints)
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
