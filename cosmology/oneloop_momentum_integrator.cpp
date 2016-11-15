//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include <iostream>
#include <fstream>

#include "oneloop_momentum_integrator.h"

#include "cuba.h"

#include "boost/timer/timer.hpp"


namespace oneloop_momentum_impl
  {

    constexpr unsigned int dimensions            = 2;   // no point doing integrals over phi, because the integrands don't depend on it
    constexpr unsigned int components            = 1;
    constexpr unsigned int points_per_invocation = 1;

    constexpr unsigned int verbosity_none        = 0;
    constexpr unsigned int verbosity_reasonable  = 1;
    constexpr unsigned int verbosity_progress    = 2;
    constexpr unsigned int verbosity_subregions  = 3;

    constexpr unsigned int samples_all           = 0;
    constexpr unsigned int samples_last          = 4;

    constexpr unsigned int min_eval              = 0;
    constexpr unsigned int max_eval              = 20000000;

    constexpr unsigned int ngiven                = 0;
    constexpr unsigned int ldxgiven              = 0;
    constexpr unsigned int nextra                = 0;

    constexpr unsigned int pcores                = 10000;   // matches default Cuba value

    constexpr unsigned int cuhre_key             = 13;      // degree-13 only available in 2-dimensions
    constexpr unsigned int divonne_key1          = 47;
    constexpr unsigned int divonne_key2          = 13;      // degree-13 only available in 2-dimensions
    constexpr unsigned int divonne_key3          = 1;
    constexpr unsigned int divonne_maxpass       = 5;
    constexpr unsigned int divonne_border        = 0;
    constexpr double       divonne_maxchisq      = 10.0;
    constexpr double       divonne_minchisq      = 0.25;


    class integrand_data
      {

      public:

        integrand_data(const FRW_model& m, const Mpc_units::energy& _k, const Mpc_units::energy& UV, const Mpc_units::energy& IR,
                       const tree_power_spectrum& _Pk)
          : model(m),
            k(_k),
            UV_cutoff(UV),
            IR_cutoff(IR),
            Pk(_Pk),
            jacobian(2.0 * 2.0*M_PI * (UV_cutoff-IR_cutoff)), // Jacobian in angular directions is 2 * 2 * pi = 4pi;
                                                              // the integral over phi isn't done (the integrand doesn't depend on it), but this accounts for its contribution;
                                                              // the theta integral has been switched to d(cos theta) which runs from -1 to +2 and contributes a factor 2
            q_range(UV_cutoff - IR_cutoff),
            k_sq(k*k)
          {
          }

        const FRW_model& model;
        const Mpc_units::energy& k;
        const Mpc_units::energy& UV_cutoff;
        const Mpc_units::energy& IR_cutoff;
        const tree_power_spectrum& Pk;

        Mpc_units::energy  jacobian;
        Mpc_units::energy  q_range;
        Mpc_units::energy2 k_sq;
      };


    static int AA_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;

        Mpc_units::energy2 k_dot_q       = z * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;

        Mpc_units::energy  k_minus_q     = std::sqrt(k_minus_q_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;

        // integral is P(q) P(k-q) alpha(q,k-q)^2
        Mpc_units::inverse_energy3 Pq         = data->Pk(q);
        Mpc_units::inverse_energy3 Pk_minus_q = k_minus_q > data->IR_cutoff && k_minus_q < data->UV_cutoff ? data->Pk(k_minus_q) : Mpc_units::inverse_energy3(0);

        Mpc_units::inverse_energy  qqPq       = q*q * Pq;
        Mpc_units::inverse_energy4 PP_prod    = qqPq * Pk_minus_q;

        double alpha1 = (k_minus_q_sq*k_dot_q + q*q*data->k_sq - q*q*k_dot_q) / (2.0 * q*q * k_minus_q_sq);

        f[0] = (2.0 * data->jacobian * PP_prod * alpha1*alpha1) / Mpc_units::Mpc3;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }


    static int AB_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;

        Mpc_units::energy2 k_dot_q       = z * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;

        Mpc_units::energy  k_minus_q     = std::sqrt(k_minus_q_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;

        // integral is P(q) P(k-q) alpha(q,k-q) gamma(q,k-q)
        Mpc_units::inverse_energy3 Pq         = data->Pk(q);
        Mpc_units::inverse_energy3 Pk_minus_q = k_minus_q > data->IR_cutoff && k_minus_q < data->UV_cutoff ? data->Pk(k_minus_q) : Mpc_units::inverse_energy3(0);

        // integral is P(q) P(k-q) alpha(q,k-q) gamma(q,k-q)
        Mpc_units::inverse_energy  qqPq       = q*q * Pq;
        Mpc_units::inverse_energy4 PP_prod    = qqPq * Pk_minus_q;

        double alpha1 = (k_minus_q_sq*k_dot_q + q*q*data->k_sq - q*q*k_dot_q) / (2.0 * q*q * k_minus_q_sq);
        double gamma1 = (k_minus_q_sq*k_dot_q - q*q*k_dot_q + data->k_sq*k_dot_q) / (2.0 * q*q * k_minus_q_sq);

        f[0] = (4.0 * data->jacobian * PP_prod * alpha1*gamma1) / Mpc_units::Mpc3;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }


    static int BB_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;

        Mpc_units::energy2 k_dot_q       = z * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;

        Mpc_units::energy  k_minus_q     = std::sqrt(k_minus_q_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;

        // integral is P(q) P(k-q) gamma(q,k-q)^2
        Mpc_units::inverse_energy3 Pq         = data->Pk(q);
        Mpc_units::inverse_energy3 Pk_minus_q = k_minus_q > data->IR_cutoff && k_minus_q < data->UV_cutoff ? data->Pk(k_minus_q) : Mpc_units::inverse_energy3(0);

        // integral is P(q) P(k-q) gamma(q,k-q)^2
        Mpc_units::inverse_energy  qqPq       = q*q * Pq;
        Mpc_units::inverse_energy4 PP_prod    = qqPq * Pk_minus_q;

        double gamma1 = (k_minus_q_sq*k_dot_q - q*q*k_dot_q + data->k_sq*k_dot_q) / (2.0 * q*q * k_minus_q_sq);

        f[0] = (2.0 * data->jacobian * PP_prod * gamma1*gamma1) / Mpc_units::Mpc3;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }


    static int D_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;

        Mpc_units::energy2 k_dot_q       = z * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;

        Mpc_units::inverse_energy3 Pq   = data->Pk(q);
        Mpc_units::inverse_energy  qqPq = q*q * Pq;

        // integral is P(q) gamma(k-r,r) alpha(k,-r)
        double gamma1 = (k_minus_q_sq*k_dot_q - q*q*k_dot_q + data->k_sq*k_dot_q) / (2.0 * q*q * k_minus_q_sq);
        double alpha2 = (2.0*data->k_sq*q*q - k_dot_q*(data->k_sq + q*q)) / (2.0 * q*q * data->k_sq);

        f[0] = 8.0 * data->jacobian * qqPq * gamma1*alpha2;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }


    static int E_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;

        Mpc_units::energy2 k_dot_q       = z * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;

        Mpc_units::inverse_energy3 Pq   = data->Pk(q);
        Mpc_units::inverse_energy  qqPq = q*q * Pq;

        // integral is P(q) gamma(k-r,r) gamma(k,-r)
        double gamma1 = (k_minus_q_sq*k_dot_q - q*q*k_dot_q + data->k_sq*k_dot_q) / (2.0 * q*q * k_minus_q_sq);
        double gamma2 = (2.0*data->k_sq*q*q - k_dot_q*(2.0*data->k_sq + 2.0*q*q - 2.0*k_dot_q)) / (2.0 * q*q * data->k_sq);

        f[0] = 8.0 * data->jacobian * qqPq * gamma1*gamma2;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }


    static int F_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;

        Mpc_units::energy2 k_dot_q       = z * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;

        Mpc_units::inverse_energy3 Pq   = data->Pk(q);
        Mpc_units::inverse_energy  qqPq = q*q * Pq;

        // integral is P(q) alpha(k-r,r) alpha(k,-r)
        double alpha1 = (k_minus_q_sq*k_dot_q + q*q*data->k_sq - q*q*k_dot_q) / (2.0 * q*q * k_minus_q_sq);
        double alpha2 = (2.0*data->k_sq*q*q - k_dot_q*(data->k_sq + q*q)) / (2.0 * q*q * data->k_sq);

        f[0] = 8.0 * data->jacobian * qqPq * alpha1*alpha2;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }


    static int G_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;

        Mpc_units::energy2 k_dot_q       = z * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;

        Mpc_units::inverse_energy3 Pq   = data->Pk(q);
        Mpc_units::inverse_energy  qqPq = q*q * Pq;

        // integral is P(q) alpha(k-r,r) gamma(k,-r)
        double alpha1 = (k_minus_q_sq*k_dot_q + q*q*data->k_sq - q*q*k_dot_q) / (2.0 * q*q * k_minus_q_sq);
        double gamma2 = (2.0*data->k_sq*q*q - k_dot_q*(2.0*data->k_sq + 2.0*q*q - 2.0*k_dot_q)) / (2.0 * q*q * data->k_sq);

        f[0] = 8.0 * data->jacobian * qqPq * alpha1*gamma2;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }


    static int J1_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;

        Mpc_units::energy2 k_dot_q       = z * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;

        Mpc_units::inverse_energy3 Pq   = data->Pk(q);
        Mpc_units::inverse_energy  qqPq = q*q * Pq;

        // integral is P(q) alpha_asym(k-r,r) alpha(k,-r)
        double alpha_asym = (data->k_sq - k_dot_q) / k_minus_q_sq;
        double alpha2 = (2.0*data->k_sq*q*q - k_dot_q*(data->k_sq + q*q)) / (2.0 * q*q * data->k_sq);
        
        f[0] = 8.0 * data->jacobian * qqPq * alpha_asym * alpha2;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }
    
    
    static int J2_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;
        
        Mpc_units::energy2 k_dot_q       = z * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;
        
        Mpc_units::inverse_energy3 Pq   = data->Pk(q);
        Mpc_units::inverse_energy  qqPq = q*q * Pq;
        
        // integral is P(q) alpha_asym(k-r,r) gamma(k,-r)
        double alpha_asym = (data->k_sq - k_dot_q) / k_minus_q_sq;
        double gamma2 = (2.0*data->k_sq*q*q - k_dot_q*(2.0*data->k_sq + 2.0*q*q - 2.0*k_dot_q)) / (2.0 * q*q * data->k_sq);
        
        f[0] = 8.0 * data->jacobian * qqPq * alpha_asym * gamma2;
        
        return(0);  // return value irrelevant unless = -999, which means stop integration
      }

  }


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

    bool failAA = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::AA_integrand, delta22.get_AA(), loop_integral_type::P22);
    bool failAB = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::AB_integrand, delta22.get_AB(), loop_integral_type::P22);
    bool failBB = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::BB_integrand, delta22.get_BB(), loop_integral_type::P22);
    
    if(failAA || failAB || failBB) delta22.mark_failed();

    bool failD  = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::D_integrand, delta13.get_D(), loop_integral_type::P13);
    bool failE  = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::E_integrand, delta13.get_E(), loop_integral_type::P13);
    bool failF  = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::F_integrand, delta13.get_F(), loop_integral_type::P13);
    bool failG  = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::G_integrand, delta13.get_G(), loop_integral_type::P13);
    bool failJ1 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::J1_integrand, delta13.get_J1(), loop_integral_type::P13);
    bool failJ2 = this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::J2_integrand, delta13.get_J2(), loop_integral_type::P13);
    
    if(failD || failE || failF || failG || failJ1 || failJ2) delta13.mark_failed();
    
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
