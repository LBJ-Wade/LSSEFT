//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include "oneloop_momentum_integrator.h"

#include "cuba.h"

#include "boost/timer/timer.hpp"


namespace oneloop_momentum_impl
  {

    constexpr unsigned int dimensions            = 3;
    constexpr unsigned int components            = 1;
    constexpr unsigned int points_per_invocation = 1;

    constexpr unsigned int verbosity_none        = 0;
    constexpr unsigned int verbosity_reasonable  = 1;
    constexpr unsigned int verbosity_progress    = 2;
    constexpr unsigned int verbosity_subregions  = 3;

    constexpr unsigned int samples_all           = 0;
    constexpr unsigned int samples_last          = 4;

    constexpr unsigned int min_eval              = 0;
    constexpr unsigned int max_eval              = 10000000;

    constexpr unsigned int ngiven                = 0;
    constexpr unsigned int ldxgiven              = 0;
    constexpr unsigned int nextra                = 0;


    class integrand_data
      {

      public:

        integrand_data(const FRW_model& m, const Mpc_units::energy& _k, const Mpc_units::energy& UV, const Mpc_units::energy& IR,
                       const std::shared_ptr<tree_power_spectrum>& _Pk)
          : model(m),
            k(_k),
            UV_cutoff(UV),
            IR_cutoff(IR),
            Pk(_Pk),
            jacobian(2.0*M_PI_2*(UV_cutoff-IR_cutoff)),    // Jacobian in angular directions in 2pi * pi = 2pi^2
            q_range(UV_cutoff - IR_cutoff),
            k_sq(k*k)
          {
          }

        const FRW_model& model;
        const Mpc_units::energy& k;
        const Mpc_units::energy& UV_cutoff;
        const Mpc_units::energy& IR_cutoff;
        const std::shared_ptr<tree_power_spectrum>& Pk;

        Mpc_units::energy  jacobian;
        Mpc_units::energy  q_range;
        Mpc_units::energy2 k_sq;
      };


    static int AA_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double theta        = 2.0 * M_PI * x[1];
        double phi          = M_PI * x[2];

        Mpc_units::energy2 k_dot_q       = std::cos(theta) * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;

        Mpc_units::energy  k_minus_q(std::sqrt(k_minus_q_sq * Mpc_units::Mpc2));

        // integral is P(q) P(k-q) alpha(q,k-q)^2
        Mpc_units::inverse_energy3 Pq         = (*(data->Pk))(q);
        Mpc_units::inverse_energy3 Pk_minus_q = k_minus_q > data->IR_cutoff ? (*(data->Pk))(k_minus_q) : Mpc_units::inverse_energy3(0);

        Mpc_units::inverse_energy  qqPq       = std::sin(theta) * q*q * Pq;
        Mpc_units::inverse_energy4 PP_prod    = qqPq * Pk_minus_q;

        double alpha1 = (k_minus_q_sq*k_dot_q + q*q*data->k_sq - q*q*k_dot_q) / (2.0 * q*q * k_minus_q_sq);

        f[0] = (2.0 * data->jacobian * PP_prod * alpha1*alpha1) / Mpc_units::Mpc3;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }


    static int AB_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double theta        = 2.0 * M_PI * x[1];
        double phi          = M_PI * x[2];

        Mpc_units::energy2 k_dot_q       = std::cos(theta) * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;

        Mpc_units::energy  k_minus_q(std::sqrt(k_minus_q_sq * Mpc_units::Mpc2));

        // integral is P(q) P(k-q) alpha(q,k-q) gamma(q,k-q)
        Mpc_units::inverse_energy3 Pq         = (*(data->Pk))(q);
        Mpc_units::inverse_energy3 Pk_minus_q = k_minus_q > data->IR_cutoff ? (*(data->Pk))(k_minus_q) : Mpc_units::inverse_energy3(0);

        // integral is P(q) P(k-q) alpha(q,k-q) gamma(q,k-q)
        Mpc_units::inverse_energy  qqPq       = std::sin(theta) * q*q * Pq;
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
        double theta        = 2.0 * M_PI * x[1];
        double phi          = M_PI * x[2];

        Mpc_units::energy2 k_dot_q       = std::cos(theta) * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;

        Mpc_units::energy  k_minus_q(std::sqrt(k_minus_q_sq * Mpc_units::Mpc2));

        // integral is P(q) P(k-q) alpha(q,k-q) gamma(q,k-q)
        Mpc_units::inverse_energy3 Pq         = (*(data->Pk))(q);
        Mpc_units::inverse_energy3 Pk_minus_q = k_minus_q > data->IR_cutoff ? (*(data->Pk))(k_minus_q) : Mpc_units::inverse_energy3(0);

        // integral is P(q) P(k-q) gamma(q,k-q)^2
        Mpc_units::inverse_energy  qqPq       = std::sin(theta) * q*q * Pq;
        Mpc_units::inverse_energy4 PP_prod    = qqPq * Pk_minus_q;

        double gamma1 = (k_minus_q_sq*k_dot_q - q*q*k_dot_q + data->k_sq*k_dot_q) / (2.0 * q*q * k_minus_q_sq);

        f[0] = (2.0 * data->jacobian * PP_prod * gamma1*gamma1) / Mpc_units::Mpc3;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }


    static int D_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double theta        = 2.0 * M_PI * x[1];
        double phi          = M_PI * x[2];

        Mpc_units::energy2 k_dot_q       = std::cos(theta) * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;

        Mpc_units::energy  k_minus_q(std::sqrt(k_minus_q_sq * Mpc_units::Mpc2));

        Mpc_units::inverse_energy3 Pq   = (*(data->Pk))(q);
        Mpc_units::inverse_energy  qqPq = std::sin(theta) * q*q * Pq;

        // integral is P(q) gamma(k-r,r) alpha(k,-r)
        double gamma1 = (k_minus_q_sq*k_dot_q - q*q*k_dot_q + data->k_sq*k_dot_q) / (2.0 * q*q * k_minus_q_sq);
        double alpha2 = (2.0*data->k_sq*q*q - k_dot_q*(data->k_sq + q*q)) / (2.0 * q*q * data->k_sq);

        f[0] = 4.0 * data->jacobian * qqPq * gamma1*alpha2;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }


    static int E_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double theta        = 2.0 * M_PI * x[1];
        double phi          = M_PI * x[2];

        Mpc_units::energy2 k_dot_q       = std::cos(theta) * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;

        Mpc_units::energy  k_minus_q(std::sqrt(k_minus_q_sq * Mpc_units::Mpc2));

        Mpc_units::inverse_energy3 Pq   = (*(data->Pk))(q);
        Mpc_units::inverse_energy  qqPq = std::sin(theta) * q*q * Pq;

        // integral is P(q) gamma(k-r,r) gamma(k,-r)
        double gamma1 = (k_minus_q_sq*k_dot_q - q*q*k_dot_q + data->k_sq*k_dot_q) / (2.0 * q*q * k_minus_q_sq);
        double gamma2 = (2.0*data->k_sq*q*q - k_dot_q*(2.0*data->k_sq + 2.0*q*q - 2.0*k_dot_q)) / (2.0 * q*q * data->k_sq);

        f[0] = 4.0 * data->jacobian * qqPq * gamma1*gamma2;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }


    static int F_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double theta        = 2.0 * M_PI * x[1];
        double phi          = M_PI * x[2];

        Mpc_units::energy2 k_dot_q       = std::cos(theta) * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;

        Mpc_units::energy  k_minus_q(std::sqrt(k_minus_q_sq * Mpc_units::Mpc2));

        Mpc_units::inverse_energy3 Pq   = (*(data->Pk))(q);
        Mpc_units::inverse_energy  qqPq = std::sin(theta) * q*q * Pq;

        // integral is P(q) alpha(k-r,r) alpha(k,-r)
        double alpha1 = (k_minus_q_sq*k_dot_q + q*q*data->k_sq - q*q*k_dot_q) / (2.0 * q*q * k_minus_q_sq);
        double alpha2 = (2.0*data->k_sq*q*q - k_dot_q*(data->k_sq + q*q)) / (2.0 * q*q * data->k_sq);

        f[0] = 4.0 * data->jacobian * qqPq * alpha1*alpha2;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }


    static int G_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double theta        = 2.0 * M_PI * x[1];
        double phi          = M_PI * x[2];

        Mpc_units::energy2 k_dot_q       = std::cos(theta) * data->k * q;
        Mpc_units::energy2 k_minus_q_sq  = q*q + data->k_sq - 2.0*k_dot_q;

        Mpc_units::energy  k_minus_q(std::sqrt(k_minus_q_sq * Mpc_units::Mpc2));

        Mpc_units::inverse_energy3 Pq   = (*(data->Pk))(q);
        Mpc_units::inverse_energy  qqPq = std::sin(theta) * q*q * Pq;

        // integral is P(q) alpha(k-r,r) gamma(k,-r)
        double alpha1 = (k_minus_q_sq*k_dot_q + q*q*data->k_sq - q*q*k_dot_q) / (2.0 * q*q * k_minus_q_sq);
        double gamma2 = (2.0*data->k_sq*q*q - k_dot_q*(2.0*data->k_sq + 2.0*q*q - 2.0*k_dot_q)) / (2.0 * q*q * data->k_sq);

        f[0] = 4.0 * data->jacobian * qqPq * alpha1*gamma2;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }

  }


oneloop_momentum_integrator::oneloop_momentum_integrator(double a, double r)
  : abs_err(std::fabs(a)),
    rel_err(std::fabs(r))
  {
    // seed random number generator
    mersenne_twister.seed(random_device());
  }


loop_integral oneloop_momentum_integrator::integrate(const FRW_model& model, const Mpc_units::energy& k,
                                                     const k_token& k_tok, const Mpc_units::energy& UV_cutoff,
                                                     const UV_token& UV_tok, const Mpc_units::energy& IR_cutoff,
                                                     const IR_token& IR_tok, std::shared_ptr<tree_power_spectrum>& Pk)
  {
    boost::timer::cpu_timer timer;

    inverse_energy3_kernel AA;
    inverse_energy3_kernel AB;
    inverse_energy3_kernel BB;
    dimless_kernel         D;
    dimless_kernel         E;
    dimless_kernel         F;
    dimless_kernel         G;

    bool fail = false;

    fail = fail || this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::AA_integrand, AA);
    fail = fail || this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::AB_integrand, AB);
    fail = fail || this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::BB_integrand, BB);
    fail = fail || this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::D_integrand, D);
    fail = fail || this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::E_integrand, E);
    fail = fail || this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::F_integrand, F);
    fail = fail || this->kernel_integral(model, k, UV_cutoff, IR_cutoff, Pk, &oneloop_momentum_impl::G_integrand, G);

    timer.stop();

    loop_integral container(k, k_tok, UV_cutoff, UV_tok, IR_cutoff, IR_tok, fail, AA, AB, BB, D, E, F, G);
    container.set_integration_metadata(timer.elapsed().wall);

    return container;
  }


template <typename KernelRecord>
bool oneloop_momentum_integrator::kernel_integral(const FRW_model& model, const Mpc_units::energy& k,
                                                  const Mpc_units::energy& UV_cutoff, const Mpc_units::energy& IR_cutoff,
                                                  std::shared_ptr<tree_power_spectrum>& Pk, integrand_t integrand, KernelRecord& result)
  {
    cubareal integral[oneloop_momentum_impl::dimensions];
    cubareal error[oneloop_momentum_impl::dimensions];
    cubareal prob[oneloop_momentum_impl::dimensions];

    int regions;
    int evaluations;
    int fail;

    std::unique_ptr<oneloop_momentum_impl::integrand_data> data = std::make_unique<oneloop_momentum_impl::integrand_data>(model, k, UV_cutoff, IR_cutoff, Pk);

    Divonne(oneloop_momentum_impl::dimensions, oneloop_momentum_impl::components,
            integrand, data.get(),
            oneloop_momentum_impl::points_per_invocation,
            this->rel_err, this->abs_err,
            oneloop_momentum_impl::verbosity_none | oneloop_momentum_impl::samples_last,
            this->mersenne_twister(),                                                          // seed for internal Cuba random number generator
            oneloop_momentum_impl::min_eval, oneloop_momentum_impl::max_eval,
            LSSEFT_DEFAULT_DIVONNE_KEY1, LSSEFT_DEFAULT_DIVONNE_KEY2, LSSEFT_DEFAULT_DIVONNE_KEY3,
            LSSEFT_DEFAULT_DIVONNE_MAXPASS,
            LSSEFT_DEFAULT_DIVONNE_BORDER, LSSEFT_DEFAULT_DIVONNE_MAXCHISQ, LSSEFT_DEFAULT_DIVONNE_MINCHISQ,
            oneloop_momentum_impl::ngiven, oneloop_momentum_impl::ldxgiven, nullptr,
            oneloop_momentum_impl::nextra, nullptr,
            nullptr, nullptr,
            &regions, &evaluations, &fail,
            integral, error, prob);

    if(fail != 0)
      {
        std::cerr << "Divonne result: regions = " << regions << ", evaluations = " << evaluations << ", fail = " << fail << ", value = " << integral[0] << ", error = " << error[0] << ", probability = " << prob[0] << '\n';

      }

    // an overall factor 1 / (2pi)^3 is taken out of the integrand, so remember to put it back here
    result.value       = typename KernelRecord::value_type(integral[0] / (16.0 * M_PI_2 * M_PI));
    result.regions     = regions;
    result.evaluations = evaluations;
    result.error       = error[0];

    return(fail != 0);
  }
