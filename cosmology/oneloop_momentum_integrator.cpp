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
            jacobian(2.0*M_PI_2*(UV_cutoff*Mpc_units::Mpc - IR_cutoff*Mpc_units::Mpc)),    // Jacobian in angular directions in 2pi * pi = 2pi^2
            q_min(IR_cutoff*Mpc_units::Mpc),                                               // note guaranteed to be constructed after IR_cutoff, UV_cutoff
            q_max(UV_cutoff*Mpc_units::Mpc),
            q_range(q_max - q_min),
            k_value(k*Mpc_units::Mpc)
          {
          }

        const FRW_model& model;
        const Mpc_units::energy& k;
        const Mpc_units::energy& UV_cutoff;
        const Mpc_units::energy& IR_cutoff;
        const std::shared_ptr<tree_power_spectrum>& Pk;

        double jacobian;
        double q_min;
        double q_max;
        double q_range;

        double k_value;
      };


    static int AA_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        double q     = data->q_min + x[0] * data->q_range;
        double theta = 2.0 * M_PI * x[1];
        double phi   = M_PI * x[2];

        double k_dot_q   = data->k_value*q*std::cos(theta);
        double k_minus_q = std::sqrt(data->k_value*data->k_value + q*q - 2.0*k_dot_q);

        Mpc_units::energy q_in_Mpc(q);
        Mpc_units::energy k_minus_q_in_eV(k_minus_q);

        // integral is P(q) P(k-q) alpha(q,k-q)^2
        Mpc_units::inverse_energy3 Pq         = (*(data->Pk))(q_in_Mpc);
        Mpc_units::inverse_energy3 Pk_minus_q = k_minus_q > data->q_min ? (*(data->Pk))(k_minus_q_in_eV) : Mpc_units::inverse_energy3(0);

        Mpc_units::inverse_energy  qqPq    = std::sin(theta) * q_in_Mpc * q_in_Mpc * Pq;
        Mpc_units::inverse_energy4 PP_prod = qqPq * Pk_minus_q;

        double alpha = (k_minus_q*k_minus_q*k_dot_q + q*q*data->k_value*data->k_value - q*q*k_dot_q) / (2.0 * q*q * k_minus_q*k_minus_q);

        f[0] = data->jacobian * 2.0 * (PP_prod/Mpc_units::Mpc4) * alpha*alpha;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }


    static int AB_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        double q     = data->q_min + x[0] * data->q_range;
        double theta = 2.0 * M_PI * x[1];
        double phi   = M_PI * x[2];

        double k_dot_q   = data->k_value*q*std::cos(theta);
        double k_minus_q = std::sqrt(data->k_value*data->k_value + q*q - 2.0*k_dot_q);

        Mpc_units::energy q_in_Mpc(q);
        Mpc_units::energy k_minus_q_in_eV(k_minus_q);

        // integral is P(q) P(k-q) alpha(q,k-q) gamma(q,k-q)
        Mpc_units::inverse_energy3 Pq         = (*(data->Pk))(q_in_Mpc);
        Mpc_units::inverse_energy3 Pk_minus_q = k_minus_q > data->q_min ? (*(data->Pk))(k_minus_q_in_eV) : Mpc_units::inverse_energy3(0);

        Mpc_units::inverse_energy  qqPq    = std::sin(theta) * q_in_Mpc * q_in_Mpc * Pq;
        Mpc_units::inverse_energy4 PP_prod = qqPq * Pk_minus_q;

        double alpha = (k_minus_q*k_minus_q*k_dot_q + q*q*data->k_value*data->k_value - q*q*k_dot_q) / (2.0 * q*q * k_minus_q*k_minus_q);
        double gamma = (k_minus_q*k_minus_q*k_dot_q - q*q*k_dot_q + data->k_value*data->k_value*k_dot_q) / (2.0 * q*q * k_minus_q*k_minus_q);

        f[0] = data->jacobian * 4.0 * (PP_prod/Mpc_units::Mpc4) * alpha*gamma;

        return(0);
      }


    static int BB_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        double q     = data->q_min + x[0] * data->q_range;
        double theta = 2.0 * M_PI * x[1];
        double phi   = M_PI * x[2];

        double k_dot_q   = data->k_value*q*std::cos(theta);
        double k_minus_q = std::sqrt(data->k_value*data->k_value + q*q - 2.0*k_dot_q);

        Mpc_units::energy q_in_Mpc(q);
        Mpc_units::energy k_minus_q_in_eV(k_minus_q);

        // integral is P(q) P(k-q) gamma(q,k-q)^2
        Mpc_units::inverse_energy3 Pq         = (*(data->Pk))(q_in_Mpc);
        Mpc_units::inverse_energy3 Pk_minus_q = k_minus_q > data->q_min ? (*(data->Pk))(k_minus_q_in_eV) : Mpc_units::inverse_energy3(0);

        Mpc_units::inverse_energy  qqPq    = std::sin(theta) * q_in_Mpc * q_in_Mpc * Pq;
        Mpc_units::inverse_energy4 PP_prod = qqPq * Pk_minus_q;

        double gamma = (k_minus_q*k_minus_q*k_dot_q - q*q*k_dot_q + data->k_value*data->k_value*k_dot_q) / (2.0 * q*q * k_minus_q*k_minus_q);

        f[0] = data->jacobian * 2.0 * (PP_prod/Mpc_units::Mpc4) * gamma*gamma;

        return(0);
      }


    static int D_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        double q     = data->q_min + x[0] * data->q_range;
        double theta = 2.0 * M_PI * x[1];
        double phi   = M_PI * x[2];

        double k_dot_q   = data->k_value*q*std::cos(theta);
        double k_minus_q = std::sqrt(data->k_value*data->k_value + q*q - 2.0*k_dot_q);

        Mpc_units::energy q_in_Mpc(q);
        Mpc_units::inverse_energy3 Pq  = (*(data->Pk))(q_in_Mpc);
        Mpc_units::inverse_energy qqPq = std::sin(theta) * q_in_Mpc * q_in_Mpc * Pq;

        // integral is P(q) gamma(k-r,r) alpha(k,-r)
        double gamma1 = (k_minus_q*k_minus_q*k_dot_q - q*q*k_dot_q + data->k_value*data->k_value*k_dot_q) / (2.0 * q*q * k_minus_q*k_minus_q);
        double alpha2 = (2.0*data->k_value*data->k_value*q*q - k_dot_q*(data->k_value*data->k_value + q*q)) / (2.0 * q*q * data->k_value*data->k_value);

        f[0] = data->jacobian * 3.0 * (qqPq/Mpc_units::Mpc) * gamma1*alpha2;

        return(0);
      }


    static int E_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        double q     = data->q_min + x[0] * data->q_range;
        double theta = 2.0 * M_PI * x[1];
        double phi   = M_PI * x[2];

        double k_dot_q   = data->k_value*q*std::cos(theta);
        double k_minus_q = std::sqrt(data->k_value*data->k_value + q*q - 2.0*k_dot_q);

        Mpc_units::energy q_in_Mpc(q);
        Mpc_units::inverse_energy3 Pq  = (*(data->Pk))(q_in_Mpc);
        Mpc_units::inverse_energy qqPq = std::sin(theta) * q_in_Mpc * q_in_Mpc * Pq;

        // integral is P(q) gamma(k-r,r) gamma(k,-r)
        double gamma1 = (k_minus_q*k_minus_q*k_dot_q - q*q*k_dot_q + data->k_value*data->k_value*k_dot_q) / (2.0 * q*q * k_minus_q*k_minus_q);
        double gamma2 = (2.0*data->k_value*data->k_value - k_dot_q*(2.0*data->k_value*data->k_value + 2.0*q*q - 2.0*k_dot_q)) / (2.0 * q*q * data->k_value*data->k_value);

        f[0] = data->jacobian * 3.0 * (qqPq/Mpc_units::Mpc) * gamma1*gamma2;

        return(0);
      }


    static int F_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        double q     = data->q_min + x[0] * data->q_range;
        double theta = 2.0 * M_PI * x[1];
        double phi   = M_PI * x[2];

        double k_dot_q   = data->k_value*q*std::cos(theta);
        double k_minus_q = std::sqrt(data->k_value*data->k_value + q*q - 2.0*k_dot_q);

        Mpc_units::energy q_in_Mpc(q);
        Mpc_units::inverse_energy3 Pq  = (*(data->Pk))(q_in_Mpc);
        Mpc_units::inverse_energy qqPq = std::sin(theta) * q_in_Mpc * q_in_Mpc * Pq;

        // integral is P(q) alpha(k-r,r) alpha(k,-r)
        double alpha1 = (k_minus_q*k_minus_q*k_dot_q + q*q*data->k_value*data->k_value - q*q*k_dot_q) / (2.0 * q*q * k_minus_q*k_minus_q);
        double alpha2 = (2.0*data->k_value*data->k_value*q*q - k_dot_q*(data->k_value*data->k_value + q*q)) / (2.0 * q*q * data->k_value*data->k_value);

        f[0] = data->jacobian * 3.0 * (qqPq/Mpc_units::Mpc) * alpha1*alpha2;

        return(0);
      }


    static int G_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        double q     = data->q_min + x[0] * data->q_range;
        double theta = 2.0 * M_PI * x[1];
        double phi   = M_PI * x[2];

        double k_dot_q   = data->k_value*q*std::cos(theta);
        double k_minus_q = std::sqrt(data->k_value*data->k_value + q*q - 2.0*k_dot_q);

        Mpc_units::energy q_in_Mpc(q);
        Mpc_units::inverse_energy3 Pq  = (*(data->Pk))(q_in_Mpc);
        Mpc_units::inverse_energy qqPq = std::sin(theta) * q_in_Mpc * q_in_Mpc * Pq;

        // integral is P(q) alpha(k-r,r) gamma(k,-r)
        double alpha1 = (k_minus_q*k_minus_q*k_dot_q + q*q*data->k_value*data->k_value - q*q*k_dot_q) / (2.0 * q*q * k_minus_q*k_minus_q);
        double gamma2 = (2.0*data->k_value*data->k_value - k_dot_q*(2.0*data->k_value*data->k_value + 2.0*q*q - 2.0*k_dot_q)) / (2.0 * q*q * data->k_value*data->k_value);

        f[0] = data->jacobian * 3.0 * (qqPq/Mpc_units::Mpc) * alpha1*gamma2;

        return(0);
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
