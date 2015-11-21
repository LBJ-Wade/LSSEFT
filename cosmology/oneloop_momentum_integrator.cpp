//
// Created by David Seery on 21/11/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "oneloop_momentum_integrator.h"

#include "cuba.h"


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
    constexpr unsigned int max_eval              = 50000;


    class integrand_data
      {

      public:

        integrand_data(const FRW_model& m, const eV_units::energy& _k, const eV_units::energy& UV, const eV_units::energy& IR,
                       const std::shared_ptr<tree_power_spectrum>& _Pk)
          : model(m),
            k(_k),
            UV_cutoff(UV),
            IR_cutoff(IR),
            Pk(_Pk),
            jacobian(2.0*M_PI_2*(static_cast<double>(UV_cutoff) - static_cast<double>(IR_cutoff))),    // Jacobian in angular directions in 2pi * pi = 2pi^2
            q_min(static_cast<double>(IR_cutoff)),                                                     // note guaranteed to be constructed after IR_cutoff, UV_cutoff
            q_max(static_cast<double>(UV_cutoff)),
            q_range(q_max - q_min),
            k_value(static_cast<double>(k))
          {
          }

        const FRW_model& model;
        const eV_units::energy& k;
        const eV_units::energy& UV_cutoff;
        const eV_units::energy& IR_cutoff;
        const std::shared_ptr<tree_power_spectrum>& Pk;

        double jacobian;
        double q_min;
        double q_max;
        double q_range;

        double k_value;
      };


    static int A_integrand(const int *ndim, const cubareal x[], const int *ncomp, cubareal f[], void *userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);

        double q     = data->q_min + x[0] * data->q_range;
        double theta = 2.0 * M_PI * x[1];
        double phi   = M_PI * x[2];

        double k_dot_q   = data->k_value*q*std::cos(theta);
        double k_minus_q = std::sqrt(data->k_value*data->k_value + q*q - 2.0*k_dot_q);

        eV_units::energy q_in_eV(q);
        eV_units::energy k_minus_q_in_eV(k_minus_q);

        // integral is P(q) P(k-q) alpha(q,k-q)
        eV_units::inverse_energy3 Pq         = (*(data->Pk))(q_in_eV);
        eV_units::inverse_energy3 Pk_minus_q = k_minus_q > data->q_min ? (*(data->Pk))(k_minus_q_in_eV) : eV_units::inverse_energy3(0);

        eV_units::inverse_energy  qqPq = q_in_eV*q_in_eV*Pq;
        eV_units::inverse_energy4 PP_prod = qqPq*Pk_minus_q;

        double alpha = (k_minus_q*k_minus_q*k_dot_q + q*q*data->k_value*data->k_value - q*q*k_dot_q) / (2.0 * q*q * k_minus_q*k_minus_q);

        f[0] = data->jacobian * static_cast<double>(PP_prod) * alpha;

        return(0);  // return value irrelevant unless = -999, which means stop integration
      }

  }


oneloop_momentum_integrator::oneloop_momentum_integrator(double a, double r)
  : abs_err(std::fabs(a)),
    rel_err(std::fabs(r))
  {
  }


loop_integral oneloop_momentum_integrator::integrate(const FRW_model& model, const eV_units::energy& k,
                                                     const k_token& k_tok, const eV_units::energy& UV_cutoff,
                                                     const UV_token& UV_tok, const eV_units::energy& IR_cutoff,
                                                     const IR_token& IR_tok, std::shared_ptr<tree_power_spectrum>& Pk)
  {
    eV_units::inverse_energy3 A = this->A_integral(model, k, UV_cutoff, IR_cutoff, Pk);

    return loop_integral(k, k_tok, UV_cutoff, UV_tok, IR_cutoff, IR_tok,
                         A, eV_units::inverse_energy3(0.0), 0.0, 0.0, 0.0, 0.0);
  }


eV_units::inverse_energy3 oneloop_momentum_integrator::A_integral(const FRW_model& model, const eV_units::energy& k,
                                                                  const eV_units::energy& UV_cutoff,
                                                                  const eV_units::energy& IR_cutoff,
                                                                  std::shared_ptr<tree_power_spectrum>& Pk)
  {
    cubareal integral[oneloop_momentum_impl::dimensions];
    cubareal error[oneloop_momentum_impl::dimensions];
    cubareal prob[oneloop_momentum_impl::dimensions];

    int regions;
    int evaluations;
    int fail;

    std::unique_ptr<oneloop_momentum_impl::integrand_data> data = std::make_unique<oneloop_momentum_impl::integrand_data>(model, k, UV_cutoff, IR_cutoff, Pk);

    Cuhre(oneloop_momentum_impl::dimensions, oneloop_momentum_impl::components,
          &oneloop_momentum_impl::A_integrand, data.get(),
          oneloop_momentum_impl::points_per_invocation,
          this->rel_err, this->abs_err,
          oneloop_momentum_impl::verbosity_none | oneloop_momentum_impl::samples_last,
          oneloop_momentum_impl::min_eval, oneloop_momentum_impl::max_eval,
          LSSEFT_DEFAULT_CUHRE_KEY, nullptr, nullptr,
          &regions, &evaluations, &fail,
          integral, error, prob);

    std::cout << "Cuhre result: regions = " << regions << ", evaluations = " << evaluations << ", fail = " << fail << ", value = " << integral[0] << ", error = " << error[0] << ", probability = " << prob[0] << '\n';

    // an overall factor 2 / (2pi)^3 is taken out of the integrand, so remember to put it back here
    return eV_units::inverse_energy3(integral[0] / (8.0 * M_PI_2*M_PI));
  }
