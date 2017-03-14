//
// Created by David Seery on 07/12/2016.
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

#include <cmath>

#include "Pk_filter.h"


namespace Pk_filter_impl
  {
    
    constexpr unsigned int dimensions = 2;       // Cuhre only works in >= 2 dimensions, not 1
    constexpr unsigned int components = 1;
    constexpr unsigned int points_per_invocation = 1;
    
    constexpr unsigned int verbosity_none = 0;
    constexpr unsigned int verbosity_reasonable = 1;
    constexpr unsigned int verbosity_progress = 2;
    constexpr unsigned int verbosity_subregions = 3;
    
    constexpr unsigned int samples_all = 0;
    constexpr unsigned int samples_last = 4;
    
    constexpr unsigned int min_eval = 0;
    constexpr unsigned int max_eval = 20000000;
    
    constexpr unsigned int ngiven = 0;
    constexpr unsigned int ldxgiven = 0;
    constexpr unsigned int nextra = 0;
    
    constexpr unsigned int pcores = 10000;   // matches default Cuba value
    
    constexpr unsigned int cuhre_key = 13;      // degree-13 only available in 2-dimensions
    constexpr unsigned int divonne_key1 = 47;
    constexpr unsigned int divonne_key2 = 13;      // degree-13 only available in 2-dimensions
    constexpr unsigned int divonne_key3 = 1;
    constexpr unsigned int divonne_maxpass = 5;
    constexpr unsigned int divonne_border = 0;
    constexpr double divonne_maxchisq = 10.0;
    constexpr double divonne_minchisq = 0.25;
    
    
    class integrand_data
      {
      
      public:
        
        integrand_data(double sl_min, double sl_max, double kl, double lm, const filterable_Pk& _Pk, const approx_Pk& _Pka)
          : slog_min(sl_min),
            slog_max(sl_max),
            klog(kl),
            lambda(lm),
            Pk(_Pk),
            Pk_approx(_Pka),
            jacobian(slog_max - slog_min),
            slog_range(slog_max - slog_min)
          {
          }

        const double slog_min;
        const double slog_max;
        const double klog;
        const double lambda;

        const filterable_Pk& Pk;
        const approx_Pk& Pk_approx;

        const double jacobian;
        const double slog_range;
      };
    
    
    static int filter_integrand(const int* ndim, const cubareal* x, const int* ncomp, cubareal* f, void* userdata)
      {
        Pk_filter_impl::integrand_data* data = static_cast<Pk_filter_impl::integrand_data*>(userdata);
    
        double slog = data->slog_min + data->slog_range*x[0];
        Mpc_units::energy s = std::pow(10.0, slog) / Mpc_units::Mpc;
    
        f[0] = data->jacobian * (data->Pk(s) / data->Pk_approx(s)) * std::exp(-(data->klog-slog)*(data->klog-slog) / (2.0*data->lambda*data->lambda));
    
        return(0);  // return value irrelevant unless = -999, which means stop integration
      }
    
    
    static int window_integrand(const int* ndim, const cubareal* x, const int* ncomp, cubareal* f, void* userdata)
      {
        Pk_filter_impl::integrand_data* data = static_cast<Pk_filter_impl::integrand_data*>(userdata);
        
        double slog = data->slog_min + data->slog_range*x[0];
        
        f[0] = data->jacobian * std::exp(-(data->klog-slog)*(data->klog-slog) / (2.0*data->lambda*data->lambda));
        
        return(0);  // return value irrelevant unless = -999, which means stop integration
      }
    
  }   // namespace Pk_filter_impl


std::pair< Mpc_units::inverse_energy3, Mpc_units::inverse_energy3 >
Pk_filter::operator()(const FRW_model& model, const filterable_Pk& Pk_lin, const Mpc_units::energy& k)
  {
    // build reference Eisenstein & Hu power spectrum
    std::unique_ptr<approx_Pk> Papprox = this->eisenstein_hu(model, Pk_lin);
    
    // get maximum available scale from linear power spectrum
    const Mpc_units::energy k_min = SPLINE_PK_DEFAULT_BOTTOM_CLEARANCE * Pk_lin.get_db().get_k_min();
    const Mpc_units::energy k_max = SPLINE_PK_DEFAULT_TOP_CLEARANCE * Pk_lin.get_db().get_k_max();
    
    const double klog = std::log10(k * Mpc_units::Mpc);
    
    const double slog_max = std::log10(k_max * Mpc_units::Mpc);
    const double slog_min = std::log10(k_min * Mpc_units::Mpc);

    const Mpc_units::energy pivot = this->params.get_pivot();
    const double amplitude = this->params.get_amplitude();
    const double index = this->params.get_index();
    
    const double lambda = amplitude * std::pow(k/pivot, index);
    
    try
      {
        const double filtered_Pk = this->integrate(slog_min, slog_max, klog, lambda, Pk_lin, *Papprox, Pk_filter_impl::filter_integrand);
        const double volume      = this->integrate(slog_min, slog_max, klog, lambda, Pk_lin, *Papprox, Pk_filter_impl::window_integrand);
    
        const Mpc_units::inverse_energy3 P_nw = (*Papprox)(k) * filtered_Pk / volume;
    
        return std::make_pair(P_nw, (*Papprox)(k));
      }
    catch(runtime_exception& xe)
      {
      }
    
    std::ostringstream msg;
    msg << LSSEFT_PK_FILTER_FAIL << " k = " << k * Mpc_units::Mpc << " h/Mpc";
    throw runtime_exception(exception_type::filter_failure, msg.str());
  }


double Pk_filter::integrate(const double slog_min, const double slog_max, const double klog, const double lambda,
                            const filterable_Pk& Pk_lin, const approx_Pk& Papprox, integrand_t integrand)
  {
    cubareal integral[Pk_filter_impl::dimensions];
    cubareal error[Pk_filter_impl::dimensions];
    cubareal prob[Pk_filter_impl::dimensions];
    
    int regions;
    int evaluations;
    int fail;

    std::unique_ptr<Pk_filter_impl::integrand_data> data =
      std::make_unique<Pk_filter_impl::integrand_data>(slog_min, slog_max, klog, lambda, Pk_lin, Papprox);
    
    Cuhre(Pk_filter_impl::dimensions,
          Pk_filter_impl::components,
          integrand, data.get(),
          Pk_filter_impl::points_per_invocation,
          this->params.get_relerr(), this->params.get_abserr(),
          Pk_filter_impl::verbosity_none | Pk_filter_impl::samples_last,
          Pk_filter_impl::min_eval, Pk_filter_impl::max_eval,
          Pk_filter_impl::cuhre_key,
          nullptr, nullptr,
          &regions, &evaluations, &fail,
          integral, error, prob);
    
    if(fail != 0 || !std::isfinite(integral[0]) || !std::isfinite(error[0]))
      throw runtime_exception(exception_type::filter_failure);
    
    return integral[0];
  }


std::unique_ptr<approx_Pk> Pk_filter::eisenstein_hu(const FRW_model& model, const filterable_Pk& Pk_lin)
  {
    // extract database of power spectrum sample points
    const tree_Pk::database_type& db = Pk_lin.get_db();
    
    // compute Eisenstein & Hu parameters
    double omega_m          = model.get_omega_m();
    double omega_cc         = model.get_omega_cc();
    double h                = model.get_h();
    double Neff             = model.get_Neff();
    Mpc_units::energy T_CMB = model.get_T_CMB();
    double f_baryon         = model.get_f_baryon();
    double z_star           = model.get_z_star();
    double z_drag           = model.get_z_drag();
    double z_eq             = model.get_z_eq();
    double A_curv           = model.get_A_curv();
    double ns               = model.get_ns();
    Mpc_units::energy k_piv = model.get_k_piv();
    
    constexpr double omega0 = 1.0;
    
    double Theta27 = T_CMB / Mpc_units::Kelvin / 2.7;

    double omega_b = f_baryon * omega_m;
    double omega_c = omega_m - omega_b;
    
    // formulae for Eisenstein & Hu fitting function from astro-ph/9710252
    
    // y and s defined on p. 4 of Eisenstein & Hu paper
    double y = (1.0 + z_eq) / (1.0 + z_drag);
    Mpc_units::inverse_energy s = 44.5 * std::log(9.83 / (omega_m * h * h))
                                  / std::sqrt(1.0 + 10.0 * std::pow(omega_b * h * h, 3.0 / 4.0)) * Mpc_units::Mpc;
    
    // fs defined at top of Sec 3.1 on p.4 of Eisenstein & Hu paper
    double fb  = omega_b / omega_m;
    double fc  = omega_c / omega_m;
    double fcb = fc + fb;
    
    // ps defined in Eq. (11) on p.5 of Eisenstein & Hu paper
    double pc  = (1.0/4.0) * (5.0 - std::sqrt(1.0 + 24.0*fc));
    double pb  = (1.0/4.0) * (5.0 - std::sqrt(1.0 + 24.0*fb));
    double pcb = (1.0/4.0) * (5.0 - std::sqrt(1.0 + 24.0*fcb));
    
    // alpha_nu defined in Eq. (15) on p.6 of Eisenstein & Hu paper
    double anu = (fc / fcb) * (5.0 - 2.0 * (pc + pcb)) / (5.0 - 4.0 * pcb) * std::pow(1.0 + y, pcb - pc) *
                 (1.0 + (pc - pcb) / 2.0 * (1.0 + 1.0 / ((3.0 - 4.0 * pc) * (7.0 - 4.0 * pcb))) / (1 + y));
    
    approx_Pk::database_type approx_db;
    
    // normalization is adjusted so that we match the input power spectrum on the largest scale
    boost::optional<double> normalization = boost::none;

    for(tree_Pk::database_type::const_record_iterator t = db.record_cbegin(); t != db.record_cend(); ++t)
      {
        const Mpc_units::energy& k = t->get_wavenumber();

        // Gamma_eff defined in Eq. (16) on p.6 of Eisenstein & Hu paper
        double sqrt_anu = std::sqrt(anu);
        double Gamma_eff = omega_m*h*h * (sqrt_anu + (1.0 - sqrt_anu) / (1.0 + std::pow(0.43*k*s, 4)));
        
        // q_eff defined in Eq. (5) p.4 of Eisenstein & Hu paper
        double q_eff = (k*h * Mpc_units::Mpc) * Theta27*Theta27 / Gamma_eff;
        
        // beta_c = 1 if there are no neutrinos
        double L = std::log(std::exp(1.0) + 1.84*sqrt_anu*q_eff);
        double C = 14.4 + 325.0 / (1.0 + 60.5*std::pow(q_eff, 1.08));
        
        double T = L / (L + C*q_eff*q_eff);
        
        // power spectrum is proportional to k^4 T^2 P_Phi(k), but we don't have to get the constant
        // of proportionality correct since it will be adjusted below
        Mpc_units::inverse_energy3 Pk = T*T * std::pow(k*h / k_piv, ns-1.0) * k * Mpc_units::Mpc4;
        
        if(normalization)
          {
            // normalization has already been calculated
            Pk = *normalization * Pk;
          }
        else
          {
            // adjust normalization to match input power spectrum
            // (accounts for linear growth factor between Eisenstein & Hu formula which is valid at z=0,
            // and the input power spectrum which will be defined at z >> 1)
            normalization = t->get_Pk() / Pk;
            Pk = t->get_Pk();
          }
        
        approx_db.add_record(k, Pk);
      }
    
    return std::make_unique<approx_Pk>(approx_db);
  }
