//
// Created by David Seery on 15/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_DELTA_13_H
#define LSSEFT_DELTA_13_H


#include "shared.h"


namespace oneloop_momentum_impl
  {
    
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
        
        f[0] = 8.0 * data->jacobian_2d * qqPq * gamma1*alpha2;
        
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
        
        f[0] = 8.0 * data->jacobian_2d * qqPq * gamma1*gamma2;
        
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
        
        f[0] = 8.0 * data->jacobian_2d * qqPq * alpha1*alpha2;
        
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
        
        f[0] = 8.0 * data->jacobian_2d * qqPq * alpha1*gamma2;
        
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
        
        f[0] = 8.0 * data->jacobian_2d * qqPq * alpha_asym * alpha2;
        
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
        
        f[0] = 8.0 * data->jacobian_2d * qqPq * alpha_asym * gamma2;
        
        return(0);  // return value irrelevant unless = -999, which means stop integration
      }
    
  }   // namespace oneloop_momentum_impl

#endif //LSSEFT_DELTA_13_H
