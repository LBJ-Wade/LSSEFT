//
// Created by David Seery on 15/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_DELTA_RSD_13_H
#define LSSEFT_DELTA_RSD_13_H


#include "shared.h"


namespace oneloop_momentum_impl
  {
    
    // all of these kernels are just 1D integrals, but we treat them as 2d integrals with a trivial
    // z direction. The Cuhre integrator does not seem to work for just 1 dimension, as can
    // be confirmed eg. by using the Mathematica interface
    
    static int RSD13_a_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
    
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
    
        double atanh = (q > data->k) ? std::atanh(data->k / q) : std::atanh(q / data->k);
        Mpc_units::energy2 intg = (data->k / q) * (data->k / q) * (data->k / q) * data->k * (q - atanh * data->k);
    
        f[0] = data->jacobian_1d * intg * data->Pk(q);
    
        return 0;
      }
    
    
    static int RSD13_b_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        
        double atanh = (q > data->k) ? std::atanh(data->k / q) : std::atanh(q / data->k);
        Mpc_units::energy2 intg = (data->k / q) * data->k * data->k * atanh;
        
        f[0] = data->jacobian_1d * intg * data->Pk(q);
        
        return 0;
      }
    
    
    static int RSD13_c_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        
        Mpc_units::energy2 intg = data->k * data->k;
        
        f[0] = data->jacobian_1d * intg * data->Pk(q);
        
        return 0;
      }
    
    
    static int RSD13_d_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        
        double atanh = (q > data->k) ? std::atanh(data->k / q) : std::atanh(q / data->k);
        Mpc_units::energy2 intg = data->k * q * atanh;
        
        f[0] = data->jacobian_1d * intg * data->Pk(q);
        
        return 0;
      }
    
    
    static int RSD13_e_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        
        Mpc_units::energy2 intg = q * q;
        
        f[0] = data->jacobian_1d * intg * data->Pk(q);
        
        return 0;
      }
    
    
    static int RSD13_f_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        
        double atanh = (q > data->k) ? std::atanh(data->k / q) : std::atanh(q / data->k);
        Mpc_units::energy2 intg = (q / data->k) * q * q * atanh;
        
        f[0] = data->jacobian_1d * intg * data->Pk(q);
        
        return 0;
      }
    
    
    static int RSD13_g_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        
        double atanh = (q > data->k) ? std::atanh(data->k / q) : std::atanh(q / data->k);
        Mpc_units::energy2 intg = (q / data->k) * (q / data->k) * (q / data->k) * q * (data->k - q * atanh);
        
        f[0] = data->jacobian_1d * intg * data->Pk(q);
        
        return 0;
      }
    
  }   // namespace oneloop_momentum_impl

#endif //LSSEFT_DELTA_RSD_13_H
