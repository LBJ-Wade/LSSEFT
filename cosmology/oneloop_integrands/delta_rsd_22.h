//
// Created by David Seery on 15/11/2016.
// Copyright (c) 2016 University of Sussex. All rights reserved.
//

#ifndef LSSEFT_DELTA_RSD_22_H
#define LSSEFT_DELTA_RSD_22_H


#include "shared.h"


namespace oneloop_momentum_impl
  {
    
    static int RSD22_A1_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
    
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;

        Mpc_units::energy2 s_sq = data->k*data->k + q*q - 2*data->k*q*z;
        Mpc_units::energy  s    = std::sqrt(s_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;
        
        Mpc_units::inverse_energy3 Pq = data->Pk(q);
        Mpc_units::inverse_energy3 Ps = s > data->IR_cutoff && s < data->UV_cutoff ? data->Pk(s) : Mpc_units::inverse_energy3(0);
        
        Mpc_units::energy2 intg = (data->k*data->k*data->k*data->k/(s_sq*s_sq)) * (data->k - 2.0*q*z) * (data->k - 2.0*q*z) * (z*z - 1.0);
        
        f[0] = (data->jacobian_2d * intg * Pq * Ps) / Mpc_units::Mpc3;

        return 0;
      }
    
    
    static int RSD22_A2_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;
        
        Mpc_units::energy2 s_sq = data->k*data->k + q*q - 2*data->k*q*z;
        Mpc_units::energy  s    = std::sqrt(s_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;
        
        Mpc_units::inverse_energy3 Pq = data->Pk(q);
        Mpc_units::inverse_energy3 Ps = s > data->IR_cutoff && s < data->UV_cutoff ? data->Pk(s) : Mpc_units::inverse_energy3(0);
        
        Mpc_units::energy2 intg = (data->k*data->k*data->k*data->k/(s_sq*s_sq)) * q * (q + data->k*z - 2.0*q*z*z) * (z*z - 1.0);
        
        f[0] = (data->jacobian_2d * intg * Pq * Ps) / Mpc_units::Mpc3;
        
        return 0;
      }
    
    
    static int RSD22_A3_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;
        
        Mpc_units::energy2 s_sq = data->k*data->k + q*q - 2*data->k*q*z;
        Mpc_units::energy  s    = std::sqrt(s_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;
        
        Mpc_units::inverse_energy3 Pq = data->Pk(q);
        Mpc_units::inverse_energy3 Ps = s > data->IR_cutoff && s < data->UV_cutoff ? data->Pk(s) : Mpc_units::inverse_energy3(0);
        
        Mpc_units::energy2 intg = (data->k*data->k*data->k*data->k/(s_sq*s_sq)) * q * z * (data->k - q*z) * (z*z - 1.0);
        
        f[0] = (data->jacobian_2d * intg * Pq * Ps) / Mpc_units::Mpc3;
        
        return 0;
      }
    
    
    static int RSD22_A5_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;
        
        Mpc_units::energy2 s_sq = data->k*data->k + q*q - 2*data->k*q*z;
        Mpc_units::energy  s    = std::sqrt(s_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;
        
        Mpc_units::inverse_energy3 Pq = data->Pk(q);
        Mpc_units::inverse_energy3 Ps = s > data->IR_cutoff && s < data->UV_cutoff ? data->Pk(s) : Mpc_units::inverse_energy3(0);
        
        Mpc_units::energy2 intg = (data->k*data->k*data->k*data->k/(s_sq*s_sq)) * z * (data->k*data->k*z + q*q*z*(2.0*z*z-1.0) + data->k*(q - 3.0*q*z*z));
        
        f[0] = (data->jacobian_2d * intg * Pq * Ps) / Mpc_units::Mpc3;
        
        return 0;
      }
    
    
    static int RSD22_B2_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;
        
        Mpc_units::energy2 s_sq = data->k*data->k + q*q - 2*data->k*q*z;
        Mpc_units::energy  s    = std::sqrt(s_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;
        
        Mpc_units::inverse_energy3 Pq = data->Pk(q);
        Mpc_units::inverse_energy3 Ps = s > data->IR_cutoff && s < data->UV_cutoff ? data->Pk(s) : Mpc_units::inverse_energy3(0);
        
        Mpc_units::energy2 intg = (data->k*data->k*data->k*data->k/(s_sq*s_sq)) * q * q * (z*z - 1.0);
        
        f[0] = (data->jacobian_2d * intg * Pq * Ps) / Mpc_units::Mpc3;
        
        return 0;
      }
    
    
    static int RSD22_B3_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;
        
        Mpc_units::energy2 s_sq = data->k*data->k + q*q - 2*data->k*q*z;
        Mpc_units::energy  s    = std::sqrt(s_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;
        
        Mpc_units::inverse_energy3 Pq = data->Pk(q);
        Mpc_units::inverse_energy3 Ps = s > data->IR_cutoff && s < data->UV_cutoff ? data->Pk(s) : Mpc_units::inverse_energy3(0);
        
        Mpc_units::energy2 intg = (data->k*data->k*data->k*data->k/(s_sq*s_sq)) * (2.0*q*q + (data->k - 2.0*q*z)*(data->k + 2.0*q*z + data->k*z*z - 2.0*q*z*z*z));
        
        f[0] = (data->jacobian_2d * intg * Pq * Ps) / Mpc_units::Mpc3;
        
        return 0;
      }
    
    
    static int RSD22_B6_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;
        
        Mpc_units::energy2 s_sq = data->k*data->k + q*q - 2*data->k*q*z;
        Mpc_units::energy  s    = std::sqrt(s_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;
        
        Mpc_units::inverse_energy3 Pq = data->Pk(q);
        Mpc_units::inverse_energy3 Ps = s > data->IR_cutoff && s < data->UV_cutoff ? data->Pk(s) : Mpc_units::inverse_energy3(0);
        
        Mpc_units::energy2 intg = (data->k*data->k*data->k*data->k/(s_sq*s_sq)) * z * (data->k - q*z) * (q + data->k*z - 2.0*q*z*z);
        
        f[0] = (data->jacobian_2d * intg * Pq * Ps) / Mpc_units::Mpc3;
        
        return 0;
      }
    
    
    static int RSD22_B8_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;
        
        Mpc_units::energy2 s_sq = data->k*data->k + q*q - 2*data->k*q*z;
        Mpc_units::energy  s    = std::sqrt(s_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;
        
        Mpc_units::inverse_energy3 Pq = data->Pk(q);
        Mpc_units::inverse_energy3 Ps = s > data->IR_cutoff && s < data->UV_cutoff ? data->Pk(s) : Mpc_units::inverse_energy3(0);
        
        Mpc_units::energy2 intg = (data->k*data->k*data->k*data->k/(s_sq*s_sq)) * (2.0*data->k*data->k*z*z + data->k*q*z*(3.0-7.0*z*z) + q*q*(1.0-5.0*z*z+6.0*z*z*z*z));
        
        f[0] = (data->jacobian_2d * intg * Pq * Ps) / Mpc_units::Mpc3;
        
        return 0;
      }
    
    
    static int RSD22_B9_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;
        
        Mpc_units::energy2 s_sq = data->k*data->k + q*q - 2*data->k*q*z;
        Mpc_units::energy  s    = std::sqrt(s_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;
        
        Mpc_units::inverse_energy3 Pq = data->Pk(q);
        Mpc_units::inverse_energy3 Ps = s > data->IR_cutoff && s < data->UV_cutoff ? data->Pk(s) : Mpc_units::inverse_energy3(0);
        
        Mpc_units::energy2 intg = (data->k*data->k*data->k*data->k/(s_sq*s_sq)) * z * (2.0*data->k*data->k*z*z + q*q*z*(-1.0+3.0*z*z) + data->k*(q - 5.0*q*z*z));
        
        f[0] = (data->jacobian_2d * intg * Pq * Ps) / Mpc_units::Mpc3;
        
        return 0;
      }
    
    
    static int RSD22_C1_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;
        
        Mpc_units::energy2 s_sq = data->k*data->k + q*q - 2*data->k*q*z;
        Mpc_units::energy  s    = std::sqrt(s_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;
        
        Mpc_units::inverse_energy3 Pq = data->Pk(q);
        Mpc_units::inverse_energy3 Ps = s > data->IR_cutoff && s < data->UV_cutoff ? data->Pk(s) : Mpc_units::inverse_energy3(0);
        
        Mpc_units::energy2 intg = (data->k*data->k*data->k*data->k/(s_sq*s_sq)) * (q + 2.0*data->k*z - 3.0*q*z*z) * (q + data->k*z - 2.0*q*z*z);
        
        f[0] = (data->jacobian_2d * intg * Pq * Ps) / Mpc_units::Mpc3;
        
        return 0;
      }
    
    
    static int RSD22_C2_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;
        
        Mpc_units::energy2 s_sq = data->k*data->k + q*q - 2*data->k*q*z;
        Mpc_units::energy  s    = std::sqrt(s_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;
        
        Mpc_units::inverse_energy3 Pq = data->Pk(q);
        Mpc_units::inverse_energy3 Ps = s > data->IR_cutoff && s < data->UV_cutoff ? data->Pk(s) : Mpc_units::inverse_energy3(0);
        
        Mpc_units::energy2 intg = (data->k*data->k*data->k*data->k/(s_sq*s_sq)) * z * (data->k - q*z) * (q + 2.0*data->k*z - 3.0*q*z*z);
        
        f[0] = (data->jacobian_2d * intg * Pq * Ps) / Mpc_units::Mpc3;
        
        return 0;
      }
    
    
    static int RSD22_C4_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;
        
        Mpc_units::energy2 s_sq = data->k*data->k + q*q - 2*data->k*q*z;
        Mpc_units::energy  s    = std::sqrt(s_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;
        
        Mpc_units::inverse_energy3 Pq = data->Pk(q);
        Mpc_units::inverse_energy3 Ps = s > data->IR_cutoff && s < data->UV_cutoff ? data->Pk(s) : Mpc_units::inverse_energy3(0);
        
        Mpc_units::energy2 intg = (data->k*data->k*data->k*data->k/(s_sq*s_sq)) * (z*z - 1.0) * (2.0*data->k*data->k - 12.0*data->k*q*z + 3.0*q*q*(-1.0 + 5.0*z*z));
        
        f[0] = (data->jacobian_2d * intg * Pq * Ps) / Mpc_units::Mpc3;
        
        return 0;
      }
    
    
    static int RSD22_D1_integrand(const int* ndim, const cubareal x[], const int* ncomp, cubareal f[], void* userdata)
      {
        oneloop_momentum_impl::integrand_data* data = static_cast<integrand_data*>(userdata);
        
        Mpc_units::energy q = data->IR_cutoff + x[0] * data->q_range;
        double z            = 2.0*x[1] - 1.0;
        
        Mpc_units::energy2 s_sq = data->k*data->k + q*q - 2*data->k*q*z;
        Mpc_units::energy  s    = std::sqrt(s_sq * Mpc_units::Mpc2) / Mpc_units::Mpc;
        
        Mpc_units::inverse_energy3 Pq = data->Pk(q);
        Mpc_units::inverse_energy3 Ps = s > data->IR_cutoff && s < data->UV_cutoff ? data->Pk(s) : Mpc_units::inverse_energy3(0);
        
        Mpc_units::energy2 intg = (data->k*data->k*data->k*data->k/(s_sq*s_sq)) * (8.0*data->k*q*z*(3.0 - 5.0*z*z) + 4.0*data->k*data->k*(-1.0 + 3.0*z*z) + q*q*(3.0 - 30.0*z*z + 35.0*z*z*z*z));
        
        f[0] = (data->jacobian_2d * intg * Pq * Ps) / Mpc_units::Mpc3;
        
        return 0;
      }
    
  }   // namespace oneloop_momentum_impl

#endif //LSSEFT_DELTA_RSD_22_H
