//
// Created by David Seery on 17/08/2015.
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

#include <utility>

#include "oneloop_growth_integrator.h"
#include "constants.h"
#include "EdS_growth.h"

#include "units/Mpc_units.h"

#include "boost/numeric/odeint.hpp"


typedef std::vector<double> state_vector;

constexpr unsigned int RHO_M        = 0;
constexpr unsigned int RHO_R        = 1;
constexpr unsigned int ELEMENT_Dlin    = 2;
constexpr unsigned int ELEMENT_A    = 3;
constexpr unsigned int ELEMENT_B    = 4;
constexpr unsigned int ELEMENT_D    = 5;
constexpr unsigned int ELEMENT_E    = 6;
constexpr unsigned int ELEMENT_F    = 7;
constexpr unsigned int ELEMENT_G    = 8;
constexpr unsigned int ELEMENT_J    = 9;
constexpr unsigned int ELEMENT_dDlindz = 10;
constexpr unsigned int ELEMENT_dAdz = 11;
constexpr unsigned int ELEMENT_dBdz = 12;
constexpr unsigned int ELEMENT_dDdz = 13;
constexpr unsigned int ELEMENT_dEdz = 14;
constexpr unsigned int ELEMENT_dFdz = 15;
constexpr unsigned int ELEMENT_dGdz = 16;
constexpr unsigned int ELEMENT_dJdz = 17;

constexpr unsigned int STATE_SIZE   = 18;


class oneloop_functor
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    oneloop_functor(const FRW_model& m, const growth_params& p);

    //! destructor is default
    ~oneloop_functor() = default;


    // INTERFACE

  public:

    //! compute RHS of ODE system
    void operator()(const state_vector& x, state_vector& dxdz, double z);

    //! compute ics for ODE system
    void ics(state_vector& x, double z);


    // INTERNAL DATA

  private:

    //! reference to FRW model
    const FRW_model& model;
    
    //! reference to parameter block
    const growth_params& params;

    //! cache H0
    Mpc_units::energy H0;

    //! cache rho_cc
    Mpc_units::energy4 rho_cc;

  };


class oneloop_observer
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    oneloop_observer(oneloop_growth& c, const growth_params& p);

    //! destructor is default
    ~oneloop_observer() = default;


    // INTERFACE

  public:

    //! store
    void operator()(const state_vector& x, double z);


    // TIMING FUNCTIONS

  public:

    //! start integration timer
    void start_timer();

    //! stop integration timer
    void stop_timer();

    //! read integration timer
    boost::timer::nanosecond_type read_timer();


    // INTERNAL DATA

  private:

    //! reference to oneloop_growth container
    oneloop_growth& container;
    
    //! capture parameters used for this integration
    const growth_params& params;
    
    //! capture whether we are uisng EdS mode
    bool use_EdS;

    //! integration timer
    boost::timer::cpu_timer timer;

  };


// ONELOOP_INTEGRATOR METHODS


oneloop_growth_integrator::oneloop_growth_integrator(const growth_params& p, const growth_params_token& t)
  : params(p),
    token(t)
  {

  }


growth_integrator_data
oneloop_growth_integrator::integrate(const FRW_model& model, z_database& z_db)
  {
    // set up empty oneloop_growth container
    std::unique_ptr<oneloop_growth> ctr = std::make_unique<oneloop_growth>(this->token, z_db);

    // set up a functor for the ODE system
    oneloop_functor rhs(model, this->params);

    // set up an observer
    oneloop_observer obs(*ctr, this->params);

    // set up a state vector
    state_vector x(STATE_SIZE);

    // get initial time -- note use of reverse iterator to get last z!
    z_database::reverse_value_iterator max_z = z_db.value_rbegin();

    double init_z = *max_z;
    rhs.ics(x, init_z);

    // set up vector of sample times (needed until the Boost version of odeint catchs up to the github version, which includes my patch)
    std::vector<double> z_sample;
    std::copy(z_db.value_rbegin(), z_db.value_rend(), std::back_inserter(z_sample));

    // set up stepper
    auto stepper = boost::numeric::odeint::make_dense_output< boost::numeric::odeint::runge_kutta_dopri5<state_vector> >(this->params.get_abserr(), this->params.get_relerr());

    // run the integration
    // depending whether EdS mode is set in the growth parameter block, this either writes the EdS approximations
    // or the full one-loop result into the data container ctr
    obs.start_timer();
    size_t steps = boost::numeric::odeint::integrate_times(stepper, rhs, x, z_sample.begin(), z_sample.end(), -1E-3, obs);
    obs.stop_timer();
    
    return growth_integrator_data(std::move(ctr), obs.read_timer(), steps);
  }


// ONELOOP_FUNCTOR METHODS


oneloop_functor::oneloop_functor(const FRW_model& m, const growth_params& p)
  : model(m),
    params(p),
    H0(model.get_h() * 100.0 * Mpc_units::Kilometre / Mpc_units::Second / Mpc_units::Mpc),
    rho_cc(3.0 * H0*H0 * Mpc_units::PlanckMass*Mpc_units::PlanckMass * model.get_omega_cc())
  {
  }


void oneloop_functor::operator()(const state_vector& x, state_vector& dxdz, double z)
  {
    // compute rho, remembering that we use 1/Mpc^4 as units for background densities
    // during the integration
    state_vector::value_type rho = x[RHO_M] + x[RHO_R] + this->rho_cc * Mpc_units::Mpc4;

    // compute H in the same units
    state_vector::value_type Mp_in_Mpc_units = Mpc_units::PlanckMass * Mpc_units::Mpc;
    state_vector::value_type H               = std::sqrt(rho / (3.0 * Mp_in_Mpc_units*Mp_in_Mpc_units));

    state_vector::value_type Hdot            = -(3.0*x[RHO_M] + 4.0*x[RHO_R]) / (6.0 * Mp_in_Mpc_units*Mp_in_Mpc_units);
    state_vector::value_type epsilon         = -Hdot/(H*H);

    state_vector::value_type Omega_m         = x[RHO_M]/rho;
    state_vector::value_type Omega_r         = x[RHO_R]/rho;

    state_vector::value_type one_plus_z      = 1.0+z;
    state_vector::value_type one_plus_z_sq   = (1.0+z)*(1.0+z);

    // evolve background
    dxdz[RHO_M] = 3.0*x[RHO_M] / one_plus_z;
    dxdz[RHO_R] = 4.0*x[RHO_R] / one_plus_z;

    // evolve linear growth factor
    dxdz[ELEMENT_Dlin]    = x[ELEMENT_dDlindz];
    dxdz[ELEMENT_dDlindz] = (1.0 - epsilon) * x[ELEMENT_dDlindz] / one_plus_z
                            + (3.0 * Omega_m / 2.0) * x[ELEMENT_Dlin] / one_plus_z_sq;

    // evolve one-loop kernels -- derivatives
    dxdz[ELEMENT_A] = x[ELEMENT_dAdz];
    dxdz[ELEMENT_B] = x[ELEMENT_dBdz];
    dxdz[ELEMENT_D] = x[ELEMENT_dDdz];
    dxdz[ELEMENT_E] = x[ELEMENT_dEdz];
    dxdz[ELEMENT_F] = x[ELEMENT_dFdz];
    dxdz[ELEMENT_G] = x[ELEMENT_dGdz];
    dxdz[ELEMENT_J] = x[ELEMENT_dJdz];

    // evolve one-loop kernels -- second derivatives
    dxdz[ELEMENT_dAdz] = (3.0 * Omega_m / 2.0) * x[ELEMENT_Dlin] * x[ELEMENT_Dlin] / one_plus_z_sq
                         + (1.0 - epsilon) * x[ELEMENT_dAdz] / one_plus_z
                         + (3.0 * Omega_m / 2.0) * x[ELEMENT_A] / one_plus_z_sq;

    dxdz[ELEMENT_dBdz] = x[ELEMENT_dDlindz] * x[ELEMENT_dDlindz]
                         + (1.0 - epsilon) * x[ELEMENT_dBdz] / one_plus_z
                         + (3.0 * Omega_m / 2.0) * x[ELEMENT_B] / one_plus_z_sq;

    dxdz[ELEMENT_dDdz] = x[ELEMENT_dDlindz] * x[ELEMENT_dAdz]
                         + (1.0 - epsilon) * x[ELEMENT_dDdz] / one_plus_z
                         + (3.0 * Omega_m / 2.0) * x[ELEMENT_D] / one_plus_z_sq;

    dxdz[ELEMENT_dEdz] = x[ELEMENT_dDlindz] * x[ELEMENT_dBdz]
                         + (1.0 - epsilon) * x[ELEMENT_dEdz] / one_plus_z
                         + (3.0 * Omega_m / 2.0) * x[ELEMENT_E] / one_plus_z_sq;

    dxdz[ELEMENT_dFdz] = (3.0 * Omega_m / 2.0) * x[ELEMENT_Dlin] * x[ELEMENT_A] / one_plus_z_sq
                         + (1.0 - epsilon) * x[ELEMENT_dFdz] / one_plus_z
                         + (3.0 * Omega_m / 2.0) * x[ELEMENT_F] / one_plus_z_sq;

    dxdz[ELEMENT_dGdz] = (3.0 * Omega_m / 2.0) * x[ELEMENT_Dlin] * x[ELEMENT_B] / one_plus_z_sq
                         + (1.0 - epsilon) * x[ELEMENT_dGdz] / one_plus_z
                         + (3.0 * Omega_m / 2.0) * x[ELEMENT_G] / one_plus_z_sq;

    dxdz[ELEMENT_dJdz] = x[ELEMENT_dDlindz] * x[ELEMENT_dDlindz] * x[ELEMENT_Dlin]
                         + (1.0 - epsilon) * x[ELEMENT_dJdz] / one_plus_z
                         + (3.0 * Omega_m / 2.0) * x[ELEMENT_J] / one_plus_z_sq;
  }


void oneloop_functor::ics(state_vector& x, double z)
  {
    // compute (a0/a)^3 and (a0/a)^4
    state_vector::value_type a_three = (1.0+z)*(1.0+z)*(1.0+z);
    state_vector::value_type a_four  = (1.0+z)*a_three;

    // compute matter and radiation densities today
    // for radiation, we need the Stefan-Boltzman law and the present day CMB temperature
    Mpc_units::energy4 rho_m0 = this->model.get_omega_m() * (3.0 * this->H0*this->H0 * Mpc_units::PlanckMass*Mpc_units::PlanckMass);

    Mpc_units::energy T_CMB = this->model.get_T_CMB();

    Mpc_units::energy4 rho_r0 = g_star * radiation_constant * T_CMB*T_CMB*T_CMB*T_CMB * (1.0 + this->model.get_Neff() * (7.0/8.0) * std::pow(4.0/11.0, 4.0/3.0));

    Mpc_units::energy4 rho_m = rho_m0 * a_three;
    Mpc_units::energy4 rho_r = rho_r0 * a_four;

    // initial conditions for background; we need to decide which units are going to be used to
    // represent these during the integration.
    // We used 1/Mpc^4
    x[RHO_M] = rho_m * Mpc_units::Mpc4;
    x[RHO_R] = rho_r * Mpc_units::Mpc4;

    // initial conditions for linear growth factor
    // to get the derviative, assume that we are early in matter domination and D(z) grows like a(z).
    // so D(z) = a(z)/a_init(z). Then dD/dz = -da/dz / a_init = (a_0/a_init) / (1+z)^2 = 1/(1+z).
    state_vector::value_type D_lin = 1.0;
    state_vector::value_type f_lin = 1.0;
    
    x[ELEMENT_Dlin] = 1.0;
    x[ELEMENT_dDlindz] = -1.0 * D_lin * f_lin / (1.0+z);

    // initial conditions for one-loop growth functions/factors
    bool EdS_ics = this->params.use_EdS_ics();
    EdS_growth<state_vector::value_type> EdS(D_lin, f_lin);

    x[ELEMENT_A]    = EdS_ics ? EdS.DA() : 0.0;
    x[ELEMENT_B]    = EdS_ics ? EdS.DB() : 0.0;
    x[ELEMENT_D]    = EdS_ics ? EdS.DD() : 0.0;
    x[ELEMENT_E]    = EdS_ics ? EdS.DE() : 0.0;
    x[ELEMENT_F]    = EdS_ics ? EdS.DF() : 0.0;
    x[ELEMENT_G]    = EdS_ics ? EdS.DG() : 0.0;
    x[ELEMENT_J]    = EdS_ics ? EdS.DJ() : 0.0;
    x[ELEMENT_dAdz] = EdS_ics ? -1.0 * EdS.DA() * EdS.fA() / (1.0+z) : 0.0;
    x[ELEMENT_dBdz] = EdS_ics ? -1.0 * EdS.DB() * EdS.fB() / (1.0+z) : 0.0;
    x[ELEMENT_dDdz] = EdS_ics ? -1.0 * EdS.DD() * EdS.fD() / (1.0+z) : 0.0;
    x[ELEMENT_dEdz] = EdS_ics ? -1.0 * EdS.DE() * EdS.fE() / (1.0+z) : 0.0;
    x[ELEMENT_dFdz] = EdS_ics ? -1.0 * EdS.DF() * EdS.fF() / (1.0+z) : 0.0;
    x[ELEMENT_dGdz] = EdS_ics ? -1.0 * EdS.DG() * EdS.fG() / (1.0+z) : 0.0;
    x[ELEMENT_dJdz] = EdS_ics ? -1.0 * EdS.DJ() * EdS.fJ() / (1.0+z) : 0.0;
  }


// ONELOOP_OBSERVER METHODS


oneloop_observer::oneloop_observer(oneloop_growth& c, const growth_params& p)
  : container(c),
    params(p),
    use_EdS(p.use_EdS())
  {
    // stop timer
    timer.stop();
  }


void oneloop_observer::start_timer()
  {
    this->timer.start();
  }


void oneloop_observer::stop_timer()
  {
    this->timer.stop();
  }


boost::timer::nanosecond_type oneloop_observer::read_timer()
  {
    return(this->timer.elapsed().wall);
  }


void oneloop_observer::operator()(const state_vector& x, double z)
  {
    state_vector::value_type D_lin = x[ELEMENT_Dlin];
    state_vector::value_type f_lin = - (1.0+z) * x[ELEMENT_dDlindz] / D_lin;
    
    // set up EdS growth function calculator with these linear values
    EdS_growth<state_vector::value_type> EdS(D_lin, f_lin);
    
    state_vector::value_type A  = this->use_EdS ? EdS.DA() : x[ELEMENT_A];
    state_vector::value_type B  = this->use_EdS ? EdS.DB() : x[ELEMENT_B];
    state_vector::value_type D  = this->use_EdS ? EdS.DD() : x[ELEMENT_D];
    state_vector::value_type E  = this->use_EdS ? EdS.DE() : x[ELEMENT_E];
    state_vector::value_type F  = this->use_EdS ? EdS.DF() : x[ELEMENT_F];
    state_vector::value_type G  = this->use_EdS ? EdS.DG() : x[ELEMENT_G];
    state_vector::value_type J  = this->use_EdS ? EdS.DJ() : x[ELEMENT_J];

    state_vector::value_type fA = this->use_EdS ? EdS.fA() : - (1.0+z) * x[ELEMENT_dAdz] / (std::fabs(A) > 0.0 ? A : 1.0);
    state_vector::value_type fB = this->use_EdS ? EdS.fB() : - (1.0+z) * x[ELEMENT_dBdz] / (std::fabs(B) > 0.0 ? B : 1.0);
    state_vector::value_type fD = this->use_EdS ? EdS.fD() : - (1.0+z) * x[ELEMENT_dDdz] / (std::fabs(D) > 0.0 ? D : 1.0);
    state_vector::value_type fE = this->use_EdS ? EdS.fE() : - (1.0+z) * x[ELEMENT_dEdz] / (std::fabs(E) > 0.0 ? E : 1.0);
    state_vector::value_type fF = this->use_EdS ? EdS.fF() : - (1.0+z) * x[ELEMENT_dFdz] / (std::fabs(F) > 0.0 ? F : 1.0);
    state_vector::value_type fG = this->use_EdS ? EdS.fG() : - (1.0+z) * x[ELEMENT_dGdz] / (std::fabs(G) > 0.0 ? G : 1.0);
    state_vector::value_type fJ = this->use_EdS ? EdS.fJ() : - (1.0+z) * x[ELEMENT_dJdz] / (std::fabs(J) > 0.0 ? J : 1.0);
    
    this->container.push_back(D_lin, A, B, D, E, F, G, J,
                              f_lin, fA, fB, fD, fE, fF, fG, fJ);
  }
