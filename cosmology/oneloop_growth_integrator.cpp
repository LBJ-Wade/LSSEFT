//
// Created by David Seery on 17/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include <utility>

#include "oneloop_growth_integrator.h"
#include "constants.h"

#include "units/Mpc_units.h"

#include "boost/numeric/odeint.hpp"


typedef std::vector<double> state_vector;

constexpr unsigned int RHO_M        = 0;
constexpr unsigned int RHO_R        = 1;
constexpr unsigned int ELEMENT_g    = 2;
constexpr unsigned int ELEMENT_A    = 3;
constexpr unsigned int ELEMENT_B    = 4;
constexpr unsigned int ELEMENT_D    = 5;
constexpr unsigned int ELEMENT_E    = 6;
constexpr unsigned int ELEMENT_F    = 7;
constexpr unsigned int ELEMENT_G    = 8;
constexpr unsigned int ELEMENT_J    = 9;
constexpr unsigned int ELEMENT_dgdz = 10;
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
    oneloop_functor(const FRW_model& m);

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
    oneloop_observer(oneloop_growth& c);

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

    //! integration timer
    boost::timer::cpu_timer timer;

  };


// ONELOOP_INTEGRATOR METHODS


oneloop_growth_integrator::oneloop_growth_integrator(double a, double r)
  : abs_err(a),
    rel_err(r)
  {

  }


growth_integrator_data oneloop_growth_integrator::integrate(const FRW_model& model, z_database& z_db)
  {
    // set up empty oneloop_growth container
    std::unique_ptr<oneloop_growth> ctr = std::make_unique<oneloop_growth>(z_db);

    // set up a functor for the ODE system
    oneloop_functor rhs(model);

    // set up an observer
    oneloop_observer obs(*ctr);

    // set up a state vector
    state_vector x(STATE_SIZE);

    // get initial time -- note use of reverse iterator to get last z!
    z_database::reverse_value_iterator max_z = z_db.value_rbegin();

    double init_z = *max_z;
    rhs.ics(x, init_z);

    // set up vector of sample times (needed until the Boost version of odeint catchs up to the github version, which includes my patch)
    std::vector<double> z_sample;
    std::copy(z_db.value_rbegin(), z_db.value_rend(), std::back_inserter(z_sample));

    auto stepper = boost::numeric::odeint::make_dense_output< boost::numeric::odeint::runge_kutta_dopri5<state_vector> >(this->abs_err, this->rel_err);

    obs.start_timer();
    size_t steps = boost::numeric::odeint::integrate_times(stepper, rhs, x, z_sample.begin(), z_sample.end(), -1E-3, obs);
    obs.stop_timer();
    
    return growth_integrator_data(std::move(ctr), obs.read_timer(), steps);
  }


// ONELOOP_FUNCTOR METHODS


oneloop_functor::oneloop_functor(const FRW_model& m)
  : model(m),
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
    dxdz[ELEMENT_g]    = x[ELEMENT_dgdz];
    dxdz[ELEMENT_dgdz] = (1.0 - epsilon) * x[ELEMENT_dgdz] / one_plus_z
                         + (3.0 * Omega_m / 2.0) * x[ELEMENT_g] / one_plus_z_sq;

    // evolve one-loop kernels -- derivatives
    dxdz[ELEMENT_A] = x[ELEMENT_dAdz];
    dxdz[ELEMENT_B] = x[ELEMENT_dBdz];
    dxdz[ELEMENT_D] = x[ELEMENT_dDdz];
    dxdz[ELEMENT_E] = x[ELEMENT_dEdz];
    dxdz[ELEMENT_F] = x[ELEMENT_dFdz];
    dxdz[ELEMENT_G] = x[ELEMENT_dGdz];
    dxdz[ELEMENT_J] = x[ELEMENT_dJdz];

    // evolve one-loop kernels -- second derivatives
    dxdz[ELEMENT_dAdz] = (3.0 * Omega_m / 2.0) * x[ELEMENT_g] * x[ELEMENT_g] / one_plus_z_sq
                         + (1.0 - epsilon) * x[ELEMENT_dAdz] / one_plus_z
                         + (3.0 * Omega_m / 2.0) * x[ELEMENT_A] / one_plus_z_sq;

    dxdz[ELEMENT_dBdz] = x[ELEMENT_dgdz] * x[ELEMENT_dgdz]
                         + (1.0 - epsilon) * x[ELEMENT_dBdz] / one_plus_z
                         + (3.0 * Omega_m / 2.0) * x[ELEMENT_B] / one_plus_z_sq;

    dxdz[ELEMENT_dDdz] = x[ELEMENT_dgdz] * x[ELEMENT_dAdz]
                         + (1.0 - epsilon) * x[ELEMENT_dDdz] / one_plus_z
                         + (3.0 * Omega_m / 2.0) * x[ELEMENT_D] / one_plus_z_sq;

    dxdz[ELEMENT_dEdz] = x[ELEMENT_dgdz] * x[ELEMENT_dBdz]
                         + (1.0 - epsilon) * x[ELEMENT_dEdz] / one_plus_z
                         + (3.0 * Omega_m / 2.0) * x[ELEMENT_E] / one_plus_z_sq;

    dxdz[ELEMENT_dFdz] = (3.0 * Omega_m / 2.0) * x[ELEMENT_g] * x[ELEMENT_A] / one_plus_z_sq
                         + (1.0 - epsilon) * x[ELEMENT_dFdz] / one_plus_z
                         + (3.0 * Omega_m / 2.0) * x[ELEMENT_F] / one_plus_z_sq;

    dxdz[ELEMENT_dGdz] = (3.0 * Omega_m / 2.0) * x[ELEMENT_g] * x[ELEMENT_B] / one_plus_z_sq
                         + (1.0 - epsilon) * x[ELEMENT_dGdz] / one_plus_z
                         + (3.0 * Omega_m / 2.0) * x[ELEMENT_G] / one_plus_z_sq;

    dxdz[ELEMENT_dJdz] = x[ELEMENT_dgdz] * x[ELEMENT_dgdz] * x[ELEMENT_g]
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
    // to get the derviative, assume that we are early in matter domination and g(z) grows like a(z).
    // so g(z) = a(z)/a_init(z). Then dg/dz = -da/dz / a_init = (a_0/a_init) / (1+z)^2 = 1/(1+z).
    x[ELEMENT_g] = 1.0;
    x[ELEMENT_dgdz] = -1.0 / (1.0+z);

    // initial conditions for one-loop kernel
    x[ELEMENT_A] = x[ELEMENT_dAdz] = 0.0;
    x[ELEMENT_B] = x[ELEMENT_dBdz] = 0.0;
    x[ELEMENT_D] = x[ELEMENT_dDdz] = 0.0;
    x[ELEMENT_E] = x[ELEMENT_dEdz] = 0.0;
    x[ELEMENT_F] = x[ELEMENT_dFdz] = 0.0;
    x[ELEMENT_G] = x[ELEMENT_dGdz] = 0.0;
    x[ELEMENT_J] = x[ELEMENT_dJdz] = 0.0;
  }


// ONELOOP_OBSERVER METHODS


oneloop_observer::oneloop_observer(oneloop_growth& c)
  : container(c)
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
    state_vector::value_type g = x[ELEMENT_g];
    state_vector::value_type A = x[ELEMENT_A];
    state_vector::value_type B = x[ELEMENT_B];
    state_vector::value_type D = x[ELEMENT_D];
    state_vector::value_type E = x[ELEMENT_E];
    state_vector::value_type F = x[ELEMENT_F];
    state_vector::value_type G = x[ELEMENT_G];
    state_vector::value_type J = x[ELEMENT_J];
    
    state_vector::value_type f = - (1.0+z) * x[ELEMENT_dgdz] / g;
    state_vector::value_type fA = - (1.0+z) * x[ELEMENT_dAdz] / (std::fabs(A) > 0.0 ? A : 1.0);
    state_vector::value_type fB = - (1.0+z) * x[ELEMENT_dBdz] / (std::fabs(B) > 0.0 ? B : 1.0);
    state_vector::value_type fD = - (1.0+z) * x[ELEMENT_dDdz] / (std::fabs(D) > 0.0 ? D : 1.0);
    state_vector::value_type fE = - (1.0+z) * x[ELEMENT_dEdz] / (std::fabs(E) > 0.0 ? E : 1.0);
    state_vector::value_type fF = - (1.0+z) * x[ELEMENT_dFdz] / (std::fabs(F) > 0.0 ? F : 1.0);
    state_vector::value_type fG = - (1.0+z) * x[ELEMENT_dGdz] / (std::fabs(G) > 0.0 ? G : 1.0);
    state_vector::value_type fJ = - (1.0+z) * x[ELEMENT_dJdz] / (std::fabs(J) > 0.0 ? J : 1.0);
    
    this->container.push_back(g, A, B, D, E, F, G, J,
                              f, fA, fB, fD, fE, fF, fG, fJ);
  }
