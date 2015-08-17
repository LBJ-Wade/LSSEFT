//
// Created by David Seery on 17/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include <utility>

#include "oneloop_integrator.h"
#include "constants.h"

#include "units/eV_units.h"

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
constexpr unsigned int ELEMENT_dgdz = 9;
constexpr unsigned int ELEMENT_dAdz = 10;
constexpr unsigned int ELEMENT_dBdz = 11;
constexpr unsigned int ELEMENT_dDdz = 12;
constexpr unsigned int ELEMENT_dEdz = 13;
constexpr unsigned int ELEMENT_dFdz = 14;
constexpr unsigned int ELEMENT_dGdz = 15;

constexpr unsigned int STATE_SIZE = 16;


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

    //! cache rho_cc in eV
    double rho_cc;

    //! cache H0 in eV
    double H0;

  };


class oneloop_observer
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    oneloop_observer(oneloop& c);

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

    //! reference to oneloop container
    oneloop& container;

    //! integration timer
    boost::timer::cpu_timer timer;

  };


// ONELOOP_INTEGRATOR METHODS


oneloop_integrator::oneloop_integrator(double a, double r)
  : abs_err(a),
    rel_err(r)
  {

  }


oneloop oneloop_integrator::integrate(const FRW_model& model, std::shared_ptr<redshift_database>& z_db)
  {
    // set up empty oneloop container
    oneloop ctr(z_db);

    // set up a functor for the ODE system
    oneloop_functor rhs(model);

    // set up an observer
    oneloop_observer obs(ctr);

    // set up a state vector
    state_vector x(STATE_SIZE);

    // get initial time -- note use of reverse iterator to get last z!
    redshift_database::reverse_value_iterator max_z = z_db->value_rbegin();

    double init_z = *max_z;
    rhs.ics(x, init_z);

    // set up vector of sample times (needed until the Boost version of odeint catchs up to the github version, which includes my patch)
    std::vector<double> z_sample;
    std::copy(z_db->value_rbegin(), z_db->value_rend(), std::back_inserter(z_sample));

    auto stepper = boost::numeric::odeint::make_dense_output< boost::numeric::odeint::runge_kutta_dopri5<state_vector> >(this->abs_err, this->rel_err);

    obs.start_timer();
    size_t steps = boost::numeric::odeint::integrate_times(stepper, rhs, x, z_sample.begin(), z_sample.end(), -1E-3, obs);
    obs.stop_timer();

    ctr.set_integration_metadata(obs.read_timer(), steps);

    return(ctr);
  }


// ONELOOP_FUNCTOR METHODS


oneloop_functor::oneloop_functor(const FRW_model& m)
  : model(m)
  {
    // cache value of H0 and rho_cc
    constexpr double Mp = static_cast<double>(eV_units::PlanckMass);
    eV_units::energy H0_value = model.get_h() * 100.0 * eV_units::Kilometre / eV_units::Second / eV_units::Mpc;

    H0 = static_cast<double>(H0_value);
    rho_cc = 3.0 * H0*H0 * Mp*Mp * model.get_omega_cc();
  }


void oneloop_functor::operator()(const state_vector& x, state_vector& dxdz, double z)
  {
    constexpr double Mp = static_cast<double>(eV_units::PlanckMass);

    double rho = x[RHO_M] + x[RHO_R] + this->rho_cc;
    double H   = std::sqrt(rho / (3.0*Mp*Mp));

    double Hdot    = -(3.0*x[RHO_M] + 4.0*x[RHO_R]) / (6.0*Mp*Mp);
    double epsilon = -Hdot/(H*H);

    double Omega_m   = x[RHO_M]/rho;
    double Omega_r   = x[RHO_R]/rho;

    double one_plus_z = 1.0+z;
    double one_plus_z_sq = (1.0+z)*(1.0+z);

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
  }


void oneloop_functor::ics(state_vector& x, double z)
  {
    constexpr double Mp = static_cast<double>(eV_units::PlanckMass);

    // compute (a0/a)^3 and (a0/a)^4
    double a_three = (1.0+z)*(1.0+z)*(1.0+z);
    double a_four  = (1.0+z)*a_three;

    // compute matter and radiation densities today
    // for radiation, we need the Stefan-Boltzman law and the present day CMB temperature
    double rho_m0 = this->model.get_omega_m() * (3.0*this->H0*this->H0*Mp*Mp);

    eV_units::energy T_CMB = this->model.get_T_CMB();
    double T_CMB_in_eV = static_cast<double>(T_CMB);
    double rho_r0 = g_star * radiation_constant * T_CMB_in_eV*T_CMB_in_eV*T_CMB_in_eV*T_CMB_in_eV;

    double rho_m = rho_m0 * a_three;
    double rho_r = rho_r0 * a_four;

    // initial conditions for background
    x[RHO_M] = rho_m;
    x[RHO_R] = rho_r;

    // initial conditions for linear growth factor
    x[ELEMENT_g] = 1.0;
    x[ELEMENT_dgdz] = 0.0;

    // initial conditions for one-loop kernel
    x[ELEMENT_A] = x[ELEMENT_dAdz] = 0.0;
    x[ELEMENT_B] = x[ELEMENT_dBdz] = 0.0;
    x[ELEMENT_D] = x[ELEMENT_dDdz] = 0.0;
    x[ELEMENT_E] = x[ELEMENT_dEdz] = 0.0;
    x[ELEMENT_F] = x[ELEMENT_dFdz] = 0.0;
    x[ELEMENT_G] = x[ELEMENT_dGdz] = 0.0;
  }


// ONELOOP_OBSERVER METHODS


oneloop_observer::oneloop_observer(oneloop& c)
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
    this->container.push_back(x[ELEMENT_g], x[ELEMENT_A], x[ELEMENT_B], x[ELEMENT_D], x[ELEMENT_E], x[ELEMENT_F], x[ELEMENT_G]);
  }
