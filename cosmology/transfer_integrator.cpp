//
// Created by David Seery on 15/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//


#include <utility>

#include "transfer_integrator.h"
#include "constants.h"

#include "utilities/formatter.h"

#include "boost/numeric/odeint.hpp"


typedef std::vector<double> state_vector;

constexpr unsigned int RHO_M   = 0;
constexpr unsigned int RHO_R   = 1;
constexpr unsigned int DELTA_M = 2;
constexpr unsigned int DELTA_R = 3;
constexpr unsigned int THETA_M = 4;
constexpr unsigned int THETA_R = 5;
constexpr unsigned int PHI     = 6;

constexpr unsigned int STATE_SIZE = 7;


class transfer_functor
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    transfer_functor(const FRW_model& m, const eV_units::energy& _k);

    //! destructor is default
    ~transfer_functor() = default;


    // INTERFACE

  public:

    //! compute RHS of ODE system
    void operator()(const state_vector& x, state_vector& dxdz, double z);

    //! compute ics for ODE system
    void ics(state_vector& x, double z);

    //! estimate a suitable initial redshift for the integration
    double find_init_z();


    // INTERNAL DATA

  private:

    //! reference to FRW model
    const FRW_model& model;

    //! reference to wavenumber object
    const eV_units::energy& k;

    //! cache rho_cc in eV
    double rho_cc;

    //! cache H0 in eV
    double H0;

  };


class transfer_observer
  {

    // CONSTRUCTOR, DESTRUCTOR

  public:

    //! constructor
    transfer_observer(transfer_function& c);

    //! destructor is default
    ~transfer_observer() = default;


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

    //! reference to transfer_function container
    transfer_function& container;

    //! keep track of first invocation -- we don't want to store anything at the initial time
    bool first_call;

    //! integration timer
    boost::timer::cpu_timer timer;

  };


// TRANSFER_FUNCTOR METHODS


transfer_functor::transfer_functor(const FRW_model& m, const eV_units::energy& _k)
  : model(m),
    k(_k)
  {
    // cache value of H0 and rho_cc
    constexpr double Mp = static_cast<double>(eV_units::PlanckMass);
    eV_units::energy H0_value = model.get_h() * 100.0 * eV_units::Kilometre / eV_units::Second / eV_units::Mpc;

    H0 = static_cast<double>(H0);
    rho_cc = 3.0 * H0*H0 * Mp*Mp * model.get_omega_cc();
  }


void transfer_functor::operator()(const state_vector& x, state_vector& dxdz, double z)
  {
    constexpr double Mp = static_cast<double>(eV_units::PlanckMass);

    double rho = x[RHO_M] + x[RHO_R] + this->rho_cc;
    double H   = std::sqrt(rho / (3.0*Mp*Mp));

    double Hdot    = -(3.0*x[RHO_M] + 4.0*x[RHO_R]) / (6.0*Mp*Mp);
    double epsilon = -Hdot/(H*H);

    double Omega_m   = x[RHO_M]/rho;
    double Omega_r   = x[RHO_R]/rho;

    // k measured in eV here
    double aH                = H/(1.0+z);
    double k_over_aH         = static_cast<double>(this->k) / aH;
    double k_over_aH_squared = k_over_aH * k_over_aH;

    // evolve background
    dxdz[RHO_M] = 3.0*x[RHO_M] / (1.0+z);
    dxdz[RHO_R] = 4.0*x[RHO_R] / (1.0+z);

    // TRANSFER FUNCTIONS
    // note: Phi and delta are individually dimensionless, so the Phi and delta
    // transfer functions are dimensionless also
    // the velocity v is dimensionless, so theta = div v has dimensions of energy
    // then if we choose to measure theta in Hubble units, corresponding to
    // tilde{theta} = theta/H, then tilde{theta} is again dimensionless.
    // What we actually evolve is tilde{theta}.
    // Then, the conventional 'growth factor' is f = -theta/delta.

    // evolve Phi transfer function
    dxdz[PHI] = (1.0/3.0)*k_over_aH_squared*x[PHI]/(1.0+z)
                - (1.0/2.0)*(Omega_m*x[DELTA_M] + Omega_r*x[DELTA_R]) / (1.0+z)
                + (Omega_m + Omega_r)*x[PHI]/(1.0+z);

    // evolve velocity transfer functions
    dxdz[THETA_M] = ((2.0-epsilon)*x[THETA_M] + k_over_aH_squared*x[PHI]) / (1.0+z);
    dxdz[THETA_R] = ((1.0-epsilon)*x[THETA_R] + k_over_aH_squared*x[PHI] - (1.0/4.0)*k_over_aH_squared*x[DELTA_R]) / (1.0+z);

    // evolve density transfer functions
    dxdz[DELTA_M] = x[THETA_M]/(1.0+z) - 3.0*dxdz[PHI];
    dxdz[DELTA_R] = (4.0/3.0)*x[THETA_R]/(1.0+z) - 4.0*dxdz[PHI];
  }


void transfer_functor::ics(state_vector& x, double z)
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

    double rho = rho_m + rho_r + this->rho_cc;
    double H   = std::sqrt(rho/(3.0*Mp*Mp));

    double aH = H / (1.0 + z);
    double k_over_aH = static_cast<double>(this->k) / aH;

    // initial conditions for background objects
    x[RHO_M] = rho_m;
    x[RHO_R] = rho_r;

    // initial conditions for transfer functions
    x[DELTA_M] = 3.0 / 2.0;
    x[DELTA_R] = 2.0;
    x[THETA_M] = -k_over_aH*k_over_aH/2.0;
    x[THETA_R] = -k_over_aH*k_over_aH/2.0;
    x[PHI]     = 1.0;
  }


double transfer_functor::find_init_z()
  {
    // try to get the k-mode 10 e-folds outside the horizon
    // superhorizon_factor = exp(10)
    constexpr double superhorizon_factor = 22026.4657948;

    double Mp_over_T_CMB = eV_units::PlanckMass / this->model.get_T_CMB();
    double k_over_T_CMB  = this->k / this->model.get_T_CMB();

    double z_init = std::pow((3.0/(g_star*radiation_constant)) * Mp_over_T_CMB * k_over_T_CMB * superhorizon_factor, 1.0/3.0) - 1.0;

    if(z_init < z_eq) z_init = z_eq * superhorizon_factor;

    return(z_init);
  }


// TRANSFER_OBSERVER METHODS


transfer_observer::transfer_observer(transfer_function& c)
  : container(c),
    first_call(true)
  {
    // stop timer
    timer.stop();
  }


void transfer_observer::start_timer()
  {
    this->timer.start();
  }


void transfer_observer::stop_timer()
  {
    this->timer.stop();
  }


boost::timer::nanosecond_type transfer_observer::read_timer()
  {
    return(this->timer.elapsed().wall);
  }


void transfer_observer::operator()(const state_vector& x, double z)
  {
    if(this->first_call)
      {
        this->first_call = false;
      }
    else
      {
        this->container.push_back(x[DELTA_M], x[DELTA_R], x[THETA_M], x[THETA_R], x[PHI]);
//        std::cout << "z = " << z << ", delta_m = " << x[DELTA_M] << ", delta_r = " << x[DELTA_R] << ", theta_m = " << x[THETA_M] << ", theta_r = " << x[THETA_R] << ", Phi = " << x[PHI] << '\n';
      }
  }


// TRANSFER_INTEGRATOR METHODS


transfer_integrator::transfer_integrator(double a, double r)
  : abs_err(a),
    rel_err(r)
  {
  }


transfer_function transfer_integrator::integrate(const FRW_model& model, const eV_units::energy& k,
                                                 const wavenumber_token& tok, std::shared_ptr<redshift_database>& z_db)
  {
    // set up an empty transfer_function container
    transfer_function ctr(k, tok, z_db);

    // set up a functor for the ODE system
    transfer_functor rhs(model, k);

    // set up an observer
    transfer_observer obs(ctr);

    // set up a state vector
    state_vector x(STATE_SIZE);

    // find initial time for integation; typically guessed by asking that the k-mode
    // is sufficiently superhorizon
    double init_z = rhs.find_init_z();

    // set up initial conditions
    rhs.ics(x, init_z);

    // set up vector of sample times
    std::vector<double> z_sample{ init_z };
    std::copy(z_db->value_rbegin(), z_db->value_rend(), std::back_inserter(z_sample));

    auto stepper = boost::numeric::odeint::make_dense_output< boost::numeric::odeint::runge_kutta_dopri5<state_vector> >(this->abs_err, this->rel_err);

    obs.start_timer();
    size_t steps = boost::numeric::odeint::integrate_times(stepper, rhs, x, z_sample.begin(), z_sample.end(), -1E-3, obs);
    obs.stop_timer();

    ctr.set_integration_metadata(obs.read_timer(), steps);

    return(ctr);
  }
