//
// Created by David Seery on 17/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "oneloop_growth.h"


oneloop_growth::oneloop_growth(z_database& z)
  : z_db(std::make_unique<z_database>(z))
  {
    g_linear = std::make_unique< std::vector<double> >();
    A = std::make_unique< std::vector<double> >();
    B = std::make_unique< std::vector<double> >();
    D = std::make_unique< std::vector<double> >();
    E = std::make_unique< std::vector<double> >();
    F = std::make_unique< std::vector<double> >();
    G = std::make_unique< std::vector<double> >();
    J = std::make_unique< std::vector<double> >();
    
    f_linear = std::make_unique< std::vector<double> >();
    fA = std::make_unique< std::vector<double> >();
    fB = std::make_unique< std::vector<double> >();
    fD = std::make_unique< std::vector<double> >();
    fE = std::make_unique< std::vector<double> >();
    fF = std::make_unique< std::vector<double> >();
    fG = std::make_unique< std::vector<double> >();
    fJ = std::make_unique< std::vector<double> >();
  }


oneloop_growth::oneloop_growth()
  {
  }


void oneloop_growth::push_back(double g, double A, double B, double D, double E, double F, double G, double J,
                               double f, double fA, double fB, double fD, double fE, double fF, double fG, double fJ)
  {
    this->g_linear->push_back(g);
    this->A->push_back(A);
    this->B->push_back(B);
    this->D->push_back(D);
    this->E->push_back(E);
    this->F->push_back(F);
    this->G->push_back(G);
    this->J->push_back(J);
    
    this->f_linear->push_back(f);
    this->fA->push_back(fA);
    this->fB->push_back(fB);
    this->fD->push_back(fD);
    this->fE->push_back(fE);
    this->fF->push_back(fF);
    this->fG->push_back(fG);
    this->fJ->push_back(fJ);
  }


void oneloop_growth::set_integration_metadata(boost::timer::nanosecond_type t, size_t s)
  {
    this->integration_time = t;
    this->steps = s;
  }


oneloop_growth::oneloop_growth(oneloop_growth&& obj)
  : z_db(std::make_unique<z_database>(*obj.z_db)),
    g_linear(std::move(obj.g_linear)),
    A(std::move(obj.A)),
    B(std::move(obj.B)),
    D(std::move(obj.D)),
    E(std::move(obj.E)),
    F(std::move(obj.F)),
    G(std::move(obj.G)),
    J(std::move(obj.J)),
    f_linear(std::move(obj.f_linear)),
    fA(std::move(obj.fA)),
    fB(std::move(obj.fB)),
    fD(std::move(obj.fD)),
    fE(std::move(obj.fE)),
    fF(std::move(obj.fF)),
    fG(std::move(obj.fG)),
    fJ(std::move(obj.fJ)),
    integration_time(obj.integration_time),
    steps(obj.steps)
  {
    obj.g_linear = std::make_unique< std::vector<double> >();
    obj.A = std::make_unique< std::vector<double> >();
    obj.B = std::make_unique< std::vector<double> >();
    obj.D = std::make_unique< std::vector<double> >();
    obj.E = std::make_unique< std::vector<double> >();
    obj.F = std::make_unique< std::vector<double> >();
    obj.G = std::make_unique< std::vector<double> >();
    obj.J = std::make_unique< std::vector<double> >();
    
    obj.f_linear = std::make_unique< std::vector<double> >();
    obj.fA = std::make_unique< std::vector<double> >();
    obj.fB = std::make_unique< std::vector<double> >();
    obj.fD = std::make_unique< std::vector<double> >();
    obj.fE = std::make_unique< std::vector<double> >();
    obj.fF = std::make_unique< std::vector<double> >();
    obj.fG = std::make_unique< std::vector<double> >();
    obj.fJ = std::make_unique< std::vector<double> >();
  }
