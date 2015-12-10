//
// Created by David Seery on 17/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "oneloop_growth.h"


oneloop_growth::oneloop_growth(z_database& z)
  : z_db(z)
  {
    g_linear = std::make_unique< std::vector<double> >();
    A = std::make_unique< std::vector<double> >();
    B = std::make_unique< std::vector<double> >();
    D = std::make_unique< std::vector<double> >();
    E = std::make_unique< std::vector<double> >();
    F = std::make_unique< std::vector<double> >();
    G = std::make_unique< std::vector<double> >();
    J = std::make_unique< std::vector<double> >();
  }


void oneloop_growth::push_back(double g, double A, double B, double D, double E, double F, double G, double J)
  {
    this->g_linear->push_back(g);
    this->A->push_back(A);
    this->B->push_back(B);
    this->D->push_back(D);
    this->E->push_back(E);
    this->F->push_back(F);
    this->G->push_back(G);
    this->J->push_back(J);
  }


void oneloop_growth::set_integration_metadata(boost::timer::nanosecond_type t, size_t s)
  {
    this->integration_time = t;
    this->steps = s;
  }
