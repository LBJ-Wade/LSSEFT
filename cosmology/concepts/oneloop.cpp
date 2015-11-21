//
// Created by David Seery on 17/08/2015.
// Copyright (c) 2015 University of Sussex. All rights reserved.
//

#include "oneloop.h"


oneloop::oneloop(z_database& z)
  : z_db(z)
  {
    g_linear.reset(new std::vector<double>);
    A.reset(new std::vector<double>);
    B.reset(new std::vector<double>);
    D.reset(new std::vector<double>);
    E.reset(new std::vector<double>);
    F.reset(new std::vector<double>);
    G.reset(new std::vector<double>);
  }


void oneloop::push_back(double g, double A, double B, double D, double E, double F, double G)
  {
    this->g_linear->push_back(g);
    this->A->push_back(A);
    this->B->push_back(B);
    this->D->push_back(D);
    this->E->push_back(E);
    this->F->push_back(F);
    this->G->push_back(G);
  }


void oneloop::set_integration_metadata(boost::timer::nanosecond_type t, size_t s)
  {
    this->integration_time = t;
    this->steps = s;
  }
