/* Numerically solve the transport equation for the Van Vleck determinant
 * along a geodesic.
 *
 * Copyright (C) 2009 Barry Wardell
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This software is distributed on an "AS IS" basis, WITHOUT WARRANTY OF ANY
 * KIND, either express or implied.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "SpacetimeTensors.h"

#define DIMENSION 4
#define NUM_EQS (5+16+1)

enum GeodesicType {
  Timelike  = -1,
  Null      = 0,
  Spacelike = 1
};

void setComponents(Tensor &t, const double* v)
{
  double* components = t.getComponents();
  for (int i = 0; i < pow(DIMENSION, t.getRank()); i++) {
    components[i] = v[i];
  }
}

void getComponents(Tensor &t, double* v)
{
  double* components = t.getComponents();
  for (int i = 0; i < pow(DIMENSION, t.getRank()); i++) {
    v[i] = components[i];
  }
}

/* RHS of our system of ODEs */
int rhs (double tau, const double y[], double f[], void *params)
{
  struct geodesic_params p = *((struct geodesic_params*) params);
  Tensor::IndexType up = Tensor::UP;
  Tensor::IndexType down = Tensor::DOWN;
  double tauinv;

  if(tau == 0)
    tauinv = 0;
  else
    tauinv = 1./tau;

  Schwarzschild s;
  
  double r = y[1];
  s.r = y[1];
  s.theta = y[2];
  s.M = 1.0;
  double rp = y[4];
  s.calc_all();
  
  Tensor  u(1, up);
  Tensor du(1, up);
  Tensor  Q(2, up, down);
  Tensor VV(0);
  
  setComponents(u,  &y[0]);
  setComponents(Q,  &y[5]);
  setComponents(VV, &y[21]);
  
  du.get(0) = r/(r-2*s.M)*p.e;
  du.get(1) = rp;
  du.get(2) = 0.0;
  du.get(3) = p.l/gsl_pow_2(r);
  
  f[4] = (gsl_pow_2(p.l)*(r-3*s.M)+s.M*gsl_pow_2(r)*p.type)/gsl_pow_4(r);
  
  Tensor G = *s.Gudd;
  Tensor R = *s.Ruddd;
  Tensor dQ = (Q('a','d')*G('d','c','b') - G('a','c','d')*Q('c','b'))*du('c')
    - tauinv*(Q('a','b') + Q('a','c')*Q('c','b')) - tau*R('a','c','b','d')*du('c')*du('d');
  Tensor dVV = -tauinv*0.5 * VV * Q('a', 'a');

  getComponents(du, &f[0]);
  getComponents(dQ, &f[5]);
  getComponents(dVV,&f[21]);

  return GSL_SUCCESS;
}

int main (int argc, char * argv[])
{
  int i;

  /* Use a Runge-Kutta integrator with adaptive step-size */
  const gsl_odeiv_step_type * T = gsl_odeiv_step_rkf45;
  gsl_odeiv_step * s = gsl_odeiv_step_alloc (T, NUM_EQS);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new (1e-10, 1e-10);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (NUM_EQS);


  /* Time-like geodesic starting at r=10M and going in to r=4M */
  struct geodesic_params params = {0.94868329805051379960,3.5355339059327376220,Timelike};

  gsl_odeiv_system sys = {rhs, NULL, NUM_EQS, &params};

  double tau = 0.0, tau1 = 100.0;
  double h = 1e-6;
  double r0 = 10.0;

  /* Initial Condidions */
  double y[NUM_EQS] = {
    0, r0, M_PI_2, 0.0, 0.0, /* t, r, theta, phi, r' */
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, /* Q^a'_b' */
    1.0 /* Delta^1/2 */
  };

  while (tau < tau1)
  {
    int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &tau, tau1, &h, y);

    if (status != GSL_SUCCESS)
      break;

    /* Output the results */
    printf ("%.5f", tau);
    for(i=0; i<NUM_EQS; i++)
    {
      printf (", %.5f", y[i]);
    }
    printf("\n");

    /* Exit if step size get smaller than 10^-12 */
    if (h < 1e-13)
    {
      fprintf(stderr,"Error: step size %e less than 1e-13 is not allowed.\n",h);
      break;
    }
  }

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (s);

  return 0;
}
