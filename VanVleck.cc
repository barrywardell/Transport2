/* Numerically solve the transport equation for the Van Vleck determinant
 * along a geodesic.
 *
 * Copyright (C) 2011 Barry Wardell
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
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "SpacetimeTensors.h"

#define NUM_EQS (5+16+1)

#define SQR(x) ((x)*(x))

/* RHS of our system of ODEs */
int rhs (double tau, const double y[], double f[], void *params)
{
  struct geodesic_params p = *((struct geodesic_params*) params);

  /* Avoid dividing by zero initially - use l'Hopital's rule */
  double tauinv;
  if(tau == 0.)
    tauinv = 0;
  else
    tauinv = 1./tau;

  /* Unpack the data into tensor objects */
  Tensor x("^a"), Q("^a_b"), VV(0);
  Tensor u("^a"), dQ("^a_b"), dVV(0);
  double dr = y[4];
  
  x.setComponents(&y[0]);
  Q.setComponents(&y[5]);
  VV.setComponents(&y[21]);
  u.setComponents(&f[0]);
  
  /* Compute metric, Christoffel and Riemann at the current point */
  Schwarzschild s;
  s.t   = x(0);
  s.r   = x(1);
  s.theta = x(2);
  s.phi   = x(3);
  s.M   = 1.0;
  s.calc_all();
  Tensor& G = s.Gudd;
  Tensor& R = s.Ruddd;

  /* Geodesic equations */
  u(0) = s.r/(s.r-2*s.M)*p.e;
  u(1) = dr;
  u(2) = 0.0;
  u(3) = p.l/SQR(s.r);
  f[4] = (SQR(p.l)*(s.r-3*s.M)+s.M*SQR(s.r)*p.type)/SQR(SQR(s.r));
  
  /* Transport equation for Q */
  dQ = (Q["ac"]*G["cbd"] - G["acd"]*Q["cb"])*u["d"]
     - tauinv*(Q["ab"] + Q["ac"]*Q["cb"]) - tau*R["acbd"]*u["c"]*u["d"];

  /* Q(theta,theta) blows up as theta*cot(theta) and makes the numerical
   * integration break down. Since we know the analytic form, don't compute it
   * numerically */
  dQ(2,2) = 0.0;

  /* Transport equation for the Van Vleck determinant */
  dVV = -tauinv*0.5 * VV * Q["aa"];

  /* Pack the data back into an array for GSL */
  u.getComponents(&f[0]);
  dQ.getComponents(&f[5]);
  dVV.getComponents(&f[21]);

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
