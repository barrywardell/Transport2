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

#define SQR(x) ((x)*(x))

/* RHS of our system of ODEs */
int rhs (double tau, const double y[], double f[], void *params)
{
  struct geodesic_params p = *((struct geodesic_params*) params);
  const double M = p.s.M;

  /* Avoid dividing by zero initially - use l'Hopital's rule */
  double tauinv;
  if(tau == 0.)
    tauinv = 0;
  else
    tauinv = 1./tau;

  /* Evolved tensors */
  TensorList T(p.T);
  T.setComponents(y);

  Tensor& x  = T["x"];
  Tensor& ur = T["ur"];
  Tensor& Q  = T["Q"];
  Tensor& VV = T["VV"];

  /* Right hand side tensors */
  TensorList dT(p.T);
  dT.setComponents(f);

  Tensor& dx  = dT["x"];
  Tensor& dur = dT["ur"];
  Tensor& dQ  = dT["Q"];
  Tensor& dVV = dT["VV"];

  /* Compute metric, Christoffel and Riemann at the current point */
  Schwarzschild& s = p.s;
  s.setPoint(x);
  s.calc_all();
  Tensor& G = s["Gudd"];
  Tensor& R = s["Ruddd"];

  /* Geodesic equations */
  const double& r = x(1);
  dx(0) = r/(r-2*M)*p.e;
  dx(1) = ur();
  dx(2) = 0.0;
  dx(3) = p.l/SQR(r);
  dur() = (SQR(p.l)*(r-3*M)+M*SQR(r)*p.type)/SQR(SQR(r));
  
  /* Transport equation for Q */
  dQ = (Q["ac"]*G["cbd"] - G["acd"]*Q["cb"])*dx["d"]
     - tauinv*(Q["ab"] + Q["ac"]*Q["cb"]) - tau*R["acbd"]*dx["c"]*dx["d"];

  /* Q(theta,theta) blows up as theta*cot(theta) and makes the numerical
   * integration break down. Since we know the analytic form, don't compute it
   * numerically */
  dQ(2,2) = 0.0;

  /* Transport equation for the Van Vleck determinant */
  dVV = -tauinv*0.5 * VV * Q["aa"];

  /* Pack the data back into an array for GSL */
  dT.getComponents(f);

  return GSL_SUCCESS;
}

int main (int argc, char * argv[])
{
  /* Evolved tensors */
  TensorList T;
  T.append("x", "^a");   Tensor& x  = T["x"];  /* Position */
  T.append("ur");        Tensor& ur = T["ur"]; /* Radial 4-velocity */
  T.append("Q", "^a_b"); Tensor& Q  = T["Q"];  /* \sigma^a'_b' - \delta^a'_b' */
  T.append("VV");        Tensor& VV = T["VV"]; /* Van Vleck determinant */
  const int numEqs = T.getNumComponents();

  /* Spacetime */
  Schwarzschild schw(1.0);

  /* Use a Runge-Kutta integrator with adaptive step-size */
  const gsl_odeiv_step_type * t = gsl_odeiv_step_rkf45;
  gsl_odeiv_step * s = gsl_odeiv_step_alloc (t, numEqs);
  gsl_odeiv_control * c = gsl_odeiv_control_y_new (1e-10, 1e-10);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (numEqs);

  /* Time-like geodesic starting at r=10M and going in to r=4M */
  struct geodesic_params params =
    {0.94868329805051379960, 3.5355339059327376220, Timelike, T, schw};

  gsl_odeiv_system sys = {rhs, NULL, numEqs, &params};

  /* Initial Conditions */
  double tau = 0.0, tau1 = 100.0;
  double h = 1e-6;
  double r0 = 10.0;

  x(1) = r0;
  x(2) = M_PI_2;
  VV() = 1.0;

  /* Convert tensors to a flat C array which GSL understands */
  double y[numEqs];
  T.getComponents(y);

  while (tau < tau1)
  {
    int status = gsl_odeiv_evolve_apply (e, c, s, &sys, &tau, tau1, &h, y);

    if (status != GSL_SUCCESS)
      break;

    /* Output the results */
    T.setComponents(y);
    printf ("%.5f", tau);
    for(int i=0; i<numEqs; i++)
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
