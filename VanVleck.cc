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
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "SpacetimeTensors.h"

#define SQR(x) ((x)*(x))

/* RHS of our system of ODEs */
int rhs (double tau, const double y[], double f[], void *params)
{
  struct geodesic_params p = *((struct geodesic_params*) params);
  const double M = p.s.M;

  /* Evolved tensors */
  TensorList T(p.T);
  T.setComponents(y);

  Tensor &x = T["x"], &ur = T["ur"], &Q = T["Q"], &VV = T["VV"], &eta = T["eta"],
         &s12 = T["s12"], &s3 = T["s3"], &s13 = T["s13"], &s4 = T["s4"],
         &I = T["I"], &I1 = T["I1"], &I2 = T["I2"], &V0 = T["V0"];

  /* Right hand side tensors */
  TensorList dT(p.T);
  dT.setComponents(f);

  Tensor &dx = dT["x"], &dur = dT["ur"], &dQ = dT["Q"], &dVV = dT["VV"],
         &deta = dT["eta"], &ds12 = dT["s12"], &ds3 = dT["s3"],
         &ds13 = dT["s13"], &ds4 = dT["s4"],
         &dI = dT["I"], &dI1 = dT["I1"], &dI2 = dT["I2"], &dV0 = dT["V0"];

  /* Compute metric, Christoffel and Riemann at the current point */
  Schwarzschild& s = p.s;
  s.setPoint(x);
  s.calc_all();
  Tensor &G = s["Gudd"], &R = s["Ruddd"], &dR = s["Ruddd;d"],
         &ddR = s["Ruddd;dd"], &gu = s["guu"];

  /* Geodesic equations */
  const double& r = x(1);
  dx(0) = r/(r-2*M)*p.e;
  dx(1) = ur();
  dx(2) = 0.0;
  dx(3) = p.l/SQR(r);
  dur() = (SQR(p.l)*(r-3*M)+M*SQR(r)*p.type)/SQR(SQR(r));
  
  /* I^-1 */
  Tensor Iinv("_a^b");

  gsl_matrix * Imat = gsl_matrix_calloc(4,4);
  I.getComponents(Imat->data);

  int signum=0;
  gsl_permutation * perm = gsl_permutation_calloc (4);
  gsl_matrix * I_inv = gsl_matrix_calloc(4,4);
  gsl_matrix * lu = gsl_matrix_calloc(4,4);
  gsl_matrix_memcpy(lu, Imat);
  gsl_linalg_LU_decomp (lu, perm, &signum);
  gsl_linalg_LU_invert (lu, perm, I_inv);

  Iinv.setComponents(I_inv->data);
  gsl_matrix_free(I_inv);
  gsl_matrix_free(lu);
  gsl_matrix_free(Imat);
  gsl_permutation_free(perm);

  /* gamma */
  Tensor gamma("^a_b");

  gsl_matrix * etaMat = gsl_matrix_calloc(4,4);
  eta.getComponents(etaMat->data);

  int signum2=0;
  gsl_permutation * perm2 = gsl_permutation_calloc (4);
  gsl_matrix * gamma_mat = gsl_matrix_calloc(4,4);
  gsl_matrix * lu2 = gsl_matrix_calloc(4,4);
  gsl_matrix_memcpy(lu2, etaMat);
  gsl_linalg_LU_decomp (lu2, perm2, &signum2);
  gsl_linalg_LU_invert (lu2, perm2, gamma_mat);

  gamma.setComponents(gamma_mat->data);
  gsl_matrix_free(gamma_mat);
  gsl_matrix_free(lu2);
  gsl_matrix_free(etaMat);
  gsl_permutation_free(perm2);

  /* Transport equations, avoid dividing by zero initially - l'Hopital's rule */
  if (tau == 0.) {
    /* Transport equation for Q = xi-1 */
    dQ = (Q["ac"]*G["cbd"] - G["acd"]*Q["cb"])*dx["d"]
       - tau*R["acbd"]*dx["c"]*dx["d"];

    /* Q(theta,theta) blows up as theta*cot(theta) and makes the numerical
     * integration break down. Since we know the analytic form, don't compute it
     * numerically */
    //dQ(2,2) = 0.0;

    Tensor kd("^a_b");
    kd(0,0) = 1; kd(1,1) = 1; kd(2,2) = 1; kd(3,3) = 1;
    Tensor xi("^a_b");
    xi = Q["ab"] + kd["ab"];

    /* Transport equation for the Van Vleck determinant */
    dVV = 0;

    /* Transport equation for eta = sigma^a_b' */
    deta = eta["ac"]*G["cbd"]*dx["d"];

    /* Transport equation for sigma^a_{b'c'} */
    ds12 = s12["aec"]*G["ebd"]*dx["d"] + s12["abe"]*G["ecd"]*dx["d"]
          - 0.5*R["adbc"]*dx["d"] - 1./6.*(R["abdc"] + R["acdb"])*dx["d"]
          - eta["ad"]*R["dbec"]*dx["e"];

    /* Transport equation for sigma^a'_{b'c'} */
    ds3 = - 1./3.*(R["abdc"]+R["acdb"])*dx["d"]
         - xi["ad"]*R["dbec"]*dx["e"]
         + xi["db"]*R["adec"]*dx["e"] + xi["dc"]*R["adeb"]*dx["e"]
         - tau*dR["adbec"]*dx["d"]*dx["e"]
         + s3["aec"]*G["ebd"]*dx["d"] + s3["abe"]*G["ecd"]*dx["d"]
         - s3["ebc"]*G["aed"]*dx["d"];

    /* Transport equation for sigma^a'_{b'c' d'} */
    ds4 = 0*s4["abcd"]
        - tau*(ddR["aebfcd"]*dx["e"]*dx["f"]
        - R["afed"]*R["ebgc"]*dx["f"]*dx["g"]
        - R["afec"]*R["ebgd"]*dx["f"]*dx["g"])

        + s4["aecd"]*G["ebf"]*dx["f"] + s4["abed"]*G["ecf"]*dx["f"]
        + s4["abce"]*G["edf"]*dx["f"] - s4["ebcd"]*G["aef"]*dx["f"]

        - dR["fbecd"]*dx["e"]*xi["af"] + dR["afecd"]*dx["e"]*xi["fb"]
        + dR["afebd"]*dx["e"]*xi["fc"] + dR["afebc"]*dx["e"]*xi["fd"]

        - R["ecfd"]*dx["f"]*s3["abe"]
        - R["ebfd"]*dx["f"]*s3["ace"]
        - R["ebfc"]*dx["f"]*s3["ade"]
        + R["agfc"]*dx["f"]*s3["gbd"]
        + R["agfb"]*dx["f"]*s3["gcd"]
        + R["agfd"]*dx["f"]*s3["gbc"]

        - 1./4.*(dR["abecd"]+dR["abedc"]+dR["acebd"]+dR["acedb"]
            +dR["adebc"]+dR["adecb"])*dx["e"];

    ds13 = 0*s13["abcd"]
        + s13["aecd"]*G["ebf"]*dx["f"] + s13["abed"]*G["ecf"]*dx["f"]
        + s13["abce"]*G["edf"]*dx["f"]

        - eta["ae"]*dR["ebfcd"]*dx["f"]

        - s12["ade"]*R["ebfc"]*dx["f"]
        - s12["ace"]*R["ebfd"]*dx["f"]
        - s12["abe"]*R["ecfd"]*dx["f"]

        + 1./4.*(dR["acedb"] + dR["adecb"])*dx["e"]
        - 5./12.*(dR["abecd"] + dR["abedc"])*dx["e"]
        + 7./12.*(dR["adebc"] + dR["acebd"])*dx["e"];

    dI = 0*I["ab"] + I["cb"]*G["cad"]*dx["d"];

    dI1 = 0*I1["abc"] +
        I1["abe"]*G["ecd"]*dx["d"] - I1["aec"]*G["bed"]*dx["d"]
        + R["bedc"]*dx["d"]*Iinv["ae"]
        - 1./2.*R["badc"]*dx["d"];

    dI2 = 0*I2["abcd"]
        - I2["aecd"]*G["bef"]*dx["f"] + I2["abed"]*G["ecf"]*dx["f"]
        + I2["abce"]*G["edf"]*dx["f"]

        + R["befd"]*dx["f"]*I1["aec"] + R["befc"]*dx["f"]*I1["aed"]
        - R["ecfd"]*dx["f"]*I1["abe"] + dR["befcd"]*dx["f"]*Iinv["ae"]

        - 1./3.*(dR["baecd"]+dR["baedc"])*dx["e"];

    dV0 = 0; //FIXME missing potential
  } else {
    double tauinv = 1./tau;
    /* Transport equation for Q = xi-1 */
    dQ = (Q["ac"]*G["cbd"] - G["acd"]*Q["cb"])*dx["d"]
       - tauinv*(Q["ab"] + Q["ac"]*Q["cb"]) - tau*R["acbd"]*dx["c"]*dx["d"];

    /* Q(theta,theta) blows up as theta*cot(theta) and makes the numerical
     * integration break down. Since we know the analytic form, don't compute it
     * numerically */
    //dQ(2,2) = 0.0;

    Tensor kd("^a_b");
    kd(0,0) = 1; kd(1,1) = 1; kd(2,2) = 1; kd(3,3) = 1;
    Tensor xi("^a_b");
    xi = Q["ab"] + kd["ab"];

    /* Transport equation for the Van Vleck determinant */
    dVV = -tauinv*0.5 * VV * Q["aa"];

    /* Transport equation for eta = sigma^a_b' */
    deta = eta["ac"]*G["cbd"]*dx["d"] - tauinv*(eta["ac"]*xi["cb"] - eta["ab"]);

    /* Transport equation for sigma^a_{b'c'} */
    ds12 = s12["aec"]*G["ebd"]*dx["d"] + s12["abe"]*G["ecd"]*dx["d"]
          - tauinv*(xi["db"]*s12["adc"] + xi["dc"]*s12["adb"]
          + eta["ad"]*s3["dbc"] - s12["abc"]) - eta["ad"]*R["dbec"]*dx["e"];

    /* Transport equation for sigma^a'_{b'c'} */
    ds3 = tauinv*s3["abc"]
         - tauinv*(xi["db"]*s3["adc"] + xi["dc"]*s3["adb"]
         + xi["ad"]*s3["dbc"])
         - xi["ad"]*R["dbec"]*dx["e"]
         + xi["db"]*R["adec"]*dx["e"] + xi["dc"]*R["adeb"]*dx["e"]
         - tau*dR["adbec"]*dx["d"]*dx["e"]
         + s3["aec"]*G["ebd"]*dx["d"] + s3["abe"]*G["ecd"]*dx["d"]
         - s3["ebc"]*G["aed"]*dx["d"];

    /* Transport equation for sigma^a'_{b'c' d'} */
    ds4 = 0*s4["abcd"]
        - tau*(ddR["aebfcd"]*dx["e"]*dx["f"]
                - R["afed"]*R["ebgc"]*dx["f"]*dx["g"]
                - R["afec"]*R["ebgd"]*dx["f"]*dx["g"])

        + s4["aecd"]*G["ebf"]*dx["f"] + s4["abed"]*G["ecf"]*dx["f"]
        + s4["abce"]*G["edf"]*dx["f"] - s4["ebcd"]*G["aef"]*dx["f"]

        - dR["fbecd"]*dx["e"]*xi["af"] + dR["afecd"]*dx["e"]*xi["fb"]
        + dR["afebd"]*dx["e"]*xi["fc"] + dR["afebc"]*dx["e"]*xi["fd"]

        - R["ecfd"]*dx["f"]*s3["abe"]
        - R["ebfd"]*dx["f"]*s3["ace"]
        - R["ebfc"]*dx["f"]*s3["ade"]
        + R["agfc"]*dx["f"]*s3["gbd"]
        + R["agfb"]*dx["f"]*s3["gcd"]
        + R["agfd"]*dx["f"]*s3["gbc"]

        - tauinv*(
          s3["ecd"]*s3["aeb"] + s3["ebd"]*s3["aec"] + s3["ebc"]*s3["aed"]
          - s4["abcd"]
          + xi["ed"]*s4["aebc"] + xi["ec"]*s4["aebd"]
          + xi["eb"]*s4["aecd"] + xi["ae"]*s4["ebcd"]);

    ds13 = 0*s13["abcd"]
        + s13["aecd"]*G["ebf"]*dx["f"] + s13["abed"]*G["ecf"]*dx["f"]
        + s13["abce"]*G["edf"]*dx["f"]

        - eta["ae"]*dR["ebfcd"]*dx["f"]

        - s12["ade"]*R["ebfc"]*dx["f"]
        - s12["ace"]*R["ebfd"]*dx["f"]
        - s12["abe"]*R["ecfd"]*dx["f"]

        - tauinv*(
            s13["aebc"]*xi["ed"] + s13["aebd"]*xi["ec"] + s13["aecd"]*xi["eb"]
          + s12["aeb"]*s3["ecd"] + s12["aec"]*s3["ebd"] + s12["aed"]*s3["ebc"]
          + eta["ae"]*s4["ebcd"]
          - s13["abcd"]
        );

    dI = 0*I["ab"] + I["cb"]*G["cad"]*dx["d"];

    dI1 = 0*I1["abc"] +
        I1["abe"]*G["ecd"]*dx["d"] - I1["aec"]*G["bed"]*dx["d"]
        + R["bedc"]*dx["d"]*Iinv["ae"]
        - tauinv*xi["dc"]*I1["abd"];

    dI2 = 0*I2["abcd"]
        - I2["aecd"]*G["bef"]*dx["f"] + I2["abed"]*G["ecf"]*dx["f"]
        + I2["abce"]*G["edf"]*dx["f"]

        + R["befd"]*dx["f"]*I1["aec"] + R["befc"]*dx["f"]*I1["aed"]
        - R["ecfd"]*dx["f"]*I1["abe"] + dR["befcd"]*dx["f"]*Iinv["ae"]

        - tauinv*(I2["abed"]*xi["ec"] + I2["abec"]*xi["ed"] + I1["abe"]*s3["ecd"]);

    Tensor BoxSqrtDelta(0);
    BoxSqrtDelta = 0.5*VV*(
        0.5*(I["ab"]*I1["bac"]+gamma["ab"]*s12["bac"])*
            (I["de"]*I1["edf"]+gamma["de"]*s12["edf"])*gu["cf"]
        - I["ab"]*I1["bcd"]*I["ce"]*I1["eaf"]*gu["df"]
        - gamma["ab"]*s12["bcd"]*gamma["ce"]*s12["eaf"]*gu["df"]
        + I["ab"]*I2["bacd"]*gu["cd"] + gamma["ab"]*s13["bacd"]*gu["cd"]);

    Tensor temp(0);
    temp = -0.5*V0*Q["aa"];
    dV0(0) = tauinv*( temp(0) - V0(0) - 0.5*BoxSqrtDelta(0)); //FIXME missing potential
  }

  /* Pack the data back into an array for GSL */
  dT.getComponents(f);

  return GSL_SUCCESS;
}

int main (int argc, char * argv[])
{
  /* Evolved tensors */
  TensorList T;
  T.append("x", "^a");         Tensor& x  = T["x"];  /* Position */
  T.append("ur");              Tensor& ur = T["ur"]; /* Radial 4-velocity */
  T.append("Q", "^a_b");       Tensor& Q  = T["Q"];  /* \sigma^a'_b' - \delta^a'_b' */
  T.append("VV");              Tensor& VV = T["VV"]; /* Van Vleck determinant */
  T.append("eta", "^a_b");     Tensor& eta  = T["eta"];  /* \sigma^a'_b' */
  T.append("s12", "^a_b_c");   Tensor& s12 = T["s12"]; /* \sigma^a_{b'c'} */
  T.append("s3",  "^a_b_c");   Tensor& s3  = T["s3"];  /* \sigma^a'_{b'c'} */
  T.append("s13", "^a_b_c_d"); Tensor& s13 = T["s13"]; /* \sigma^a_{b'c'd'} */
  T.append("s4",  "^a_b_c_d"); Tensor& s4  = T["s4"];  /* \sigma^a'_{b'c'd'} */
  T.append("I",  "_a^b");      Tensor& I   = T["I"];   /* g_a'^b */
  T.append("I1", "_a^b_c");    Tensor& I1  = T["I1"];  /* g_a^b'_;c' */
  T.append("I2", "_a^b_c_d");  Tensor& I2  = T["I2"];  /* g_a^b'_{;c' d'} */
  T.append("V0");              Tensor& V0 = T["V0"]; /* V_0 */

  const int numEqs = T.getNumComponents();

  /* Spacetime */
  Schwarzschild schw(1.0);
  Tensor &R  = schw["Ruddd"];
  Tensor &gu = schw["guu"];
  Tensor &RicciScalar = schw["R"];

  /* Use a Runge-Kutta integrator with adaptive step-size */
  const gsl_odeiv_step_type * t = gsl_odeiv_step_rkf45;
  gsl_odeiv_step * s = gsl_odeiv_step_alloc (t, numEqs);
  gsl_odeiv_control * c = gsl_odeiv_control_standard_new (1e-6, 1e-6, 1.0, 1.0);
  gsl_odeiv_evolve * e = gsl_odeiv_evolve_alloc (numEqs);

  /* Time-like circular geodesic at r=10M */
  struct geodesic_params params =
    {0.956183, 3.77964, Timelike, T, schw};

  gsl_odeiv_system sys = {rhs, NULL, numEqs, &params};

  /* Initial Conditions */
  double tau = 0.0, tau1 = 100.0;
  double h = 1e-6;
  double r0 = 10.0;

  x = 0.*x["a"];
  ur() = 0.;
  Q = 0.*Q["ab"];
  VV() = 0.;
  eta = 0.*eta["ab"];
  s12 = 0.*s12["abc"];
  s3 = 0.*s3["abc"];
  s13 = 0.*s13["abcd"];
  s4 = 0.*s4["abcd"];
  I = 0*I["ab"];
  I1 = 0*I1["abc"];
  I2 = 0*I2["abcd"];
  double mass = 0.0;
  double xi = 0.0;
  V0() = 0.5*(mass*mass + (xi-1./6.)*RicciScalar());

  x(1) = r0;
  x(2) = M_PI_2;
  VV() = 1.0;
  eta(0,0) = -1.0;
  eta(1,1) = -1.0;
  eta(2,2) = -1.0;
  eta(3,3) = -1.0;
  I(0,0) = 1.0;
  I(1,1) = 1.0;
  I(2,2) = 1.0;
  I(3,3) = 1.0;

  schw.setPoint(x);
  schw.calc_all();
  s4 = s4["abcd"] - 1./3.*(R["acbd"] + R["adbc"]);
  s13 = s13["abcd"] - 1./6.*(R["acbd"] + R["adbc"]) - 1./2.*R["abcd"];
  I2 = I2["abcd"] - 1./2.*R["bacd"];

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
    printf ("%.5g", tau);
    for(int i=0; i<numEqs; i++)
    {
      printf (", %.5g", y[i]);
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
