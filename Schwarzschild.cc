/* Components of the Christoffel symbols, Riemann tensor and its derivatives
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

#include "SpacetimeTensors.h"
#include <cassert>
#include <math.h>

Spacetime::Spacetime()
{
  append("guu", "^a^b");              /* Contravariant metric */
  append("gdd", "_a_b");              /* Covariant metric */
  append("Gudd", "^a_b_c");           /* Christoffel symbol */
  append("Ruddd", "^a_b_c_d");        /* Riemann */
  append("Ruddd;d", "^a_b_c_d_e");    /* Covariant derivative of Riemann */
  append("Ruddd;dd", "^a_b_c_d_e_f"); /* Double covariant derivative of Riemann */
  append("R");                        /* Ricci scalar */
  append("R;d", "_a");                /* Covariant derivative of Ricci scalar */
}

Schwarzschild::Schwarzschild(double mass)
  : M(mass)
{
}

void Schwarzschild::calc_all()
{
  calc_guu();
  calc_gdd();
  calc_Gudd();
  calc_Ruddd();
  calc_Rudddcd();
  calc_Rudddcdd();
  calc_R();
  calc_Rcd();
}

void Schwarzschild::calc_guu()
{
  Tensor& guu = (*this)["guu"];
  double sintheta2 = sintheta*sintheta;
  double r2 = r*r;

  guu(0, 0) = -r/(r-2.0*M);
  guu(1, 1) = (r-2.0*M)/r;
  guu(2, 2) = 1/r2;
  guu(3, 3) = 1/(r2*sintheta2);
}

void Schwarzschild::calc_gdd()
{
  Tensor& gdd = (*this)["gdd"];
  double sintheta2 = sintheta*sintheta;
  double r2 = r*r;

  gdd(0, 0) = -(r-2.0*M)/r;
  gdd(1, 1) = r/(r-2.0*M);
  gdd(2, 2) = r2;
  gdd(3, 3) = r2*sintheta2;
}

void Schwarzschild::calc_Gudd()
{
  Tensor& Gudd = (*this)["Gudd"];

  Gudd(0, 0, 1) = -(M/((2*M - r)*r));
  Gudd(0, 1, 0) = -(M/((2*M - r)*r));
  Gudd(1, 0, 0) = (M*(-2*M + r))/pow(r,3);
  Gudd(1, 1, 1) = M/((2*M - r)*r);
  Gudd(1, 2, 2) = 2*M - r;
  Gudd(1, 3, 3) = (2*M - r)*pow(sin(theta),2);
  Gudd(2, 1, 2) = 1/r;
  Gudd(2, 2, 1) = 1/r;
  Gudd(2, 3, 3) = -(cos(theta)*sin(theta));
  Gudd(3, 1, 3) = 1/r;
  Gudd(3, 2, 3) = cos(theta)/sin(theta);
  Gudd(3, 3, 1) = 1/r;
  Gudd(3, 3, 2) = cos(theta)/sin(theta);
}

void Schwarzschild::calc_Ruddd()
{
  Tensor& Ruddd = (*this)["Ruddd"];

  Ruddd(0, 1, 0, 1) = (-2*M)/((2*M - r)*pow(r,2));
  Ruddd(0, 1, 1, 0) = (2*M)/((2*M - r)*pow(r,2));
  Ruddd(0, 2, 0, 2) = (M*(-2*M + r))/((2*M - r)*r);
  Ruddd(0, 2, 2, 0) = -((M*(-2*M + r))/((2*M - r)*r));
  Ruddd(0, 3, 0, 3) = -((M*pow(sin(theta),2))/r);
  Ruddd(0, 3, 3, 0) = (M*pow(sin(theta),2))/r;
  Ruddd(1, 0, 0, 1) = (2*M*(1 - (2*M)/r))/pow(r,3);
  Ruddd(1, 0, 1, 0) = (-2*M*(1 - (2*M)/r))/pow(r,3);
  Ruddd(1, 2, 1, 2) = (M*(1 - (2*M)/r))/(2*M - r);
  Ruddd(1, 2, 2, 1) = -((M*(1 - (2*M)/r))/(2*M - r));
  Ruddd(1, 3, 1, 3) = (M*(1 - (2*M)/r)*pow(sin(theta),2))/(2*M - r);
  Ruddd(1, 3, 3, 1) = -((M*(1 - (2*M)/r)*pow(sin(theta),2))/(2*M - r));
  Ruddd(2, 0, 0, 2) = -((M*(-2*M + r))/pow(r,4));
  Ruddd(2, 0, 2, 0) = (M*(-2*M + r))/pow(r,4);
  Ruddd(2, 1, 1, 2) = -(M/((2*M - r)*pow(r,2)));
  Ruddd(2, 1, 2, 1) = M/((2*M - r)*pow(r,2));
  Ruddd(2, 3, 2, 3) = (2*M*pow(sin(theta),2))/r;
  Ruddd(2, 3, 3, 2) = (-2*M*pow(sin(theta),2))/r;
  Ruddd(3, 0, 0, 3) = (M*(2*M - r))/pow(r,4);
  Ruddd(3, 0, 3, 0) = -((M*(2*M - r))/pow(r,4));
  Ruddd(3, 1, 1, 3) = -(M/((2*M - r)*pow(r,2)));
  Ruddd(3, 1, 3, 1) = M/((2*M - r)*pow(r,2));
  Ruddd(3, 2, 2, 3) = (-2*M)/r;
  Ruddd(3, 2, 3, 2) = (2*M)/r;
}

void Schwarzschild::calc_Rudddcd()
{
  Tensor& Rudddcd = (*this)["Ruddd;d"];

  Rudddcd(0, 1, 0, 1, 1) = (6*M)/((2*M - r)*pow(r,3));
  Rudddcd(0, 1, 0, 2, 2) = (3*M)/pow(r,2);
  Rudddcd(0, 1, 0, 3, 3) = (3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(0, 1, 1, 0, 1) = (-6*M)/((2*M - r)*pow(r,3));
  Rudddcd(0, 1, 2, 0, 2) = (-3*M)/pow(r,2);
  Rudddcd(0, 1, 3, 0, 3) = (-3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(0, 2, 0, 1, 2) = (3*M)/pow(r,2);
  Rudddcd(0, 2, 0, 2, 1) = (3*M)/pow(r,2);
  Rudddcd(0, 2, 1, 0, 2) = (-3*M)/pow(r,2);
  Rudddcd(0, 2, 2, 0, 1) = (-3*M)/pow(r,2);
  Rudddcd(0, 3, 0, 1, 3) = (3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(0, 3, 0, 3, 1) = (3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(0, 3, 1, 0, 3) = (-3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(0, 3, 3, 0, 1) = (-3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(1, 0, 0, 1, 1) = (6*M*(2*M - r))/pow(r,5);
  Rudddcd(1, 0, 0, 2, 2) = (3*M*pow(-2*M + r,2))/pow(r,4);
  Rudddcd(1, 0, 0, 3, 3) = (3*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,4);
  Rudddcd(1, 0, 1, 0, 1) = (6*M*(-2*M + r))/pow(r,5);
  Rudddcd(1, 0, 2, 0, 2) = (-3*M*pow(-2*M + r,2))/pow(r,4);
  Rudddcd(1, 0, 3, 0, 3) = (-3*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,4);
  Rudddcd(1, 2, 1, 2, 1) = (3*M)/pow(r,2);
  Rudddcd(1, 2, 2, 1, 1) = (-3*M)/pow(r,2);
  Rudddcd(1, 2, 2, 3, 3) = (-3*M*(2*M - r)*pow(sin(theta),2))/r;
  Rudddcd(1, 2, 3, 2, 3) = (3*M*(2*M - r)*pow(sin(theta),2))/r;
  Rudddcd(1, 3, 1, 3, 1) = (3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(1, 3, 2, 3, 2) = (3*M*(2*M - r)*pow(sin(theta),2))/r;
  Rudddcd(1, 3, 3, 1, 1) = (-3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(1, 3, 3, 2, 2) = (-3*M*(2*M - r)*pow(sin(theta),2))/r;
  Rudddcd(2, 0, 0, 1, 2) = (3*M*(-2*M + r))/pow(r,5);
  Rudddcd(2, 0, 0, 2, 1) = (3*M*(-2*M + r))/pow(r,5);
  Rudddcd(2, 0, 1, 0, 2) = (3*M*(2*M - r))/pow(r,5);
  Rudddcd(2, 0, 2, 0, 1) = (3*M*(2*M - r))/pow(r,5);
  Rudddcd(2, 1, 1, 2, 1) = (3*M)/((2*M - r)*pow(r,3));
  Rudddcd(2, 1, 2, 1, 1) = (-3*M)/((2*M - r)*pow(r,3));
  Rudddcd(2, 1, 2, 3, 3) = (-3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(2, 1, 3, 2, 3) = (3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(2, 3, 1, 2, 3) = (3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(2, 3, 1, 3, 2) = (-3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(2, 3, 2, 1, 3) = (-3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(2, 3, 2, 3, 1) = (-6*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(2, 3, 3, 1, 2) = (3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(2, 3, 3, 2, 1) = (6*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd(3, 0, 0, 1, 3) = (3*M*(-2*M + r))/pow(r,5);
  Rudddcd(3, 0, 0, 3, 1) = (3*M*(-2*M + r))/pow(r,5);
  Rudddcd(3, 0, 1, 0, 3) = (3*M*(2*M - r))/pow(r,5);
  Rudddcd(3, 0, 3, 0, 1) = (3*M*(2*M - r))/pow(r,5);
  Rudddcd(3, 1, 1, 3, 1) = (3*M)/((2*M - r)*pow(r,3));
  Rudddcd(3, 1, 2, 3, 2) = (3*M)/pow(r,2);
  Rudddcd(3, 1, 3, 1, 1) = (-3*M)/((2*M - r)*pow(r,3));
  Rudddcd(3, 1, 3, 2, 2) = (-3*M)/pow(r,2);
  Rudddcd(3, 2, 1, 2, 3) = (-3*M)/pow(r,2);
  Rudddcd(3, 2, 1, 3, 2) = (3*M)/pow(r,2);
  Rudddcd(3, 2, 2, 1, 3) = (3*M)/pow(r,2);
  Rudddcd(3, 2, 2, 3, 1) = (6*M)/pow(r,2);
  Rudddcd(3, 2, 3, 1, 2) = (-3*M)/pow(r,2);
  Rudddcd(3, 2, 3, 2, 1) = (-6*M)/pow(r,2);
}

void Schwarzschild::calc_Rudddcdd()
{
  Tensor& Rudddcdd = (*this)["Ruddd;dd"];

  Rudddcdd(0, 1, 0, 1, 0, 0) = (6*pow(M,2))/pow(r,6);
  Rudddcdd(0, 1, 0, 1, 1, 1) = (6*M*(-9*M + 4*r))/(pow(r,4)*pow(-2*M + r,2));
  Rudddcdd(0, 1, 0, 1, 2, 2) = (-12*M)/pow(r,3);
  Rudddcdd(0, 1, 0, 1, 3, 3) = (-12*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(0, 1, 0, 2, 1, 2) = (-12*M)/pow(r,3);
  Rudddcdd(0, 1, 0, 2, 2, 1) = (3*M*(-9*M + 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 1, 0, 3, 1, 3) = (-12*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(0, 1, 0, 3, 3, 1) = (-3*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 1, 1, 0, 0, 0) = (-6*pow(M,2))/pow(r,6);
  Rudddcdd(0, 1, 1, 0, 1, 1) = (6*M*(9*M - 4*r))/(pow(r,4)*pow(-2*M + r,2));
  Rudddcdd(0, 1, 1, 0, 2, 2) = (12*M)/pow(r,3);
  Rudddcdd(0, 1, 1, 0, 3, 3) = (12*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(0, 1, 1, 2, 2, 0) = (3*pow(M,2))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 1, 1, 3, 3, 0) = (3*pow(M,2)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 1, 2, 0, 1, 2) = (12*M)/pow(r,3);
  Rudddcdd(0, 1, 2, 0, 2, 1) = (3*M*(9*M - 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 1, 2, 1, 2, 0) = (-3*pow(M,2))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 1, 3, 0, 1, 3) = (12*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(0, 1, 3, 0, 3, 1) = (3*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 1, 3, 1, 3, 0) = (-3*pow(M,2)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 2, 0, 1, 1, 2) = (-12*M)/pow(r,3);
  Rudddcdd(0, 2, 0, 1, 2, 1) = (3*M*(-9*M + 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 2, 0, 2, 0, 0) = (3*pow(M,2)*(2*M - r))/pow(r,5);
  Rudddcdd(0, 2, 0, 2, 1, 1) = (3*M*(-9*M + 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 2, 0, 2, 2, 2) = (-9*M*(2*M - r))/pow(r,2);
  Rudddcdd(0, 2, 0, 2, 3, 3) = (-3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(0, 2, 0, 3, 2, 3) = (-3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(0, 2, 0, 3, 3, 2) = (-3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(0, 2, 1, 0, 1, 2) = (12*M)/pow(r,3);
  Rudddcdd(0, 2, 1, 0, 2, 1) = (3*M*(9*M - 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 2, 2, 0, 0, 0) = (3*pow(M,2)*(-2*M + r))/pow(r,5);
  Rudddcdd(0, 2, 2, 0, 1, 1) = (3*M*(9*M - 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 2, 2, 0, 2, 2) = (9*M*(2*M - r))/pow(r,2);
  Rudddcdd(0, 2, 2, 0, 3, 3) = (3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(0, 2, 2, 3, 3, 0) = (3*pow(M,2)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(0, 2, 3, 0, 2, 3) = (3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(0, 2, 3, 0, 3, 2) = (3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(0, 2, 3, 2, 3, 0) = (-3*pow(M,2)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(0, 3, 0, 1, 1, 3) = (-12*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(0, 3, 0, 1, 3, 1) = (-3*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 3, 0, 2, 2, 3) = (-3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(0, 3, 0, 2, 3, 2) = (-3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(0, 3, 0, 3, 0, 0) = (3*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(0, 3, 0, 3, 1, 1) = (-3*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 3, 0, 3, 2, 2) = (-3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(0, 3, 0, 3, 3, 3) = (-9*M*(2*M - r)*pow(sin(theta),4))/pow(r,2);
  Rudddcdd(0, 3, 1, 0, 1, 3) = (12*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(0, 3, 1, 0, 3, 1) = (3*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 3, 2, 0, 2, 3) = (3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(0, 3, 2, 0, 3, 2) = (3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(0, 3, 2, 3, 2, 0) = (-3*pow(M,2)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(0, 3, 3, 0, 0, 0) = (-3*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(0, 3, 3, 0, 1, 1) = (3*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(0, 3, 3, 0, 2, 2) = (3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(0, 3, 3, 0, 3, 3) = (9*M*(2*M - r)*pow(sin(theta),4))/pow(r,2);
  Rudddcdd(0, 3, 3, 2, 2, 0) = (3*pow(M,2)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 0, 0, 1, 0, 0) = (6*pow(M,2)*pow(-2*M + r,2))/pow(r,8);
  Rudddcdd(1, 0, 0, 1, 1, 1) = (6*M*(-9*M + 4*r))/pow(r,6);
  Rudddcdd(1, 0, 0, 1, 2, 2) = (-12*M*pow(-2*M + r,2))/pow(r,5);
  Rudddcdd(1, 0, 0, 1, 3, 3) = (-12*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(1, 0, 0, 2, 1, 2) = (-12*M*pow(-2*M + r,2))/pow(r,5);
  Rudddcdd(1, 0, 0, 2, 2, 1) = (-3*M*(9*M - 4*r)*(2*M - r))/pow(r,5);
  Rudddcdd(1, 0, 0, 3, 1, 3) = (-12*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(1, 0, 0, 3, 3, 1) = (-3*M*(9*M - 4*r)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(1, 0, 1, 0, 0, 0) = (-6*pow(M,2)*pow(-2*M + r,2))/pow(r,8);
  Rudddcdd(1, 0, 1, 0, 1, 1) = (6*M*(9*M - 4*r))/pow(r,6);
  Rudddcdd(1, 0, 1, 0, 2, 2) = (12*M*pow(-2*M + r,2))/pow(r,5);
  Rudddcdd(1, 0, 1, 0, 3, 3) = (12*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(1, 0, 1, 2, 2, 0) = (3*pow(M,2)*(2*M - r))/pow(r,5);
  Rudddcdd(1, 0, 1, 3, 3, 0) = (3*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(1, 0, 2, 0, 1, 2) = (12*M*pow(-2*M + r,2))/pow(r,5);
  Rudddcdd(1, 0, 2, 0, 2, 1) = (3*M*(9*M - 4*r)*(2*M - r))/pow(r,5);
  Rudddcdd(1, 0, 2, 1, 2, 0) = (3*pow(M,2)*(-2*M + r))/pow(r,5);
  Rudddcdd(1, 0, 3, 0, 1, 3) = (12*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(1, 0, 3, 0, 3, 1) = (3*M*(9*M - 4*r)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(1, 0, 3, 1, 3, 0) = (-3*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(1, 2, 0, 1, 2, 0) = (3*pow(M,2)*(-2*M + r))/pow(r,5);
  Rudddcdd(1, 2, 1, 0, 2, 0) = (3*pow(M,2)*(2*M - r))/pow(r,5);
  Rudddcdd(1, 2, 1, 2, 0, 0) = (3*pow(M,2)*(2*M - r))/pow(r,5);
  Rudddcdd(1, 2, 1, 2, 1, 1) = (3*M*(-9*M + 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(1, 2, 1, 2, 2, 2) = (3*M*(-2*M + r))/pow(r,2);
  Rudddcdd(1, 2, 1, 2, 3, 3) = (-9*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 2, 1, 3, 2, 3) = (3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 2, 1, 3, 3, 2) = (3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 2, 2, 1, 0, 0) = (3*pow(M,2)*(-2*M + r))/pow(r,5);
  Rudddcdd(1, 2, 2, 1, 1, 1) = (3*M*(9*M - 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(1, 2, 2, 1, 2, 2) = (3*M*(2*M - r))/pow(r,2);
  Rudddcdd(1, 2, 2, 1, 3, 3) = (9*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 2, 2, 3, 1, 3) = (12*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 2, 2, 3, 3, 1) = (3*M*(9*M - 4*r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 2, 3, 1, 2, 3) = (-3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 2, 3, 1, 3, 2) = (-3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 2, 3, 2, 1, 3) = (-12*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 2, 3, 2, 3, 1) = (-3*M*(9*M - 4*r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 3, 0, 1, 3, 0) = (-3*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(1, 3, 1, 0, 3, 0) = (3*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(1, 3, 1, 2, 2, 3) = (3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 3, 1, 2, 3, 2) = (3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 3, 1, 3, 0, 0) = (3*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(1, 3, 1, 3, 1, 1) = (-3*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(1, 3, 1, 3, 2, 2) = (-9*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 3, 1, 3, 3, 3) = (-3*M*(2*M - r)*pow(sin(theta),4))/pow(r,2);
  Rudddcdd(1, 3, 2, 1, 2, 3) = (-3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 3, 2, 1, 3, 2) = (-3*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 3, 2, 3, 1, 2) = (-12*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 3, 2, 3, 2, 1) = (-3*M*(9*M - 4*r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 3, 3, 1, 0, 0) = (-3*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(1, 3, 3, 1, 1, 1) = (3*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(1, 3, 3, 1, 2, 2) = (9*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 3, 3, 1, 3, 3) = (3*M*(2*M - r)*pow(sin(theta),4))/pow(r,2);
  Rudddcdd(1, 3, 3, 2, 1, 2) = (12*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(1, 3, 3, 2, 2, 1) = (3*M*(9*M - 4*r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(2, 0, 0, 1, 1, 2) = (12*M*(2*M - r))/pow(r,6);
  Rudddcdd(2, 0, 0, 1, 2, 1) = (3*M*(9*M - 4*r))/pow(r,6);
  Rudddcdd(2, 0, 0, 2, 0, 0) = (-3*pow(M,2)*pow(-2*M + r,2))/pow(r,8);
  Rudddcdd(2, 0, 0, 2, 1, 1) = (3*M*(9*M - 4*r))/pow(r,6);
  Rudddcdd(2, 0, 0, 2, 2, 2) = (9*M*pow(-2*M + r,2))/pow(r,5);
  Rudddcdd(2, 0, 0, 2, 3, 3) = (3*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(2, 0, 0, 3, 2, 3) = (3*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(2, 0, 0, 3, 3, 2) = (3*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(2, 0, 1, 0, 1, 2) = (12*M*(-2*M + r))/pow(r,6);
  Rudddcdd(2, 0, 1, 0, 2, 1) = (3*M*(-9*M + 4*r))/pow(r,6);
  Rudddcdd(2, 0, 2, 0, 0, 0) = (3*pow(M,2)*pow(-2*M + r,2))/pow(r,8);
  Rudddcdd(2, 0, 2, 0, 1, 1) = (3*M*(-9*M + 4*r))/pow(r,6);
  Rudddcdd(2, 0, 2, 0, 2, 2) = (-9*M*pow(-2*M + r,2))/pow(r,5);
  Rudddcdd(2, 0, 2, 0, 3, 3) = (-3*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(2, 0, 2, 3, 3, 0) = (-3*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(2, 0, 3, 0, 2, 3) = (-3*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(2, 0, 3, 0, 3, 2) = (-3*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(2, 0, 3, 2, 3, 0) = (3*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(2, 1, 0, 1, 2, 0) = (-3*pow(M,2))/pow(r,6);
  Rudddcdd(2, 1, 1, 0, 2, 0) = (3*pow(M,2))/pow(r,6);
  Rudddcdd(2, 1, 1, 2, 0, 0) = (3*pow(M,2))/pow(r,6);
  Rudddcdd(2, 1, 1, 2, 1, 1) = (3*M*(-9*M + 4*r))/(pow(r,4)*pow(-2*M + r,2));
  Rudddcdd(2, 1, 1, 2, 2, 2) = (-3*M)/pow(r,3);
  Rudddcdd(2, 1, 1, 2, 3, 3) = (-9*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(2, 1, 1, 3, 2, 3) = (3*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(2, 1, 1, 3, 3, 2) = (3*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(2, 1, 2, 1, 0, 0) = (-3*pow(M,2))/pow(r,6);
  Rudddcdd(2, 1, 2, 1, 1, 1) = (3*M*(9*M - 4*r))/(pow(r,4)*pow(-2*M + r,2));
  Rudddcdd(2, 1, 2, 1, 2, 2) = (3*M)/pow(r,3);
  Rudddcdd(2, 1, 2, 1, 3, 3) = (9*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(2, 1, 2, 3, 1, 3) = (12*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(2, 1, 2, 3, 3, 1) = (3*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(2, 1, 3, 1, 2, 3) = (-3*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(2, 1, 3, 1, 3, 2) = (-3*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(2, 1, 3, 2, 1, 3) = (-12*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(2, 1, 3, 2, 3, 1) = (-3*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(2, 3, 0, 2, 3, 0) = (3*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(2, 3, 0, 3, 2, 0) = (-3*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(2, 3, 1, 2, 1, 3) = (-12*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(2, 3, 1, 2, 3, 1) = (-3*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(2, 3, 1, 3, 1, 2) = (12*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(2, 3, 1, 3, 2, 1) = (3*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(2, 3, 2, 0, 3, 0) = (-3*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(2, 3, 2, 1, 1, 3) = (12*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(2, 3, 2, 1, 3, 1) = (3*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(2, 3, 2, 3, 0, 0) = (-6*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(2, 3, 2, 3, 1, 1) = (6*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(2, 3, 2, 3, 2, 2) = (12*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(2, 3, 2, 3, 3, 3) = (12*M*(2*M - r)*pow(sin(theta),4))/pow(r,2);
  Rudddcdd(2, 3, 3, 0, 2, 0) = (3*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(2, 3, 3, 1, 1, 2) = (-12*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(2, 3, 3, 1, 2, 1) = (-3*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(2, 3, 3, 2, 0, 0) = (6*pow(M,2)*(2*M - r)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(2, 3, 3, 2, 1, 1) = (-6*M*(9*M - 4*r)*pow(sin(theta),2))/((2*M - r)*pow(r,3));
  Rudddcdd(2, 3, 3, 2, 2, 2) = (-12*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(2, 3, 3, 2, 3, 3) = (-12*M*(2*M - r)*pow(sin(theta),4))/pow(r,2);
  Rudddcdd(3, 0, 0, 1, 1, 3) = (12*M*(2*M - r))/pow(r,6);
  Rudddcdd(3, 0, 0, 1, 3, 1) = (3*M*(9*M - 4*r))/pow(r,6);
  Rudddcdd(3, 0, 0, 2, 2, 3) = (3*M*pow(-2*M + r,2))/pow(r,5);
  Rudddcdd(3, 0, 0, 2, 3, 2) = (3*M*pow(-2*M + r,2))/pow(r,5);
  Rudddcdd(3, 0, 0, 3, 0, 0) = (-3*pow(M,2)*pow(-2*M + r,2))/pow(r,8);
  Rudddcdd(3, 0, 0, 3, 1, 1) = (3*M*(9*M - 4*r))/pow(r,6);
  Rudddcdd(3, 0, 0, 3, 2, 2) = (3*M*pow(-2*M + r,2))/pow(r,5);
  Rudddcdd(3, 0, 0, 3, 3, 3) = (9*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(3, 0, 1, 0, 1, 3) = (12*M*(-2*M + r))/pow(r,6);
  Rudddcdd(3, 0, 1, 0, 3, 1) = (3*M*(-9*M + 4*r))/pow(r,6);
  Rudddcdd(3, 0, 2, 0, 2, 3) = (-3*M*pow(-2*M + r,2))/pow(r,5);
  Rudddcdd(3, 0, 2, 0, 3, 2) = (-3*M*pow(-2*M + r,2))/pow(r,5);
  Rudddcdd(3, 0, 2, 3, 2, 0) = (3*pow(M,2)*(2*M - r))/pow(r,5);
  Rudddcdd(3, 0, 3, 0, 0, 0) = (3*pow(M,2)*pow(-2*M + r,2))/pow(r,8);
  Rudddcdd(3, 0, 3, 0, 1, 1) = (3*M*(-9*M + 4*r))/pow(r,6);
  Rudddcdd(3, 0, 3, 0, 2, 2) = (-3*M*pow(-2*M + r,2))/pow(r,5);
  Rudddcdd(3, 0, 3, 0, 3, 3) = (-9*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,5);
  Rudddcdd(3, 0, 3, 2, 2, 0) = (3*pow(M,2)*(-2*M + r))/pow(r,5);
  Rudddcdd(3, 1, 0, 1, 3, 0) = (-3*pow(M,2))/pow(r,6);
  Rudddcdd(3, 1, 1, 0, 3, 0) = (3*pow(M,2))/pow(r,6);
  Rudddcdd(3, 1, 1, 2, 2, 3) = (3*M)/pow(r,3);
  Rudddcdd(3, 1, 1, 2, 3, 2) = (3*M)/pow(r,3);
  Rudddcdd(3, 1, 1, 3, 0, 0) = (3*pow(M,2))/pow(r,6);
  Rudddcdd(3, 1, 1, 3, 1, 1) = (3*M*(-9*M + 4*r))/(pow(r,4)*pow(-2*M + r,2));
  Rudddcdd(3, 1, 1, 3, 2, 2) = (-9*M)/pow(r,3);
  Rudddcdd(3, 1, 1, 3, 3, 3) = (-3*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(3, 1, 2, 1, 2, 3) = (-3*M)/pow(r,3);
  Rudddcdd(3, 1, 2, 1, 3, 2) = (-3*M)/pow(r,3);
  Rudddcdd(3, 1, 2, 3, 1, 2) = (-12*M)/pow(r,3);
  Rudddcdd(3, 1, 2, 3, 2, 1) = (3*M*(-9*M + 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(3, 1, 3, 1, 0, 0) = (-3*pow(M,2))/pow(r,6);
  Rudddcdd(3, 1, 3, 1, 1, 1) = (3*M*(9*M - 4*r))/(pow(r,4)*pow(-2*M + r,2));
  Rudddcdd(3, 1, 3, 1, 2, 2) = (9*M)/pow(r,3);
  Rudddcdd(3, 1, 3, 1, 3, 3) = (3*M*pow(sin(theta),2))/pow(r,3);
  Rudddcdd(3, 1, 3, 2, 1, 2) = (12*M)/pow(r,3);
  Rudddcdd(3, 1, 3, 2, 2, 1) = (3*M*(9*M - 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(3, 2, 0, 2, 3, 0) = (3*pow(M,2)*(-2*M + r))/pow(r,5);
  Rudddcdd(3, 2, 0, 3, 2, 0) = (3*pow(M,2)*(2*M - r))/pow(r,5);
  Rudddcdd(3, 2, 1, 2, 1, 3) = (12*M)/pow(r,3);
  Rudddcdd(3, 2, 1, 2, 3, 1) = (3*M*(9*M - 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(3, 2, 1, 3, 1, 2) = (-12*M)/pow(r,3);
  Rudddcdd(3, 2, 1, 3, 2, 1) = (3*M*(-9*M + 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(3, 2, 2, 0, 3, 0) = (3*pow(M,2)*(2*M - r))/pow(r,5);
  Rudddcdd(3, 2, 2, 1, 1, 3) = (-12*M)/pow(r,3);
  Rudddcdd(3, 2, 2, 1, 3, 1) = (3*M*(-9*M + 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(3, 2, 2, 3, 0, 0) = (6*pow(M,2)*(2*M - r))/pow(r,5);
  Rudddcdd(3, 2, 2, 3, 1, 1) = (6*M*(-9*M + 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(3, 2, 2, 3, 2, 2) = (-12*M*(2*M - r))/pow(r,2);
  Rudddcdd(3, 2, 2, 3, 3, 3) = (-12*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
  Rudddcdd(3, 2, 3, 0, 2, 0) = (3*pow(M,2)*(-2*M + r))/pow(r,5);
  Rudddcdd(3, 2, 3, 1, 1, 2) = (12*M)/pow(r,3);
  Rudddcdd(3, 2, 3, 1, 2, 1) = (3*M*(9*M - 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(3, 2, 3, 2, 0, 0) = (-6*pow(M,2)*(2*M - r))/pow(r,5);
  Rudddcdd(3, 2, 3, 2, 1, 1) = (6*M*(9*M - 4*r))/((2*M - r)*pow(r,3));
  Rudddcdd(3, 2, 3, 2, 2, 2) = (12*M*(2*M - r))/pow(r,2);
  Rudddcdd(3, 2, 3, 2, 3, 3) = (12*M*(2*M - r)*pow(sin(theta),2))/pow(r,2);
}

void Schwarzschild::calc_R()
{
  Tensor& R = (*this)["R"];

  R() = 0;
}

/* Covariant derivative of Ricci scalar contracted with 4-velocity */
void Schwarzschild::calc_Rcd()
{
  Tensor& Rcd = (*this)["R;d"];

  Rcd(0) = 0;
  Rcd(1) = 0;
  Rcd(2) = 0;
  Rcd(3) = 0;
}

void Schwarzschild::setPoint(const Tensor& x)
{
  assert(x.getRank() == 1);
  t = x(0);
  r = x(1);
  theta = x(2);
  phi = x(3);

  sintheta = sin(theta);
  costheta = cos(theta);
}
