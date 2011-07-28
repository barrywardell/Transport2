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
#include <math.h>

Schwarzschild::Schwarzschild()
{
  Tensor::IndexType u = Tensor::UP;
  Tensor::IndexType d = Tensor::DOWN;

  guu     = new Tensor(2, u, u);
  gdd     = new Tensor(2, d, d);
  Gudd    = new Tensor(3, u, d, d);
  Ruddd   = new Tensor(4, u, d, d, d);
  Rudddcd = new Tensor(5, u, d, d, d, d);
  R       = new Tensor(0);
  Rcd     = new Tensor(1, d);
}

void Schwarzschild::calc_all()
{
  calc_guu();
  calc_gdd();
  calc_Gudd();
  calc_Ruddd();
  calc_Rudddcd();
  calc_R();
  calc_Rcd();
}

void Schwarzschild::calc_guu()
{
  double sintheta2 = sintheta*sintheta;
  double r2 = r*r;

  guu->get(0, 0) = -r/(r-2.0*M);
  guu->get(1, 1) = (r-2.0*M)/r;
  guu->get(2, 2) = 1/r2;
  guu->get(3, 3) = 1/(r2*sintheta2);
}

void Schwarzschild::calc_gdd()
{
  double sintheta2 = sintheta*sintheta;
  double r2 = r*r;

  gdd->get(0, 0) = -(r-2.0*M)/r;
  gdd->get(1, 1) = r/(r-2.0*M);
  gdd->get(2, 2) = r2;
  gdd->get(3, 3) = r2*sintheta2;
}

void Schwarzschild::calc_Gudd()
{
  Gudd->get(0, 0, 1) = -(M/((2*M - r)*r));
  Gudd->get(0, 1, 0) = -(M/((2*M - r)*r));
  Gudd->get(1, 0, 0) = (M*(-2*M + r))/pow(r,3);
  Gudd->get(1, 1, 1) = M/((2*M - r)*r);
  Gudd->get(1, 2, 2) = 2*M - r;
  Gudd->get(1, 3, 3) = (2*M - r)*pow(sin(theta),2);
  Gudd->get(2, 1, 2) = 1/r;
  Gudd->get(2, 2, 1) = 1/r;
  Gudd->get(2, 3, 3) = -(cos(theta)*sin(theta));
  Gudd->get(3, 1, 3) = 1/r;
  Gudd->get(3, 2, 3) = cos(theta)/sin(theta);
  Gudd->get(3, 3, 1) = 1/r;
  Gudd->get(3, 3, 2) = cos(theta)/sin(theta);
}

void Schwarzschild::calc_Ruddd()
{
  Ruddd->get(0, 1, 0, 1) = (-2*M)/((2*M - r)*pow(r,2));
  Ruddd->get(0, 1, 1, 0) = (2*M)/((2*M - r)*pow(r,2));
  Ruddd->get(0, 2, 0, 2) = (M*(-2*M + r))/((2*M - r)*r);
  Ruddd->get(0, 2, 2, 0) = -((M*(-2*M + r))/((2*M - r)*r));
  Ruddd->get(0, 3, 0, 3) = -((M*pow(sin(theta),2))/r);
  Ruddd->get(0, 3, 3, 0) = (M*pow(sin(theta),2))/r;
  Ruddd->get(1, 0, 0, 1) = (2*M*(1 - (2*M)/r))/pow(r,3);
  Ruddd->get(1, 0, 1, 0) = (-2*M*(1 - (2*M)/r))/pow(r,3);
  Ruddd->get(1, 2, 1, 2) = (M*(1 - (2*M)/r))/(2*M - r);
  Ruddd->get(1, 2, 2, 1) = -((M*(1 - (2*M)/r))/(2*M - r));
  Ruddd->get(1, 3, 1, 3) = (M*(1 - (2*M)/r)*pow(sin(theta),2))/(2*M - r);
  Ruddd->get(1, 3, 3, 1) = -((M*(1 - (2*M)/r)*pow(sin(theta),2))/(2*M - r));
  Ruddd->get(2, 0, 0, 2) = -((M*(-2*M + r))/pow(r,4));
  Ruddd->get(2, 0, 2, 0) = (M*(-2*M + r))/pow(r,4);
  Ruddd->get(2, 1, 1, 2) = -(M/((2*M - r)*pow(r,2)));
  Ruddd->get(2, 1, 2, 1) = M/((2*M - r)*pow(r,2));
  Ruddd->get(2, 3, 2, 3) = (2*M*pow(sin(theta),2))/r;
  Ruddd->get(2, 3, 3, 2) = (-2*M*pow(sin(theta),2))/r;
  Ruddd->get(3, 0, 0, 3) = (M*(2*M - r))/pow(r,4);
  Ruddd->get(3, 0, 3, 0) = -((M*(2*M - r))/pow(r,4));
  Ruddd->get(3, 1, 1, 3) = -(M/((2*M - r)*pow(r,2)));
  Ruddd->get(3, 1, 3, 1) = M/((2*M - r)*pow(r,2));
  Ruddd->get(3, 2, 2, 3) = (-2*M)/r;
  Ruddd->get(3, 2, 3, 2) = (2*M)/r;
}

void Schwarzschild::calc_Rudddcd()
{
  Rudddcd->get(0, 1, 0, 1, 1) = (6*M)/((2*M - r)*pow(r,3));
  Rudddcd->get(0, 1, 0, 2, 2) = (3*M)/pow(r,2);
  Rudddcd->get(0, 1, 0, 3, 3) = (3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(0, 1, 1, 0, 1) = (-6*M)/((2*M - r)*pow(r,3));
  Rudddcd->get(0, 1, 2, 0, 2) = (-3*M)/pow(r,2);
  Rudddcd->get(0, 1, 3, 0, 3) = (-3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(0, 2, 0, 1, 2) = (3*M)/pow(r,2);
  Rudddcd->get(0, 2, 0, 2, 1) = (3*M)/pow(r,2);
  Rudddcd->get(0, 2, 1, 0, 2) = (-3*M)/pow(r,2);
  Rudddcd->get(0, 2, 2, 0, 1) = (-3*M)/pow(r,2);
  Rudddcd->get(0, 3, 0, 1, 3) = (3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(0, 3, 0, 3, 1) = (3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(0, 3, 1, 0, 3) = (-3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(0, 3, 3, 0, 1) = (-3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(1, 0, 0, 1, 1) = (6*M*(2*M - r))/pow(r,5);
  Rudddcd->get(1, 0, 0, 2, 2) = (3*M*pow(-2*M + r,2))/pow(r,4);
  Rudddcd->get(1, 0, 0, 3, 3) = (3*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,4);
  Rudddcd->get(1, 0, 1, 0, 1) = (6*M*(-2*M + r))/pow(r,5);
  Rudddcd->get(1, 0, 2, 0, 2) = (-3*M*pow(-2*M + r,2))/pow(r,4);
  Rudddcd->get(1, 0, 3, 0, 3) = (-3*M*pow(-2*M + r,2)*pow(sin(theta),2))/pow(r,4);
  Rudddcd->get(1, 2, 1, 2, 1) = (3*M)/pow(r,2);
  Rudddcd->get(1, 2, 2, 1, 1) = (-3*M)/pow(r,2);
  Rudddcd->get(1, 2, 2, 3, 3) = (-3*M*(2*M - r)*pow(sin(theta),2))/r;
  Rudddcd->get(1, 2, 3, 2, 3) = (3*M*(2*M - r)*pow(sin(theta),2))/r;
  Rudddcd->get(1, 3, 1, 3, 1) = (3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(1, 3, 2, 3, 2) = (3*M*(2*M - r)*pow(sin(theta),2))/r;
  Rudddcd->get(1, 3, 3, 1, 1) = (-3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(1, 3, 3, 2, 2) = (-3*M*(2*M - r)*pow(sin(theta),2))/r;
  Rudddcd->get(2, 0, 0, 1, 2) = (3*M*(-2*M + r))/pow(r,5);
  Rudddcd->get(2, 0, 0, 2, 1) = (3*M*(-2*M + r))/pow(r,5);
  Rudddcd->get(2, 0, 1, 0, 2) = (3*M*(2*M - r))/pow(r,5);
  Rudddcd->get(2, 0, 2, 0, 1) = (3*M*(2*M - r))/pow(r,5);
  Rudddcd->get(2, 1, 1, 2, 1) = (3*M)/((2*M - r)*pow(r,3));
  Rudddcd->get(2, 1, 2, 1, 1) = (-3*M)/((2*M - r)*pow(r,3));
  Rudddcd->get(2, 1, 2, 3, 3) = (-3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(2, 1, 3, 2, 3) = (3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(2, 3, 1, 2, 3) = (3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(2, 3, 1, 3, 2) = (-3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(2, 3, 2, 1, 3) = (-3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(2, 3, 2, 3, 1) = (-6*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(2, 3, 3, 1, 2) = (3*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(2, 3, 3, 2, 1) = (6*M*pow(sin(theta),2))/pow(r,2);
  Rudddcd->get(3, 0, 0, 1, 3) = (3*M*(-2*M + r))/pow(r,5);
  Rudddcd->get(3, 0, 0, 3, 1) = (3*M*(-2*M + r))/pow(r,5);
  Rudddcd->get(3, 0, 1, 0, 3) = (3*M*(2*M - r))/pow(r,5);
  Rudddcd->get(3, 0, 3, 0, 1) = (3*M*(2*M - r))/pow(r,5);
  Rudddcd->get(3, 1, 1, 3, 1) = (3*M)/((2*M - r)*pow(r,3));
  Rudddcd->get(3, 1, 2, 3, 2) = (3*M)/pow(r,2);
  Rudddcd->get(3, 1, 3, 1, 1) = (-3*M)/((2*M - r)*pow(r,3));
  Rudddcd->get(3, 1, 3, 2, 2) = (-3*M)/pow(r,2);
  Rudddcd->get(3, 2, 1, 2, 3) = (-3*M)/pow(r,2);
  Rudddcd->get(3, 2, 1, 3, 2) = (3*M)/pow(r,2);
  Rudddcd->get(3, 2, 2, 1, 3) = (3*M)/pow(r,2);
  Rudddcd->get(3, 2, 2, 3, 1) = (6*M)/pow(r,2);
  Rudddcd->get(3, 2, 3, 1, 2) = (-3*M)/pow(r,2);
  Rudddcd->get(3, 2, 3, 2, 1) = (-6*M)/pow(r,2);
}

void Schwarzschild::calc_R()
{
  R->get(0) = 0;
}

/* Covariant derivative of Ricci scalar contracted with 4-velocity */
void Schwarzschild::calc_Rcd()
{
  Rcd->get(0) = 0;
  Rcd->get(1) = 0;
  Rcd->get(2) = 0;
  Rcd->get(3) = 0;
}

