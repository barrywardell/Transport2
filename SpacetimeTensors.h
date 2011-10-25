/* A class containing useful tensors for a spacetime.
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
#ifndef SPACETIMETENSORS_H
#define SPACETIMETENSORS_H

#include "Tensor.h"
#include "TensorList.h"

using namespace Mosquito;

enum GeodesicType {
  Timelike  = -1,
  Null      = 0,
  Spacelike = 1
};

class Spacetime : public TensorList {
  public:
    virtual void calc_all() {};
    virtual void calc_guu() {};
    virtual void calc_gdd() {};
    virtual void calc_Gudd() {};
    virtual void calc_Ruddd() {};
    virtual void calc_Rudddcd() {};
    virtual void calc_R() {};
    virtual void calc_Rcd() {};

    virtual void setPoint(const Tensor& t) {};

    Spacetime();
    virtual ~Spacetime() {};
};


class Schwarzschild : public Spacetime {
  public:
    Schwarzschild(double mass = 1.0);

    /* Functions to compute tensor components */
    void calc_all();
    void calc_guu();
    void calc_gdd();
    void calc_Gudd();
    void calc_Ruddd();
    void calc_Rudddcd();
    void calc_Rudddcdd();
    void calc_R();
    void calc_Rcd();

    void setPoint(const Tensor& t);

    /* Spacetime parameters */
    const double M; /* Black Hole Mass */

  private:
    /* Spacetime coordinates */
    double t, r, theta, phi;

    double sintheta;
    double costheta;

};

struct geodesic_params {
  double e;           /* "Energy" constant of motion */
  double l;           /* "Angular momentum" constant of motion */
  int type;           /* Type of geodesic. 0=null, -1=time-like */
  TensorList& T;      /* List of evolved tensors */
  Schwarzschild& s;   /* Spacetime */
};

#endif
