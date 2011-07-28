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

struct geodesic_params {
  double e; /* "Energy" constant of motion */
  double l; /* "Angular momentum" constant of motion */
  int type; /* Type of geodesic. 0=null, -1=time-like */
};

class Spacetime {
  public:
    virtual void calc_all() {};
    virtual void calc_guu() {};
    virtual void calc_gdd() {};
    virtual void calc_Gudd() {};
    virtual void calc_Ruddd() {};
    virtual void calc_Rudddcd() {};
    virtual void calc_R() {};
    virtual void calc_Rcd() {};
    virtual ~Spacetime() {};
    
    Tensor* guu;     /* Contravariant metric */
    Tensor* gdd;     /* Covariant metric */
    Tensor* Gudd;    /* Christoffel symbol */
    Tensor* Ruddd;   /* Riemann */
    Tensor* Rudddcd; /* Covariant derivative of Riemann */
    Tensor* R;       /* Ricci scalar */
    Tensor* Rcd;     /* Covariant derivative of Ricci scalar */
};


class Schwarzschild : public Spacetime {
  public:
    /* Spacetime parameters */
    double M; /* Black Hole Mass */
    
    /* Spacetime coordinates */
    double t, r, theta, phi;
    double sintheta;
    double costheta;

    /* Functions to compute tensor components */
    void calc_all();
    void calc_guu();
    void calc_gdd();
    void calc_Gudd();
    void calc_Ruddd();
    void calc_Rudddcd();
    void calc_R();
    void calc_Rcd();

    Schwarzschild();
    virtual ~Schwarzschild() {};
};

#endif
