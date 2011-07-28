#include "SpacetimeTensors.h"
#include "Tensor.h"
#include <math.h>
#include <stdio.h>

int main(int argc, char* argv[])
{
  Tensor::IndexType u = Tensor::UP;
  Tensor::IndexType d = Tensor::DOWN;

  Schwarzschild s;
  s.M     = 1.0;
  s.t     = 0;
  s.r     = 10.0;
  s.theta = M_PI_2;
  s.phi   = 0;
  s.calc_all();

  Tensor Q(2, u, d);
  Tensor u1(1, u);

  u1.get(0) = sqrt(s.r/(s.r-3*s.M));
  u1.get(1) = 0;
  u1.get(2) = 0;
  u1.get(3) = sqrt(s.M/(s.r-3*s.M))/s.r;

  Tensor G = *s.Gudd;
  Tensor R = *s.Ruddd;
  Tensor Qrhs = (Q('a', 'd')*G('d', 'c', 'b') - G('a', 'c', 'd')*Q('c', 'b'))*u1('c') 
  - Q('a', 'b') - Q('a','c')*Q('c','b') - R('a','c','b','d')*u1('c')*u1('d');
  printf("%g\n", Qrhs.get(1,1));
}