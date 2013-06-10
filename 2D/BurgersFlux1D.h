#ifndef BURGERSFLUX1D_H
#define BURGERSFLUX1D_H

#include "Flux.h"

class Burgers1D : public Flux
{
 public:

 Burgers1D() : Flux(1,1){}

  inline VectorField evaluate(VectorField u, int d)
  {
    m_evaluated_flux_d[0] = 0.5*u[0]*u[0];
    return m_evaluated_flux_d;
  }

  inline SField get_max_eigenvalue(VectorField u, int d)
  {
    return u[0];
  }
};

#endif
