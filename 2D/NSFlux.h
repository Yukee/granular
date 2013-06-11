#ifndef NSFLUX_H
#define NSFLUX_H

#include "Flux.h"

/* Incompressible Navier-Stokes flux, used in the incompressible NS equations system. As many solved dimensions as space dimensions. */

class NSFlux : public Flux
{
 public:

 NSFlux(int i) : Flux(i,i){}

  inline VectorField evaluate(const VectorField & u, const int d)
  {
    for(int i=0;i<m_solved_dimensions;i++) m_evaluated_flux_d[i] = u[i]*u[d];
    return m_evaluated_flux_d;
  }

  inline SField get_max_eigenvalue(const VectorField & u, const int d)
  {
    // max eigenvalue is 2*u, we want the module.
    m_max_eigenvalue = 2*u[d].module();
    return m_max_eigenvalue;
  }

  inline int delta(int i,int j)
  {
    int delta = 0;
    if(i == j) delta = 1;
    return delta;
  }

};

#endif
