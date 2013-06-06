#ifndef NSFLUX_H
#define NSFLUX_H

#include "Flux.h"

/* Incompressible Navier-Stokes flux, used in the incompressible NS equations system. As many solved dimensions as space dimensions. */

class NSFlux : public Flux
{
 public:

 NSFlux(int i) : Flux(i,i){}

  inline TensorField evaluate(VectorField u)
  {
    for(int i=0;i<m_solved_dimensions;i++) for(int d=0;d<m_space_dimensions;d++) m_evaluated_flux[d][i] = u[i]*u[d];
    return m_evaluated_flux;
  }

  inline VectorField evaluate(VectorField u, int d)
  {
    for(int i=0;i<m_solved_dimensions;i++) m_evaluated_flux_d[i] = u[i]*u[d];
    return m_evaluated_flux_d;
  }

  inline Vector<TensorField> evaluate_jacobian(VectorField u)
  {
    for(int d=0;d<m_space_dimensions;d++) for(int i=0;i<m_solved_dimensions;i++) for(int j=0;j<m_solved_dimensions;j++) m_evaluated_flux_jacobian[d][i][j] = delta(i,j)*u[d] + delta(d,j)*u[i];
    return m_evaluated_flux_jacobian;
  }

  inline PeriodicField get_max_eigenvalue(VectorField u, int d)
  {
    // max eigenvalue is 2*u, we want the module.
    return 2*u[d].module();
  }

  inline int delta(int i,int j)
  {
    int delta = 0;
    if(i == j) delta = 1;
    return delta;
  }

  inline bool has_exact_jacobian()
  {
    return true;
  }
};

#endif
