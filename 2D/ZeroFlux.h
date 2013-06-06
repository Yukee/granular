#ifndef ZEROFLUX_H
#define ZEROFLUX_H

#include "Flux.h"

class ZeroFlux : public Flux
{
 public:

 ZeroFlux(int d, int i) : Flux(d,i){}

  inline TensorField evaluate(VectorField u)
  {
    for(int d=0;d<m_space_dimensions;d++) for(int i=0;i<m_solved_dimensions;i++) m_evaluated_flux[d][i] = 0*u[i];
    return m_evaluated_flux;
  }

  inline VectorField evaluate(VectorField u, int d)
  {
    for(int i=0;i<m_solved_dimensions;i++) m_evaluated_flux_d[i] = 0*u[i];
    return m_evaluated_flux_d;
  }

  inline ScalarField get_max_eigenvalue(VectorField u, int d)
  {
    return 0*u[0];
  }

  inline Vector<TensorField> evaluate_jacobian(VectorField u)
  {
    for(int d=0;d<m_space_dimensions;d++) for(int i=0;i<m_solved_dimensions;i++) for(int j=0;j<m_solved_dimensions;j++) m_evaluated_flux_jacobian[d][i][j] = 0*u[i];
    return m_evaluated_flux_jacobian;
  }

  inline bool has_exact_jacobian()
  {
    return true;
  }
};

#endif
