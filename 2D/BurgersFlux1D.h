#ifndef BURGERSFLUX1D_H
#define BURGERSFLUX1D_H

#include "Flux.h"

class Burgers1D : public Flux
{
 public:

 Burgers1D() : Flux(1,1){}

  inline TensorField evaluate(VectorField u)
  {
    m_evaluated_flux[0][0] = 0.5*u[0]*u[0];
    return m_evaluated_flux;
  }

  inline Vector<TensorField> evaluate_jacobian(VectorField u)
  {
    m_evaluated_flux_jacobian[0][0][0] = u[0];
    return m_evaluated_flux_jacobian;
  }

  inline bool has_exact_jacobian()
  {
    return true;
  }
};

#endif
