#ifndef LINNULLFLUX_H
#define LINNULLFLUX_H

#include "Flux.h"

/* Linear flux, used to test the travelling wave solution. Special case to test the behaviour of the solver when one flux is zero. */

class LinNullFlux : public Flux 
{
 public:

 LinNullFlux() : Flux(1,2){ 
    //sets speed of the travelling wave to 1 by default
    m_param = 1;}

  inline VectorField evaluate(const VectorField & u, const int d)
  {
    m_evaluated_flux_d[0] = m_param*u[0];
    m_evaluated_flux_d[1] = 0*u[1];
    return m_evaluated_flux_d;
  }

  inline SField get_max_eigenvalue(const VectorField & u, const int d)
  {
    // max eigenvalue is 1
    m_max_eigenvalue.resize_field(u[d].get_range());
    m_max_eigenvalue = 1;
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
