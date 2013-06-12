#ifndef LINFLUX_H
#define LINFLUX_H

#include "Flux.h"

/* Linear flux, used to test the travelling wave solution */

class LinFlux : public Flux
{
 public:

 LinFlux(int d, int i) : Flux(d,i){ 
    //sets speed of the travelling wave to 1 by default
    m_param = 1;}

  inline VectorField evaluate(const VectorField & u, const int d)
  {
    for(int i=0;i<m_solved_dimensions;i++) m_evaluated_flux_d[i] = m_param*u[i];
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
