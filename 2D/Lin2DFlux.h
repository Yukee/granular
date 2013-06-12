#ifndef LIN2DFLUX_H
#define LIN2DFLUX_H

#include "Flux.h"

/* Linear flux, used to test the travelling wave solution in 2D case. */

class Lin2DFlux: public Flux
{
 private:
  SField m_xspeed;
  SField m_yspeed;

 public:

 Lin2DFlux() : Flux(2,1){ 
    //sets speed of the travelling wave to 1 by default
    m_param = 1;}

  inline VectorField evaluate(const VectorField & u, const int d)
  {
    switch(d){
    case 0:
      m_evaluated_flux_d[0] = m_xspeed*u[0];
      break;

    case 1:
      m_evaluated_flux_d[0] = m_yspeed*u[0];
      break;

    default:
      throw std::invalid_argument("In Lin2DFlux::evaluate: d must be 0 or 1.");
    }
    return m_evaluated_flux_d;
  }

  inline SField get_max_eigenvalue(const VectorField & u, const int d)
  {
    switch(d){
    case 0:
      m_max_eigenvalue = m_xspeed;
      break;

    case 1:
      m_max_eigenvalue = m_yspeed;
      break;

    default:
      throw std::invalid_argument("In Lin2DFlux::get_max_eigenvalue: d must be 0 or 1.");
    }
    return m_max_eigenvalue;
  }

  inline void set_parameter(const SField & v)
  {
    m_xspeed = v;
    m_yspeed = 2*v;
  }

};

#endif
