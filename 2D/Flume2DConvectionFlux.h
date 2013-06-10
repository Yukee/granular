#ifndef FLUMEFLUX_H
#define FLUMEFLUX_H

#include "Flux.h"

class Flume2DConvectionFlux: public Flux
{
 private:

  SField m_segregation_rate;

 public:

 Flume2DConvectionFlux(): Flux(2,3) {m_segregation_rate = 0;}

  inline void set_segregation_rate(SField sr)
  {
    m_segregation_rate = sr;
  }

  inline virtual VectorField evaluate(VectorField u, int d)
  {
    switch(d){
    case 0:
      m_evaluated_flux_d[0] = u[0]*u[1];
      m_evaluated_flux_d[1] = 0;
      m_evaluated_flux_d[2] = 0;
      break;

    case 1:
      m_evaluated_flux_d[0] = u[0]*u[2] - m_segregation_rate*u[0]*(1-u[0]);
      m_evaluated_flux_d[1] = 0;
      m_evaluated_flux_d[2] = 0;
      break;

    default:
      throw std::invalid_argument("In FlumeFlux::evaluate dimension must either be 0 or 1");
    }
    return m_evaluated_flux_d;
  }

  inline virtual SField get_max_eigenvalue(VectorField u, int d)
  {
    switch(d){
    case 0:
      m_max_eigenvalue = u[1].module();
      break;

    case 1:
      m_max_eigenvalue = u[2] - m_segregation_rate*(1-2*u[0]);
      break;

    default:
      throw std::invalid_argument("In FlumeFlux::get_max_eigenvalue dimension must either be 0 or 1");
    }
    return m_max_eigenvalue;
  }

};

#endif
