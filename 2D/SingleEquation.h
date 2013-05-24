#ifndef SINGLEEQ_H
#define SINGLEEQ_H

#include "Equation.h"

class SingleEquation : public Equation
{
  public:
  SingleEquation(convectionFlux, convectionFluxJacobian);
  inline Vector<ScalarField> get_convectionFlux(ScalarField u)
  {
      return m_f(u);
  }
  inline Vector<ScalarField> get_convectionFluxJacobian(ScalarField u)
  {
      return m_df(u);
  }
  inline Vector<ScalarField> get_speed_field()
  {
      std::cout << "no speed field in this configuration" << std::endl;
      return Vector<ScalarField>(0);
  }
  inline void set_speed_field(Vector<ScalarField>, ScalarField)
  {
      std::cout << "no speed field in this configuration" << std::endl;
  }
};

#endif
