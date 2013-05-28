#include "Equation.h"

Equation::Equation()
{
}

Vector<VectorField> Equation::get_convectionFlux(VectorField u)
{
  return m_conv.evaluate(u);
}

Vector< Vector<VectorField> > Equation::get_convectionFluxJacobian(VectorField u)
{
  if(!m_conv.has_exact_jacobian()) std::cout << "In Equation::get_convectionFluxJacobian: the exact jacobian of the convetion flux is unknown." << std::endl;
  return m_conv_jac.evaluate_jacobian(u);
}

Vector<VectorField> Equation::get_diffusionFlux(VectorField u)
{
  return m_diff.evaluate(u);
}
