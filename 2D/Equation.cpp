#include "Equation.h"

Equation::~Equation()
{
  if(m_conv) delete[] m_conv;
  if(m_diff) delete[] m_diff;
}

TensorField Equation::get_convectionFlux(VectorField u)
{
  return m_conv->evaluate(u);
}

Vector<TensorField> Equation::get_convectionFluxJacobian(VectorField u)
{
  if(!m_conv->has_exact_jacobian()) std::cout << "In Equation::get_convectionFluxJacobian: the exact jacobian of the convetion flux is unknown." << std::endl;
  return m_conv->evaluate_jacobian(u);
}

TensorField Equation::get_diffusionFlux(VectorField u)
{
  return m_diff->evaluate(u);
}
