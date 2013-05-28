#ifndef EQUATION_H
#define	EQUATION_H

#include "Vector.h"
#include "ScalarField.h"
#include "Flux.h"

typedef Vector<ScalarField> VectorField;
typedef Vector<VectorField> TensorField;

class Equation
{
protected:
  Flux m_conv; //convection flux
  Flux m_diff; //diffusion flux

public:
  Equation();
 Equation(Flux f, Flux diff) : m_conv(f), m_diff(diff){}
  TensorField get_convectionFlux(VectorField);
  Vector<TensorField> get_convectionFluxJacobian(VectorField);
  TensorField get_diffusionFlux(VectorField);
};

#endif


