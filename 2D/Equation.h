#ifndef EQUATION_H
#define	EQUATION_H

#include "Vector.h"
#include "ScalarField.h"
#include "Flux.h"
#include "ZeroFlux.h"

typedef Vector<ScalarField> VectorField;
typedef Vector<VectorField> TensorField;

class Equation
{
protected:
  Flux *m_conv; //convection flux
  Flux *m_diff; //diffusion flux

public:
  Equation() {m_conv =  new ZeroFlux(1,1); m_diff = new ZeroFlux(1,1);}
 Equation(Flux *f) : m_conv(f) {m_diff = new ZeroFlux(f->get_space_dimensions(),f->get_solved_dimensions());}
  Equation(Flux *f, Flux *diff) : m_conv(f), m_diff(diff){}
  ~Equation();
  TensorField get_convectionFlux(VectorField);
  Vector<TensorField> get_convectionFluxJacobian(VectorField);
  TensorField get_diffusionFlux(VectorField);
};

#endif


