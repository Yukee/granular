#ifndef FLUX_H
#define FLUX_H

#include "Vector.h"
#include "ScalarField.h"

typedef Vector<ScalarField> VectorField;
typedef Vector<VectorField> TensorField;

class Flux
{
 protected:
  int m_space_dimensions;
  int m_solved_dimensions;
  TensorField m_evaluated_flux;
  Vector<TensorField> m_evaluated_flux_jacobian;
  //The positions on the numerical grid your flux is evaluated on; useless if the flux is not explicitely dependant on the position (which is the case if your equation is quasilinear)
  VectorField m_positions;

 public:
  Flux();
  Flux(int space_dimensions, int solved_dimensions);
  virtual TensorField evaluate(VectorField)=0; 
  virtual bool has_exact_jacobian()=0;
  virtual Vector<TensorField> evaluate_jacobian(VectorField)=0;
  int get_space_dimensions();
  int get_solved_dimensions();
};

#endif
