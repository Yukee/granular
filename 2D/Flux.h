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
  VectorField m_evaluated_flux_d; // m_evaluated_flux[d]
  Vector<TensorField> m_evaluated_flux_jacobian;

  // no need to know each eigenvalue for our pb, but perhaps useful in a more general case!
  TensorField m_eigenvalues;

  // The positions on the numerical grid your flux is evaluated on; useless if the flux is not explicitely dependant on the position (which is the case if your equation is quasilinear)
  VectorField m_positions;

 public:

  // by default space and solved dimensions are 1
  Flux();
  Flux(int space_dimensions, int solved_dimensions);

  // returns the flux evaluated at gridpoints
  virtual TensorField evaluate(VectorField)=0;

  //returns the restriction to 1 space dimension of the flux evaluated at gridpoints 
  virtual VectorField evaluate(VectorField, int)=0;

  // tells if the exact jacobian of the flux is specified or not
  virtual bool has_exact_jacobian()=0;

  // returns the jac evaluated at gridpoints
  virtual Vector<TensorField> evaluate_jacobian(VectorField)=0;

  // returns the maximal eigenvalue (in module) of the restriction to 1 space dimension of the jacobian
  virtual ScalarField get_max_eigenvalue(VectorField, int)=0;

  int get_space_dimensions();
  int get_solved_dimensions();
};

#endif
