#ifndef FLUX_H
#define FLUX_H

#include "Vector.h"
#include "ScalarField.h"
#include "PeriodicField.h"
#include "NullField.h"
#include "PrescribedField.h"

typedef PrescribedField SField;
typedef Vector<SField> VectorField;
typedef Vector<VectorField> TensorField;

class Flux
{
 protected:
  int m_space_dimensions;
  int m_solved_dimensions;
  VectorField m_evaluated_flux_d; // m_evaluated_flux[d]
  SField m_max_eigenvalue;

  // The positions on the numerical grid your flux is evaluated on; useless if the flux is not explicitely dependant on the position (which is the case if your equation is quasilinear)
  VectorField m_positions;

 public:

  // by default space and solved dimensions are 1
  Flux();
  Flux(int space_dimensions, int solved_dimensions);

  //returns the restriction to 1 space dimension of the flux evaluated at gridpoints 
  virtual VectorField evaluate(VectorField, int)=0;

  // returns the maximal eigenvalue (in module) of the restriction to 1 space dimension of the jacobian
  virtual SField get_max_eigenvalue(VectorField, int)=0;

  int get_space_dimensions();
  int get_solved_dimensions();
};

#endif
