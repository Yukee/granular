#include <iostream>
#include <fstream>
#include <math.h>
#include <limits.h>
#include "Vector.h"
#include "Flux.h"
#include "ZeroFlux.h"
#include "BurgersFlux1D.h"
#include "NSFlux.h"
#include "Equation.h"
//#include "FD1Solver.h"
//#include "timeSolver.h"
//#include "EulerSolver.h"
//#include "RK2Solver.h"
//#include "RK3Solver.h"
#include "ScalarField.h"

using namespace std;
int main()
{
  // tests in a 42D space, with a system of 3 equations (ie a 3D vector field)
  Vector<int> R(42,1);
  ScalarField s(R);
  VectorField v(3,s);
  Flux *zero = new ZeroFlux(42,3);
  for(int d=0;d<42;d++) cout << zero->evaluate(v)[0][d] << endl;
  
  // Navier-Stokes in a 3D space
  int d = 3;
  Vector<int> x(d,1);
  VectorField u(d, ScalarField (x));
  for(int i=0;i<d;i++) u[i](Vector<int> (d,0)) = 2;
  u[2](Vector<int> (d,0)) = 1;
  Flux *f = new NSFlux(d);
  for(int i=0;i<d;i++) for(int j=0;j<d;j++) cout << f->evaluate(u)[i][j] << endl;
  for(int dim=0;dim<d;dim++) 
    {
      cout << "Jacobian of the flux " << dim << ":" << endl;
      for(int i=0;i<d;i++) for(int j=0;j<d;j++) cout << f->evaluate_jacobian(u)[dim][i][j] << endl;
    }

  // Navier-Stokes equations without forces and pressure gradient (= conservation of momentum)
  Equation NSEq(f);
  for(int dim=0;dim<d;dim++) 
    {
      cout << "Jacobian of the flux " << dim << ":" << endl;
      for(int i=0;i<d;i++) for(int j=0;j<d;j++) cout <<NSEq.get_diffusionFlux(u)[i][j]<< endl;
    }

  return 0;
}
