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
#include "FD1Solver.h"
#include "timeSolver.h"
#include "EulerSolver.h"
//#include "RK2Solver.h"
//#include "RK3Solver.h"
#include "PeriodicField.h"
#include "PeriodicField.h"
#include "NullField.h"

using namespace std;
int main()
{
  int dim = 1;
  Vector<double> dx(dim,0.1); Vector<double> xI(dim,2*M_PI); Vector<double> lfc(dim,0);
  Flux *bf = new NSFlux(dim);
  Equation *eq = new Equation(bf);
  FD1Solver sol(dx, xI, eq, lfc);

  double dt = 0.01; double T = 10;
  Vector<int> xr (dim); for(int d=0;d<dim;d++) xr[d] = xI[d]/dx[d];
  VectorField pos = sol.get_position();
  VectorField u0 (dim, PeriodicField (xr) );
  for(int d=0;d<dim;d++) for(int it=0;it<u0[d].get_size();++it)
  			   u0[d][it] = 0.5 + cos(pos[0][it]);
   EulerSolver tse(dt, T, &sol, u0);

   tse.get_solution("burgers",0.01);

  return 0;
}
