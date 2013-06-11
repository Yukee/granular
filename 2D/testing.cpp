#include <iostream>
#include <fstream>
#include <math.h>
#include <limits.h>
#include "Vector.h"
#include "Flux.h"
#include "ZeroFlux.h"
#include "BurgersFlux1D.h"
#include "NSFlux.h"
#include "Flume2DConvectionFlux.h"
#include "Equation.h"
#include "FD1Solver.h"
#include "timeSolver.h"
#include "EulerSolver.h"
//#include "RK2Solver.h"
//#include "RK3Solver.h"
#include "ScalarField.h"
#include "PeriodicField.h"
#include "NullField.h"
#include "PrescribedField.h"
#include "WriteVectorField.h"

// Flume pb
//#include "Flume3D.h"
#include "Flume2D.h"

using namespace std;
int main()
{
  int dim = 1;
  Vector<double> dx(dim,0.1); Vector<double> xI(dim); xI[0]=2*M_PI; Vector<double> llc(dim,0);
  Flux *ptrf = new NSFlux(dim);
  Equation *eq = new Equation(ptrf);
  FD1Solver sol(dx, xI, eq, llc);

  Vector<int> xr (dim); for(int d=0;d<dim;d++) xr[d] = xI[d]/dx[d];
  VectorField pos = sol.get_position();

  VectorField u0(dim, SField(xr));
  for(int it=0;it<u0[0].get_size();++it) u0[0][it] = 0.5 + sin( pos[0][it] );

  double dt = 0.01; double T = 10;
  EulerSolver ts(dt, T, &sol, u0);
  ts.get_solution("burgers",0.01);

  return 0;
}
