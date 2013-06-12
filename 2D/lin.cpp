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
#include "LinFlux.h"
#include "LinNullFlux.h"
#include "Lin2DFlux.h"

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
#include "ImmersedField.h"
#include "WriteVectorField.h"

// Flume pb
//#include "Flume3D.h"
//#include "Flume2D.h"

using namespace std;
int main()
{
  int dim = 2;
  int syst = 1;
  Vector<double> dx(dim,0.1); Vector<double> xI(dim,2*M_PI); Vector<double> llc(dim,0);
  Flux *ptrf = new Lin2DFlux();
  Equation *eq = new Equation(ptrf);
  FD1Solver sol(dx, xI, eq, llc);

  Vector<int> xr (dim); for(int d=0;d<dim;d++) xr[d] = xI[d]/dx[d];
  VectorField pos = sol.get_position();

  // sets the propagation speed of the travelling wave
  SField xpropag_speed(xr); xpropag_speed = 1;
  ptrf->set_parameter(xpropag_speed);
  
  // initial condition
  VectorField u0(syst, SField(xr));
  for(int it=0;it<u0[0].get_size();++it) {
    u0[0][it] = sin( pos[0][it] ) + sin( pos[1][it] );
  }

  cout << "space dim:" << eq->get_space_dimensions() << endl;
  cout << "solved dim:" << eq->get_solved_dimensions() << endl;

  double dt = 0.005; double T = 10;
  EulerSolver ts(dt, T, &sol, u0);
  ts.get_solution("linear2D",1);

  return 0;
}
