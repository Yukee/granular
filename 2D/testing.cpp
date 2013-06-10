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

#include <vector>

using namespace std;
int main()
{
  int dim = 2;
  Vector<double> dx(dim,0.01); Vector<double> xI(dim); xI[0]=3; xI[1]=1; Vector<double> lfc(dim,0); lfc[0]=-xI[0];
  Flux *ptrf = new Flume2DConvectionFlux();
  Equation *eq = new Equation(ptrf); //don't forget to set the sr somewhere!!
  FD1Solver sol(dx, xI, eq, lfc);

  // double dt = 0.01; double T = 10;
  Vector<int> xr (dim); for(int d=0;d<dim;d++) xr[d] = xI[d]/dx[d];
  VectorField pos = sol.get_position();
  // VectorField u0 (dim, SField (xr) );
  // cout << "u0 bounds : " << u0[0].get_bounds()[0].get_range() << endl;

  // for(int d=0;d<dim;d++) for(int it=0;it<u0[d].get_size();++it)
  // 			   u0[d][it] = 0.5 + cos(pos[0][it]);
  // cout << "u0 bounds : " << u0[0].get_bounds()[0].get_range() << endl;
  // EulerSolver tse(dt, T, &sol, u0);

  // tse.get_solution("NS2D",0.01);

  // Sets initial fields and writes the data in files in the Results/Flume2D_initial/ folder
  SField phi(xr);
  SField bound(xr);
  VectorField u0(dim, SField (xr));
  for(int it=0;it<phi.get_size();++it)
    {
      phi[it] = phi0(pos[0][it], pos[1][it]);
      bound[it] = boundary(pos[0][it], pos[1][it]);
      u0[0][it] = u(pos[0][it], pos[1][it]);
      u0[1][it] = w(pos[0][it], pos[1][it]);
    }
  // the velocity field is null outside the domain, as well as the sr
  u0 = bound*u0;
  SField sr = 1.3*bound; (*ptrf).set_segregation_rate(sr);
  // the actual solved field is a vector comtaining the velocity and the concentration fields
  VectorField uInit (3); uInit[0] = phi; uInit[1] = u0[0]; uInit[2] = u0[1];
  
  fstream init;
  init.open("Results/Flume2D_initial/phi0.tsv",ios::out);
  phi.write_in_file_matrixform(init);
  fstream speed;
  speed.open("Results/Flume2D_initial/speed.tsv",ios::out);
  write_VectorField(u0, pos, speed);
  fstream b;
  b.open("Results/Flume2D_initial/boundary.tsv",ios::out);
  bound.write_in_file_matrixform(b);

  fstream f;
  f.open("Results/Flume2D_initial/convection_flux.tsv",ios::out);
  VectorField jac(2);
  jac[0] = ptrf->get_max_eigenvalue(uInit, 0);
  jac[1] = ptrf->get_max_eigenvalue(uInit, 1);
  write_VectorField(jac, pos, f);

  return 0;
}
