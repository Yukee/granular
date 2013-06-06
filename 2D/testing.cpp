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
  // tests in a 42D space, with a system of 3 equations (ie a 3D vector field)
  Vector<int> R(42,1);
  PeriodicField s(R);
  VectorField v(3,s);
  Flux *zero = new ZeroFlux(42,3);
  for(int d=0;d<42;d++) cout << (zero->evaluate(v)[d] == zero->evaluate(v,d) ) << endl;

  // Navier-Stokes in a 3D space
  int d = 3;
  Vector<int> x(d,1);
  VectorField u(d, PeriodicField (x));
  for(int i=0;i<d;i++) u[i](Vector<int> (d,0)) = 2;
  u[2](Vector<int> (d,0)) = 1;
  Flux *f = new NSFlux(d);
  for(int i=0;i<d;i++) for(int j=0;j<d;j++) cout << ( f->evaluate(u)[i][j] == f->evaluate(u,i)[j] ) << endl;
  for(int dim=0;dim<d;dim++) 
    {
      cout << "Jacobian of the flux " << dim << ":" << endl;
      for(int i=0;i<d;i++) for(int j=0;j<d;j++) cout << f->evaluate_jacobian(u)[dim][i][j] << endl;
    }

  // Navier-Stokes equations without forces and pressure gradient (= conservation of momentum)
  Equation NSEq(f);
  for(int dim=0;dim<d;dim++) 
    {
      for(int i=0;i<d;i++) for(int j=0;j<d;j++) cout <<NSEq.get_diffusionFlux(u,i)[j]<< endl;
    }

  PeriodicField a(Vector<int> (2,2));
  Vector<int> y(2); y[0] = 1; y[1] = 0;
  Vector<int> z(2); z[0] = 0; z[1] = 1;
  for(int i=0;i<2;i++) for(int j=0;j<2;j++) a(i*y+j*z) = i+j;
  cout << a;
  for(int d=0;d<a.get_size();d++)
    {
      a[d] = a.get_pos(d)[0] + a.get_pos(d)[1];
    }
  cout << a;

  cout << "********Periodic and null fields test area************" << endl;
  ScalarField *ptr = 0;

  PeriodicField pf(Vector<int> (3,2));
  ptr = &pf;
  for(int it=0;it<ptr->get_size();++it) (*ptr)[it] = ptr->get_pos(it)[0] + ptr->get_pos(it)[1] + ptr->get_pos(it)[2];
  cout << (*ptr) << endl;
  Vector< Vector<int> > b = (Vector<int> (3)).get_base_vectors(1,0);
  int n = 3;
  for(int i=0;i<n;i++) for(int j=0;j<n;j++) for(int k=0;k<n;k++){
	cout << i << j << k << (*ptr)(i*b[0]+j*b[1]+k*b[2]) << endl; // no runtime error: the periodic field re-ranges the position
      }

  NullField nf(Vector<int> (2,2));
  ptr = &nf;
  for(int it=0;it<ptr->get_size();++it) (*ptr)[it] = 1;
  cout << (*ptr)(Vector<int> (2,2)); // no runtime error: the null field is 0 outside its range

  PeriodicField sf(Vector<int> (2,2));
  ptr = &sf;
  //cout << (*ptr)(Vector<int> (2,2)) << endl; // runtime error: you are outside the range of the field

  cout << "******** Initialization test area ****************" << endl;
  int m_m = 3;
  int m_n = 4;
  VectorField vect(m_m, PeriodicField (Vector<int> (3,2)));
  cout << vect.size() << "\t" << vect[2].get_range() << endl;

  TensorField tens(m_m, VectorField (m_n, PeriodicField (Vector<int> (4,2))));
  cout << tens.size() << "\t" << tens[0].size() << "\t" << tens[0][0].get_range() << endl;

  VectorField vect2;
  vect2 =  VectorField (m_m, PeriodicField (Vector<int> (3,2)));
  cout << vect.size() << "\t" << vect[2].get_range() << endl;


  cout << "****************** Solver test area *******************" << endl;
  int dim = 2;
  Vector<double> dx(dim,0.1); Vector<double> xI(dim,2*M_PI); Vector<double> lfc(dim,0);
  Flux *bf = new NSFlux(dim);
  Equation *eq = new Equation(bf);
  FD1Solver sol(dx, xI, eq, lfc);

  double dt = 0.001; double T = 2;
  Vector<int> xr (dim); for(int d=0;d<dim;d++) xr[d] = xI[d]/dx[d];
  VectorField pos = sol.get_position();
  VectorField u0 (dim, PeriodicField (xr) );
  for(int d=0;d<dim;d++) for(int it=0;it<u0[d].get_size();++it)
  			   u0[d][it] = cos(pos[d][it]);
   EulerSolver tse(dt, T, &sol, u0);

   tse.get_solution("caca1D");

  return 0;
}
