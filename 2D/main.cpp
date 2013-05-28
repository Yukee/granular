#include <iostream>
#include <fstream>
#include <math.h>
#include <limits.h>
#include "Vector.h"
#include "Flux.h"
#include "BurgersFlux1D.h"
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
  Flux *f = new Burgers1D();
  return 0;
}
