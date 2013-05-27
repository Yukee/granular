#include <iostream>
#include <fstream>
#include <math.h>
#include <limits.h>
#include "Vector.h"
#include "Equation.h"
#include "FD1Solver.h"
#include "timeSolver.h"
#include "EulerSolver.h"
//#include "RK2Solver.h"
#include "RK3Solver.h"
#include "Vector.h"
#include "ScalarField.h"
#include "SegregationEquation.h"
#include "SingleEquation.h"

using namespace std;

double H=1;double W=1.;double U=2.3;int n=4;int m=2;double alpha=0.5;double uF=2.;double h=0.05*H;double sr=0.14;double shift=0.01;
double yShape(double x)
{
    return W*sqrt(tanh(-x/H));
}
double flowDepth(double x)//flow depth in the y=0 plane
{
    x += shift;//BAD bcs you don't shift the other x dependant quantities as well.
    return H*yShape(x)/W;
}

double csch(double x)
{
    return 1/sinh(x);
}

double sech(double x)
{
    return 1/cosh(x);
}

double coth(double x)
{
    return 1/tanh(x);
}

ScalarField flowSurface(VectorField x)//flow surface in the y=0 plane
{
    Vector<int> r = x[0].get_range();
    ScalarField surf(r);
    Vector<int> xp(2,1);Vector<int> zp(2,1);
    xp[1] = 0; zp[0] = 0;
    for(int j=0;j<r[0];j++)
    {
        for(int k=0;k<r[1];k++)
        {
            if( x[1](j*xp+k*zp)<H*flowDepth(x[0](j*xp+k*zp)) ) surf(j*xp+k*zp)=1;
            else surf(j*xp+k*zp)=0;
        }
    }

    return surf;
}

VectorField speed(VectorField x)//here the position x is resized to fit in the boundaries
{
    Vector<int> r = x[0].get_range();
    VectorField speed(x[0].get_space_dimension());
    for(unsigned int dir=0;dir<x[0].get_space_dimension();dir++) speed[dir].resize_field(r);
    Vector<int> xp(2,1);Vector<int> zp(2,1); xp[1] = 0; zp[0] = 0;

    double xcurrent;
    double zcurrent;
    for(int j=0;j<r[0];j++)
    {
        for(int k=0;k<r[1];k++)
        {
            xcurrent = x[0](j*xp+k*zp);
            zcurrent = x[1](j*xp+k*zp);

            speed[0](j*xp+k*zp)=-2. + (2. - 0.31846153846153846*pow(coth(xcurrent),3)*pow(-1.*tanh(xcurrent),3.5))*
            (0.5 + 1.*zcurrent*pow(coth(xcurrent),4)*pow(-1.*tanh(xcurrent),3.5));

            speed[1](j*xp+k*zp)=-1.*pow(zcurrent,2)*coth(xcurrent)*((4.*pow(sech(xcurrent),2)*pow(tanh(xcurrent),3))/pow(-1.*tanh(xcurrent),3.5) +
            (3.5*pow(sech(xcurrent),2)*(-1.*pow(zcurrent,8)*
            pow(0.5 + 1.*zcurrent*pow(coth(xcurrent),4)*pow(-1.*tanh(xcurrent),3.5),8)*
            pow(0. + 1.1146153846153846*coth(xcurrent)*pow(csch(xcurrent),2)*pow(-1.*tanh(xcurrent),2.5) +
            1.2738461538461539*pow(coth(xcurrent),2)*pow(csch(xcurrent),2)*pow(-1.*tanh(xcurrent),3.5),8) +
            pow(tanh(xcurrent),4)))/pow(-1.*tanh(xcurrent),4.5));

        }
    }

    return speed;
}

ScalarField dv(VectorField x)
{
  Vector<int> r = x[0].get_range();
  ScalarField dv(r);

  Vector<int> xp(2,1);Vector<int> zp(2,1); xp[1] = 0; zp[0] = 0;
  double xcurrent;
  double zcurrent;
  for(int j=0;j<r[0];j++)
    {
      for(int k=0;k<r[1];k++)
	{
	  xcurrent = x[0](j*xp+k*zp);
	  zcurrent = x[1](j*xp+k*zp);
	  dv(j*xp+k*zp) = 0.31846153846153846*pow(coth(xcurrent),2)*pow(csch(xcurrent),2)*
	    (0.5 + zcurrent*pow(coth(xcurrent),4)*pow(tanh(xcurrent),7));
	}
    }

  return dv;
}

VectorField f(ScalarField u)
{
    return VectorField (2,0.5*u*u);
}

VectorField df(ScalarField u)
{
    return VectorField (2,u);
}

VectorField linf(ScalarField u)
{
    return VectorField (2,u);
}

VectorField lindf(ScalarField u)
{
    return VectorField (2,1+0*u); //because f*** you
}

int main()
{
    Vector<int> x(2,1);
    Vector<int> y(2,1);
    x[1] = 0; y[0] = 0;

    Vector<double> deltaX(2);
    double deltaT = 0.01;
    double tInterval = 10;
    Vector<double> xInterval(2);
    Vector<double> lowerLeftCorner(2);

    deltaX[0] = 0.05;
    deltaX[1] = 0.01;
    xInterval[0] = 6;// -3<x<0
    xInterval[1] = H;// 0<z<1
    lowerLeftCorner[0] = -xInterval[0]; lowerLeftCorner[1] = 0;
    Vector<int> nxSteps(2);
    nxSteps[0] = xInterval[0]/deltaX[0];
    nxSteps[1] = xInterval[1]/deltaX[1];
    for(int i=0;i<2;i++) if(nxSteps[i]==INT_MAX) cout << "Number of space steps exceeds the largest possible value " << INT_MAX << endl;
    int ntSteps = tInterval/deltaT;
    if(ntSteps==INT_MAX) cout << "Number of time steps exceeds the largest possible value " << INT_MAX << endl;

    /***********************Granular flow 2D*****************************************/
    Equation *flow2DEq = new SegregationEquation(&speed, &dv, sr, 2);
    FD1Solver *flow2Dsolver = new FD1Solver(2, deltaX, xInterval, flow2DEq, prescribedWestAndSouth, &flowSurface, lowerLeftCorner);

    ScalarField psi0(nxSteps);//initial conditions, non extended outside the boundaries
    VectorField pos = flow2Dsolver->get_resized_position();
    for(int j=0;j<nxSteps[0];j++)
    {
        for(int k=0;k<nxSteps[1];k++)
        {
            if( pos[1](j*x+k*y)<(H-h)*flowDepth(pos[0](j*x+k*y)) ) psi0(j*x+k*y)=1;
            else psi0(j*x+k*y)=0;
        }
    }
    Vector<double> uWest(nxSteps[1]); Vector<double> uSouth(nxSteps[0]);
    for(int k=0;k<nxSteps[1];k++)
    {
        if( pos[1](k*y)<(H-h)*flowDepth(pos[0](k*y) - nxSteps[0]) ) uWest[k] = 1;
        else uWest[k] = 0;
    }
    flow2Dsolver->set_uWest(uWest);
    for(int j=0;j<nxSteps[0];j++)
    {
        uSouth[j] = 1;
    }
    flow2Dsolver->set_uSouth(uSouth);

    timeSolver *flow2DRK3 = new RK3Solver(deltaT, tInterval, flow2Dsolver, psi0);
    flow2DRK3->get_solution("concentration_with_source");

    /************************Burgers 2D***************************************************/
   // Equation *burgersEq = new SingleEquation(&f, &df);
   // FD1Solver *burgersSolver = new FD1Solver(2, deltaX, xInterval, burgersEq, periodic, &flowSurface, lowerLeftCorner);
   // ScalarField u0(nxSteps);
   // VectorField pos = burgersSolver->get_resized_position();
   // for(int j=0;j<nxSteps[0];j++)
   // {
   //     for(int k=0;k<nxSteps[1];k++)
   //     {
   //         u0(j*x+k*y) = 0.5 + sin(pos[0](j*x+k*y));
   //     }
   // }
   // timeSolver *burgersRK3 = new RK3Solver(deltaT, tInterval, burgersSolver, u0);
   // cout << currentT << endl;
   // burgersRK3->get_solution("test_burgers_1D_new");

    /**********************Wave 2D*******************************************************/
//    Equation *linEq = new SingleEquation(&f, &df);
//    FD1Solver *linSolver = new FD1Solver(2, deltaX, xInterval, linEq, periodic, &flowSurface, lowerLeftCorner);
//    ScalarField u0(nxSteps);
//    VectorField pos = linSolver->get_resized_position();
//    for(int j=0;j<nxSteps[0];j++)
//    {
//        for(int k=0;k<nxSteps[1];k++)
//        {
//            u0(j*x+k*y) = 0.5 + sin(pos[0](j*x+k*y));
//        }
//    }
//    timeSolver *linRK3 = new RK3Solver(deltaT, tInterval, linSolver, u0);
//    linRK3->get_solution("test_burgers_1D_propagating_shock");

    return 0;
}
