#include "FD1Solver.h"
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

FD1Solver::FD1Solver(int spaceDimension, Vector<double> deltaX, Vector<double> xInterval, Equation *eq, boundaryConditions bc, boundarySurface bs, Vector<double> lowerLeftCorner) :
   m_spaceDimension(spaceDimension), m_deltaX(deltaX), m_xInterval(xInterval), m_eq(eq), m_bc(bc), m_bs(bs), m_lowerLeftCorner(lowerLeftCorner)
{
    x.resize(m_spaceDimension); y.resize(m_spaceDimension);
    x[0] = 1; x[1] = 0; y[0] = 0; y[1] = 1;

    for(int i=0;i<m_spaceDimension;i++) m_nxSteps[i] = m_xInterval[i]/m_deltaX[i];
    m_position.resize(m_spaceDimension);
    if(lowerLeftCorner.size() != m_spaceDimension) throw invalid_argument("lowerLeftCorner dimension is not equal to the space dimension");
    for(int dir=0;dir<m_spaceDimension;dir++)
    {
        m_position[dir].resize_field(m_nxSteps+2*(x+y));
        //assuming 2D space
        for(int j=0;j<m_nxSteps[0]+2;j++) for(int k=0;k<m_nxSteps[1]+2;k++) {m_position[dir](j*x+k*y) = m_lowerLeftCorner[dir] + m_deltaX[dir]*((1-dir)*(j-1) + dir*(k-1));}//shift of j and k since we use extended BC;m_position(j,k) returns position at the point j-1,k-1 of the grid
    }

    cout << "begining calculating boundaries" << endl;
    m_test_boundaries = m_bs(m_position);
    cout << "finished calculating boundaries" << endl;

    resize_pos();

    if(bc == prescribedWestAndSouth)
    {
        cout << "begining calculating speed field" << endl;
        m_eq->set_speed_field(m_resizedPos, m_bs(m_resizedPos));
        cout << "finished calculating speed field" << endl;
    }

    upper_right_intermediate_un_values.resize(m_spaceDimension);
    lower_right_intermediate_un_values.resize(m_spaceDimension);
    upper_left_intermediate_un_values.resize(m_spaceDimension);
    lower_left_intermediate_un_values.resize(m_spaceDimension);

    right_localSpeed.resize(m_spaceDimension);
    left_localSpeed.resize(m_spaceDimension);

    right_convection_flux.resize(m_spaceDimension);
    left_convection_flux.resize(m_spaceDimension);
}

FD1Solver::~FD1Solver()
{
    if(m_eq) delete[] m_eq;
}

double FD1Solver::check_CFL(double deltaT)
{
    double newDeltaT = deltaT;
    Vector<ScalarField> waveSpeed(m_spaceDimension);
    Vector<double> maxFrequency(m_spaceDimension);
    double overallMaxFreq=0;
    for(int dir=0;dir<m_spaceDimension;dir++)
    {
        waveSpeed[dir] = right_localSpeed[dir].max_field(left_localSpeed[dir]);
        maxFrequency[dir] = waveSpeed[dir].get_max()/m_deltaX[dir];
        if(maxFrequency[dir]>overallMaxFreq) overallMaxFreq = maxFrequency[dir];
    }
    if(newDeltaT>1/(8*overallMaxFreq)) newDeltaT = 1/(9*overallMaxFreq);
    return newDeltaT;
}

double FD1Solver::minmod(double a, double b)
{
  double minmod;
  if(a*b<=0) minmod = 0;
  else
    {
      if(a>0) minmod = min(a,b);
      else minmod = max(a,b);
    }
  return minmod;
}

double FD1Solver::three_pts_derivative(Vector<int> j, int dir)
{
  double deriv;
    Vector<int> jp;
    Vector<int> jm;
    jp = j + (1-dir)*x + dir*y; //translation along the x axis if dir=0, along y axis if dir=1
    jm = j + (dir-1)*x + (-dir)*y;

    //bc TODO replace that by a second layer of un in set_un, and don't forget to shift by 2*(x+y) un(Vector<int>)
    if(m_bc == periodic)
    {
        if(jm[0]==-2) jm[0]+=m_nxSteps[0];
        if(jm[1]==-2) jm[1]+=m_nxSteps[1];
        if(jp[0]==m_nxSteps[0]+1) jp[0]-=m_nxSteps[0];
        if(jp[1]==m_nxSteps[1]+1) jp[1]-=m_nxSteps[1];
    }

    else
    {
        if(jm[0]==-2) return (un(jp) - un(j))/m_deltaX[dir];
        if(jm[1]==-2) return (un(jp) - un(j))/m_deltaX[dir];
        if(jp[0]==m_nxSteps[0]+1) return (un(j) - un(jm))/m_deltaX[dir];
        if(jp[1]==m_nxSteps[1]+1) return (un(j) - un(jm))/m_deltaX[dir];
    }

    deriv = minmod( (un(j) - un(jm))/m_deltaX[dir], (un(jp) - un(j))/m_deltaX[dir] );
    return deriv;
}

// double FD1Solver::intermediate_uxt_values(Vector<int> j, grid_position p, bound b, int dir)
// {
//     Vector<int> jp;
//     Vector<int> jm;
//     jp = j + (1-dir)*x + dir*y;
//     jm = j + (dir-1)*x + (-dir)*y;

//     if(p == lft)
//     {
//         if(b == lower) return un(jm) + m_deltaX[dir]*un_derivative[dir](jm)*0.5;
//         if(b == upper) return un(j) - m_deltaX[dir]*un_derivative[dir](j)*0.5;
//     }
//     if(p == rght)
//     {
//         if(b == lower) return un(j) + m_deltaX[dir]*un_derivative[dir](j)*0.5;
//         if(b == upper) return un(jp) - m_deltaX[dir]*un_derivative[dir](jp)*0.5;
//     }

//     return -1;
// }

Vector<double> FD1Solver::u_values_at_cell_edges(Vector<int> j)
{
  Vector<double> values(4*m_spaceDimension);
  Vector<int> jp(m_spaceDimension);
  Vector<int> jm(m_spaceDimension);

  for(int d=0;d<m_spaceDimension;d++)
    {
      jp = j + (1-d)*x + d*y;
      jm = j + (d-1)*x + (-d)*y;

      values[4*d] = un(jm) + m_deltaX[d]*three_pts_derivative(jm,d)*0.5;//un-1/2-
      values[4*d+1] = un(j) - m_deltaX[d]*three_pts_derivative(j,d)*0.5;//un-1/2+
      values[4*d+2] = un(j) + m_deltaX[d]*three_pts_derivative(j,d)*0.5;//un+1/2-
      values[4*d+3] = un(jp) - m_deltaX[d]*three_pts_derivative(jp,d)*0.5;;//un+1/2+
    }
  return values;
}

void FD1Solver::compute_intermediate_un_values()
{
    for(int dir=0;dir<m_spaceDimension;dir++)
    {
        upper_right_intermediate_un_values[dir].resize_field(m_nxSteps);
        lower_right_intermediate_un_values[dir].resize_field(m_nxSteps);
        upper_left_intermediate_un_values[dir].resize_field(m_nxSteps);
        lower_left_intermediate_un_values[dir].resize_field(m_nxSteps);
    }

    Vector<int> p;
    Vector<double> cellEdges;
        for(int j=0;j<m_nxSteps[0];j++)
        {
            for(int k=0;k<m_nxSteps[1];k++)
            {
                p[0] = j; p[1] = k;
		cellEdges = u_values_at_cell_edges(p);
		for(int d=0;d<m_spaceDimension;d++)
		  {
		    lower_left_intermediate_un_values[d](p) = cellEdges[4*d];
		    upper_left_intermediate_un_values[d](p) = cellEdges[4*d+1];
		    lower_right_intermediate_un_values[d](p) = cellEdges[4*d+2];
		    upper_right_intermediate_un_values[d](p) = cellEdges[4*d+3];
		  }
            }
        }
}

void FD1Solver::compute_localSpeed()
{
    ScalarField lowerSpeed;
    ScalarField upperSpeed;
    //Modify here if you want to use implicit definition of the jacobian

    for(int dir=0;dir<m_spaceDimension;dir++)
    {
        upperSpeed = m_eq->get_convectionFluxJacobian(upper_right_intermediate_un_values[dir])[dir];
        lowerSpeed = m_eq->get_convectionFluxJacobian(lower_right_intermediate_un_values[dir])[dir];
        right_localSpeed[dir] = upperSpeed.max_field(lowerSpeed);

        upperSpeed = m_eq->get_convectionFluxJacobian(upper_left_intermediate_un_values[dir])[dir];
        lowerSpeed = m_eq->get_convectionFluxJacobian(lower_left_intermediate_un_values[dir])[dir];
        left_localSpeed[dir] = upperSpeed.max_field(lowerSpeed);
    }
}

void FD1Solver::compute_numerical_convection_flux()
{
    //no particular bc on the flux?
    ScalarField upperFlux;
    ScalarField lowerFlux;

    for(int dir=0;dir<m_spaceDimension;dir++)
    {
        upperFlux = m_eq->get_convectionFlux(upper_right_intermediate_un_values[dir])[dir];
        lowerFlux = m_eq->get_convectionFlux(lower_right_intermediate_un_values[dir])[dir];
        right_convection_flux[dir] = 0.5*(upperFlux + lowerFlux) + (-0.5)*right_localSpeed[dir]*(upper_right_intermediate_un_values[dir] + (-1)*lower_right_intermediate_un_values[dir]);

        upperFlux = m_eq->get_convectionFlux(upper_left_intermediate_un_values[dir])[dir];
        lowerFlux = m_eq->get_convectionFlux(lower_left_intermediate_un_values[dir])[dir];
        left_convection_flux[dir] = 0.5*(upperFlux + lowerFlux) + (-0.5)*left_localSpeed[dir]*(upper_left_intermediate_un_values[dir] + (-1)*lower_left_intermediate_un_values[dir]);
    }
}

void FD1Solver::set_un(ScalarField Un)
{
    Vector<int> rBc; rBc = m_nxSteps + 2*(x+y);
    m_un.resize_field(rBc);

    if(m_bc == null)
    {
        //assuming 2D
        for(int j=0;j<m_nxSteps[0]+2;j++)
        {
            for(int k=0;k<m_nxSteps[1]+2;k++)
            {
                //assuming that test_boundaries is written with supplementary borders (same as m_un)
                if(m_test_boundaries(j*x+k*y)==0) m_un(j*x+k*y) = 0;
                else
                {
                    if(j==0 || k==0) throw invalid_argument("In FD1Solver::set_un: boundary surface is on the rectangular integration domain; it should be strictly inside");
                    m_un(j*x+k*y) = Un((j-1)*x + (k-1)*y);
                }
            }
        }
    }

    if(m_bc == prescribedWestAndSouth)
    {
        //assuming 2D
        for(int j=1;j<m_nxSteps[0]+2;j++)
        {
            for(int k=1;k<m_nxSteps[1]+2;k++)
            {
                //assuming that test_boundaries is written with supplementary borders (same as m_un)
                if(m_test_boundaries(j*x+k*y)==0) m_un(j*x+k*y) = 0;
                else
                {
                    if(j==0 || k==0) throw invalid_argument("In FD1Solver::set_un: boundary surface is on the rectangular integration domain; it should be strictly inside, except on the west and south sides of the rectangle");
                    m_un(j*x+k*y) = Un((j-1)*x + (k-1)*y);
                }
            }
        }

        for(int k=1;k<m_nxSteps[1]+1;k++)
        {
            if(uWest.size()!=m_nxSteps[1]) throw invalid_argument("In FD1Solver::set_un: field value on the west boundary should have as many elements as the number of vertical steps. Perhaps have you forgot to set uWest?");
            m_un(k*y) = uWest[k-1];
        }
        m_un(0*y) = 0; m_un((m_nxSteps[1]+1)*y) = 0;//Actually we don't even have to set this values since they won't be used

        for(int j=1;j<m_nxSteps[0]+1;j++)
        {
            if(uSouth.size()!=m_nxSteps[0]) throw invalid_argument("In FD1Solver::set_un: field value on the south boundary should have as many elements as the number of horizontal steps. Perhaps have you forgot to set uSouth?");
            m_un(j*x) = uSouth[j-1];
        }
        m_un((m_nxSteps[0]+1)*x) = 0;
    }

    if(m_bc == periodic)
    {
        Vector<int> p1; Vector<int> p2;
        for(int dir=0;dir<2;dir++)
        {
            p1=0*p1;
            p2 = m_nxSteps + (-1)*(x+y);
            for(int i=1;i<m_nxSteps[dir]+1;i++)
            {
                p1[dir] = i; p2[dir] = i-1;
                m_un(p1) = Un(p2);
            }
        }

        for(int j=1;j<m_nxSteps[0]+1;j++) for(int k=1;k<m_nxSteps[1]+1;k++)
        {
            m_un(j*x + k*y) = Un((j-1)*x + (k-1)*y);
        }

        for(int dir=0;dir<2;dir++)
        {
            p1 = m_nxSteps + (x+y);
            p2=0*p2;
            for(int i=1;i<m_nxSteps[dir]+1;i++)
            {
                p1[dir] = i; p2[dir]=i-1;
                m_un(p1) = Un(p2);
            }
        }
    }
}

//the flux gradient may be infinite, if you take a too large time step.
ScalarField FD1Solver::get_numerical_flux_gradient(ScalarField un)
{
    set_un(un);
    compute_intermediate_un_values();
    compute_localSpeed();
    compute_numerical_convection_flux();
    ScalarField flux_gradient(m_nxSteps);
    for(int j=0;j<m_nxSteps[0];j++) for(int k=0;k<m_nxSteps[1];k++) flux_gradient(j*x+k*y) = 0;
    for(int dir=0;dir<m_spaceDimension;dir++)
    {
        flux_gradient = flux_gradient + (1./m_deltaX[dir])*( right_convection_flux[dir] + (-1)*left_convection_flux[dir] );
    }
    return flux_gradient;
}

double FD1Solver::un(Vector<int> j)
{
  return m_un(j+x+y);
}


void FD1Solver::resize_pos()
{
    m_resizedPos.resize(m_spaceDimension);
    for(int dir=0;dir<m_spaceDimension;dir++) m_resizedPos[dir].resize_field(m_nxSteps);
    for(int dir=0;dir<m_spaceDimension;dir++)
    {
    for(int j=1;j<m_nxSteps[0]+1;j++)
    {
        for(int k=1;k<m_nxSteps[1]+1;k++)
        {
            m_resizedPos[dir]((j-1)*x+(k-1)*y)=m_position[dir](j*x+k*y);
        }
    }
    }
}
