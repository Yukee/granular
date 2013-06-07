#include "EulerSolver.h"
#include "FD1Solver.h"
#include <iostream>
#include <boost/lexical_cast.hpp>

using namespace std;

EulerSolver::EulerSolver(double deltaT, double T, FD1Solver *spatialSolver, VectorField initial_conditions) :
  timeSolver(deltaT, T, spatialSolver, initial_conditions)
{
  m_ntSteps = m_T/m_deltaT;
}

void EulerSolver::get_solution(string name)
{
  Vector<double> deltaX = m_spatialSolver->get_deltaX();
  Vector<double> lowerLeftCorner = m_spatialSolver->get_lowerLeftCorner();
  fstream data;
  string path;

  double testDeltaT;
  double newDeltaT = m_deltaT;
  double currenttime = 0;
  int writingCounter = 0;

  VectorField un1 = m_un;
  VectorField df;

  PeriodicField unity = m_spatialSolver->get_unity();

  for(currenttime=0;currenttime<=m_T;currenttime += newDeltaT)
    {
      df = m_spatialSolver->get_numerical_flux_gradient(m_un);

      testDeltaT = m_spatialSolver->check_CFL(newDeltaT);
      if(testDeltaT!=newDeltaT) newDeltaT = testDeltaT;

      un1 = m_un  - (newDeltaT*unity)*df;

      if( (int)(1*currenttime) == writingCounter )
        {
	  cout << "writing file number " << writingCounter << endl;
	  //path = "Results/" + name + "_" + to_string(writingCounter) + ".tsv"; Only in C++ 11
	  path =  "Results/" + name + "_" + boost::lexical_cast<string>(writingCounter) + ".tsv";
	  data.open(path.c_str(), ios::out);
	  un1[0].write_in_file(data, deltaX, lowerLeftCorner);
	  data.close();
	  writingCounter++;
        }

      m_un = un1;
    }
}
