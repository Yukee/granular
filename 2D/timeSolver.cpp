#include "timeSolver.h"
#include <iostream>
#include <fstream>

using namespace std;

timeSolver::timeSolver(double deltaT, double T, FD1Solver *spatialSolver, VectorField initial_conditions)
: m_deltaT(deltaT), m_T(T), m_spatialSolver(spatialSolver)
{
    m_un = initial_conditions; //Why can't I use the initialization list?
    m_ntSteps = m_T/m_deltaT;
}
