#include "EulerSolver.h"
#include "FD1Solver.h"
#include <iostream>

using namespace std;

EulerSolver::EulerSolver(double deltaT, double T, FD1Solver *spatialSolver, ScalarField initial_conditions) :
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
    long long int writingCounter = 0;

    ScalarField un1 = m_un;
    ScalarField df;

    for(currenttime=0;currenttime<=m_T;currenttime += newDeltaT)
    {
        df = m_spatialSolver->get_numerical_flux_gradient(m_un);

        testDeltaT = m_spatialSolver->check_CFL(newDeltaT);
        if(testDeltaT!=newDeltaT) newDeltaT = testDeltaT;

        un1 = m_un + (-1)*newDeltaT*df;

        if( (int)(10*currenttime) == writingCounter )
        {
            cout << "writing file number " << writingCounter << endl;
            path = "Results/" + name + "_" + to_string(writingCounter) + ".tsv";
            data.open(path.c_str(), ios::out);
            un1.write_in_file(data, deltaX, lowerLeftCorner);
            data.close();
            writingCounter++;
        }

        m_un = un1;
    }
}
