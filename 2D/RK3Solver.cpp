#include "RK3Solver.h"
#include "FD1Solver.h"
#include <string>
#include <boost/lexical_cast.hpp>

using namespace std;

RK3Solver::RK3Solver(double deltaT, double T, FD1Solver *spatialSolver, ScalarField initial_conditions) :
    timeSolver(deltaT, T, spatialSolver, initial_conditions)
{
    m_ntSteps = m_T/m_deltaT;
}

void RK3Solver::get_solution(string name)
{
    Vector<double> deltaX = m_spatialSolver->get_deltaX();
    Vector<double> lowerLeftCorner = m_spatialSolver->get_lowerLeftCorner();
    fstream data;
    string path;

    ScalarField u;
    ScalarField df;
    double testDeltaT;
    double newDeltaT = m_deltaT;
    double currenttime = 0;
    int writingCounter = 0;

    for(currenttime=0;currenttime<=m_T;currenttime += newDeltaT)
    {

        newDeltaT = m_deltaT;
        df = m_spatialSolver->get_numerical_flux_gradient(m_un);
	testDeltaT = m_spatialSolver->check_CFL(newDeltaT);
        if(testDeltaT!=newDeltaT) newDeltaT = testDeltaT;

        u = m_un + (-1)*newDeltaT*df;

        df = m_spatialSolver->get_numerical_flux_gradient(u);
        u = (3/4.)*m_un + (1/4.)*(u + (-1)*newDeltaT*df);

        df = m_spatialSolver->get_numerical_flux_gradient(u);
        u = (1/3.)*m_un + (2/3.)*(u + (-1)*newDeltaT*df);

        if( (int)(currenttime/10) == writingCounter )
        {
            cout << "writing file number " << writingCounter << endl;
            //path = "Results/" + name + "_" + to_string(writingCounter) + ".tsv"; only in C++ 11
	    path =  "Results/" + name + "_" + boost::lexical_cast<string>(writingCounter) + ".tsv";
            data.open(path.c_str(), ios::out);
            u.write_in_file(data, deltaX, lowerLeftCorner);
            data.close();
            writingCounter++;
        }
        m_un = u;
    }
}


