#include "RK2Solver.h"
#include "FD1Solver.h"

using namespace std;

RK2Solver::RK2Solver(double deltaT, double T, FD1Solver *spatialSolver, std::vector< std::vector<double> > initial_conditions) :
    timeSolver(deltaT, T, spatialSolver, initial_conditions)
{
    m_ntSteps = m_T/m_deltaT;
}

void RK2Solver::get_solution(const char *filename)
{
    fstream data(filename, ios::out);
    vector<double> deltaX = m_spatialSolver->get_deltaX();
    vector<int> nxSteps = m_spatialSolver->get_nxSteps();

    // initial data
    for(int j=0;j<nxSteps[0];j++)
    {
        for(int k=0;k<nxSteps[1];k++)
        {
            data << j*deltaX[0] << " " << k*deltaX[1] << " " << m_un[j][k] << endl;
        }
    }
    data << endl << endl;

    vector< vector<double> > u = m_un;
    vector< vector <double> > df;

    for(int n=0;n<m_ntSteps;n++)
    {
        df = m_spatialSolver->get_numerical_flux_gradient(m_un);
        for(int j=0;j<nxSteps[0];j++)
        {
            for(int k=0;k<nxSteps[1];k++)
            {
                u[j][k] = m_un[j][k] - m_deltaT*df[j][k];
            }
        }

        df = m_spatialSolver->get_numerical_flux_gradient(u);
        for(int j=0;j<nxSteps[0];j++)
        {
            for(int k=0;k<nxSteps[1];k++)
            {
                u[j][k] = 0.5*m_un[j][k] + 0.5*(u[j][k] - m_deltaT*df[j][k]);
                data << j*deltaX[0] << " " << k*deltaX[1] << " " << u[j][k] << endl;
            }
        }
        data << endl << endl;
        m_un=u;
    }
}

