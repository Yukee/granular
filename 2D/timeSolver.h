#ifndef TIMESOLVER_H
#define TIMESOLVER_H

#include <vector>
#include <string>
#include <fstream>
#include "FD1Solver.h"
#include "ScalarField.h"

class timeSolver
{
    public:
        timeSolver(double deltaT, double T, FD1Solver *spatialSolver, ScalarField  initial_conditions);

        virtual void get_solution(std::string path) = 0;

    protected:
        double m_deltaT;
        double m_T;
        int m_ntSteps;
        FD1Solver *m_spatialSolver;
        ScalarField m_un;
};

#endif // TIMESOLVER_H
