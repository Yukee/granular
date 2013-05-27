#ifndef EULERSOLVER_H
#define EULERSOLVER_H

#include "timeSolver.h"

class EulerSolver : public timeSolver
{
    public:
        EulerSolver(double deltaT, double T, FD1Solver *spatialSolver, ScalarField initial_conditions);

        void get_solution(std::string name);
};

#endif // EULERSOLVER_H
