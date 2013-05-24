#ifndef RK3SOLVER_H
#define RK3SOLVER_H

#include "ScalarField.h"
#include "timeSolver.h"

class RK3Solver : public timeSolver
{
    public:
        RK3Solver(double deltaT, double T, FD1Solver *spatialSolver, ScalarField initial_conditions);

        void get_solution(std::string );
};

#endif


