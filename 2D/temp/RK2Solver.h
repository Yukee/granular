#ifndef RK2SOLVER_H
#define RK2SOLVER_H

#include "timeSolver.h"

class RK2Solver : public timeSolver
{
    public:
        RK2Solver(double deltaT, double T, FD1Solver *spatialSolver, std::vector< std::vector<double> > initial_conditions);

        void get_solution(const char *filename);
};

#endif

