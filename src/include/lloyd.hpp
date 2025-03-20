#include "problem.hpp"
#include "solver.hpp"

class Lloyd : public Solver
{
public:
    SolveStatus solve(Problem* prob);
};