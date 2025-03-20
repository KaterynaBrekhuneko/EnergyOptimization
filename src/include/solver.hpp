#pragma once

#include "problem.hpp"
#include "miscellaneous.hpp"

class Solver {

public:
    virtual ~Solver() {};

    virtual SolveStatus solve(Problem* prob) = 0;
};