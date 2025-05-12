#pragma once

#include "problem.hpp"
#include "miscellaneous.hpp"

class Solver {

public:
    virtual ~Solver() {};

    virtual void solve(Problem* prob) = 0;
};