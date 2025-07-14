#pragma once

#include "problem.hpp"
#include "global_optimization.hpp"

void perform_edge_flips(Problem *problem, bool ignore);

void locally_optimize_solution(Problem *problem);
void locally_optimize_obtuse(Problem *problem);
void locally_optimize_constraint(Problem *problem);

void globally_optimize_solution(Problem *problem);
void globally_optimize_obtuse(Problem *problem);

void fix_boundary(Problem *problem);