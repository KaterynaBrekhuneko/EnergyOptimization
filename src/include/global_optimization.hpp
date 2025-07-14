#pragma once

#include "problem.hpp"
#include "local_optimization.hpp"


std::vector<Point> globally_optimize_position(std::vector<Point>& steiner_points, std::vector<Polygon>& triangles, Problem *problem, bool debug);