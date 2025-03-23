#include "problem.hpp"
#include "solver.hpp"

std::vector<Polygon> find_neighborhood(const Point& steiner, std::vector<Polygon>& triangles);
bool is_interior_vertex(const Point& steiner, const Polygon& boundary);
bool is_in_the_neighborhood(const Point& p, const std::vector<Polygon>& neighborhood);
int find_point_index(const Point& s, const Polygon& polygon);
void update_neighborhood(std::vector<Polygon>& neighborhood, const Point& s, const Point& new_s);
double norm(const Point& gradient);
Segment find_boundary_segment(const Point& p, const Polygon& polygon);
Point project_onto_segment(const Point p, const Segment& s);
bool is_on_constraint(const Point& steiner, Problem* problem, Segment* constraint);
bool is_interior_vertex_with_tolerance(const Point& p, const Polygon& polygon);
bool is_on_constraint_with_tolerance(const Point& steiner, Problem* problem, Segment* constraint);

std::vector<Point> line_search(std::vector<Polygon>& triangles, const std::vector<Point>& steiner, const std::vector<Point>& gradient, double s_max = 1.0, double shrink = 0.8, double max_iters = 64, double armijo_const = 1e-4);

double value_ln(const Point& old_s, const Point& new_s, const std::vector<Polygon> neighborhood);
double value_sigmoid(const Point& old_s, const Point& new_s, const std::vector<Polygon> neighborhood);
double value_refined_sigmoid(const Point& old_s, const Point& new_s, const std::vector<Polygon> neighborhood);

Point calculate_gradient_equilateral(Point& s, std::vector<Polygon>& triangles);
Point calculate_gradient_ln(Point& s, std::vector<Polygon>& triangles);
Point calculate_gradient_sigmoid(Point& s, std::vector<Polygon>& triangles);
Point calculate_gradient_refined_sigmoid(Point& s, std::vector<Polygon>& triangles);

Point locally_optimize_position_constraint(Point steiner, std::vector<Polygon>& triangles, Problem *problem);
Point locally_optimize_position(Point steiner, std::vector<Polygon>& triangles, Problem *problem);