#include "problem.hpp"
#include "tinyAD_optimization.hpp"

void refine(Problem* problem);
void step_by_step_mesh(Problem* problem);
void classic_delaunay_refinement(Problem* problem);

Point project_onto_segment_refinement(const Point p, const Segment& s);
bool is_on_constraint_refinement(const Point& steiner, Problem* problem);
bool is_on_boundary_refinement(Point& p, Polygon& boundary);
bool is_boundary_segment(Polygon& boundary, Point& a, Point& b);
std::pair<Point, bool> find_boundary_point(Problem* problem, Point& a, Point& b, Point& c);
bool point_within_segment(Point& a, Point& b, Point& p);

void fix_boundary_refinement(Problem *problem);