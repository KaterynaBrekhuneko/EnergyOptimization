#include "problem.hpp"
#include "solver.hpp"

bool is_on_constraint_lloyd(const Point& steiner, Problem* problem);
bool is_on_boundary_lloyd(Point& p, Polygon& boundary);
Point compute_voronoi_centroid(const Vertex_handle& v, CDT& cdt);
void iterate_lloyd(CDT& cdt, Problem *problem, std::vector<Segment>& constraints, int iterations);