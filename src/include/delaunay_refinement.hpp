#include "problem.hpp"
#include "tinyAD_optimization.hpp"
#include "mesh_statistics.hpp"
#include "gradient_descent_optimization.hpp"

Mesh_Statistics refine(Problem* problem);
void step_by_step_mesh(Problem* problem);
void classic_delaunay_refinement(Problem* problem);
Mesh_Statistics offcenter_delaunay_refinement(Problem* problem);
void pentagon_refinement(Problem* problem);

void locally_optimize_obtuse_refinement(Problem *problem);

Point project_onto_segment_refinement(const Point p, const Segment& s);
bool is_on_constraint_refinement(const Point& steiner, Problem* problem);
bool is_constrained_segment(Problem* problem, Point& a, Point& b);
bool is_on_boundary_refinement(Point& p, Polygon& boundary);
bool is_boundary_segment(Polygon& boundary, Point& a, Point& b);
std::pair<Point, bool> find_boundary_point(Problem* problem, Point& a, Point& b, Point& c);
bool point_within_segment(Point& a, Point& b, Point& p);

void fix_boundary_refinement(Problem *problem);
void fix_constraint_refinement(Problem *problem);
//void fix_rhombus(Problem* problem, CDT& cdt);

Mesh_Statistics uniform_mesh(Problem* problem);

void mesh_equilateral_single(Problem* problem);
void mesh_cgal(Problem* problem);
Mesh_Statistics mesh_equilateral(Problem* problem);

double mean_absolute_deviation(Problem* problem);
double mean_aspect_ratio(Problem* problem);
std::vector<Polygon> save_min_max_angle(Problem* problem, Mesh_Statistics* stats);

void save_angle_stats_for_plot(Problem* problem, std::string path);
void save_aspect_ratios_for_plot(Problem* problem, std::string path);