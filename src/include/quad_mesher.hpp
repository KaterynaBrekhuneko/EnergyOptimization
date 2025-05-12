#include "problem.hpp"

void build_quad_mesh_medians(Problem* problem);
void build_quad_mesh_gmsh(Problem* problem);

void to_IPE_quad(std::string path, std::vector<Point> points, std::vector<Segment> constraints, Polygon boundary, std::vector<Point> steiner, std::vector<Polygon> quads);
void to_SVG_quad(std::string path, std::vector<Point> points, std::vector<Segment> constraints, Polygon boundary, std::vector<Point> steiner, std::vector<Polygon> quads);

void add_text_and_overwrite(Problem* problem, const std::string& path);