#pragma once

#include <fstream>
#include <format>
#include <json.hpp>
#include <vector>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Boolean_set_operations_2/oriented_side.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_2.h>

#define BLUE    "\033[34m"
#define GREEN   "\033[32m"
#define RED     "\033[31m" 
#define RESET   "\033[0m"

using json = nlohmann::json;
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K_i;
typedef CGAL::Polygon_2<K> Polygon;
typedef Polygon::Vertex_iterator VertexIterator;
typedef Polygon::Edge_const_iterator EdgeIterator;
typedef CGAL::Point_2<K> Point;
typedef K::Vector_2 Vector;
typedef K::Line_2 Line;
typedef K::Segment_2 Segment;
typedef K::FT FT;

typedef CGAL::Delaunay_mesh_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CDT::Face_handle Face_handle;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Vertex_circulator Vertex_circulator;

typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;

class Problem
{
private:
    std::string name;

    Polygon boundary;
    std::vector<int> boundary_indices;

    std::vector<Point> points;
    std::vector<Point> steiner;

    std::vector<Segment> constraints;
    std::vector<std::tuple<int, int>> constraints_indices;

    std::vector<Polygon> triangulation;

    Problem* best_mesh;
    int num_obtuse;

public:
    Problem(std::string file_name);
    //Problem(char *file_name);

    Problem(Problem *problem){
        name = problem->name;
        points = problem->points;
        boundary = problem->boundary;
        constraints = problem->constraints;
        num_obtuse = -1;
    }

    std::string get_name() { return name; };
    Polygon get_boundary() { return boundary; };
    std::vector<int> get_boundary_indices() { return boundary_indices; };
    std::vector<Point> get_points() { return points; };
    std::vector<Segment> get_constraints() { return constraints; };
    std::vector<std::tuple<int, int>> get_constraints_indices() { return constraints_indices; };

    std::vector<Point> get_steiner() { return steiner; };
    void set_steiner(std::vector<Point> new_steiner) { steiner = new_steiner; };
    void add_steiner(Point s) { steiner.push_back(s); };
    void update_steiner(Point s, Point new_s);

    std::vector<Polygon> get_triangulation(){return triangulation; };
    void set_triangulation(std::vector<Polygon> t){triangulation = t; };
    void add_triangle(Polygon triangle) { triangulation.push_back(triangle); };

    void set_num_obtuse(int num) { num_obtuse = num; };
    int get_num_obtuse() { return num_obtuse; };

    void set_best_mesh(Problem* result) { best_mesh = result; };
    Problem* set_best_mesh() { return best_mesh; };

    void clear_solution() { triangulation.clear(); steiner.clear(); };

    void update_triangulation(Point s, Point new_s);
    void remove_triangle(Polygon t);

    void load_solution();
    void visualize_solution(std::vector<Polygon> problematic_triangles);
    void visualize_solution_voronoi(std::vector<std::pair<Point, Point>> voronoi_edges);

    void save_intermidiate_result(std::string path);
    void write_problem_to_json(const std::string& output_file);

    CDT generate_CDT() {
        std::vector<Segment> all_constraints = constraints;
        for (int i = 0; i < boundary.size(); i++) all_constraints.push_back(Segment(boundary[i], boundary[(i+1)%boundary.size()]));

        CDT cdt;
        cdt.insert(points.begin(), points.end());
        cdt.insert(steiner.begin(), steiner.end());
        cdt.insert_constraints(all_constraints.begin(), all_constraints.end());
        return cdt;
    }

    bool triangle_is_inside(CDT& cdt, const Face_handle& triangle){
        Vector a = Vector(0, 0);
        for (int i = 0; i < 3; i++) {
            Point p = cdt.point(triangle->vertex(i));
            a += Vector(p.x(), p.y());
        }
        Point c = Point(a.x() / 3, a.y() / 3);
        return CGAL::oriented_side(c, boundary) == CGAL::POSITIVE;
    }

    std::vector<std::vector<Point>> grab_triangulation(CDT& cdt) {
        std::vector<std::vector<Point>> triangulation;
        for (const auto& face : cdt.finite_face_handles()) {
            Face_handle fh = face;
            if(triangle_is_inside(cdt, face)){
                std::vector<Point> triangle;
                for (int i = 0; i < 3; i++) {
                    Point p = face->vertex(i)->point();
                    triangle.push_back(p);
                }
                triangulation.push_back(triangle);
            }
        }
        return triangulation;
    }

    void update_problem(CDT& cdt, std::set<Point>& point_set){
        clear_solution();
    
        for (const auto& p : cdt.finite_vertex_handles()) {
            if (point_set.find(p->point()) == point_set.end()) {
                add_steiner(p->point());
            }
        }
        for (auto triangle : grab_triangulation(cdt)) {
            add_triangle(Polygon(triangle.begin(), triangle.end()));
        }
    }

    void generate_uniform_mesh(CDT& cdt, double max_size){
        std::set<Point> point_set(points.begin(), points.end());

        Mesher mesher(cdt);
        mesher.set_criteria(Criteria(0.125, max_size/2));
        mesher.refine_mesh();

        update_problem(cdt, point_set);
    }
};

int find_obtuse_angle(Polygon& triangle);
bool is_obtuse_triangle(Polygon& triangle);
int count_obtuse_triangles(Problem* problem);

double squared_distance(Point& a, Point& b);
double distance(Point& a, Point& b);

void to_IPE(std::string path, std::vector<Point> points, std::vector<Segment> constraints, std::vector<Segment> newConstraints, std::vector<Point> boundary, std::vector<Point> steiner, std::vector<Segment> triangulation, std::vector<Polygon> obtuseTriangles, std::vector<Polygon> problematic_triangles);
void to_SVG(std::string path, std::vector<Point> points, std::vector<Segment> constraints, Polygon boundary, std::vector<Point> steiner, std::vector<Polygon> triangles, std::vector<Polygon> obtuse_triangles);

void to_IPE_voronoi(std::string path, std::vector<Point> points, std::vector<Segment> constraints, std::vector<Segment> newConstraints, std::vector<Point> boundary, std::vector<Point> steiner, std::vector<Segment> triangulation, std::vector<Polygon> obtuseTriangles, std::vector<Polygon> problematic_triangles, std::vector<std::pair<Point, Point>> voronoi_edges);