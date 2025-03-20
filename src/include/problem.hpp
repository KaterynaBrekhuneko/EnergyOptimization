#pragma once

#include <fstream>
#include <format>
#include <json.hpp>
#include <vector>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_2.h>

using json = nlohmann::json;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Exact_predicates_inexact_constructions_kernel K_i;
typedef CGAL::Polygon_2<K> Polygon;
typedef Polygon::Vertex_iterator VertexIterator;
typedef Polygon::Edge_const_iterator EdgeIterator;
typedef CGAL::Point_2<K> Point;
typedef K::Vector_2 Vector;
typedef K::Line_2 Line;
typedef K::Segment_2 Segment;
typedef K::FT FT;

class Problem
{
private:
    std::string name;

    Polygon boundary;

    std::vector<Point> points;
    std::vector<Point> steiner;

    std::vector<Segment> constraints;
    std::vector<Polygon> triangulation;

public:
    Problem(std::string file_name);
    //Problem(char *file_name);

    // ! delete this constructor?
    Problem(Problem *problem){
        name = problem->name;
        points = problem->points;
        boundary = problem->boundary;
        constraints = problem->constraints;
    }

    std::string get_name() { return name; };
    Polygon get_boundary() { return boundary; };
    std::vector<Point> get_points() { return points; };
    std::vector<Segment> get_constraints() { return constraints; };

    std::vector<Point> get_steiner() { return steiner; };
    void add_steiner(Point s) { steiner.push_back(s); };

    std::vector<Polygon> get_triangulation(){return triangulation; };
    void set_triangulation(std::vector<Polygon> t){triangulation = t; };
    void add_triangle(Polygon triangle) { triangulation.push_back(triangle); };

    void clear_solution() { triangulation.clear(); steiner.clear(); };

    void update_triangulation(Point s, Point new_s);
    void remove_triangle(Polygon t);

    void load_solution();
    void visualize_solution();

    template <typename CDT>
    CDT generate_CDT() {
        std::vector<Segment> all_constraints = constraints;
        for (int i = 0; i < boundary.size(); i++) all_constraints.push_back(Segment(boundary[i], boundary[(i+1)%boundary.size()]));

        CDT cdt;
        cdt.insert(points.begin(), points.end());
        cdt.insert(steiner.begin(), steiner.end());
        cdt.insert_constraints(all_constraints.begin(), all_constraints.end());
        return cdt;
    }
};

int find_obtuse_angle(Polygon& triangle);
bool is_obtuse_triangle(Polygon& triangle);
int count_obtuse_triangles(Problem* problem);

double squared_distance(Point& a, Point& b);
double distance(Point& a, Point& b);

void to_IPE(std::string path, std::vector<Point> points, std::vector<Segment> constraints, std::vector<Segment> newConstraints, std::vector<Point> boundary, std::vector<Point> steiner, std::vector<Segment> triangulation, std::vector<Polygon> obtuseTriangles);