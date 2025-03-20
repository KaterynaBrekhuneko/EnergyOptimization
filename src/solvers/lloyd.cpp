#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Boolean_set_operations_2/oriented_side.h>

#include "lloyd.hpp"

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CDT::Vertex_handle Vertex_handle;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;

std::vector<std::vector<Point>> grab_triangulation(CDT *cdt, const Polygon& hull) {
    std::vector<std::vector<Point>> triangulation;

    for (auto fit : cdt->finite_face_handles()) {
        std::vector<Point> triangle;
        Vector a = Vector(0, 0);
        for (int i = 0; i < 3; i++) {
            Point p = fit->vertex(i)->point();
            a += Vector(p.x(), p.y());
            triangle.push_back(Point(p.x(), p.y()));
        }

        Point c = Point(a.x() / 3, a.y() / 3);
        if (CGAL::oriented_side(c, hull) == CGAL::POSITIVE) {
            triangulation.push_back(triangle);
        }
    }

    return triangulation;
}

bool is_on_constraint(const Point steiner, Problem* problem){
    std::vector<Segment> constraints = problem->get_constraints();
    for(Segment c : constraints){
        if(c.has_on(steiner)){
            return true;
        }
    }

    return false;
}

Point compute_voronoi_centroid(const Vertex_handle& v, CDT& cdt) {
    std::vector<Point> neighbors;
    CDT::Vertex_circulator vc = cdt.incident_vertices(v);
    if (vc != 0) {
        do {
            if (!cdt.is_infinite(vc)) {
                neighbors.push_back(vc->point());
            }
        } while (++vc != cdt.incident_vertices(v));
    }

    if (neighbors.empty()) return v->point();

    double sum_x = 0, sum_y = 0;
    for (const auto& neighbor : neighbors) {
        sum_x += CGAL::to_double(neighbor.x());
        sum_y += CGAL::to_double(neighbor.y());
    }

    return Point(sum_x / neighbors.size(), sum_y / neighbors.size());
}


void iterate_lloyd(CDT& cdt, Problem *problem, std::vector<Segment>& constraints, int iterations) {
    std::vector<Point> points = problem->get_points();
    Polygon boundary = problem->get_boundary();

    for (int iter = 0; iter < iterations; ++iter) {
        std::vector<Point> new_positions;

        for (auto v = cdt.finite_vertices_begin(); v != cdt.finite_vertices_end(); ++v) {
            Point current = v->point();

            if(std::find(points.begin(), points.end(), current) == points.end() && !boundary.has_on_boundary(current) && !is_on_constraint(current, problem)){
                new_positions.push_back(compute_voronoi_centroid(v, cdt));
            } else {
                new_positions.push_back(v->point());
            }

        }

        // Clear and reinsert points while preserving constraints
        cdt.clear();
        cdt.insert(new_positions.begin(), new_positions.end());
        cdt.insert_constraints(constraints.begin(), constraints.end());
    }
}

SolveStatus Lloyd::solve(Problem *problem) {
    std::vector<Point> points = problem->get_points();
    std::set<Point> point_set(points.begin(), points.end());
    std::vector<Point> boundary = problem->get_boundary().vertices();
    std::vector<Segment> constraints = problem->get_constraints();
    for (int i = 0; i < boundary.size(); i++) constraints.push_back(Segment(boundary[i], boundary[(i+1)%boundary.size()]));
    CGAL::Polygon_2<K> polygon(boundary.begin(), boundary.end());

    std::set<Point> steiner_points = point_set;

    auto cdt = problem->generate_CDT<CDT>();

    auto bbox = polygon.bbox();
    double max_size = std::sqrt(pow(bbox.x_span(),2) + pow(bbox.y_span(),2)) * 0.1;

    std::cout << "Max_size: " << max_size << std::endl;

    Mesher mesher(cdt);
    mesher.set_criteria(Criteria(0, max_size/2));
    mesher.refine_mesh();

    for (auto p : cdt.finite_vertex_handles()) {
        //if (!point_set.contains(p->point())) problem->add_steiner(p->point());
    }

    for (auto triangle : grab_triangulation(&cdt, polygon)) {
        problem->add_triangle(Polygon(triangle.begin(), triangle.end()));
    }

    problem->visualize_solution();

    iterate_lloyd(cdt, problem, constraints, 10);

    problem->clear_solution();

    for (auto p : cdt.finite_vertex_handles()) {
        //if (!point_set.contains(p->point())) problem->add_steiner(p->point());
    }

    for (auto triangle : grab_triangulation(&cdt, polygon)) {
        problem->add_triangle(Polygon(triangle.begin(), triangle.end()));
    }

    problem->visualize_solution();
    
    return SolveStatus::Feasible;
}