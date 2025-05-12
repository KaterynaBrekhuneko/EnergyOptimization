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
typedef CDT::Face_handle Face_handle;
typedef CDT::Vertex_circulator Vertex_circulator;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;

bool is_on_constraint_lloyd(const Point& steiner, Problem* problem){
    std::vector<Segment> constraints = problem->get_constraints();
    for(Segment c : constraints){
        if(c.has_on(steiner)){
            return true;
        }
    }

    return false;
}

bool is_on_boundary_lloyd(Point& p, Polygon& boundary){
    double tolerance = 1e-12;
    for (auto edge = boundary.edges_begin(); edge != boundary.edges_end(); ++edge) {
        double distance = CGAL::to_double(CGAL::squared_distance(*edge, p));
        if (distance <= tolerance || edge->source() == p || edge->target() == p) {
            return true;
        }
    }
    return false;
}

void Lloyd::solve(Problem *problem) {
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

    for (auto triangle : problem->grab_triangulation<CDT, Face_handle>(cdt)) {
        problem->add_triangle(Polygon(triangle.begin(), triangle.end()));
    }

    problem->visualize_solution();

    iterate_lloyd<CDT, Vertex_handle, Vertex_circulator>(cdt, problem, constraints, 10);

    problem->clear_solution();

    for (auto p : cdt.finite_vertex_handles()) {
        //if (!point_set.contains(p->point())) problem->add_steiner(p->point());
    }

    for (auto triangle : problem->grab_triangulation<CDT, Face_handle>(cdt)) {
        problem->add_triangle(Polygon(triangle.begin(), triangle.end()));
    }

    problem->visualize_solution();
    
    //return SolveStatus::Feasible;
}