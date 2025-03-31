#include "delaunay_refinement.hpp"

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Boolean_set_operations_2/oriented_side.h>
#include <CGAL/lloyd_optimize_mesh_2.h>
#include <CGAL/Meshes/Double_map_container.h>

#include <unordered_map>
#include <limits>

#define BLUE    "\033[34m"
#define GREEN   "\033[32m"
#define RED     "\033[31m" 
#define RESET   "\033[0m"

typedef CGAL::Delaunay_mesh_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;

typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>            Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria>              Mesher;

typedef CDT::Face_handle Face_handle;
typedef CDT::Edge Edge;

CDT generate_triangulation(std::set<Point>& points, std::set<Point>& steiner_points, std::vector<Segment>& constraints) {
    CDT cdt;
    cdt.insert(points.begin(), points.end());
    cdt.insert(steiner_points.begin(), steiner_points.end());
    cdt.insert_constraints(constraints.begin(), constraints.end());
    return cdt;
}

bool triangle_is_inside(CDT* cdt, const Polygon& boundary, const Face_handle& triangle){
    Vector a = Vector(0, 0);
    for (int i = 0; i < 3; i++) {
        CDT::Point p = cdt->point(triangle->vertex(i));
        a += Vector(p.x(), p.y());
    }
    Point c = Point(a.x() / 3, a.y() / 3);
    return CGAL::oriented_side(c, boundary) == CGAL::POSITIVE;
}

std::vector<std::vector<Point>> grab_triangulation(CDT* cdt, const Polygon& boundary) {
    std::vector<std::vector<Point>> triangulation;
    for (const auto& face : cdt->finite_face_handles()) {
        if(triangle_is_inside(cdt, boundary, face)){
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

double get_sizing(Problem* problem) {
    Polygon boundary = problem->get_boundary();
    auto bbox = boundary.bbox();
    double max_size = std::sqrt(std::pow(bbox.x_span(), 2) + std::pow(bbox.y_span(), 2)) * 0.1;

    std::cout << RED << "xmin: " << bbox.xmin() << " scale: " << std::max(bbox.xmax() - bbox.xmin(), bbox.ymax() - bbox.ymin()) << RESET << std::endl;
    
    return max_size / 2;
}

double get_sizing_ratio_triangle(const Face_handle& triangle, double sizing){
    Point a = triangle->vertex(0)->point();
    Point b = triangle->vertex(1)->point();
    Point c = triangle->vertex(2)->point();
    double longest_edge_length = std::sqrt(CGAL::to_double(CGAL::max(CGAL::squared_distance(a, b), CGAL::max(CGAL::squared_distance(a, c), CGAL::squared_distance(b, c)))));

    return longest_edge_length/sizing;
}

double get_sizing_ratio_edge(const Edge& edge, double sizing){
    Point p1 = edge.first->vertex((edge.second + 1) % 3)->point();
    Point p2 = edge.first->vertex((edge.second + 2) % 3)->point();
    double edge_length = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p2)));
    return edge_length/sizing;
}

double compute_max_sizing_ratio(Problem* problem, CDT& cdt, double sizing){
    double max_sizing_ratio = 0;

    for(const auto& c : cdt.constrained_edges()){
        double sizing_ratio_constraint = get_sizing_ratio_edge(c, sizing);
        max_sizing_ratio = std::max(max_sizing_ratio, sizing_ratio_constraint);
    }   
    for(const auto& t : cdt.finite_face_handles()){
        if(triangle_is_inside(&cdt, problem->get_boundary(), t)){
            double sizing_ratio_triangle = get_sizing_ratio_triangle(t, sizing);
            max_sizing_ratio = std::max(max_sizing_ratio, get_sizing_ratio_triangle(t, sizing));
        }
    }
    return max_sizing_ratio;
}

void get_constraints_and_faces(Problem* problem, CDT& cdt, std::vector<Edge>& constraints, std::vector<Face_handle>& faces) {
    constraints.clear();
    for (const auto& constraint : cdt.constrained_edges()) { 
        constraints.push_back(constraint);
    }
    faces.clear();
    for (const auto& face : cdt.finite_face_handles()) {
        if(triangle_is_inside(&cdt, problem->get_boundary(), face)){
            faces.push_back(face);
        }
    }
}

// Custom hash function for CGAL::Edge
struct EdgeHash {
    std::size_t operator()(const Edge& e) const noexcept {
        std::hash<Face_handle> face_hash;
        std::hash<int> int_hash;
        return face_hash(e.first) ^ int_hash(e.second);
    }
};

// Custom equality function for CGAL::Edge
struct EdgeEqual {
    bool operator()(const Edge& e1, const Edge& e2) const noexcept {
        return e1.first == e2.first && e1.second == e2.second;
    }
};

std::unordered_map<Edge, double, EdgeHash, EdgeEqual> create_edge_map(CDT& cdt, double sizing){
    std::unordered_map<Edge, double, EdgeHash, EdgeEqual> edge_to_ratio_map;
    for(const auto& edge : cdt.constrained_edges()){
        edge_to_ratio_map[edge] = get_sizing_ratio_edge(edge, sizing);
    }

    return edge_to_ratio_map;
}

std::unordered_map<Face_handle, double> create_face_map(Problem* problem, CDT& cdt, double sizing){
    std::unordered_map<Face_handle, double> face_to_ratio_map;
    for(const auto& face : cdt.finite_face_handles()){
        if(triangle_is_inside(&cdt, problem->get_boundary(), face)){
            face_to_ratio_map[face] = get_sizing_ratio_triangle(face, sizing);
        }
    }

    return face_to_ratio_map;
}

void update_problem(Problem* problem, CDT& cdt, std::set<Point>& point_set){
    auto boundary = problem->get_boundary();
    problem->clear_solution();

    for (const auto& p : cdt.finite_vertex_handles()) {
        if (point_set.find(p->point()) == point_set.end()) {
            problem->add_steiner(p->point());
        }
    }
    for (auto triangle : grab_triangulation(&cdt, boundary)) {
        problem->add_triangle(Polygon(triangle.begin(), triangle.end()));
    }
}

void mesh(Problem* problem, CDT& cdt, double sizing, double target_sizing){

    std::unordered_map<Edge, double, EdgeHash, EdgeEqual> edge_to_ratio_map = create_edge_map(cdt, sizing);
    std::unordered_map<Face_handle, double> face_to_ratio_map = create_face_map(problem, cdt, sizing);

    std::cout << RED << "target ratio: " << target_sizing << RESET << std::endl;

    std::vector<Point> points_to_insert;

    for(Edge c : cdt.constrained_edges()){
        double current_ratio_constraint = edge_to_ratio_map[c];
        if(current_ratio_constraint > target_sizing){
            // Calculate new steiner point
            Point midpoint = CGAL::midpoint(c.first->vertex((c.second + 1) % 3)->point(), c.first->vertex((c.second + 2) % 3)->point());
            points_to_insert.push_back(midpoint);

            // Update sizing ratios
            Face_handle t1 = c.first;
            Face_handle t2 = c.first->neighbor(c.second);
            face_to_ratio_map[t1] += 0.5*(current_ratio_constraint/std::sqrt(3) - target_sizing);
            face_to_ratio_map[t2] += 0.5*(current_ratio_constraint/std::sqrt(3) - target_sizing);
        }
    }   
    for(Face_handle t : cdt.finite_face_handles()){
        if(!triangle_is_inside(&cdt, problem->get_boundary(), t)){
            continue;
        }

        double current_ratio_face = face_to_ratio_map[t];
        if(current_ratio_face > target_sizing){
            // Calculate new steiner point
            Point circumsenter = CGAL::circumcenter(t->vertex(0)->point(), t->vertex(1)->point(), t->vertex(2)->point());

            auto side = problem->get_boundary().oriented_side(circumsenter);
            if(side == CGAL::POSITIVE || side == CGAL::ZERO){
                points_to_insert.push_back(circumsenter);

                // Update sizing ratios
                Face_handle t1 = t->neighbor(0);
                Face_handle t2 = t->neighbor(1);
                Face_handle t3 = t->neighbor(2);
                face_to_ratio_map[t1] += (1/3)*(current_ratio_face/std::sqrt(3) - target_sizing);
                face_to_ratio_map[t2] += (1/3)*(current_ratio_face/std::sqrt(3) - target_sizing);
                face_to_ratio_map[t3] += (1/3)*(current_ratio_face/std::sqrt(3) - target_sizing);
            }
        }
    }
    cdt.insert(points_to_insert.begin(), points_to_insert.end());
}

void refine(Problem* problem){
    // Preprocess the input instance, compute cdt
    std::vector<Point> points = problem->get_points();
    std::set<Point> point_set(points.begin(), points.end());
    std::set<Point> steiner_point_set;
    std::vector<Point> boundary_points = problem->get_boundary().vertices();
    std::vector<Segment> constraints = problem->get_constraints();
    Polygon boundary = problem->get_boundary();
    for (size_t i = 0; i < boundary_points.size(); i++) {
        constraints.push_back(Segment(boundary_points[i], boundary_points[(i + 1) % boundary_points.size()]));
    }

    auto box = problem->get_boundary().bbox();
    auto scale = std::max(box.xmax() - box.xmin(), box.ymax() - box.ymin());

    std::cout << "scale: " << scale << " min: " << box.xmin() << "\n" << std::endl;

    CDT cdt = generate_triangulation(point_set, steiner_point_set, constraints);
    std::cout  << "Num of points in cdt before: " << cdt.number_of_vertices() << std::endl;

    double sizing, max_sizing, target_sizing;

    for(int i = 0; i<10; i++){
        sizing = get_sizing(problem);
        max_sizing = compute_max_sizing_ratio(problem, cdt, sizing);
        target_sizing = std::max(max_sizing/std::sqrt(3), 1.0);

        if(target_sizing == 1.0){
            break;
        }

        mesh(problem, cdt, sizing, target_sizing);
        std::cout << BLUE  << "Num of points in cdt after " << i << ": " << cdt.number_of_vertices()  << RESET << std::endl;

        update_problem(problem, cdt, point_set);
        problem->visualize_solution();

        // LLoyd optimization
        /*CGAL::lloyd_optimize_mesh_2(cdt, CGAL::parameters::number_of_iterations(1000));
        update_problem(problem, cdt, point_set);
        problem->visualize_solution();*/
    }

    std::cout << "End of meshing" << std::endl;
    //optimizeTinyAD(problem);
    // LLoyd optimization
    /*CGAL::lloyd_optimize_mesh_2(cdt, CGAL::parameters::number_of_iterations(1000));
    update_problem(problem, cdt, point_set);
    problem->visualize_solution();*/

    // Regenerate CDT after optimization
    /*steiner_point_set.clear();
    for(const auto& steiner : problem->get_steiner()){
        steiner_point_set.insert(steiner);
    }
    cdt = generate_triangulation(point_set, steiner_point_set, constraints);*/
}

void step_by_step_mesh(Problem* problem){
    // Preprocess the input instance, compute and compute cdt
    std::vector<Point> points = problem->get_points();
    std::set<Point> point_set(points.begin(), points.end());
    std::set<Point> steiner_point_set;
    std::vector<Point> boundary_points = problem->get_boundary().vertices();
    std::vector<Segment> constraints = problem->get_constraints();
    Polygon boundary = problem->get_boundary();
    for (size_t i = 0; i < boundary_points.size(); i++) {
        constraints.push_back(Segment(boundary_points[i], boundary_points[(i + 1) % boundary_points.size()]));
    }

    CDT cdt = generate_triangulation(point_set, steiner_point_set, constraints);
    update_problem(problem, cdt, point_set);
    std::cout  << "Num of points in cdt before: " << cdt.number_of_vertices() << std::endl;

    // Set up meshing
    double size = get_sizing(problem);
    Criteria criteria(0.125, size);
    //Criteria criteria(0.125, std::numeric_limits<double>::max());
    std::optional<Mesher> mesher;
    mesher.emplace(cdt, criteria);
    mesher->init();
    /*Mesher mesher(cdt, Criteria(0.125, size));
    mesher.init();*/

    int iter = 1;
    int opt_count = 1;

    while (iter != 100 && !mesher->is_refinement_done()) {
        mesher->step_by_step_refine_mesh();

        //update_problem(problem, cdt, point_set);
        //problem->visualize_solution();

        if(iter%10 == 0){
            // Optimize
            update_problem(problem, cdt, point_set);
            std::cout << "Num obtuse before optimization " << opt_count << ": " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution();
            optimizeTinyAD(problem);
            std::cout << "Num obtuse after optimization " << opt_count << ": " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution();

            // Update CDT
            steiner_point_set.clear();
            for(const auto& steiner : problem->get_steiner()){
                steiner_point_set.insert(steiner);
            }
            cdt = generate_triangulation(point_set, steiner_point_set, constraints);
            std::cout << "Num obtuse after regenerating CDT after optimization " << opt_count << ": " << count_obtuse_triangles(problem) << "\n" << std::endl;
            problem->visualize_solution();
            opt_count++;

            // Reinit the mesher
            mesher.emplace(cdt, criteria);
            mesher->init();     
        }

        iter++;
    }

    std::cout << "Total num iters: " << iter << std::endl;
}