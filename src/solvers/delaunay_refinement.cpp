#include "delaunay_refinement.hpp"
#include "lloyd.hpp"

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Boolean_set_operations_2/oriented_side.h>
#include <CGAL/lloyd_optimize_mesh_2.h>
#include <CGAL/Meshes/Double_map_container.h>

#include <unordered_map>
#include <limits>

#include "global_optimization.hpp"
#include "local_optimization.hpp"

typedef CGAL::Delaunay_mesh_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;

typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>            Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria>              Mesher;

typedef CDT::Face_handle Face_handle;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Vertex_circulator Vertex_circulator;
typedef CDT::Edge Edge;

double get_sizing(Problem* problem) {
    Polygon boundary = problem->get_boundary();
    auto bbox = boundary.bbox();
    double max_size = std::sqrt(std::pow(bbox.x_span(), 2) + std::pow(bbox.y_span(), 2)) * 0.1;

    std::cout << RED << "sizing: " << max_size/2 << std::endl;
    
    return max_size/(2);
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

double get_max_edge_ength(Problem* problem){
    double max_len = 0.0;
    for(const auto& t : problem->get_triangulation()){
        Point p1 = t[0];
        Point p2 = t[1];
        Point p3 = t[2];
        double len1 = CGAL::to_double(CGAL::squared_distance(p1, p2));
        double len2 = CGAL::to_double(CGAL::squared_distance(p1, p3));
        double len3 = CGAL::to_double(CGAL::squared_distance(p2, p3));

        max_len = std::max(max_len, len1);
        max_len = std::max(max_len, len2);
        max_len = std::max(max_len, len3);
    }
    return std::sqrt(max_len);
}

double compute_max_sizing_ratio(Problem* problem, CDT& cdt, double sizing){
    double max_sizing_ratio = 0;

    for(const auto& c : cdt.constrained_edges()){
        double sizing_ratio_constraint = get_sizing_ratio_edge(c, sizing);
        max_sizing_ratio = std::max(max_sizing_ratio, sizing_ratio_constraint);
    }   
    for(const auto& t : cdt.finite_face_handles()){
        if(problem->triangle_is_inside<CDT, Face_handle>(cdt, t)){
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
        if(problem->triangle_is_inside<CDT, Face_handle>(cdt, face)){
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
        if(problem->triangle_is_inside<CDT, Face_handle>(cdt, face)){
            face_to_ratio_map[face] = get_sizing_ratio_triangle(face, sizing);
        }
    }

    return face_to_ratio_map;
}

void mesh(Problem* problem, CDT& cdt, double sizing, double target_sizing){

    std::unordered_map<Edge, double, EdgeHash, EdgeEqual> edge_to_ratio_map = create_edge_map(cdt, sizing);
    std::unordered_map<Face_handle, double> face_to_ratio_map = create_face_map(problem, cdt, sizing);

    //std::cout << RED << "target ratio: " << target_sizing << RESET << std::endl;

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
            if(cdt.is_infinite(t1)){
                face_to_ratio_map[t2] += (current_ratio_constraint/std::sqrt(3) - target_sizing);
            } else if(cdt.is_infinite(t2)){
                face_to_ratio_map[t1] += (current_ratio_constraint/std::sqrt(3) - target_sizing);
            } else {
                face_to_ratio_map[t1] += 0.5*(current_ratio_constraint/std::sqrt(3) - target_sizing);
                face_to_ratio_map[t2] += 0.5*(current_ratio_constraint/std::sqrt(3) - target_sizing);
            }   
        }
    }   
    for(Face_handle t : cdt.finite_face_handles()){
        if(!problem->triangle_is_inside<CDT, Face_handle>(cdt, t)){
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
                if(cdt.is_infinite(t1)){
                    if(cdt.is_infinite(t2)){
                        face_to_ratio_map[t3] += (current_ratio_face/std::sqrt(3) - target_sizing);
                    } else if(cdt.is_infinite(t3)){
                        face_to_ratio_map[t2] += (current_ratio_face/std::sqrt(3) - target_sizing);
                    } else {
                        face_to_ratio_map[t2] += (1/2)*(current_ratio_face/std::sqrt(3) - target_sizing);
                        face_to_ratio_map[t3] += (1/2)*(current_ratio_face/std::sqrt(3) - target_sizing);
                    }
                } else if(cdt.is_infinite(t2)){
                    if(cdt.is_infinite(t3)){
                        face_to_ratio_map[t1] += (current_ratio_face/std::sqrt(3) - target_sizing);
                    } else {
                        face_to_ratio_map[t1] += (1/2)*(current_ratio_face/std::sqrt(3) - target_sizing);
                        face_to_ratio_map[t3] += (1/2)*(current_ratio_face/std::sqrt(3) - target_sizing);
                    }
                } else if(cdt.is_infinite(t3)){
                    face_to_ratio_map[t1] += (1/2)*(current_ratio_face/std::sqrt(3) - target_sizing);
                    face_to_ratio_map[t2] += (1/2)*(current_ratio_face/std::sqrt(3) - target_sizing);
                } else {
                    face_to_ratio_map[t1] += (1/3)*(current_ratio_face/std::sqrt(3) - target_sizing);
                    face_to_ratio_map[t2] += (1/3)*(current_ratio_face/std::sqrt(3) - target_sizing);
                    face_to_ratio_map[t3] += (1/3)*(current_ratio_face/std::sqrt(3) - target_sizing);
                }
            }
        }
    }
    cdt.insert(points_to_insert.begin(), points_to_insert.end());
}

bool next_refinement_point(Problem* problem, CDT& cdt){

    auto boundary = problem->get_boundary();

    for(Face_handle t : cdt.finite_face_handles()){
        if(!problem->triangle_is_inside<CDT, Face_handle>(cdt, t)){
            continue;
        }

        Polygon triangle;
        triangle.push_back(t->vertex(0)->point());
        triangle.push_back(t->vertex(1)->point());
        triangle.push_back(t->vertex(2)->point());

        int obtuse_angle = find_obtuse_angle(triangle);
        
        if(obtuse_angle != -1){
            Point s = triangle[obtuse_angle];
            Point a = triangle[(obtuse_angle + 1)%3];
            Point b = triangle[(obtuse_angle + 2)%3];

            if(is_on_boundary_refinement(a, boundary) && is_on_boundary_refinement(b, boundary) && is_boundary_segment(boundary, a, b)){
                Segment segment(a, b);
                Point p = project_onto_segment_refinement(s, segment);
                problem->add_steiner(p);
            } else {
                // Calculate new steiner point
                Point circumsenter = CGAL::circumcenter(t->vertex(0)->point(), t->vertex(1)->point(), t->vertex(2)->point());

                auto side = problem->get_boundary().oriented_side(circumsenter);
                if(side == CGAL::POSITIVE || side == CGAL::ZERO){
                    cdt.insert(circumsenter);
                    return false;
                }
            }
        }
    }
    return true;
}

Mesh_Statistics refine(Problem* problem){
    Mesh_Statistics stats;

    std::vector<Point> points = problem->get_points();
    std::set<Point> point_set(points.begin(), points.end());
    auto constraints = problem->get_constraints();

    CDT cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);

    /*mesh_cgal(problem);
    cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);*/

    //problem->visualize_solution({});

    double sizing = get_sizing(problem);
    double max_sizing, target_sizing;
    int i = 0;

    while(true){
        max_sizing = compute_max_sizing_ratio(problem, cdt, sizing);
        target_sizing = std::max(max_sizing/std::sqrt(3), 1.0);

        if(max_sizing <= 1.0){
            break;
        }

        mesh(problem, cdt, sizing, target_sizing);
        //std::cout << BLUE  << "Num of points in cdt after " << i << ": " << cdt.number_of_vertices()  << RESET << std::endl;
        problem->update_problem<CDT, Face_handle>(cdt, point_set);
        //problem->visualize_solution({});

        // LLoyd optimization
        /*iterate_lloyd<CDT, Vertex_handle, Vertex_circulator>(cdt, problem, constraints, 1000);
        problem->update_problem<CDT, Face_handle>(cdt, point_set);*/

        // Energy
        /*optimizeTinyAD(problem);
        cdt = problem->generate_CDT<CDT>();
        problem->update_problem<CDT, Face_handle>(cdt, point_set);
        optimizeTinyAD(problem);
        cdt = problem->generate_CDT<CDT>();
        problem->update_problem<CDT, Face_handle>(cdt, point_set);
        optimizeTinyAD(problem);
        cdt = problem->generate_CDT<CDT>();
        problem->update_problem<CDT, Face_handle>(cdt, point_set);*/

        //problem->visualize_solution({});

        //!For this first need to convert to an inexact kernel
        CGAL::lloyd_optimize_mesh_2(cdt, CGAL::parameters::number_of_iterations(1000));
        problem->update_problem<CDT, Face_handle>(cdt, point_set);
        //problem->visualize_solution({});
        i++;
    }

    /*iterate_lloyd<CDT, Vertex_handle, Vertex_circulator>(cdt, problem, constraints, 1000);
    problem->update_problem<CDT, Face_handle>(cdt, point_set);*/

    stats.set_name(problem->get_name());
    stats.set_steiner_after_meshing(problem->get_steiner().size());
    stats.set_deviation(mean_absolute_deviation(problem));
    stats.set_aspect_ratio(mean_aspect_ratio(problem));
    save_min_max_angle(problem, &stats);

    stats.set_max_length(get_max_edge_ength(problem));

    //problem->visualize_solution({});

    std::cout << "num steiner: "  << problem->get_steiner().size() << std::endl;
    std::cout << "deviation: "  << mean_absolute_deviation(problem) << std::endl;
    std::cout << "max_edge_length: "  << get_max_edge_ength(problem) << RESET << std::endl;

    save_angle_stats_for_plot(problem, "../results/angle_stats_lloyd_simple-polygon_20_4bd3c2e5");
    save_aspect_ratios_for_plot(problem, "../results/aspect_ratios_lloyd_simple-polygon_20_4bd3c2e5");

    return stats;
}

void classic_delaunay_refinement(Problem* problem){
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

    CDT cdt = problem->generate_CDT<CDT>();
    std::cout  << "Num of points in cdt before: " << cdt.number_of_vertices() << std::endl;

    problem->update_problem<CDT, Face_handle>(cdt, point_set);
    problem->visualize_solution({});

    for(int i = 0; i < 30; i++){
        next_refinement_point(problem, cdt);
        problem->update_problem<CDT, Face_handle>(cdt, point_set);
    }

    fix_boundary_refinement(problem);
    cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);
    fix_boundary_refinement(problem);
    cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);
    problem->visualize_solution({});

    iterate_lloyd<CDT, Vertex_handle, Vertex_circulator>(cdt, problem, constraints, 1000);
    problem->update_problem<CDT, Face_handle>(cdt, point_set);
    optimizeTinyAD(problem);
    cdt = problem->generate_CDT<CDT>();
    problem->visualize_solution({});

    fix_boundary_refinement(problem);
    cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);
    problem->visualize_solution({});
}

//***************************CGAL step-for-step********************************************************

Point project_onto_segment_refinement(const Point p, const Segment& s) {
    Line line(s.source(), s.target()); 
    Point proj = line.projection(p); 

    // ! add sqrt?
    // Clamp projection to segment endpoints
    if (CGAL::squared_distance(proj, s.source()) + CGAL::squared_distance(proj, s.target()) 
        > CGAL::squared_distance(s.source(), s.target())) {
        throw std::runtime_error("Projected point is not on segment!");
    }
    return proj;
}

bool is_on_constraint_refinement(const Point& steiner, Problem* problem){
    std::vector<Segment> constraints = problem->get_constraints();
    for(Segment c : constraints){
        if(c.has_on(steiner)){
            return true;
        }
    }

    return false;
}

bool is_on_boundary_refinement(Point& p, Polygon& boundary){
    double tolerance = 1e-12;
    for (auto edge = boundary.edges_begin(); edge != boundary.edges_end(); ++edge) {
        double distance = CGAL::to_double(CGAL::squared_distance(*edge, p));
        if (distance <= tolerance || edge->source() == p || edge->target() == p) {
            return true;
        }
    }
    return false;
}

bool is_boundary_segment(Polygon& boundary, Point& a, Point& b){
    double tolerance = 1e-12;
    for (auto edge = boundary.edges_begin(); edge != boundary.edges_end(); ++edge) {
        // Check like this because a and b might not be input points and thus are not end points of the boundary polygon
        // which only containsinput points and is not updated i. e. split into subsegments after steiner point insertion
        double distance_a = CGAL::to_double(CGAL::squared_distance(*edge, a));
        double distance_b = CGAL::to_double(CGAL::squared_distance(*edge, b));
        if ((distance_a < tolerance || edge->has_on(a)) && (distance_b < tolerance || edge->has_on(b))) {
        //if((a == edge->source() && b == edge->target()) || (a == edge->target() || b == edge->source())){
            return true;
        }
    }
    return false;
}

// First two points are on the boundary
std::pair<Point, bool> find_boundary_point(Problem* problem, Point& a, Point& b, Point& c){
    auto triangles = problem->get_triangulation();
    auto boundary = problem->get_boundary();
    for(Polygon t : triangles){
        Point p0 = t[0];
        Point p1 = t[1];
        Point p2 = t[2];
        if(((p0 == b && p1 == c && p2 != a) || (p1 == b && p0 == c && p2 != a)) && is_on_boundary_refinement(p2, boundary)){
            return {p2, true};
        } else if(((p1 == b && p2 == c && p0 != a) || (p2 == b && p1 == c && p0 != a)) && is_on_boundary_refinement(p0, boundary)){
            return {p0, true};
        } else if(((p2 == b && p0 == c && p1 != a) || (p0 == b && p2 == c && p1 != a)) && is_on_boundary_refinement(p1, boundary)){
            return {p1, true};
        }
    }
    //std::cout << RED << "Boundary point for boundary fix was not found!" << RESET << std::endl;
    return {Point(0, 0), false};
}

bool point_within_segment(Point& a, Point& b, Point& p){
    double tolerance = 1e-10;
    Segment segment1(a, b);
    Segment segment2(b, a);
    double dist1 = CGAL::to_double(CGAL::squared_distance(segment1, p));
    double dist2 = CGAL::to_double(CGAL::squared_distance(segment2, p));
    return dist1 < tolerance || segment1.has_on(p) || dist2 < tolerance || segment2.has_on(p);
}

void fix_boundary_refinement(Problem *problem){
    std::vector<Polygon> triangulation = problem->get_triangulation();
    std::vector<Point> points = problem->get_points();
    Polygon boundary = problem->get_boundary();

    Segment constraint;

    for(Polygon t : triangulation){
        int angle = find_obtuse_angle(t);
        if(angle != -1){
            Point s = t.vertex(angle);
            Point a = t.vertex((angle + 1)%3);
            Point b = t.vertex((angle + 2)%3);

            // Case 1 (for boundary)
            if(is_on_boundary_refinement(a, boundary) && is_on_boundary_refinement(b, boundary) && is_boundary_segment(boundary, a, b)){
                Segment segment(a, b);
                try{
                    Point p = project_onto_segment_refinement(s, segment);
                    problem->add_steiner(p);

                    int new_obtuse = problem->get_num_obtuse() - 1;
                    problem->set_num_obtuse(new_obtuse);
                } catch(const std::runtime_error& e){
                    // Do not insert any new points, case could not be fixed
                    std::cout << RED << "Projected point is not on segment!"<< RESET << std::endl;
                }
            } 
            // Case 2
            else if (is_on_boundary_refinement(s, boundary)){
                //bool a_boundary = is_on_boundary_refinement(a, boundary);
                //bool b_boundary = is_on_boundary_refinement(b, boundary);

                bool a_boundary = is_boundary_segment(boundary, a, s);
                bool b_boundary = is_boundary_segment(boundary, b, s);

                Point c;
                bool is_valid;

                if(a_boundary){
                    auto [result, res_valid] = find_boundary_point(problem, a, s, b);
                    if(res_valid){
                        is_valid = res_valid;
                        c = result;
                    }
                } else if(b_boundary){
                    auto [result, res_valid] = find_boundary_point(problem, b, s, a);
                    if(res_valid){
                        is_valid = res_valid;
                        c = result;
                    }
                } else {
                    continue;
                }

                bool first_option = false;

                // if s is not constrained and can be moved
                if(std::find(points.begin(), points.end(), s) == points.end()){
                    if(a_boundary && is_valid){
                        Line line(s, a); 
                        Point proj = line.projection(b);
                        // check whether the new point position is within the corresponding segment
                        if(point_within_segment(a, c, proj)){
                            //std::cout << "aaa" << std::endl;
                            problem->update_steiner(s, proj);
                            first_option = true;

                            int new_obtuse = problem->get_num_obtuse() - 1;
                            problem->set_num_obtuse(new_obtuse);
                        }
                    } else if(b_boundary && is_valid){
                        Line line(s, b); 
                        Point proj = line.projection(a); 
                        // check whether the new point position is within the corresponding segment
                        if(point_within_segment(b, c, proj)){
                            //std::cout << "bbb" << std::endl;
                            problem->update_steiner(s, proj);
                            first_option = true;

                            int new_obtuse = problem->get_num_obtuse() - 1;
                            problem->set_num_obtuse(new_obtuse);
                        }
                    }
                } 

                // Here, do not alter the number of obtuse because we do not know whether it helps
                // Keep in mind that because of this the number of obtuse is not accurate and needs to be recomputed for statistics
                if (!first_option){
                    if(a_boundary){ // try to move b
                        if(std::find(points.begin(), points.end(), b) == points.end() && !is_on_boundary_refinement(b, boundary) && !is_on_constraint_refinement(b, problem)){
                            Line line_parallel(b, a - s);
                            Point proj = line_parallel.projection(s);
                            //problem->add_steiner(proj);
                            problem->update_steiner(b, proj);
                        }
                    } else if (b_boundary){ // try to move a
                        if(std::find(points.begin(), points.end(), a) == points.end() && !is_on_boundary_refinement(a, boundary) && !is_on_constraint_refinement(a, problem)){
                            Line line_parallel(a, b - s);
                            Point proj = line_parallel.projection(s);
                            //problem->add_steiner(proj);
                            problem->update_steiner(a, proj);
                        }
                    }
                }
            }
        }
        
    }
}

void globally_optimize_triangles(Problem *problem, bool debug){
    std::vector<Point> steiner = problem->get_steiner();
    std::vector<Polygon> triangulation = problem->get_triangulation();

    std::vector<Point> new_steiner = globally_optimize_position(steiner, triangulation, problem, debug);

    problem->clear_solution();
    for(int i = 0; i < new_steiner.size(); i++){
        problem->add_steiner(new_steiner[i]);
    }
    problem->set_triangulation(triangulation);
}

void locally_optimize_triangles(Problem *problem){
    std::vector<Polygon> triangulation = problem->get_triangulation();
    std::vector<Point> steiner;

    for(Point s : problem->get_steiner()){
        Point new_s = locally_optimize_position(s, triangulation, problem);
        steiner.push_back(new_s);
    }
    problem->clear_solution();
    for(int i = 0; i < steiner.size(); i++){
        problem->add_steiner(steiner[i]);
    }
    problem->set_triangulation(triangulation);
}

void locally_optimize_obtuse_refinement(Problem *problem){
    std::vector<Polygon> triangulation = problem->get_triangulation();
    std::vector<Point> points = problem->get_points();
    std::vector<Point> steiner;

    for(Polygon t : triangulation){
        if(is_obtuse_triangle(t)){
            Point a = t.vertex(0);
            Point b = t.vertex(1);
            Point c = t.vertex(2);

            if(std::find(points.begin(), points.end(), a) == points.end()){
                Point new_s1 = locally_optimize_position(a, triangulation, problem);
                steiner.push_back(new_s1);
            }
            if(std::find(points.begin(), points.end(), b) == points.end()){
                Point new_s2 = locally_optimize_position(b, triangulation, problem);
                steiner.push_back(new_s2);
            }
            if(std::find(points.begin(), points.end(), c) == points.end()){
                Point new_s3 = locally_optimize_position(c, triangulation, problem);
                steiner.push_back(new_s3);
            }
        }
    }
    
    /*problem->clear_solution();
    for(int i = 0; i < steiner.size(); i++){
        problem->add_steiner(steiner[i]);
    }*/
    problem->set_triangulation(triangulation);
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

    CDT cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);
    std::cout  << "Num of points in cdt before: " << cdt.number_of_vertices() << std::endl;

    // Set up meshing
    double size = get_sizing(problem);
    Criteria criteria(0, size);
    //Criteria criteria(0.125, std::numeric_limits<double>::max());
    std::optional<Mesher> mesher;
    mesher.emplace(cdt, criteria);
    mesher->init();
    /*Mesher mesher(cdt, Criteria(0.125, size));
    mesher.init();*/

    int iter = 1;
    int opt_count = 1;

    while (iter != 40 && !mesher->is_refinement_done() && (count_obtuse_triangles(problem) != 0)) {
        mesher->step_by_step_refine_mesh();

        problem->update_problem<CDT, Face_handle>(cdt, point_set);

        /*if(iter >= 20){
            problem->update_problem<CDT, Face_handle>(cdt, point_set);
            problem->visualize_solution({});
        }*/

        if(iter%10 == 0){
            // Optimize
            std::cout << "Num obtuse before optimization " << opt_count << ": " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});
            //bool debug = (opt_count == 9);
            //globally_optimize_triangles(problem, debug);
            optimizeTinyAD(problem);

            std::cout << "Num obtuse after optimization " << opt_count << ": " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});

            // Update CDT
            cdt = problem->generate_CDT<CDT>();
            std::cout << "Num obtuse after regenerating CDT after optimization " << opt_count << ": " << count_obtuse_triangles(problem) << std::endl;
            problem->update_problem<CDT, Face_handle>(cdt, point_set);
            //problem->visualize_solution({});

            // Fix boundary
            fix_boundary_refinement(problem);
            cdt = problem->generate_CDT<CDT>();
            problem->update_problem<CDT, Face_handle>(cdt, point_set);
            std::cout << "Num obtuse after fixing boundary " << opt_count << ": " << count_obtuse_triangles(problem) << "\n" << std::endl;
            problem->visualize_solution({});

            opt_count++;

            mesher.emplace(cdt, criteria);
            mesher->init();    
        }

        iter++;
    }

    std::cout << "Total num iters: " << iter << std::endl;
}

Mesh_Statistics uniform_mesh(Problem* problem){
    std::vector<Point> points = problem->get_points();
    std::set<Point> point_set(points.begin(), points.end());
    std::set<Point> steiner_point_set;
    std::vector<Point> boundary_points = problem->get_boundary().vertices();
    std::vector<Segment> constraints = problem->get_constraints();
    Polygon boundary = problem->get_boundary();
    for (size_t i = 0; i < boundary_points.size(); i++) {
        constraints.push_back(Segment(boundary_points[i], boundary_points[(i + 1) % boundary_points.size()]));
    }

    CDT cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);

    auto bbox = boundary.bbox();
    double max_size = std::sqrt(pow(bbox.x_span(),2) + pow(bbox.y_span(),2)) * 0.1;

    Mesher mesher(cdt);
    mesher.set_criteria(Criteria(0.125, max_size/2));
    mesher.refine_mesh();

    problem->update_problem<CDT, Face_handle>(cdt, point_set);

    int obtuse_meshing = count_obtuse_triangles(problem);
    problem->set_num_obtuse(obtuse_meshing);

    Mesh_Statistics statistics;
    statistics.set_name(problem->get_name());
    statistics.set_obtuse_after_meshing(obtuse_meshing);
    statistics.set_steiner_after_meshing(problem->get_steiner().size());

    //std::string path1 = "../results/MESH-" + problem->get_name() + ".svg";
    //problem->save_intermidiate_result(path1);

    optimizeTinyAD(problem);
    cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);
    optimizeTinyAD(problem);
    cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);

    int obtuse_opt = count_obtuse_triangles(problem);
    problem->set_num_obtuse(obtuse_opt);

    int obtuse_after_fix = obtuse_opt;

    int i = 0;

    //while(i != 100){
        fix_boundary_refinement(problem);
        /*int current_obtuse = problem->get_num_obtuse();
        if(current_obtuse < obtuse_after_fix){
            obtuse_after_fix = current_obtuse;
        } else {
            break;
        }
        //std::cout << RED << "---------------------------------------------" << RESET << std::endl;
        i++;
    }*/

    cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);

    std::cout << RED << "num obtuse: " << count_obtuse_triangles(problem) << std::endl;

    /*locally_optimize_obtuse_refinement(problem);
    cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);

    fix_boundary_refinement(problem);
    cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);*/

    obtuse_opt = count_obtuse_triangles(problem);
    problem->set_num_obtuse(obtuse_opt);

    statistics.set_obtuse_after_optimization(obtuse_opt);
    statistics.set_steiner_after_optimization(problem->get_steiner().size());

    std::string path2 = "../results/OPT-" + problem->get_name() + ".svg";
    problem->save_intermidiate_result(path2);

    std::cout << RED << "num obtuse: " << obtuse_opt << std::endl;
    std::cout << "num steiner: " << problem->get_steiner().size() << RESET << std::endl;

    return statistics;
}

//***********************************************Equilateral********************************************************

void mesh_cgal(Problem* problem){
    std::vector<Point> points = problem->get_points();
    std::set<Point> point_set(points.begin(), points.end());
    std::set<Point> steiner_point_set;
    std::vector<Point> boundary_points = problem->get_boundary().vertices();
    std::vector<Segment> constraints = problem->get_constraints();
    Polygon boundary = problem->get_boundary();
    for (size_t i = 0; i < boundary_points.size(); i++) {
        constraints.push_back(Segment(boundary_points[i], boundary_points[(i + 1) % boundary_points.size()]));
    }

    CDT cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);

    auto bbox = boundary.bbox();
    double max_size = std::sqrt(pow(bbox.x_span(),2) + pow(bbox.y_span(),2)) * 0.1;

    Mesher mesher(cdt);
    mesher.set_criteria(Criteria(0.125, max_size/2));
    mesher.refine_mesh();

    problem->update_problem<CDT, Face_handle>(cdt, point_set);
    //problem->visualize_solution({});
}

double angle_degrees(const Point& a, const Point& b, const Point& c) {
    Vector u = a - b;
    Vector v = c - b;
    double dot = CGAL::to_double(u * v);
    double norm_u = std::sqrt(CGAL::to_double(u.squared_length()));
    double norm_v = std::sqrt(CGAL::to_double(v.squared_length()));
    double cos_angle = std::max(-1.0, std::min(1.0, dot / (norm_u * norm_v))); // clamp
    return std::acos(cos_angle) * 180.0 / M_PI;
}

double get_triangle_area(const Point& p1, const Point& p2, const Point& p3){
    auto area = (p1.x() * (p2.y() - p3.y()) + p2.x() * (p3.y() - p1.y()) + p3.x() * (p1.y() - p2.y()))/2.0;
    return std::abs(CGAL::to_double(area));
}

double get_shortest_edge(const Point& a, const Point& b, const Point& c){
    auto ab = CGAL::squared_distance(a, b);
    auto bc = CGAL::squared_distance(b, c);
    auto ac = CGAL::squared_distance(a, c);

    if(ab <= bc && ab <= ac){
        return std::sqrt(CGAL::to_double(ab));
    } else if (bc <= ab && bc <= ac) {
        return std::sqrt(CGAL::to_double(bc));
    } else {
        return std::sqrt(CGAL::to_double(ac));
    }
}

double get_aspect_ratio(const Point& a, const Point& b, const Point& c){
    if(get_triangle_area(a, b, c) < 1e-6){
        return INFINITY;
    }
    double shortest_edge = get_shortest_edge(a, b, c);
    if(shortest_edge < 1e-6){
        return INFINITY;
    }
    double radius = std::sqrt(CGAL::to_double(CGAL::squared_radius(a, b, c)));  
    return radius/shortest_edge;
}

void save_aspect_ratios_for_plot(Problem* problem, std::string path){
    auto triangulation = problem->get_triangulation();

    std::vector<double> ratios;

    for (const Polygon& triangle : triangulation) {
        if (triangle.size() != 3) continue; // skip non-triangles

        const Point& A = triangle[0];
        const Point& B = triangle[1];
        const Point& C = triangle[2];

        ratios.push_back(get_aspect_ratio(A, B, C));
    }

    std::ofstream out(path);
    out << "ratio\n"; 
    for (double ratio : ratios) {
        out << ratio << "\n";
    }
    out.close();
}

double mean_aspect_ratio(Problem* problem){
    double sum = 0.0;

    auto triangulation = problem->get_triangulation();

    std::vector<double> ratios;

    for (const Polygon& triangle : triangulation) {
        if (triangle.size() != 3) continue;

        const Point& A = triangle[0];
        const Point& B = triangle[1];
        const Point& C = triangle[2];

        sum += get_aspect_ratio(A, B, C);
    }

    return sum/(triangulation.size());
}

void save_angle_stats_for_plot(Problem* problem, std::string path){
    auto triangulation = problem->get_triangulation();

    std::vector<double> all_angles;

    for (const Polygon& triangle : triangulation) {
        if (triangle.size() != 3) continue; // skip non-triangles

        const Point& A = triangle[0];
        const Point& B = triangle[1];
        const Point& C = triangle[2];

        all_angles.push_back(angle_degrees(C, A, B));
        all_angles.push_back(angle_degrees(A, B, C));
        all_angles.push_back(angle_degrees(B, C, A));
    }

    std::ofstream out(path);
    out << "angle\n"; 
    for (double angle : all_angles) {
        out << angle << "\n";
    }
    out.close();
}

std::vector<Polygon> save_min_max_angle(Problem* problem, Mesh_Statistics* stats){
    auto triangulation = problem->get_triangulation();

    std::vector<Polygon> problematic_triangles;

    double min_angle = 360.0;
    double max_angle = 0.0;

    for (const Polygon& triangle : triangulation) {
        if (triangle.size() != 3) continue; // skip non-triangles

        const Point& A = triangle[0];
        const Point& B = triangle[1];
        const Point& C = triangle[2];

        auto angle_a = angle_degrees(C, A, B);
        auto angle_b = angle_degrees(A, B, C);
        auto angle_c = angle_degrees(B, C, A);

        /*if(angle_a<1.0 || angle_b<1.0 || angle_c<1.0){
            problematic_triangles.push_back(triangle);
        }
        if(angle_a < min_angle && angle_a>=1.0) min_angle = angle_a;
        if(angle_b < min_angle && angle_b>=1.0) min_angle = angle_b;
        if(angle_c < min_angle && angle_c>=1.0) min_angle = angle_c;*/

        if(angle_a < min_angle) min_angle = angle_a;
        if(angle_b < min_angle) min_angle = angle_b;
        if(angle_c < min_angle) min_angle = angle_c;

        if(angle_a > max_angle) max_angle = angle_a;
        if(angle_b > max_angle) max_angle = angle_b;
        if(angle_c > max_angle) max_angle = angle_c;
    }

    stats->set_min_angle(min_angle);
    stats->set_max_angle(max_angle);
    return problematic_triangles;
}

double mean_absolute_deviation(Problem* problem){
    auto triangulation = problem->get_triangulation();
    double sum = 0.0;

    for(const Polygon& triangle : triangulation){
        if (triangle.size() != 3) continue; // skip non-triangles

        const Point& A = triangle[0];
        const Point& B = triangle[1];
        const Point& C = triangle[2];

        auto angle_a = angle_degrees(C, A, B);
        auto angle_b = angle_degrees(A, B, C);
        auto angle_c = angle_degrees(B, C, A);

        sum = sum + std::abs(angle_a - 60) + std::abs(angle_b - 60) +std::abs(angle_c - 60);
    }
    sum = sum/(3*triangulation.size());
    return sum;
}

void mesh_equilateral_single(Problem* problem){
    auto constraints = problem->get_constraints();
    std::vector<Point> points = problem->get_points();
    std::set<Point> point_set(points.begin(), points.end());
    mesh_cgal(problem);

    std::cout << "num initial: " << problem->get_points().size() << " num steiner: " << problem->get_steiner().size() << std::endl;
    std::cout << "mean absolute deviation initial: " << mean_absolute_deviation(problem) << std::endl;
    // Stats after initial meshing
    //save_angle_stats_for_plot(problem, "../results/angle_data_initial.csv");
    //save_aspect_ratios_for_plot(problem, "../results/aspect_ratios_initial.csv");
    //problem->visualize_solution({});

    // Stats after optimization meshing
    /*CDT cdt = problem->generate_CDT<CDT>();
    iterate_lloyd<CDT, Vertex_handle, Vertex_circulator>(cdt, problem, constraints, 1000);
    problem->update_problem<CDT, Face_handle>(cdt, point_set);
    /*CGAL::lloyd_optimize_mesh_2(cdt, CGAL::parameters::number_of_iterations(1000));
    problem->update_problem<CDT, Face_handle>(cdt, point_set);*/
    
    optimizeTinyAD(problem);
    CDT cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);
    optimizeTinyAD(problem);
    cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);
    optimizeTinyAD(problem);
    cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);
    std::cout << "mean absolute deviation lloyd: " << mean_absolute_deviation(problem) << std::endl;

    //save_angle_stats_for_plot(problem, "../results/angle_data_lloyd.csv");
    //save_aspect_ratios_for_plot(problem, "../results/aspect_ratios_lloyd.csv");
    problem->visualize_solution({});
}

Mesh_Statistics mesh_equilateral(Problem* problem){
    Mesh_Statistics stats;

    std::vector<Point> points = problem->get_points();
    auto constraints = problem->get_constraints();
    std::set<Point> point_set(points.begin(), points.end());

    mesh_cgal(problem);
    CDT cdt = problem->generate_CDT<CDT>();
    //refine(problem, cdt);
    //problem->visualize_solution({});

    /*optimizeTinyAD(problem);
    cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);
    optimizeTinyAD(problem);
    cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);*/
    /*optimizeTinyAD(problem);
    cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);*/

    cdt = problem->generate_CDT<CDT>();
    iterate_lloyd<CDT, Vertex_handle, Vertex_circulator>(cdt, problem, constraints, 1000);
    problem->update_problem<CDT, Face_handle>(cdt, point_set);

    stats.set_name(problem->get_name());
    stats.set_steiner_after_meshing(problem->get_steiner().size());
    stats.set_deviation(mean_absolute_deviation(problem));
    std::vector<Polygon> problematic_triangles = save_min_max_angle(problem, &stats);

    //problem->visualize_solution({});

    //std::cout << "deviation: " << stats.get_deviation() << " min angle: " << stats.get_min_angle() << " max angle: " << stats.get_max_angle() << std::endl;

    return stats;
}