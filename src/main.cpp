#include "main.hpp"
#include "local_optimization.hpp"
#include "global_optimization.hpp"
#include "tinyAD_optimization.hpp"
#include "delaunay_refinement.hpp"
#include "quad_mesher.hpp"

#include <gmsh.h>

namespace fs = std::filesystem;

Segment find_longest_edge(Polygon& triangle){
    auto ab = CGAL::squared_distance(triangle.vertex(0), triangle.vertex(1));
    auto bc = CGAL::squared_distance(triangle.vertex(1), triangle.vertex(2));
    auto ac = CGAL::squared_distance(triangle.vertex(0), triangle.vertex(2));

    if(ab >= bc && ab >= ac){
        return Segment(triangle.vertex(0), triangle.vertex(1));
    } else if (bc > ab && bc >= ac){
        return Segment(triangle.vertex(1), triangle.vertex(2));
    } else {
        return Segment(triangle.vertex(0), triangle.vertex(2));
    }
}

bool contains_edge(const Polygon& triangle, const Point& p1, const Point& p2) {
    for (size_t i = 0; i < 3; ++i) {
        size_t j = (i + 1) % 3;  
        if ((triangle[i] == p1 && triangle[j] == p2) || 
            (triangle[i] == p2 && triangle[j] == p1)) {
            return true;  
        }
    }
    return false;
}

Point find_opposite_point(Polygon& triangle, Segment& longestEdge){
    for (int i = 0; i < 3; i++) {
        if(triangle[i] != longestEdge.source() && triangle[i] != longestEdge.target()){
            return triangle[i];
        }
    }
    throw std::runtime_error("Opposite point not found!");
}

double obtuseness_factor(Polygon& triangle) {
    Point a = triangle.vertex(0);
    Point b = triangle.vertex(1);
    Point c = triangle.vertex(2);

    double ab = CGAL::to_double(CGAL::squared_distance(a, b));
    double bc = CGAL::to_double(CGAL::squared_distance(b, c));
    double ac = CGAL::to_double(CGAL::squared_distance(c, a));

    // Find the maximum squared edge length
    double max_sq = std::max({ab, bc, ac});

    // Compute obtuseness factor
    double obtuseness = max_sq - (ab + bc + ac - max_sq);

    return obtuseness; // If positive, the triangle is obtuse
}

void flip(Polygon& triangle, std::vector<Polygon>& triangulation, bool ignore){
    Segment longestEdge = find_longest_edge(triangle);

    std::vector<Polygon> removed_triangles;

    for (const auto& t : triangulation) {
        if (contains_edge(t, longestEdge.source(), longestEdge.target())) {
            removed_triangles.push_back(t);
        }
    }

    if(removed_triangles.size() == 2){
        triangulation.erase(std::remove_if(triangulation.begin(), triangulation.end(),
        [&](const Polygon& tri) { return contains_edge(tri, longestEdge.source(), longestEdge.target()); }),
        triangulation.end());

        Point a = find_opposite_point(removed_triangles[0], longestEdge);
        Point b = find_opposite_point(removed_triangles[1], longestEdge);

        double obtuseness0 = obtuseness_factor(triangle);

        Polygon t1;
        t1.push_back(a);
        t1.push_back(b);
        t1.push_back(longestEdge.source());
        
        Polygon t2;
        t2.push_back(a);
        t2.push_back(b);
        t2.push_back(longestEdge.target());

        double obtuseness1 = std::max(obtuseness_factor(t1), obtuseness_factor(t2));

        if(obtuseness1 < obtuseness0 || ignore){
            triangulation.push_back(t1);
            triangulation.push_back(t2);
        } else {
            triangulation.push_back(removed_triangles[0]);
            triangulation.push_back(removed_triangles[1]);
        }
    }
}

void perform_edge_flips(Problem *problem, bool ignore){
    std::vector<Polygon> triangulation = problem->get_triangulation();

    std::vector<Polygon> to_flip = triangulation;
    for(Polygon& triangle : to_flip){
        if(is_obtuse_triangle(triangle)){
            flip(triangle, triangulation, ignore);
        }
    }

    problem->set_triangulation(triangulation);
}

void locally_optimize_solution(Problem *problem){
    std::vector<Polygon> triangulation = problem->get_triangulation();
    std::vector<Point> steiner;

    //std::cout << "boundary size: " << (problem->get_boundary()).size() << std::endl;

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

void locally_optimize_obtuse(Problem *problem){
    std::vector<Polygon> triangulation = problem->get_triangulation();
    std::vector<Point> points = problem->get_points();

    for(Polygon t : triangulation){
        if(is_obtuse_triangle(t)){
            Point a = t.vertex(0);
            Point b = t.vertex(1);
            Point c = t.vertex(2);

            if(std::find(points.begin(), points.end(), a) == points.end()){
                locally_optimize_position(a, triangulation, problem);
            }
            if(std::find(points.begin(), points.end(), b) == points.end()){
                locally_optimize_position(b, triangulation, problem);
            }
            if(std::find(points.begin(), points.end(), c) == points.end()){
                locally_optimize_position(c, triangulation, problem);
            }
        }
    }
    problem->set_triangulation(triangulation);
}

void locally_optimize_constraint(Problem *problem){
    std::vector<Polygon> triangulation = problem->get_triangulation();
    std::vector<Point> points = problem->get_points();

    for(Polygon t : triangulation){
        if(is_obtuse_triangle(t)){
            Point a = t.vertex(0);
            Point b = t.vertex(1);
            Point c = t.vertex(2);

            if(std::find(points.begin(), points.end(), a) == points.end()){
                locally_optimize_position_constraint(a, triangulation, problem);
            }
            if(std::find(points.begin(), points.end(), b) == points.end()){
                locally_optimize_position_constraint(b, triangulation, problem);
            }
            if(std::find(points.begin(), points.end(), c) == points.end()){
                locally_optimize_position_constraint(c, triangulation, problem);
            }
        }
    }
    problem->set_triangulation(triangulation);
}

void globally_optimize_solution(Problem *problem){
    std::vector<Point> steiner = problem->get_steiner();
    std::vector<Polygon> triangulation = problem->get_triangulation();

    std::vector<Point> new_steiner = globally_optimize_position(steiner, triangulation, problem, false);

    problem->clear_solution();
    for(int i = 0; i < new_steiner.size(); i++){
        problem->add_steiner(new_steiner[i]);
    }
    problem->set_triangulation(triangulation);
}

void globally_optimize_obtuse(Problem *problem){
    std::vector<Polygon> triangulation = problem->get_triangulation();
    std::vector<Point> points = problem->get_points();

    std::set<Point> obtuse_steiner_set;

    for(Polygon t : triangulation){
        if(is_obtuse_triangle(t)){
            Point a = t.vertex(0);
            Point b = t.vertex(1);
            Point c = t.vertex(2);

            if(std::find(points.begin(), points.end(), a) == points.end()){
                obtuse_steiner_set.insert(a);
            }
            if(std::find(points.begin(), points.end(), b) == points.end()){
                obtuse_steiner_set.insert(b);
            }
            if(std::find(points.begin(), points.end(), c) == points.end()){
                obtuse_steiner_set.insert(c);
            }
        }
    }

    std::vector<Point> obtuse_steiner_points(obtuse_steiner_set.begin(), obtuse_steiner_set.end());

    //std::cout << "Num obtuse steiner: " << obtuse_steiner_points.size() << std::endl;

    std::vector<Point> new_steiner = globally_optimize_position(obtuse_steiner_points, triangulation, problem, false);
    
    problem->clear_solution();
    for(int i = 0; i < new_steiner.size(); i++){
        problem->add_steiner(new_steiner[i]);
    }
    problem->set_triangulation(triangulation);
}

void fix_constraints(Problem *problem){
    std::vector<Polygon> triangulation = problem->get_triangulation();

    for(Polygon t : triangulation){
        int angle = find_obtuse_angle(t);
        if(angle != -1){
            Point s = t.vertex(angle);
            Point a = t.vertex((angle + 1)%3);
            Point b = t.vertex((angle + 2)%3);
        }
    }
}


void fix_boundary(Problem *problem){
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

            //std::cout << "ttt" << std::endl;

            // Case 1 (for boundary)
            if(boundary.has_on_boundary(a) && boundary.has_on_boundary(b)){
            //if(is_boundary_vertex_with_tolerance(a, boundary) && is_boundary_vertex_with_tolerance(b, boundary)){
                Segment s1 = find_boundary_segment(a, boundary);
                Segment s2 = find_boundary_segment(b, boundary);
                if(s1 != s2){
                    //std::cout << "Special case while fixing boundary: segments do not match!" << std::endl;
                } else {
                    Point p = project_onto_segment(s, s1);

                    problem->remove_triangle(t);

                    Polygon t1;
                    t1.push_back(a);
                    t1.push_back(s);
                    t1.push_back(p);
                    problem->add_triangle(t1);

                    Polygon t2;
                    t2.push_back(b);
                    t2.push_back(s);
                    t2.push_back(p);
                    problem->add_triangle(t2);

                    problem->add_steiner(p);
                }
            } 
            // Case 2
            if (boundary.has_on_boundary(s)){
            //if(is_boundary_vertex_with_tolerance(s, boundary)){
                //std::cout << "222" << std::endl;
                Segment segment = find_boundary_segment(s, boundary);
                double distance_a = std::sqrt(CGAL::to_double(CGAL::squared_distance(segment, a)));
                double distance_b = std::sqrt(CGAL::to_double(CGAL::squared_distance(segment, b)));
                if(std::find(points.begin(), points.end(), s) == points.end()){
                    if(!segment.has_on(a)){
                    //if(distance_a > 1e-6){
                        Point p = project_onto_segment(a, segment);
                        problem->update_triangulation(s, p);
                    } else if (!segment.has_on(b)){
                    //} else if (distance_b > 1e-6){
                        Point p = project_onto_segment(b, segment);
                        problem->update_triangulation(s, p);
                    }
                } else {
                    if(!segment.has_on(a) && std::find(points.begin(), points.end(), a) == points.end() && !boundary.has_on_boundary(a) && !is_on_constraint(a, problem, &constraint)){
                    //if(distance_a > 1e-6 && std::find(points.begin(), points.end(), a) == points.end() && !is_boundary_vertex_with_tolerance(a, boundary) && !is_on_constraint(a, problem, &constraint)){
                        Line line_parallel(a, b - s);
                        Point p = line_parallel.projection(s);
                        problem->update_triangulation(a, p);

                    } else if (!segment.has_on(b) && std::find(points.begin(), points.end(), b) == points.end() && !boundary.has_on_boundary(b) && !is_on_constraint(b, problem, &constraint)){
                    //} else if (distance_b > 1e-6 && std::find(points.begin(), points.end(), b) == points.end() && !is_boundary_vertex_with_tolerance(b, boundary) && !is_on_constraint(b, problem, &constraint)){
                        Line line_parallel(b, a - s);
                        Point p = line_parallel.projection(s);
                        problem->update_triangulation(b, p);
                    }
                }
            }
        }
    }
}

double get_function_value(Problem* problem){
    std::vector<Polygon> triangulation = problem->get_triangulation();
    double cost = 0.0;

    for(Polygon triangle : triangulation){
        Point a = triangle.vertex(0);
        Point b = triangle.vertex(1);
        Point c = triangle.vertex(2);

        /*double ab2 = CGAL::to_double(CGAL::squared_distance(a, b));
        double bc2 = CGAL::to_double(CGAL::squared_distance(b, c));
        double ca2 = CGAL::to_double(CGAL::squared_distance(c, a));

        double w1 = (ab2 + bc2 - ca2)/(2*std::sqrt(ab2)*std::sqrt(bc2)); // cos alpha
        double w2 = (ab2 + ca2 - bc2)/(2*std::sqrt(ab2)*std::sqrt(ca2)); // cos beta
        double w3 = (bc2 + ca2 - ab2)/(2*std::sqrt(bc2)*std::sqrt(ca2)); // cos gamma*/

        Vector ab = (b - a)/((b - a).squared_length());
        Vector bc = (c - b)/((c - b).squared_length());
        Vector ca = (a - c)/((a - c).squared_length());

        auto w1 = CGAL::to_double((-ca)*(ab));
        auto w2 = CGAL::to_double((-ab)*(bc));
        auto w3 = CGAL::to_double((-bc)*(ca));

        double e1 = 1/(1 + std::exp(-1*(acos(w1) - M_PI/2)));
        double e2 = 1/(1 + std::exp(-1*(acos(w2) - M_PI/2)));
        double e3 = 1/(1 + std::exp(-1*(acos(w3) - M_PI/2)));

        cost = cost + e1 + e2 + e3;
    }

    return cost;
}

int main(int argc, char **argv)
{
    Problem *problem = new Problem("../instances_presentation/simple-polygon-exterior-20_60_53ad6d23.instance.json");
    problem->visualize_solution();
    //Problem *problem = new Problem(argv[1]);

    /*for (const auto& entry : fs::directory_iterator("../instances/challenge_instances_cgshop25")) {
        if (entry.is_regular_file() && entry.path().extension() == ".json") {
            Problem *problem = new Problem("../instances/challenge_instances_cgshop25/" + entry.path().filename().string());
            build_quad_mesh_gmsh(problem);
            //std::cout << entry.path().filename().string() << std::endl;
        }
    }*/

    //* Quad Mesh
    //build_quad_mesh_gmsh(problem);
    //build_quad_mesh_medians(problem);
    
    //* Custom Delaunay refinement
    //refine(problem);

    //* Step by step CGAL refinement
    //step_by_step_mesh(problem);
    //refine(problem);
    //classic_delaunay_refinement(problem);


    //* Just testing optimization
    /*Solver* solver = new Mesh();
    solver->solve(problem);
    std::cout << "Num obtuse before: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution();

    globally_optimize_solution(problem);
    std::cout << "Num obtuse after: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution();*/

    /*std::cout << "Function value before: " << get_function_value(problem) << std::endl;
    globally_optimize_solution(problem);
    std::cout << "Num obtuse after first optimization: " << countObtuse(problem->get_triangulation()) << std::endl;
    std::cout << "Function value after: " << get_function_value(problem) << std::endl;
    problem->visualize_solution();*/
  
    /*optimizeTinyAD(problem);
    std::cout << "Num obtuse after optimization 2: " << countObtuse(problem->get_triangulation()) << std::endl;*/


    /*parseOptions(argc, argv);

    std::string path = "../instances_presentation";

    /*for (const auto& entry : fs::directory_iterator(path)) {
        Problem *problem = new Problem(entry.path().string());
        Solver* solver = new Mesh();
        solver->solve(problem);
        problem->visualize_solution();
        problem->output();
    }*/

    /*Problem *problem = new Problem("../instances_presentation/simple-polygon_250_432b4814.instance.json");
    //problem->load_solution();

    Solver* solver = new Mesh();
    solver->solve(problem);

    //std::cout << "triangles size: " << (problem->get_triangulation()).size() << " steiner size: " << (problem->get_steiner()).size() << std::endl;

    //problem->visualize_solution();
    std::cout << "Num obtuse before: " << countObtuse(problem->get_triangulation()) << std::endl;

    /*int num = 0;
    int num_exact = 0;
    for(Point s : problem->get_steiner()){
        if(is_boundary_vertex_with_tolerance(s, problem->get_boundary())){
            num++;
        }
        if(!is_interior_vertex(s, problem->get_boundary())){
            num_exact++;
        }
    }

    std::cout << "num_exact: " << num_exact << " num: " << num << std::endl;*/

    /*locally_optimize_solution(problem);
    //globally_optimize_solution(problem);
    std::cout << "Num obtuse after first optimization: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution();

    /*globally_optimize_obtuse(problem);
    std::cout << "Num obtuse after second optimization: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution();

    globally_optimize_obtuse(problem);
    std::cout << "Num obtuse after second optimization: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution();*/

    /*perform_edge_flips(problem, true);
    std::cout << "Num obtuse after flips 1: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution();

    fix_boundary(problem);
    std::cout << "Num obtuse after fixing boundary 1: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution();

    locally_optimize_obtuse(problem);
    //globally_optimize_obtuse(problem);
    std::cout << "Num obtuse after second optimization: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution();

    fix_boundary(problem);
    std::cout << "Num obtuse after fixing boundary 2: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution();

    perform_edge_flips(problem, false);
    std::cout << "Num obtuse after flips 2: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution();

    locally_optimize_obtuse(problem);
    //globally_optimize_obtuse(problem);
    std::cout << "Num obtuse after third optimization: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution();

    fix_boundary(problem);
    std::cout << "Num obtuse after fixing boundary 3: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution();

    locally_optimize_obtuse(problem);
    //globally_optimize_obtuse(problem);
    std::cout << "Num obtuse after third optimization: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution();

    fix_boundary(problem);
    std::cout << "Num obtuse after fixing boundary 3: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution();*/


    //Problem *problem = new Problem(argv[1]);

    //std::string name = problem->get_name();
    //std::string path = "../solutions/ipe/SOLUTION-" + name + ".ipe";

    /*Solver* solver;

    switch (algorithm) {
        case 0:
            solver = new Stupid();
            break;
        case 1:
            solver = new Mesh();
            break;
        case 2:
            solver = new Manual();
            break;
    }

    if (solver) solver->solve(problem);

    problem->visualize_solution();*/


    return 0;
}
