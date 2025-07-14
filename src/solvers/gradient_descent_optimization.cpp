#include "gradient_descent_optimization.hpp"


// Return the longest edge of a triangle
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

// Returns <true> if the polygon triangle contains edge p1p2
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

// Given a triangle an its edge longest_edge returns the third point
Point find_opposite_point(Polygon& triangle, Segment& longest_edge){
    for (int i = 0; i < 3; i++) {
        if(triangle[i] != longest_edge.source() && triangle[i] != longest_edge.target()){
            return triangle[i];
        }
    }
    throw std::runtime_error("Opposite point not found!");
}

// Returns measure of how obtuse a triangle is
double obtuseness_factor(Polygon& triangle) {
    Point a = triangle.vertex(0);
    Point b = triangle.vertex(1);
    Point c = triangle.vertex(2);

    double ab = CGAL::to_double(CGAL::squared_distance(a, b));
    double bc = CGAL::to_double(CGAL::squared_distance(b, c));
    double ac = CGAL::to_double(CGAL::squared_distance(c, a));

    double max_sq = std::max({ab, bc, ac});

    double obtuseness = max_sq - (ab + bc + ac - max_sq);

    return obtuseness; // If positive, the triangle is obtuse
}

// Flips the longest edge of the triangle is the obtusness score decreases (score not checked if ignore = true)
std::pair<std::vector<Polygon>, std::vector<Polygon>> flip(Polygon& triangle, std::vector<Polygon>& triangulation, bool ignore){
    Segment longest_edge = find_longest_edge(triangle);

    std::vector<Polygon> removed_triangles;
    std::vector<Polygon> new_triangles;

    for (const auto& t : triangulation) {
        if (contains_edge(t, longest_edge.source(), longest_edge.target())) {
            removed_triangles.push_back(t);
        }
    }

    if(removed_triangles.size() == 2){
        Point a = find_opposite_point(removed_triangles[0], longest_edge);
        Point b = find_opposite_point(removed_triangles[1], longest_edge);

        double obtuseness0 = obtuseness_factor(triangle);

        Polygon t1;
        t1.push_back(a);
        t1.push_back(b);
        t1.push_back(longest_edge.source());
        
        Polygon t2;
        t2.push_back(a);
        t2.push_back(b);
        t2.push_back(longest_edge.target());

        double obtuseness1 = std::max(obtuseness_factor(t1), obtuseness_factor(t2));

        if(obtuseness1 < obtuseness0 || ignore){
            new_triangles.push_back(t1);
            new_triangles.push_back(t2);
            return {removed_triangles, new_triangles};
        } else {
            removed_triangles.clear();
            new_triangles.clear();
            return {removed_triangles, new_triangles};
        }
    }

    return {removed_triangles, new_triangles};
}

// Performs edge flips for all obtuse triangles
void perform_edge_flips(Problem *problem, bool ignore){
    std::vector<Polygon> triangulation = problem->get_triangulation();

    std::vector<Polygon> removed_triangles;
    std::vector<Polygon> new_triangles;

    std::vector<Polygon> to_flip = triangulation;
    for(Polygon& triangle : to_flip){
        if(is_obtuse_triangle(triangle)){
            auto [current_remove, current_new] = flip(triangle, triangulation, ignore);

            //add and remove if it has not been done before
            for (const Polygon& tri : current_remove) {
                auto it = std::find(removed_triangles.begin(), removed_triangles.end(), tri);
                if (it == removed_triangles.end()) {
                    removed_triangles.push_back(tri);
                }
            }
            for (const Polygon& tri : current_new) {
                auto it = std::find(new_triangles.begin(), new_triangles.end(), tri);
                if (it == new_triangles.end()) {
                    new_triangles.push_back(tri);
                }
            }
        }
    }

    triangulation.erase(
    std::remove_if(triangulation.begin(), triangulation.end(),
        [&](const Polygon& tri) {
            return std::find(removed_triangles.begin(), removed_triangles.end(), tri) != removed_triangles.end();
        }),triangulation.end());

    for(auto& triangle : new_triangles){
        triangulation.push_back(triangle);
    }

    problem->set_triangulation(triangulation);
}

// Performs local optimization for all Steiner points in the triangulation
void locally_optimize_solution(Problem *problem){
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

// Performs local optimization for vertices of obtuse triangles in the triangulation
void locally_optimize_obtuse(Problem *problem){
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

    problem->set_triangulation(triangulation);
}

// Performs local optimization for all Steiner points in the triangulation (including along constraints)
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

// Performs global optimization for all Steiner points in the triangulation
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

// Performs global optimization for vertices of obtuse triangles in the triangulation
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

    std::vector<Point> new_steiner = globally_optimize_position(obtuse_steiner_points, triangulation, problem, false);
    
    problem->set_triangulation(triangulation);
}

// Addresses boundary-adjacent obtuse triangles
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