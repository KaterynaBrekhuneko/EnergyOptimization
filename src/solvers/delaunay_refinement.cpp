#include "delaunay_refinement.hpp"

typedef CDT::Edge Edge;

void print_current_time_refinement(){
    auto now = std::chrono::system_clock::now();              
    std::time_t now_c = std::chrono::system_clock::to_time_t(now); 

    std::tm* local_time = std::localtime(&now_c);

    std::cout << "Current time: "
              << std::put_time(local_time, "%Y-%m-%d %H:%M:%S")
              << std::endl;
}

Segment get_shortest_edge(const Point& a, const Point& b, const Point& c){
    auto ab = CGAL::squared_distance(a, b);
    auto bc = CGAL::squared_distance(b, c);
    auto ac = CGAL::squared_distance(a, c);

    if(ab <= bc && ab <= ac){
        return Segment(a, b);
    } else if (bc <= ab && bc <= ac) {
        return Segment(b, c);
    } else {
        return Segment(a, c);
    }
}

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

double get_max_edge_length (Problem* problem){
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
        if(problem->triangle_is_inside(cdt, t)){
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
        if(problem->triangle_is_inside(cdt, face)){
            faces.push_back(face);
        }
    }
}

struct EdgeHash {
    std::size_t operator()(const Edge& e) const noexcept {
        std::hash<Face_handle> face_hash;
        std::hash<int> int_hash;
        return face_hash(e.first) ^ int_hash(e.second);
    }
};

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
        if(problem->triangle_is_inside(cdt, face)){
            face_to_ratio_map[face] = get_sizing_ratio_triangle(face, sizing);
        }
    }

    return face_to_ratio_map;
}

void mesh(Problem* problem, CDT& cdt, double sizing, double target_sizing){

    std::unordered_map<Edge, double, EdgeHash, EdgeEqual> edge_to_ratio_map = create_edge_map(cdt, sizing);
    std::unordered_map<Face_handle, double> face_to_ratio_map = create_face_map(problem, cdt, sizing);

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
        if(!problem->triangle_is_inside(cdt, t)){
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

void fix_rhombus(Problem* problem, CDT& cdt){
    auto triangles = problem->get_triangulation();
    auto boundary = problem->get_boundary();
    auto points = problem->get_points();

    std::set<Point> rhombus_centers;

    std::vector<Point> points_to_remove;

    for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
        Vertex_handle v = vit;
        std::vector<Face_handle> incident_faces;
        CDT::Face_circulator fc_start = cdt.incident_faces(v), fc = fc_start;

        // Collect incident faces (triangles)
        if (fc != 0) {
            do {
                if (!cdt.is_infinite(fc)) {
                    incident_faces.push_back(fc);
                }
                ++fc;
            } while (fc != fc_start);
        }

        // Rhombus pattern typically has 4 triangles around a point
        if (incident_faces.size() != 4) continue;

        bool has_obtuse = false;
        Point s;
        for (auto fh : incident_faces){
            Polygon triangle;
            triangle.push_back(fh->vertex(0)->point());
            triangle.push_back(fh->vertex(1)->point());
            triangle.push_back(fh->vertex(2)->point());

            int obtuse = find_obtuse_angle(triangle);
            if(obtuse != -1){
                has_obtuse = true;
                s = fh->vertex(obtuse)->point();
            }
        }
        if(has_obtuse && std::find(points.begin(), points.end(), s) == points.end() && std::find(points_to_remove.begin(), points_to_remove.end(), s) == points_to_remove.end() && !is_on_boundary_refinement(s, boundary)){
            points_to_remove.push_back(s);
        }
    }

    auto steiner = problem->get_steiner();
    steiner.erase(
        std::remove_if(
            steiner.begin(),
            steiner.end(),
            [&points_to_remove](const Point& p) {
                return std::find(points_to_remove.begin(), points_to_remove.end(), p) != points_to_remove.end();
            }
        ),
        steiner.end()
    );
    problem->set_steiner(steiner);
}

Point project_point_onto_boundary(Problem* problem, Point& p){
    auto boundary = problem->get_boundary();
    
    double min_dist = DBL_MAX;
    Segment e;

    for (auto edge = boundary.edges_begin(); edge != boundary.edges_end(); ++edge) {
        double current_dist = CGAL::to_double(CGAL::squared_distance(*edge, p));
        if(current_dist < min_dist){
            min_dist = current_dist;
            e = *edge;
        }
    }

    Point proj = project_onto_segment_refinement(p, e);
    return proj;
}

std::pair<bool, Point> point_on_other_side_of_constraint(Problem* problem, Point& s, Point& new_p){
    auto constraints = problem->get_constraints();

    Segment segment(s, new_p);

    for(Segment c : constraints){
        const auto result = intersection(segment, c);

        if (result) {
            if (const Segment* s_intersection = boost::get<Segment>(&*result)) {
                return {false, Point()};
            } else if (const Point* p_intersection = boost::get<Point>(&*result)){
                if(*p_intersection != s ){
                    return {true, *p_intersection};
                } else {
                    continue;
                }
            }
        }
    }
    return {true, new_p};
}

bool next_refinement_point(Problem* problem, CDT& cdt){
    auto boundary = problem->get_boundary();

    for(Face_handle t : cdt.finite_face_handles()){
        if(!problem->triangle_is_inside(cdt, t)){
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

            if (is_on_constraint_refinement(a, problem) && is_on_constraint_refinement(b, problem) && is_constrained_segment(problem, a, b)){
                Segment segment(a, b);
                Point p = project_onto_segment_refinement(s, segment);
                cdt.insert(p);
                problem->add_steiner(p);
                return false;
            }
        }
    }



    for(Face_handle t : cdt.finite_face_handles()){
        if(!problem->triangle_is_inside(cdt, t)){
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
                cdt.insert(p);
                problem->add_steiner(p);
            } else {
                // Calculate new steiner point
                Point circumsenter = CGAL::circumcenter(t->vertex(0)->point(), t->vertex(1)->point(), t->vertex(2)->point());
                auto side = problem->get_boundary().oriented_side(circumsenter);

                if(side == CGAL::POSITIVE || side == CGAL::ZERO){
                    auto [is_valid, new_point] = point_on_other_side_of_constraint(problem, s, circumsenter);
                    if(is_valid){
                        problem->add_steiner(new_point);
                        cdt.insert(new_point);
                        return false;
                    } else {
                        continue;
                    }
                } else {
                    Point proj = project_point_onto_boundary(problem, circumsenter);
                    problem->add_steiner(proj);
                    cdt.insert(proj);
                    return false;
                }
            }
        }
    }
    return true;
}

//* Batch refinement
Mesh_Statistics refine(Problem* problem){
    Mesh_Statistics stats;

    std::vector<Point> points = problem->get_points();
    std::set<Point> point_set(points.begin(), points.end());
    auto constraints = problem->get_constraints();

    CDT cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);

    /*mesh_cgal(problem);
    cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);*/

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
        problem->update_problem(cdt, point_set);
        //problem->visualize_solution({});

        // LLoyd optimization
        iterate_lloyd(cdt, problem, constraints, 1000);
        problem->update_problem(cdt, point_set);

        // Energy
        /*optimize_TinyAD(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);
        optimize_TinyAD(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);
        optimize_TinyAD(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);*/

        //problem->visualize_solution({});

        //!For this first need to convert to an inexact kernel
        /*CGAL::lloyd_optimize_mesh_2(cdt, CGAL::parameters::number_of_iterations(1000));
        problem->update_problem(cdt, point_set);*/
        //problem->visualize_solution({});
        i++;
    }

    /*iterate_lloyd(cdt, problem, constraints, 1000);
    problem->update_problem(cdt, point_set);*/

    stats.set_name(problem->get_name());
    stats.set_steiner_after_meshing(problem->get_steiner().size());
    stats.set_deviation(mean_absolute_deviation(problem));
    stats.set_aspect_ratio(mean_aspect_ratio(problem));
    save_min_max_angle(problem, &stats);

    stats.set_max_length(get_max_edge_length (problem));

    problem->visualize_solution({});

    std::cout << "num steiner: "  << problem->get_steiner().size() << std::endl;
    std::cout << "deviation: "  << mean_absolute_deviation(problem) << std::endl;
    std::cout << "max_edge_length: "  << get_max_edge_length (problem) << RESET << std::endl;

    save_angle_stats_for_plot(problem, "../results/angle_stats_lloyd_simple-polygon_20_4bd3c2e5");
    save_aspect_ratios_for_plot(problem, "../results/aspect_ratios_lloyd_simple-polygon_20_4bd3c2e5");

    return stats;
}

//* Traditional circumcenter Delaunay
void classic_delaunay_refinement(Problem* problem){
    // Preprocess the input instance, compute cdt
    std::vector<Point> points = problem->get_points();
    std::set<Point> point_set(points.begin(), points.end());
    Polygon boundary = problem->get_boundary();

    auto box = problem->get_boundary().bbox();
    auto scale = std::max(box.xmax() - box.xmin(), box.ymax() - box.ymin());

    CDT cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);
    //problem->visualize_solution({});
    std::cout  << "Num of obtuse in cdt before: " << count_obtuse_triangles(problem) << std::endl;

    int obtuse_after_fix = count_obtuse_triangles(problem);
    int j = 0;
    while(j != 10){
        fix_boundary_refinement(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);

        int current_obtuse = problem->get_num_obtuse();
        if(current_obtuse < obtuse_after_fix){
            obtuse_after_fix = current_obtuse;
        } else {
            break;
        }
        j++;
    }

    j = 0;
    while(j != 10){
        fix_constraint_refinement(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);
        j++;
    }

    //problem->visualize_solution({});

    for(int i = 0; i < 110; i++){
        bool is_refinement_done = next_refinement_point(problem, cdt);
        problem->update_problem(cdt, point_set);

        if(count_obtuse_triangles(problem) == 0){
            break;
        }

        optimize_TinyAD(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);

        j = 0;
        while(j != 10){
            fix_constraint_refinement(problem);
            cdt = problem->generate_CDT();
            problem->update_problem(cdt, point_set);
            j++;
        }

        obtuse_after_fix = count_obtuse_triangles(problem);
        j = 0;
        while(j != 10){
            fix_boundary_refinement(problem);
            cdt = problem->generate_CDT();
            problem->update_problem(cdt, point_set);

            int current_obtuse = problem->get_num_obtuse();
            if(current_obtuse < obtuse_after_fix){
                obtuse_after_fix = current_obtuse;
            } else {
                break;
            }
            j++;
        }

        optimize_TinyAD(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);

        j = 0;
        while(j != 10){
            fix_constraint_refinement(problem);
            cdt = problem->generate_CDT();
            problem->update_problem(cdt, point_set);
            j++;
        }

        obtuse_after_fix = count_obtuse_triangles(problem);
        j = 0;
        while(j != 10){
            fix_boundary_refinement(problem);
            cdt = problem->generate_CDT();
            problem->update_problem(cdt, point_set);

            int current_obtuse = problem->get_num_obtuse();
            if(current_obtuse < obtuse_after_fix){
                obtuse_after_fix = current_obtuse;
            } else {
                break;
            }
            j++;
        }
    }

    std::cout << RED << "num obtuse: " << count_obtuse_triangles(problem) << RESET << std::endl;
    std::cout << RED << "num steiner: " << problem->get_steiner().size() << RESET << std::endl;
    print_current_time_refinement();
    problem->visualize_solution({});

    std::string path = "../normal_delaunay_refinement/svg/" + problem->get_name() + ".svg";
    problem->save_intermidiate_result(path);
}

//****************************************Offcenters*****************************************************

bool is_encroached(CDT& cdt, Point& p){
    for (auto eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); ++eit) {
        Face_handle f = eit->first;
        int i = eit->second;

        Face_handle neighbor = f->neighbor(i);
        // check if this is a boundary edge
        if (cdt.is_infinite(f) || cdt.is_infinite(neighbor)) {
            // check if the edge encroaches p
            Point a = f->vertex(cdt.cw(i))->point();
            Point b = f->vertex(cdt.ccw(i))->point();
            Point midpoint = CGAL::midpoint(a, b);
            FT radius_squared = CGAL::squared_distance(a, b) / 4.0;
            FT dist_squared = CGAL::squared_distance(p, midpoint);
            if(dist_squared <= radius_squared){
                return true;
            }
        }
    }
    return false;
}

double circumradius(const Point& p, const Point& q, const Point& r) {
    Point c = CGAL::circumcenter(p, q, r);
    return std::sqrt(CGAL::to_double(CGAL::squared_distance(p, c)));
}

Point find_offcenter(const Point& P, const Point& Q, const Point& C1) {
    Point mid_PQ = CGAL::midpoint(P, Q);
    Vector dir = mid_PQ - C1;

    double target_ratio = std::sqrt(2.0);
    double pq_length = std::sqrt(CGAL::to_double(CGAL::squared_distance(P, Q)));

    // Check bounds
    auto ratio_at = [&](double t) {
        Point C = C1 + t * dir;
        double R = circumradius(P, Q, C);
        return R / pq_length;
    };

    double t_low = 0.0, t_high = 1.0, eps = 1e-8;
    Point best_point;

    double ratio = ratio_at(0.0);

    while (t_high - t_low > eps) {
        double t = 0.5 * (t_low + t_high);
        ratio = ratio_at(t);

        if (std::abs(ratio - target_ratio) < eps) {
            best_point = C1 + t * dir;
            break;
        }

        if (ratio < target_ratio)
            t_high = t;
        else
            t_low = t;
    }


    return C1 + 0.5 * (t_low + t_high) * dir; // fallback: midpoint of final interval
}

bool next_refinement_point_offcenter(Problem* problem, CDT& cdt){
    auto boundary = problem->get_boundary();
    double beta = 2;

    std::vector<Point> encroached_points;

    for(Face_handle t : cdt.finite_face_handles()){
        if(!problem->triangle_is_inside(cdt, t)){
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
                cdt.insert(p);
                problem->add_steiner(p);
            } else {
                // Calculate new steiner point which is offsenter
                Point c;

                Segment pq = get_shortest_edge(a, b, s);
                Point p = pq[0];
                Point q = pq[1];
                Point c1 = CGAL::circumcenter(t->vertex(0)->point(), t->vertex(1)->point(), t->vertex(2)->point());

                Point c2 = CGAL::circumcenter(p, q, c1);

                // Check the ratio
                auto ratio = CGAL::to_double(CGAL::squared_distance(c1, c2) / pq.squared_length());
                if((CGAL::squared_distance(c1, c2) / pq.squared_length()) <= beta){
                    c = c1;
                } else {
                    c = find_offcenter(p, q, c1);
                }

                if(is_encroached(cdt, c)){
                    encroached_points.push_back(c);
                    continue;
                }

                auto side = problem->get_boundary().oriented_side(c);

                if(side == CGAL::POSITIVE || side == CGAL::ZERO){
                    auto [is_valid, new_point] = point_on_other_side_of_constraint(problem, s, c);
                    if(is_valid){
                        problem->add_steiner(new_point);
                        cdt.insert(new_point);
                        return false;
                    } else {
                        continue;
                    }
                } else {
                    Point proj = project_point_onto_boundary(problem, c);
                    problem->add_steiner(proj);
                    cdt.insert(proj);
                    return false;
                }
            }
        }
    }

    for(Point c : encroached_points){
        auto side = problem->get_boundary().oriented_side(c);

        if(side == CGAL::POSITIVE || side == CGAL::ZERO){
            problem->add_steiner(c);
            cdt.insert(c);
            return false;
        } else {
            Point proj = project_point_onto_boundary(problem, c);
            problem->add_steiner(proj);
            cdt.insert(proj);
            return false;
        }
    }


    return true;
}

//* Off-center Delaunay
Mesh_Statistics offcenter_delaunay_refinement(Problem* problem){
    std::vector<Point> points = problem->get_points();
    std::set<Point> point_set(points.begin(), points.end());
    Polygon boundary = problem->get_boundary();

    auto box = problem->get_boundary().bbox();
    auto scale = std::max(box.xmax() - box.xmin(), box.ymax() - box.ymin());

    CDT cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);
    problem->visualize_solution({});
    std::cout  << "Num of obtuse in cdt before: " << count_obtuse_triangles(problem) << std::endl;

    int obtuse_after_fix = count_obtuse_triangles(problem);
    int j = 0;
    while(j != 10){
        fix_boundary_refinement(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);

        int current_obtuse = problem->get_num_obtuse();
        if(current_obtuse < obtuse_after_fix){
            obtuse_after_fix = current_obtuse;
        } else {
            break;
        }
        j++;
    }

    j = 0;
    while(j != 10){
        fix_constraint_refinement(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);
        j++;
    }

    //problem->visualize_solution({});

    for(int i = 0; i < 400; i++){
        bool is_refinement_done = next_refinement_point_offcenter(problem, cdt);
        problem->update_problem(cdt, point_set);

        if(count_obtuse_triangles(problem) == 0 /*|| is_refinement_done*/){
            break;
        }

        optimize_TinyAD(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);

        j = 0;
        while(j != 10){
            fix_constraint_refinement(problem);
            cdt = problem->generate_CDT();
            problem->update_problem(cdt, point_set);
            j++;
        }

        obtuse_after_fix = count_obtuse_triangles(problem);
        j = 0;
        while(j != 10){
            fix_boundary_refinement(problem);
            cdt = problem->generate_CDT();
            problem->update_problem(cdt, point_set);

            int current_obtuse = problem->get_num_obtuse();
            if(current_obtuse < obtuse_after_fix){
                obtuse_after_fix = current_obtuse;
            } else {
                break;
            }
            j++;
        }

        optimize_TinyAD(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);

        j = 0;
        while(j != 10){
            fix_constraint_refinement(problem);
            cdt = problem->generate_CDT();
            problem->update_problem(cdt, point_set);
            j++;
        }

        obtuse_after_fix = count_obtuse_triangles(problem);
        j = 0;
        while(j != 10){
            fix_boundary_refinement(problem);
            cdt = problem->generate_CDT();
            problem->update_problem(cdt, point_set);

            int current_obtuse = problem->get_num_obtuse();
            if(current_obtuse < obtuse_after_fix){
                obtuse_after_fix = current_obtuse;
            } else {
                break;
            }
            j++;
        }
    }

    int final_obtuse = count_obtuse_triangles(problem);
    int final_num_steiner = problem->get_steiner().size();

    std::cout << RED << "num obtuse: " << final_obtuse << RESET << std::endl;
    std::cout << RED << "num steiner: " << final_num_steiner << RESET << std::endl;

    Mesh_Statistics stats;
    stats.set_name(problem->get_name());
    stats.set_obtuse_after_meshing(final_obtuse);
    stats.set_steiner_after_meshing(final_num_steiner);

    print_current_time_refinement();
    problem->visualize_solution({});

    std::string path = "../offcenter_delaunay_refinement/svg/" + problem->get_name() + ".svg";
    problem->save_intermidiate_result(path);

    return stats;
}

//*******************************************************************************************************

//***************************CGAL step-for-step**********************************************************

Point project_onto_segment_refinement(const Point p, const Segment& s) {
    Line line(s.source(), s.target()); 
    Point proj = line.projection(p); 

    // Clamp projection to segment endpoints
    if (CGAL::squared_distance(proj, s.source()) + CGAL::squared_distance(proj, s.target()) 
        > CGAL::squared_distance(s.source(), s.target())) {
        throw std::runtime_error("Projected point is not on segment!");
    }
    return proj;
}

bool is_on_constraint_refinement(const Point& steiner, Problem* problem){
    std::vector<Segment> constraints = problem->get_constraints();
    double tolerance = 1e-12;
    for(Segment c : constraints){
        double distance = CGAL::to_double(CGAL::squared_distance(c, steiner));
        if(c.has_on(steiner) || distance <= tolerance){
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
        // Check like this because a and b might not be input points and thus are not end points of the boundary polygon which only containsinput points and is not updated i. e. split into subsegments after steiner point insertion
        double distance_a = CGAL::to_double(CGAL::squared_distance(*edge, a));
        double distance_b = CGAL::to_double(CGAL::squared_distance(*edge, b));
        if ((distance_a < tolerance || edge->has_on(a)) && (distance_b < tolerance || edge->has_on(b))) {
            return true;
        }
    }
    return false;
}

bool is_constrained_segment(Problem* problem, Point& a, Point& b){
    std::vector<Segment> constraints = problem->get_constraints();
    double tolerance = 1e-12;
    for(Segment c : constraints){
        // Check like this because a and b might not be input points and thus are not end points of a constraint
        // which only contain input points and are not updated i. e. split into subsegments after steiner point insertion
        double distance_a = CGAL::to_double(CGAL::squared_distance(c, a));
        double distance_b = CGAL::to_double(CGAL::squared_distance(c, b));
        if ((distance_a < tolerance || c.has_on(a)) && (distance_b < tolerance || c.has_on(b))) {
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
       if(p0 != a && p1 != a && p2!= a){
            if(p0 == b && is_boundary_segment(boundary, p0, p1)){
                return {p1, true};
            } else if (p0 == b && is_boundary_segment(boundary, p0, p2)){
                return {p2, true};
            } else if(p1 == b && is_boundary_segment(boundary, p1, p0)){
                return {p0, true};
            } else if (p1 == b && is_boundary_segment(boundary, p1, p2)){
                return {p2, true};
            } else if(p2 == b && is_boundary_segment(boundary, p2, p0)){
                return {p0, true};
            } else if (p2 == b && is_boundary_segment(boundary, p2, p1)){
                return {p1, true};
            }
       }
    }
    return {Point(0, 0), false};
}

// First two points are on a constraint
std::pair<Point, bool> find_constraint_point(Problem* problem, Point& a, Point& b, Point& c){
    auto triangles = problem->get_triangulation();
    auto constraints = problem->get_constraints();
    for(Polygon t : triangles){
        Point p0 = t[0];
        Point p1 = t[1];
        Point p2 = t[2];
       if(p0 != a && p1 != a && p2!= a){
            if(p0 == b && is_constrained_segment(problem, p0, p1)){
                return {p1, true};
            } else if (p0 == b && is_constrained_segment(problem, p0, p2)){
                return {p2, true};
            } else if(p1 == b && is_constrained_segment(problem, p1, p0)){
                return {p0, true};
            } else if (p1 == b && is_constrained_segment(problem, p1, p2)){
                return {p2, true};
            } else if(p2 == b && is_constrained_segment(problem, p2, p0)){
                return {p0, true};
            } else if (p2 == b && is_constrained_segment(problem, p2, p1)){
                return {p1, true};
            }
       }
    }
    return {Point(0, 0), false};
}

bool point_within_segment(Point& a, Point& b, Point& p){
    double tolerance = 1e-10;
    Segment segment1(a, b);
    Segment segment2(b, a);
    double dist1 = CGAL::to_double(CGAL::squared_distance(segment1, p));
    double dist2 = CGAL::to_double(CGAL::squared_distance(segment2, p));
    return (dist1 < tolerance || segment1.has_on(p) || dist2 < tolerance || segment2.has_on(p)) && (a != p && b != p);
}

std::tuple<Point, Polygon, bool> find_opposite_point(Problem* problem, Point a, Point b, Point s){
    auto triangles = problem->get_triangulation();
    for(const auto t : triangles){
        Point p0 = t[0];
        Point p1 = t[1];
        Point p2 = t[2];
        if((p0 == a && p1 == b && p2 != s) || (p1 == a && p0 == b && p2 != s) ){
            return {p2, t, true};
        } else if((p0 == a && p2 == b && p1 != s) || (p2 == a && p0 == b && p1 != s)){
            return {p1, t, true};
        } else if((p1 == a && p2 == b && p0 != s) || (p2 == a && p1 == b && p0 != s)){
            return {p0, t, true};
        }
    }
    return {Point(0,0), Polygon(), false};
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
            if (std::find(points.begin(), points.end(), s) != points.end() && is_on_boundary_refinement(a, boundary) && is_on_boundary_refinement(b, boundary) && !is_boundary_segment(boundary, a, b)){

                Point c1, c2;
                bool is_valid1, is_valid2 = false;
                
                auto [result1, res_valid1] = find_boundary_point(problem, s, a, b);
                if(res_valid1){
                    is_valid1 = res_valid1;
                    c1 = result1;
                }
                auto [result2, res_valid2] = find_boundary_point(problem, s, b, a);
                if(res_valid2){
                    is_valid2 = res_valid2;
                    c2 = result2;
                }

                auto [c3, t_neighbor, is_valid3] = find_opposite_point(problem, a, b, s);

                // move a if possible
                if(std::find(points.begin(), points.end(), a) == points.end() && is_valid1 && is_valid3){
                    Line line(s, a); 
                    Point proj = line.projection(c3);
                    // check whether the new point position is within the corresponding segment
                    if(point_within_segment(s, c1, proj) && a != proj){
                        problem->update_steiner(a, proj);
                        int new_obtuse = problem->get_num_obtuse() - 1;
                        problem->set_num_obtuse(new_obtuse);
                    }
                }

                // move b if possible
                if(std::find(points.begin(), points.end(), b) == points.end() && is_valid2 && is_valid3){
                    Line line(s, b); 
                    Point proj = line.projection(c3);
                    // check whether the new point position is within the corresponding segment
                    if(point_within_segment(s, c2, proj) && b != proj){
                        problem->update_steiner(b, proj);
                        int new_obtuse = problem->get_num_obtuse() - 1;
                        problem->set_num_obtuse(new_obtuse);
                    }
                }
            }
        }
        
    }
}

void flip_corner_triangle(Problem* problem){
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

            if (std::find(points.begin(), points.end(), s) != points.end() && (is_on_boundary_refinement(a, boundary) || is_on_constraint_refinement(a, problem)) && (is_on_boundary_refinement(b, boundary) || is_on_constraint_refinement(b, problem)) && !is_boundary_segment(boundary, a, b) && (is_boundary_segment(boundary, a, s) || is_constrained_segment(problem, a, s)) && (is_boundary_segment(boundary, b, s) || is_constrained_segment(problem, b, s))){

                auto [c, t_neighbor, is_valid] = find_opposite_point(problem, a, b, s);

                if(is_valid){
                    problem->remove_triangle(t);
                    problem->remove_triangle(t_neighbor);

                    Polygon t1;
                    t1.push_back(a);
                    t1.push_back(s);
                    t1.push_back(c);
                    problem->add_triangle(t1);

                    Polygon t2;
                    t2.push_back(b);
                    t2.push_back(s);
                    t2.push_back(c);
                    problem->add_triangle(t2);
                }
            }
        }
    }
}

void fix_constraint_refinement(Problem *problem){
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
            if(is_on_constraint_refinement(a, problem) && is_on_constraint_refinement(b, problem) && is_constrained_segment(problem, a, b)){
                Segment segment(a, b);
                try{
                    Point p = project_onto_segment_refinement(s, segment);
                    problem->add_steiner(p);
                } catch(const std::runtime_error& e){
                    // Do not insert any new points, case could not be fixed
                    std::cout << RED << "Projected point is not on segment!"<< RESET << std::endl;
                }
            } 
            // Case 2
            else if (is_on_constraint_refinement(s, problem)){
                bool a_constraint = is_constrained_segment(problem, a, s);
                bool b_constraint = is_constrained_segment(problem, b, s);

                Point c;
                bool is_valid;

                if(a_constraint){
                    auto [result, res_valid] = find_constraint_point(problem, a, s, b);
                    if(res_valid){
                        is_valid = res_valid;
                        c = result;
                    }
                } else if(b_constraint){
                    auto [result, res_valid] = find_constraint_point(problem, b, s, a);
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
                    if(a_constraint && is_valid){
                        Line line(s, a); 
                        Point proj = line.projection(b);
                        // check whether the new point position is within the corresponding segment
                        if(point_within_segment(a, c, proj)){
                            problem->update_steiner(s, proj);
                            first_option = true;
                        }
                    } else if(b_constraint && is_valid){
                        Line line(s, b); 
                        Point proj = line.projection(a); 
                        // check whether the new point position is within the corresponding segment
                        if(point_within_segment(b, c, proj)){
                            problem->update_steiner(s, proj);
                            first_option = true;
                        }
                    }
                } 

                if (!first_option){
                    if(a_constraint){ // try to move b
                        if(std::find(points.begin(), points.end(), b) == points.end() && !is_on_boundary_refinement(b, boundary) && !is_on_constraint_refinement(b, problem)){
                            Line line_parallel(b, a - s);
                            Point proj = line_parallel.projection(s);
                            //problem->add_steiner(proj);
                            problem->update_steiner(b, proj);
                        }
                    } else if (b_constraint){ // try to move a
                        if(std::find(points.begin(), points.end(), a) == points.end() && !is_on_boundary_refinement(a, boundary) && !is_on_constraint_refinement(a, problem)){
                            Line line_parallel(a, b - s);
                            Point proj = line_parallel.projection(s);
                            //problem->add_steiner(proj);
                            problem->update_steiner(a, proj);
                        }
                    }
                }
            } 
            if (std::find(points.begin(), points.end(), s) != points.end() && is_on_constraint_refinement(a, problem) && is_on_constraint_refinement(b, problem) && !is_constrained_segment(problem, a, b)){

                Point c1, c2;
                bool is_valid1, is_valid2 = false;
                
                auto [result1, res_valid1] = find_constraint_point(problem, s, a, b);
                if(res_valid1){
                    is_valid1 = res_valid1;
                    c1 = result1;
                }
                auto [result2, res_valid2] = find_constraint_point(problem, s, b, a);
                if(res_valid2){
                    is_valid2 = res_valid2;
                    c2 = result2;
                }

                auto [c3, t_neighbor, is_valid3] = find_opposite_point(problem, a, b, s);

                // move a if possible
                if(std::find(points.begin(), points.end(), a) == points.end() && is_valid1 && is_valid3){
                    Line line(s, a); 
                    Point proj = line.projection(c3);
                    // check whether the new point position is within the corresponding segment
                    if(point_within_segment(s, c1, proj) && a != proj){
                        problem->update_steiner(a, proj);
                        int new_obtuse = problem->get_num_obtuse() - 1;
                        problem->set_num_obtuse(new_obtuse);
                    }
                }

                // move b if possible
                if(std::find(points.begin(), points.end(), b) == points.end() && is_valid2 && is_valid3){
                    Line line(s, b); 
                    Point proj = line.projection(c3);
                    // check whether the new point position is within the corresponding segment
                    if(point_within_segment(s, c2, proj) && b != proj){
                        problem->update_steiner(b, proj);
                        int new_obtuse = problem->get_num_obtuse() - 1;
                        problem->set_num_obtuse(new_obtuse);
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
    problem->set_triangulation(triangulation);
}

//* Step by step triangulation using the Terminator point insertion
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

    CDT cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);
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

        problem->update_problem(cdt, point_set);
        if(iter%10 == 0){
            // Optimize
            std::cout << "Num obtuse before optimization " << opt_count << ": " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});
            optimize_TinyAD(problem);

            std::cout << "Num obtuse after optimization " << opt_count << ": " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});

            // Update CDT
            cdt = problem->generate_CDT();
            std::cout << "Num obtuse after regenerating CDT after optimization " << opt_count << ": " << count_obtuse_triangles(problem) << std::endl;
            problem->update_problem(cdt, point_set);
            //problem->visualize_solution({});

            // Fix boundary
            fix_boundary_refinement(problem);
            cdt = problem->generate_CDT();
            problem->update_problem(cdt, point_set);
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

//* Uniform mesh with Terminator + optimization after construction
Mesh_Statistics uniform_mesh(Problem* problem){
    std::vector<Point> points = problem->get_points();
    std::set<Point> point_set(points.begin(), points.end());
    Polygon boundary = problem->get_boundary();

    CDT cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);

    auto bbox = boundary.bbox();
    double max_size = std::sqrt(pow(bbox.x_span(),2) + pow(bbox.y_span(),2)) * 0.1;

    problem->generate_uniform_mesh(cdt, max_size);
    problem->visualize_solution({});

    int obtuse_meshing = count_obtuse_triangles(problem);
    problem->set_num_obtuse(obtuse_meshing);
    int final_num_obtuse = obtuse_meshing;

    Mesh_Statistics statistics;
    statistics.set_name(problem->get_name());
    statistics.set_obtuse_after_meshing(obtuse_meshing);
    statistics.set_steiner_after_meshing(problem->get_steiner().size());

    optimize_TinyAD(problem);
    cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);
    final_num_obtuse = std::min(final_num_obtuse, count_obtuse_triangles(problem));
    //problem->visualize_solution({});

    fix_rhombus(problem, cdt);
    cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);
    final_num_obtuse = std::min(final_num_obtuse, count_obtuse_triangles(problem));
    //problem->visualize_solution({});

    optimize_TinyAD(problem);
    cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);
    final_num_obtuse = std::min(final_num_obtuse, count_obtuse_triangles(problem));
    //problem->visualize_solution({});

    int obtuse_opt = final_num_obtuse;
    problem->set_num_obtuse(obtuse_opt);

    int obtuse_after_fix = obtuse_opt;
    int i = 0;
    while(i != 100){
        fix_boundary_refinement(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);

        int current_obtuse = problem->get_num_obtuse();
        if(current_obtuse < obtuse_after_fix){
            obtuse_after_fix = current_obtuse;
        } else {
            break;
        }
        i++;
    }

    final_num_obtuse = std::min(final_num_obtuse, count_obtuse_triangles(problem));

    i = 0;
    while(i != 10){
        fix_constraint_refinement(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);
        i++;
    }

    final_num_obtuse = std::min(final_num_obtuse, count_obtuse_triangles(problem));

    optimize_TinyAD(problem);
    cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);
    final_num_obtuse = std::min(final_num_obtuse, count_obtuse_triangles(problem));
    //problem->visualize_solution({});

    fix_rhombus(problem, cdt);
    cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);
    final_num_obtuse = std::min(final_num_obtuse, count_obtuse_triangles(problem));
    //problem->visualize_solution({});

    optimize_TinyAD(problem);
    cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);
    final_num_obtuse = std::min(final_num_obtuse, count_obtuse_triangles(problem));
    //problem->visualize_solution({});

    obtuse_after_fix = obtuse_opt;
    i = 0;
    while(i != 100){
        fix_boundary_refinement(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);

        int current_obtuse = problem->get_num_obtuse();
        if(current_obtuse < obtuse_after_fix){
            obtuse_after_fix = current_obtuse;
        } else {
            break;
        }
        i++;
    }
    final_num_obtuse = std::min(final_num_obtuse, count_obtuse_triangles(problem));

    i = 0;
    while(i != 10){
        fix_constraint_refinement(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);
        i++;
    }

    final_num_obtuse = std::min(final_num_obtuse, count_obtuse_triangles(problem));
    //problem->visualize_solution({});

    flip_corner_triangle(problem);
    final_num_obtuse = std::min(final_num_obtuse, count_obtuse_triangles(problem));
    //problem->visualize_solution({});

    obtuse_after_fix = obtuse_opt;
    i = 0;
    while(i != 100){
        fix_boundary_refinement(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);

        int current_obtuse = problem->get_num_obtuse();
        if(current_obtuse < obtuse_after_fix){
            obtuse_after_fix = current_obtuse;
        } else {
            break;
        }
        i++;
    }
    final_num_obtuse = std::min(final_num_obtuse, count_obtuse_triangles(problem));

    i = 0;
    while(i != 10){
        fix_constraint_refinement(problem);
        cdt = problem->generate_CDT();
        problem->update_problem(cdt, point_set);
        i++;
    }

    final_num_obtuse = std::min(final_num_obtuse, count_obtuse_triangles(problem));
    problem->visualize_solution({});

    problem->set_num_obtuse(final_num_obtuse);

    statistics.set_obtuse_after_optimization(final_num_obtuse);
    statistics.set_steiner_after_optimization(problem->get_steiner().size());

    std::string path2 = "../results/OPT-" + problem->get_name() + ".svg";
    problem->save_intermidiate_result(path2);

    std::cout << RED << "num obtuse: " << count_obtuse_triangles(problem) << RESET << std::endl;

    return statistics;
}

//******************************************************************************************************************

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

    CDT cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);

    auto bbox = boundary.bbox();
    double max_size = std::sqrt(pow(bbox.x_span(),2) + pow(bbox.y_span(),2)) * 0.1;

    Mesher mesher(cdt);
    mesher.set_criteria(Criteria(0.125, max_size/4));
    mesher.refine_mesh();

    problem->update_problem(cdt, point_set);
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

double get_aspect_ratio(const Point& a, const Point& b, const Point& c){
    if(get_triangle_area(a, b, c) < 1e-6){
        return INFINITY;
    }
    double shortest_edge = std::sqrt(CGAL::to_double(get_shortest_edge(a, b, c).squared_length()));
    if(shortest_edge < 1e-6){
        return INFINITY;
    }
    double radius = std::sqrt(CGAL::to_double(CGAL::squared_radius(a, b, c)));  
    return radius/shortest_edge;
}

void save_aspect_ratios_for_plot(Problem* problem, std::string path){
    auto triangulation = problem->get_triangulation();

    std::vector<double> ratios;

    double n1 = 0;
    double n2 = 0;

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
        if(ratio <= 0.7){
            n1 = n1+1;
        } else {
            n2 = n2+1;
        }
    }
    out.close();

    std::cout << "Percentage: " << n1/(n1+n2) << std::endl;
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
        if (triangle.size() != 3) continue;

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
    
    optimize_TinyAD(problem);
    CDT cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);
    optimize_TinyAD(problem);
    cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);

    std::cout << "mean absolute deviation: " << mean_absolute_deviation(problem) << std::endl;

    //save_angle_stats_for_plot(problem, "../results/angle_data_lloyd.csv");
    save_aspect_ratios_for_plot(problem, "../results/aspect_ratios_sq_penalty_20.csv");
    problem->visualize_solution({});
}

Mesh_Statistics mesh_equilateral(Problem* problem){
    Mesh_Statistics stats;

    std::vector<Point> points = problem->get_points();
    auto constraints = problem->get_constraints();
    std::set<Point> point_set(points.begin(), points.end());

    mesh_cgal(problem);
    CDT cdt = problem->generate_CDT();
    //refine(problem, cdt);
    //problem->visualize_solution({});

    optimize_TinyAD(problem);
    cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);
    optimize_TinyAD(problem);
    cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);
    /*optimize_TinyAD(problem);
    cdt = problem->generate_CDT();
    problem->update_problem(cdt, point_set);*/

    /*cdt = problem->generate_CDT();
    iterate_lloyd(cdt, problem, constraints, 1000);
    problem->update_problem(cdt, point_set);*/

    stats.set_name(problem->get_name());
    stats.set_steiner_after_meshing(problem->get_steiner().size());
    stats.set_deviation(mean_absolute_deviation(problem));
    std::vector<Polygon> problematic_triangles = save_min_max_angle(problem, &stats);

    //problem->visualize_solution({});

    //std::cout << "deviation: " << stats.get_deviation() << " min angle: " << stats.get_min_angle() << " max angle: " << stats.get_max_angle() << std::endl;

    return stats;
}