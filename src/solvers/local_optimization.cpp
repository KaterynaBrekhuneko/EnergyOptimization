#include "local_optimization.hpp"

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>

typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CDT::Edge_iterator Edge_iterator;

bool armijo_condition(const std::vector<Point>& gradient, double f, double f_new, double step_size, double armijo_const){
    double dot_product = 0.0;
    for(int i = 0; i < gradient.size(); i++){
        dot_product = dot_product - CGAL::to_double(gradient[i].x()*gradient[i].x() - gradient[i].y()*gradient[i].y());
    }
    return f_new <= f + armijo_const*step_size*dot_product;
}

std::vector<Point> line_search(std::vector<Polygon>& triangles, const std::vector<Point>& steiner, const std::vector<Point>& gradient, double s_max, double shrink, double max_iters, double armijo_const){
    // Check input
   assert(steiner.size() == gradient.size());

    if (s_max <= 0.0)
        std::runtime_error("Max step size not positive.");

    // Also try a step size of 1.0 (if valid)
    const bool try_one = s_max > 1.0;

    // Calculate initial function value
    std::vector<std::vector<Polygon>> neighborhoods;
    double f = 0.0;
    int i = 0;
    for(const Point& s : steiner){
        std::vector<Polygon> neighborhood = find_neighborhood(s, triangles, false);
        neighborhoods.push_back(neighborhood);
        f += value_refined_sigmoid(s, s, neighborhood);
        i++;
    }

    std::vector<Point> new_steiner = steiner;
    double step_size = s_max;
    for (int i = 0; i < max_iters; i++){
        for(int j = 0; j<steiner.size(); j++){
            Point current = new_steiner[j];
            Point current_grad = gradient[j];
            //std::cout << "grad at " << i << ": " << gradient[j] << std::endl;
            new_steiner[j] = Point(current.x() - step_size*current_grad.x(), current.y() - step_size*current_grad.y());
        }

        // Calculate new function value
        double f_new = 0.0;
        for(int k = 0; k < steiner.size(); k++){
            Point old_s = steiner[k];
            Point new_s = new_steiner[k];
            f_new += value_refined_sigmoid(old_s, new_s, neighborhoods[k]);
        }

        if (armijo_condition(gradient, f, f_new, step_size, armijo_const)){
            //std::cout << "Armijo did find improvement!" << std::endl;
            return new_steiner;
        }

        if (try_one && step_size > 1.0 && step_size * shrink < 1.0)
            step_size = 1.0;
        else
            step_size *= shrink;
    }

    std::cout << "Armijo did not find improvement!" << std::endl;

    return steiner;
}

Point project_onto_segment(const Point p, const Segment& s) {
    K::Line_2 line(s.source(), s.target()); 
    Point proj = line.projection(p); 

    // ! add sqrt?
    // Clamp projection to segment endpoints
    if (CGAL::squared_distance(proj, s.source()) + CGAL::squared_distance(proj, s.target()) 
        > CGAL::squared_distance(s.source(), s.target())) {
        std::cout << "Problems projecting a point!" << std::endl;
        return (CGAL::squared_distance(proj, s.source()) < CGAL::squared_distance(proj, s.target())) 
               ? s.source() 
               : s.target();
    }
    return proj;
}

Point move_point_on_segment(Point& p, const Segment& s, double step_size, std::vector<Polygon> neighborhood) {
    double gradient_norm_start = norm(calculate_gradient_refined_sigmoid(p, neighborhood));

    // find boundary points on the segment
    std::vector<Point> ab;
    for(const Polygon& triangle : neighborhood){
        for(Point point : triangle.vertices()){
            if(s.has_on(point) && point!=p){
                ab.push_back(point);
            }
        }
    }

    if(ab.size() != 2){
        std::cout << "ab size: " << ab.size() << std::endl;
        return p;
    }
    Point a = ab[0];
    Point b = ab[1];

    Vector segment_vector = a - b;  // Direction of segment
    Vector unit_vector = segment_vector / std::sqrt(CGAL::to_double(segment_vector.squared_length()));  // Normalize

    // Try 1. dir

    // Move point
    Point p1 = p + step_size * unit_vector; 
    p1 = project_onto_segment(p1, s);

    // ! Add sqrt?
    // Ensure the new point stays on the segment
    if (CGAL::squared_distance(p1, a) + CGAL::squared_distance(p1, b) 
        > CGAL::squared_distance(a, b)) {
        std::cout << "Point not on Segment" << std::endl;
        step_size = step_size/10;
        return p;
    }

    update_neighborhood(neighborhood, p, p1);
    double gradient_norm_opt = norm(calculate_gradient_refined_sigmoid(p1, neighborhood));

    if(gradient_norm_opt > gradient_norm_start){
        return p1;
    }

    // Try 2. dir
    Point p2 = p - step_size * unit_vector;
    p2 = project_onto_segment(p2, s);
    // ! Add sqrt?
    if (CGAL::squared_distance(p2, a) + CGAL::squared_distance(p2, b) 
        > CGAL::squared_distance(a, b)) {
        //std::cout << "Point not on Segment" << std::endl;
        step_size = step_size/10;
        return p;
    }

    update_neighborhood(neighborhood, p1, p2);
    gradient_norm_opt = norm(calculate_gradient_refined_sigmoid(p2, neighborhood));

    if(gradient_norm_opt > gradient_norm_start){
        return p2;
    }

    return p;
}

Point locally_optimize_position_constraint(Point steiner, std::vector<Polygon>& triangles, Problem *problem){
    Point s = steiner;
    Segment constraint;

    double step_size = 100;
    double TOL = 1e-7;
    const int MAX_ITER = 100;

    int iteration = 0;

    if(!is_interior_vertex(s, problem->get_boundary())){
        constraint = find_boundary_segment(s, problem->get_boundary());

    } else if (is_on_constraint(s, problem, &constraint)){
        
    } else {
        return steiner;
    }

    try{
        std::vector<Polygon> neighborhood = find_neighborhood(s, triangles, false);
        Point gradient = calculate_gradient_refined_sigmoid(s, neighborhood);
        
        while(norm(gradient) > TOL && iteration < MAX_ITER){
            Point tmp = move_point_on_segment(s, constraint, step_size, neighborhood);
            update_neighborhood(neighborhood, s, tmp);
            s = tmp;

            gradient = calculate_gradient_refined_sigmoid(s, neighborhood);
            iteration++;
        }

        update_neighborhood(triangles, steiner, s);

    }catch(const CGAL::Assertion_exception& e){
        std::cout << "CGAL Exception 2 \n" << std::endl;
    }

    return s;
}

//--------------------------------From this point code looks good------------------------------

Segment find_boundary_segment(const Point& p, const Polygon& polygon) {
    size_t n = polygon.size();
    
    for (size_t i = 0; i < n; ++i) {
        Segment edge(polygon.vertex(i), polygon.vertex((i + 1) % n)); 
        double distance = std::sqrt(CGAL::to_double(CGAL::squared_distance(edge, p)));
        if (edge.has_on(p) || distance < 1e-6) {
            return edge; 
        }
    }
    
    throw std::runtime_error("Point is not on any segment of the polygon (should not happen).");
}

std::vector<Polygon> find_neighborhood(const Point& steiner, std::vector<Polygon>& triangles, bool debug){    
    std::vector<Polygon> neighborhood;
    int count = 0;
    for (const Polygon& triangle : triangles){
        if(find_point_index(steiner, triangle, (debug && count == 9)) != -1){
            neighborhood.push_back(triangle);
        }
        count++;
    }
    return neighborhood;
}

bool is_interior_vertex(const Point& steiner, const Polygon& boundary){
    return !boundary.has_on_boundary(steiner);
}

bool is_interior_vertex_with_tolerance(const Point& p, const Polygon& polygon) {
    double tolerance = 1e-12;
    for (auto edge_it = polygon.edges_begin(); edge_it != polygon.edges_end(); ++edge_it) {
        double distance = std::sqrt(CGAL::to_double(CGAL::squared_distance(*edge_it, p)));
        //std::cout << "distance: " << distance << std::endl;
        if (distance <= tolerance) {
            return false;
        }
    }
    return true;
}

bool is_on_constraint(const Point& steiner, Problem* problem, Segment* constraint){
    std::vector<Segment> constraints = problem->get_constraints();
    for(Segment c : constraints){
        if(c.has_on(steiner)){
            constraint = &c;
            return true;
        }
    }

    return false;
}

bool is_on_constraint_with_tolerance(const Point& steiner, Problem* problem, Segment* constraint){
    double tolerance = 1e-12;
    std::vector<Segment> constraints = problem->get_constraints();
    for(Segment c : constraints){
        double distance = std::sqrt(CGAL::to_double(CGAL::squared_distance(c, steiner)));
        if(c.has_on(steiner) || distance <= tolerance){
            constraint = &c;
            return true;
        }
    }

    return false;
}

bool is_in_the_neighborhood(const Point& p, const std::vector<Polygon>& neighborhood){
    CGAL::Bounded_side side;
    for(Polygon neighbor : neighborhood){
        side = CGAL::bounded_side_2(neighbor.vertices_begin(), neighbor.vertices_end(), p, K());

        if(side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY){
            return true;
        }
    }
    return false;
}

int find_point_index(const Point& s, const Polygon& triangle, bool debug) {
    double tolerance = 1e-12;
    int index = 0;

    Point a = triangle[0];
    Point b = triangle[1];
    Point c = triangle[2];

    if(debug){
        std::cout << "a: " << a << " b: " << b << " c: " << c << " s: " << s <<std::endl;
    }

    if(std::abs(CGAL::to_double(a.x()) - CGAL::to_double(s.x())) < tolerance && std::abs(CGAL::to_double(a.y()) - CGAL::to_double(s.y())) < tolerance){
        return 0;
    }
    if(std::abs(CGAL::to_double(b.x()) - CGAL::to_double(s.x())) < tolerance && std::abs(CGAL::to_double(b.y()) - CGAL::to_double(s.y())) < tolerance){
        return 1;
    }
    if(std::abs(CGAL::to_double(c.x()) - CGAL::to_double(s.x())) < tolerance && std::abs(CGAL::to_double(c.y()) - CGAL::to_double(s.y())) < tolerance){
        return 2;
    }

    return -1; 
}

void update_neighborhood(std::vector<Polygon>& neighborhood, const Point& s, const Point& new_s) {
    double tolerance = 1e-12;
    for (auto& triangle : neighborhood) {
        for (size_t i = 0; i < triangle.size(); ++i) {
            if (std::abs(CGAL::to_double(triangle.vertex(i).x()) - CGAL::to_double(s.x())) < tolerance && std::abs(CGAL::to_double(triangle.vertex(i).y()) - CGAL::to_double(s.y())) < tolerance) {  
                triangle[i] = new_s;  
            }
        }
    }
}

double norm(const Point& gradient) {
    return std::sqrt(CGAL::to_double(gradient.x()*gradient.x() + gradient.y()*gradient.y()));
}

double value_ln(const Point& old_s, const Point& new_s, const std::vector<Polygon> neighborhood){
    double value = 0;

    if(!is_in_the_neighborhood(new_s, neighborhood)){
        return INFINITY;
    }

    for(const Polygon& triangle : neighborhood){
        int index = find_point_index(old_s, triangle, false);

        Point a = triangle.vertex((index + 1) % 3);
        Point b = triangle.vertex((index + 2) % 3);

        double ab2 = CGAL::to_double(CGAL::squared_distance(a, b));
        double sa2 = CGAL::to_double(CGAL::squared_distance(new_s, a));
        double sb2 = CGAL::to_double(CGAL::squared_distance(new_s, b));

        double w1 = std::acos((ab2 + sa2 - sb2)/(2*std::sqrt(ab2)*std::sqrt(sa2))); // cos alpha
        double w2 = std::acos((ab2 + sb2 - sa2)/(2*std::sqrt(ab2)*std::sqrt(sb2))); // cos beta
        double w3 = std::acos((sa2 + sb2 - ab2)/(2*std::sqrt(sa2)*std::sqrt(sb2))); // cos gamma
        
        if( std::isnan(w1) || std::isnan(w2) || std::isnan(w3)){
            continue;
        }

        double ln1 = std::log(1 + std::exp(w1 - M_PI/2));
        double ln2 = std::log(1 + std::exp(w2 - M_PI/2));
        double ln3 = std::log(1 + std::exp(w3 - M_PI/2));

        value = value + ln1 + ln2 + ln3;
    }

    return value;
}

double value_sigmoid(const Point& old_s, const Point& new_s, const std::vector<Polygon> neighborhood){
    double value = 0;

    double k = 15;

    if(!is_in_the_neighborhood(new_s, neighborhood)){
        return INFINITY;
    }

    for(const Polygon& triangle : neighborhood){
        int index = find_point_index(old_s, triangle, false);

        Point a = triangle.vertex((index + 1) % 3);
        Point b = triangle.vertex((index + 2) % 3);

        double ab2 = CGAL::to_double(CGAL::squared_distance(a, b));
        double sa2 = CGAL::to_double(CGAL::squared_distance(new_s, a));
        double sb2 = CGAL::to_double(CGAL::squared_distance(new_s, b));

        double w1 = std::acos((ab2 + sa2 - sb2)/(2*std::sqrt(ab2)*std::sqrt(sa2))); // cos alpha
        double w2 = std::acos((ab2 + sb2 - sa2)/(2*std::sqrt(ab2)*std::sqrt(sb2))); // cos beta
        double w3 = std::acos((sa2 + sb2 - ab2)/(2*std::sqrt(sa2)*std::sqrt(sb2))); // cos gamma

        if( std::isnan(w1) || std::isnan(w2) || std::isnan(w3)){
            continue;
        }

        double sigm1 =  1 / (1 + std::exp(-k * (w1-M_PI/2)));
        double sigm2 =  1 / (1 + std::exp(-k * (w2-M_PI/2)));
        double sigm3 =  1 / (1 + std::exp(-k * (w3-M_PI/2)));

        value = value + sigm1 + sigm2 + sigm3;
    }

    return value;
}

double value_refined_sigmoid(const Point& old_s, const Point& new_s, const std::vector<Polygon> neighborhood){
    double value = 0;

    double k = 15;
    double scale = 1.2;
    double m = 0.1;
    double k_s = 10;

    if(!is_in_the_neighborhood(new_s, neighborhood)){
        return INFINITY;
    }

    for(const Polygon& triangle : neighborhood){
        int index = find_point_index(old_s, triangle, false);

        Point a = triangle.vertex((index + 1) % 3);
        Point b = triangle.vertex((index + 2) % 3);

        double ab2 = CGAL::to_double(CGAL::squared_distance(a, b));
        double sa2 = CGAL::to_double(CGAL::squared_distance(new_s, a));
        double sb2 = CGAL::to_double(CGAL::squared_distance(new_s, b));

        double cos1 = (ab2 + sa2 - sb2)/(2*std::sqrt(ab2)*std::sqrt(sa2));
        double cos2 = (ab2 + sb2 - sa2)/(2*std::sqrt(ab2)*std::sqrt(sb2));
        double cos3 = (sa2 + sb2 - ab2)/(2*std::sqrt(sa2)*std::sqrt(sb2));

        double w1 = std::acos((ab2 + sa2 - sb2)/(2*std::sqrt(ab2)*std::sqrt(sa2))); // cos alpha
        double w2 = std::acos((ab2 + sb2 - sa2)/(2*std::sqrt(ab2)*std::sqrt(sb2))); // cos beta
        double w3 = std::acos((sa2 + sb2 - ab2)/(2*std::sqrt(sa2)*std::sqrt(sb2))); // cos gamma

        if( std::isnan(w1) || std::isnan(w2) || std::isnan(w3)){
            continue;
        }

        double sigm1 =  scale / (1 + std::exp(-k * (w1-M_PI/2))) + (m * (w1 - M_PI/2))/(1 + std::exp(-k_s * (w1-M_PI/2)));
        double sigm2 =  scale / (1 + std::exp(-k * (w2-M_PI/2))) + (m * (w2 - M_PI/2))/(1 + std::exp(-k_s * (w2-M_PI/2)));
        double sigm3 =  scale / (1 + std::exp(-k * (w3-M_PI/2))) + (m * (w3 - M_PI/2))/(1 + std::exp(-k_s * (w3-M_PI/2)));

        value = value + sigm1 + sigm2 + sigm3;
    }

    return value;
}

// f(x) = (x-60)^2
Point calculate_gradient_equilateral(Point& s, std::vector<Polygon>& triangles){
    double dx = 0.0;
    double dy = 0.0;

    for(const Polygon& triangle : triangles){
        int index = find_point_index(s, triangle, false);

        //std::cout << "index: " << index << std::endl;

        Point a = triangle.vertex((index + 1) % 3);
        Point b = triangle.vertex((index + 2) % 3);

        //std::cout << "a: " << a << " b: " << b << " steiner: " << s << std::endl;

        double ab2 = CGAL::to_double(CGAL::squared_distance(a, b));
        double sa2 = CGAL::to_double(CGAL::squared_distance(s, a));
        double sb2 = CGAL::to_double(CGAL::squared_distance(s, b));

        double w1 = (ab2 + sa2 - sb2)/(2*std::sqrt(ab2)*std::sqrt(sa2));
        double w2 = (ab2 + sb2 - sa2)/(2*std::sqrt(ab2)*std::sqrt(sb2));
        double w3 = (sa2 + sb2 - ab2)/(2*std::sqrt(sa2)*std::sqrt(sb2));

        double d1 = std::pow((acos(w1) - M_PI/3), 2);
        double d2 = std::pow((acos(w2) - M_PI/3), 2);
        double d3 = std::pow((acos(w3) - M_PI/3), 2);

        // Calculate d/dx
        double nx1 = 4*CGAL::to_double(b.x() - a.x())*std::sqrt(ab2)*std::sqrt(sa2) - (ab2 + sa2 - sb2)*std::sqrt(ab2)*2*CGAL::to_double(s.x() - a.x())/std::sqrt(sa2);
        double nx2 = 4*CGAL::to_double(a.x() - b.x())*std::sqrt(ab2)*std::sqrt(sb2) - (ab2 + sb2 - sa2)*std::sqrt(ab2)*2*CGAL::to_double(s.x() - b.x())/std::sqrt(sb2);
        double nx3 = 4*CGAL::to_double(2*s.x() - a.x() - b.x())*std::sqrt(sa2)*std::sqrt(sb2) - (sa2 + sb2 - ab2)*(std::sqrt(sb2)*2*CGAL::to_double(s.x() - a.x())/std::sqrt(sa2) + std::sqrt(sa2)*2*CGAL::to_double(s.x() - b.x())/std::sqrt(sb2));

        dx += 2*d1*(-1/std::sqrt(1 - w1*w1))*(nx1/(4*ab2*sa2));
        dx += 2*d2*(-1/std::sqrt(1 - w2*w2))*(nx2/(4*ab2*sb2));
        dx += 2*d3*(-1/std::sqrt(1 - w3*w3))*(nx3/(4*sa2*sb2));

        // Calculate d/dy
        double ny1 = 4*CGAL::to_double(b.y() - a.y())*std::sqrt(ab2)*std::sqrt(sa2) - (ab2 + sa2 - sb2)*std::sqrt(ab2)*2*CGAL::to_double(s.y() - a.y())/std::sqrt(sa2);
        double ny2 = 4*CGAL::to_double(a.y() - b.y())*std::sqrt(ab2)*std::sqrt(sb2) - (ab2 + sb2 - sa2)*std::sqrt(ab2)*2*CGAL::to_double(s.y() - b.y())/std::sqrt(sb2);
        double ny3 = 4*CGAL::to_double(2*s.y() - a.y() - b.y())*std::sqrt(sa2)*std::sqrt(sb2) - (sa2 + sb2 - ab2)*(std::sqrt(sb2)*2*CGAL::to_double(s.y() - a.y())/std::sqrt(sa2) + std::sqrt(sa2)*2*CGAL::to_double(s.y() - b.y())/std::sqrt(sb2));

        dy += 2*d1*(-1/std::sqrt(1 - w1*w1))*(ny1/(4*ab2*sa2));
        dy += 2*d2*(-1/std::sqrt(1 - w2*w2))*(ny2/(4*ab2*sb2));
        dy += 2*d3*(-1/std::sqrt(1 - w3*w3))*(ny3/(4*sa2*sb2));
    }

    return Point(dx,dy);
}

// f(x) = ln(1 + e^(theta-90))
Point calculate_gradient_ln(Point& s, std::vector<Polygon>& triangles){
    double dx = 0.0;
    double dy = 0.0;

    for(const Polygon& triangle : triangles){
        int index = find_point_index(s, triangle, false);

        Point a = triangle.vertex((index + 1) % 3);
        Point b = triangle.vertex((index + 2) % 3);

        double ab2 = CGAL::to_double(CGAL::squared_distance(a, b));
        double sa2 = CGAL::to_double(CGAL::squared_distance(s, a));
        double sb2 = CGAL::to_double(CGAL::squared_distance(s, b));

        double w1 = (ab2 + sa2 - sb2)/(2*std::sqrt(ab2)*std::sqrt(sa2)); // cos alpha
        double w2 = (ab2 + sb2 - sa2)/(2*std::sqrt(ab2)*std::sqrt(sb2)); // cos beta
        double w3 = (sa2 + sb2 - ab2)/(2*std::sqrt(sa2)*std::sqrt(sb2)); // cos gamma

        double e1 = std::exp(acos(w1) - M_PI/2);
        double e2 = std::exp(acos(w2) - M_PI/2);
        double e3 = std::exp(acos(w3) - M_PI/2);

        // Calculate d/dx
        double nx1 = 4*CGAL::to_double(b.x() - a.x())*std::sqrt(sa2)*std::sqrt(ab2) - (sa2 + ab2 - sb2)*std::sqrt(ab2)*2*CGAL::to_double(s.x() - a.x())/std::sqrt(sa2);
        double nx2 = 4*CGAL::to_double(a.x() - b.x())*std::sqrt(sb2)*std::sqrt(ab2) - (sb2 + ab2 - sa2)*std::sqrt(ab2)*2*CGAL::to_double(s.x() - b.x())/std::sqrt(sb2);
        double nx3 = 4*CGAL::to_double(2*s.x() - a.x() - b.x())*std::sqrt(sa2)*std::sqrt(sb2) - (sa2 + sb2 - ab2)*(std::sqrt(sb2)*2*CGAL::to_double(s.x() - a.x())/std::sqrt(sa2) + std::sqrt(sa2)*2*CGAL::to_double(s.x() - b.x())/std::sqrt(sb2));

        dx += (1/(1+e1))*e1*(-1/(std::sqrt(1-w1*w1)))*(nx1/(4*ab2*sa2));
        dx += (1/(1+e2))*e2*(-1/(std::sqrt(1-w2*w2)))*(nx2/(4*ab2*sb2));
        dx += (1/(1+e3))*e3*(-1/(std::sqrt(1-w3*w3)))*(nx3/(4*sa2*sb2));

        // Calculate d/dx
        double ny1 = 4*CGAL::to_double(b.y() - a.y())*std::sqrt(sa2)*std::sqrt(ab2) - (sa2 + ab2 - sb2)*std::sqrt(ab2)*2*CGAL::to_double(s.y() - a.y())/std::sqrt(sa2);
        double ny2 = 4*CGAL::to_double(a.y() - b.y())*std::sqrt(sb2)*std::sqrt(ab2) - (sb2 + ab2 - sa2)*std::sqrt(ab2)*2*CGAL::to_double(s.y() - b.y())/std::sqrt(sb2);
        double ny3 = 4*CGAL::to_double(2*s.y() - a.y() - b.y())*std::sqrt(sa2)*std::sqrt(sb2) - (sa2 + sb2 - ab2)*(std::sqrt(sb2)*2*CGAL::to_double(s.y() - a.y())/std::sqrt(sa2) + std::sqrt(sa2)*2*CGAL::to_double(s.y() - b.y())/std::sqrt(sb2));

        dy += (1/(1+e1))*e1*(-1/(std::sqrt(1-w1*w1)))*(ny1/(4*ab2*sa2));
        dy += (1/(1+e2))*e2*(-1/(std::sqrt(1-w2*w2)))*(ny2/(4*ab2*sb2));
        dy += (1/(1+e3))*e3*(-1/(std::sqrt(1-w3*w3)))*(ny3/(4*sa2*sb2));
    }

    return Point(dx,dy);
}

Point calculate_gradient_sigmoid(Point& s, std::vector<Polygon>& triangles){
    double dx = 0.0;
    double dy = 0.0;

    double k = 15;

    for(const Polygon& triangle : triangles){

        int index = find_point_index(s, triangle, false);

        Point a = triangle.vertex((index + 1) % 3);
        Point b = triangle.vertex((index + 2) % 3);

        double ab2 = CGAL::to_double(CGAL::squared_distance(a, b));
        double sa2 = CGAL::to_double(CGAL::squared_distance(s, a));
        double sb2 = CGAL::to_double(CGAL::squared_distance(s, b));;

        double w1 = (ab2 + sa2 - sb2)/(2*std::sqrt(ab2)*std::sqrt(sa2)); // cos alpha
        double w2 = (ab2 + sb2 - sa2)/(2*std::sqrt(ab2)*std::sqrt(sb2)); // cos beta
        double w3 = (sa2 + sb2 - ab2)/(2*std::sqrt(sa2)*std::sqrt(sb2)); // cos gamma

        double e1 = std::exp(-k*(acos(w1) - M_PI/2));
        double e2 = std::exp(-k*(acos(w2) - M_PI/2));
        double e3 = std::exp(-k*(acos(w3) - M_PI/2));

        // Calculate d/dx
        double nx1 = 4*CGAL::to_double(b.x() - a.x())*std::sqrt(sa2)*std::sqrt(ab2) - (sa2 + ab2 - sb2)*std::sqrt(ab2)*2*CGAL::to_double(s.x() - a.x())/std::sqrt(sa2);
        double nx2 = 4*CGAL::to_double(a.x() - b.x())*std::sqrt(sb2)*std::sqrt(ab2) - (sb2 + ab2 - sa2)*std::sqrt(ab2)*2*CGAL::to_double(s.x() - b.x())/std::sqrt(sb2);
        double nx3 = 4*CGAL::to_double(2*s.x() - a.x() - b.x())*std::sqrt(sa2)*std::sqrt(sb2) - (sa2 + sb2 - ab2)*(std::sqrt(sb2)*2*CGAL::to_double(s.x() - a.x())/std::sqrt(sa2) + std::sqrt(sa2)*2*CGAL::to_double(s.x() - b.x())/std::sqrt(sb2));

        dx += ((e1*k)/((1+e1)*(1+e1)))*(-1/(std::sqrt(1-w1*w1)))*(nx1/(4*ab2*sa2));
        dx += ((e2*k)/((1+e2)*(1+e2)))*(-1/(std::sqrt(1-w2*w2)))*(nx2/(4*ab2*sb2));
        dx += ((e3*k)/((1+e3)*(1+e3)))*(-1/(std::sqrt(1-w3*w3)))*(nx3/(4*sa2*sb2));

        // Calculate d/dx
        double ny1 = 4*CGAL::to_double(b.y() - a.y())*std::sqrt(sa2)*std::sqrt(ab2) - (sa2 + ab2 - sb2)*std::sqrt(ab2)*2*CGAL::to_double(s.y() - a.y())/std::sqrt(sa2);
        double ny2 = 4*CGAL::to_double(a.y() - b.y())*std::sqrt(sb2)*std::sqrt(ab2) - (sb2 + ab2 - sa2)*std::sqrt(ab2)*2*CGAL::to_double(s.y() - b.y())/std::sqrt(sb2);
        double ny3 = 4*CGAL::to_double(2*s.y() - a.y() - b.y())*std::sqrt(sa2)*std::sqrt(sb2) - (sa2 + sb2 - ab2)*(std::sqrt(sb2)*2*CGAL::to_double(s.y() - a.y())/std::sqrt(sa2) + std::sqrt(sa2)*2*CGAL::to_double(s.y() - b.y())/std::sqrt(sb2));

        dy += ((e1*k)/((1+e1)*(1+e1)))*(-1/(std::sqrt(1-w1*w1)))*(ny1/(4*ab2*sa2));
        dy += ((e2*k)/((1+e2)*(1+e2)))*(-1/(std::sqrt(1-w2*w2)))*(ny2/(4*ab2*sb2));
        dy += ((e3*k)/((1+e3)*(1+e3)))*(-1/(std::sqrt(1-w3*w3)))*(ny3/(4*sa2*sb2));
    }

    return Point(dx,dy);
}

Point calculate_gradient_refined_sigmoid(Point& s, std::vector<Polygon>& triangles){
    double scale = 1.2;
    double m = 0.1;
    double k_s = 10;

    Point gradSigmoid = calculate_gradient_sigmoid(s, triangles);
    double dx = scale*CGAL::to_double(gradSigmoid.x());
    double dy = scale*CGAL::to_double(gradSigmoid.y());

    for(const Polygon& triangle : triangles){
        int index = find_point_index(s, triangle, false);

        Point a = triangle.vertex((index + 1) % 3);
        Point b = triangle.vertex((index + 2) % 3);

        double ab2 = CGAL::to_double(CGAL::squared_distance(a, b));
        double sa2 = CGAL::to_double(CGAL::squared_distance(s, a));
        double sb2 = CGAL::to_double(CGAL::squared_distance(s, b));

        double w1 = (ab2 + sa2 - sb2)/(2*std::sqrt(ab2)*std::sqrt(sa2)); // cos alpha
        double w2 = (ab2 + sb2 - sa2)/(2*std::sqrt(ab2)*std::sqrt(sb2)); // cos beta
        double w3 = (sa2 + sb2 - ab2)/(2*std::sqrt(sa2)*std::sqrt(sb2)); // cos gamma

        double e1 = std::exp(-k_s*(acos(w1) - M_PI/2));
        double e2 = std::exp(-k_s*(acos(w2) - M_PI/2));
        double e3 = std::exp(-k_s*(acos(w3) - M_PI/2));

        // Calculate d/dx
        double nx1 = 4*CGAL::to_double(b.x() - a.x())*std::sqrt(sa2)*std::sqrt(ab2) - (sa2 + ab2 - sb2)*std::sqrt(ab2)*2*CGAL::to_double(s.x() - a.x())/std::sqrt(sa2);
        double nx2 = 4*CGAL::to_double(a.x() - b.x())*std::sqrt(sb2)*std::sqrt(ab2) - (sb2 + ab2 - sa2)*std::sqrt(ab2)*2*CGAL::to_double(s.x() - b.x())/std::sqrt(sb2);
        double nx3 = 4*CGAL::to_double(2*s.x() - a.x() - b.x())*std::sqrt(sa2)*std::sqrt(sb2) - (sa2 + sb2 - ab2)*(std::sqrt(sb2)*2*CGAL::to_double(s.x() - a.x())/std::sqrt(sa2) + std::sqrt(sa2)*2*CGAL::to_double(s.x() - b.x())/std::sqrt(sb2));

        double mx1 = (-1/(std::sqrt(1-w1*w1)))*(nx1/(4*ab2*sa2));
        double mx2 = (-1/(std::sqrt(1-w2*w2)))*(nx2/(4*ab2*sb2));
        double mx3 = (-1/(std::sqrt(1-w3*w3)))*(nx3/(4*sa2*sb2));

        dx += (m*(1+e1)*mx1 + m*(acos(w1) - M_PI/2)*e1*k_s*mx1)/((1+e1)*(1+e1));
        dx += (m*(1+e2)*mx2 + m*(acos(w2) - M_PI/2)*e2*k_s*mx2)/((1+e2)*(1+e2));
        dx += (m*(1+e3)*mx3 + m*(acos(w3) - M_PI/2)*e3*k_s*mx3)/((1+e3)*(1+e3));

        // Calculate d/dx
        double ny1 = 4*CGAL::to_double(b.y() - a.y())*std::sqrt(sa2)*std::sqrt(ab2) - (sa2 + ab2 - sb2)*std::sqrt(ab2)*2*CGAL::to_double(s.y() - a.y())/std::sqrt(sa2);
        double ny2 = 4*CGAL::to_double(a.y() - b.y())*std::sqrt(sb2)*std::sqrt(ab2) - (sb2 + ab2 - sa2)*std::sqrt(ab2)*2*CGAL::to_double(s.y() - b.y())/std::sqrt(sb2);
        double ny3 = 4*CGAL::to_double(2*s.y() - a.y() - b.y())*std::sqrt(sa2)*std::sqrt(sb2) - (sa2 + sb2 - ab2)*(std::sqrt(sb2)*2*CGAL::to_double(s.y() - a.y())/std::sqrt(sa2) + std::sqrt(sa2)*2*CGAL::to_double(s.y() - b.y())/std::sqrt(sb2));

        double my1 = (-1/(std::sqrt(1-w1*w1)))*(ny1/(4*ab2*sa2));
        double my2 = (-1/(std::sqrt(1-w2*w2)))*(ny2/(4*ab2*sb2));
        double my3 = (-1/(std::sqrt(1-w3*w3)))*(ny3/(4*sa2*sb2));

        dy += (m*(1+e1)*my1 + m*(acos(w1) - M_PI/2)*e1*k_s*my1)/((1+e1)*(1+e1));
        dy += (m*(1+e2)*my2 + m*(acos(w2) - M_PI/2)*e2*k_s*my2)/((1+e2)*(1+e2));
        dy += (m*(1+e3)*my3 + m*(acos(w3) - M_PI/2)*e3*k_s*my3)/((1+e3)*(1+e3));
    }

    return Point(dx,dy);
}

Point calculate_gradient_sigmoid_equilateral(Point& s, std::vector<Polygon>& triangles){
    double dx = 0.0;
    double dy = 0.0;

    double k = 1;

    for(const Polygon& triangle : triangles){

        int index = find_point_index(s, triangle, false);

        Point a = triangle.vertex((index + 1) % 3);
        Point b = triangle.vertex((index + 2) % 3);

        double ab2 = CGAL::to_double(CGAL::squared_distance(a, b));
        double sa2 = CGAL::to_double(CGAL::squared_distance(s, a));
        double sb2 = CGAL::to_double(CGAL::squared_distance(s, b));;

        double w1 = (ab2 + sa2 - sb2)/(2*std::sqrt(ab2)*std::sqrt(sa2)); // cos alpha
        double w2 = (ab2 + sb2 - sa2)/(2*std::sqrt(ab2)*std::sqrt(sb2)); // cos beta
        double w3 = (sa2 + sb2 - ab2)/(2*std::sqrt(sa2)*std::sqrt(sb2)); // cos gamma

        double e1 = std::exp(-k*(acos(w1) - M_PI/3));
        double e2 = std::exp(-k*(acos(w2) - M_PI/3));
        double e3 = std::exp(-k*(acos(w3) - M_PI/3));

        // Calculate d/dx
        double nx1 = 4*CGAL::to_double(b.x() - a.x())*std::sqrt(sa2)*std::sqrt(ab2) - (sa2 + ab2 - sb2)*std::sqrt(ab2)*2*CGAL::to_double(s.x() - a.x())/std::sqrt(sa2);
        double nx2 = 4*CGAL::to_double(a.x() - b.x())*std::sqrt(sb2)*std::sqrt(ab2) - (sb2 + ab2 - sa2)*std::sqrt(ab2)*2*CGAL::to_double(s.x() - b.x())/std::sqrt(sb2);
        double nx3 = 4*CGAL::to_double(2*s.x() - a.x() - b.x())*std::sqrt(sa2)*std::sqrt(sb2) - (sa2 + sb2 - ab2)*(std::sqrt(sb2)*2*CGAL::to_double(s.x() - a.x())/std::sqrt(sa2) + std::sqrt(sa2)*2*CGAL::to_double(s.x() - b.x())/std::sqrt(sb2));

        dx += ((e1*k)/((1+e1)*(1+e1)))*(-1/(std::sqrt(1-w1*w1)))*(nx1/(4*ab2*sa2));
        dx += ((e2*k)/((1+e2)*(1+e2)))*(-1/(std::sqrt(1-w2*w2)))*(nx2/(4*ab2*sb2));
        dx += ((e3*k)/((1+e3)*(1+e3)))*(-1/(std::sqrt(1-w3*w3)))*(nx3/(4*sa2*sb2));

        // Calculate d/dx
        double ny1 = 4*CGAL::to_double(b.y() - a.y())*std::sqrt(sa2)*std::sqrt(ab2) - (sa2 + ab2 - sb2)*std::sqrt(ab2)*2*CGAL::to_double(s.y() - a.y())/std::sqrt(sa2);
        double ny2 = 4*CGAL::to_double(a.y() - b.y())*std::sqrt(sb2)*std::sqrt(ab2) - (sb2 + ab2 - sa2)*std::sqrt(ab2)*2*CGAL::to_double(s.y() - b.y())/std::sqrt(sb2);
        double ny3 = 4*CGAL::to_double(2*s.y() - a.y() - b.y())*std::sqrt(sa2)*std::sqrt(sb2) - (sa2 + sb2 - ab2)*(std::sqrt(sb2)*2*CGAL::to_double(s.y() - a.y())/std::sqrt(sa2) + std::sqrt(sa2)*2*CGAL::to_double(s.y() - b.y())/std::sqrt(sb2));

        dy += ((e1*k)/((1+e1)*(1+e1)))*(-1/(std::sqrt(1-w1*w1)))*(ny1/(4*ab2*sa2));
        dy += ((e2*k)/((1+e2)*(1+e2)))*(-1/(std::sqrt(1-w2*w2)))*(ny2/(4*ab2*sb2));
        dy += ((e3*k)/((1+e3)*(1+e3)))*(-1/(std::sqrt(1-w3*w3)))*(ny3/(4*sa2*sb2));
    }

    return Point(dx,dy);
}

Point calculate_gradient_refined_sigmoid_equilateral(Point& s, std::vector<Polygon>& triangles){
    double scale = 1;
    double m = 0.1;
    double k_s = 1;

    Point gradSigmoid = calculate_gradient_sigmoid(s, triangles);
    double dx = scale*CGAL::to_double(gradSigmoid.x());
    double dy = scale*CGAL::to_double(gradSigmoid.y());

    for(const Polygon& triangle : triangles){
        int index = find_point_index(s, triangle, false);

        Point a = triangle.vertex((index + 1) % 3);
        Point b = triangle.vertex((index + 2) % 3);

        double ab2 = CGAL::to_double(CGAL::squared_distance(a, b));
        double sa2 = CGAL::to_double(CGAL::squared_distance(s, a));
        double sb2 = CGAL::to_double(CGAL::squared_distance(s, b));

        double w1 = (ab2 + sa2 - sb2)/(2*std::sqrt(ab2)*std::sqrt(sa2)); // cos alpha
        double w2 = (ab2 + sb2 - sa2)/(2*std::sqrt(ab2)*std::sqrt(sb2)); // cos beta
        double w3 = (sa2 + sb2 - ab2)/(2*std::sqrt(sa2)*std::sqrt(sb2)); // cos gamma

        double e1 = std::exp(-k_s*(acos(w1) - M_PI/3));
        double e2 = std::exp(-k_s*(acos(w2) - M_PI/3));
        double e3 = std::exp(-k_s*(acos(w3) - M_PI/3));

        // Calculate d/dx
        double nx1 = 4*CGAL::to_double(b.x() - a.x())*std::sqrt(sa2)*std::sqrt(ab2) - (sa2 + ab2 - sb2)*std::sqrt(ab2)*2*CGAL::to_double(s.x() - a.x())/std::sqrt(sa2);
        double nx2 = 4*CGAL::to_double(a.x() - b.x())*std::sqrt(sb2)*std::sqrt(ab2) - (sb2 + ab2 - sa2)*std::sqrt(ab2)*2*CGAL::to_double(s.x() - b.x())/std::sqrt(sb2);
        double nx3 = 4*CGAL::to_double(2*s.x() - a.x() - b.x())*std::sqrt(sa2)*std::sqrt(sb2) - (sa2 + sb2 - ab2)*(std::sqrt(sb2)*2*CGAL::to_double(s.x() - a.x())/std::sqrt(sa2) + std::sqrt(sa2)*2*CGAL::to_double(s.x() - b.x())/std::sqrt(sb2));

        double mx1 = (-1/(std::sqrt(1-w1*w1)))*(nx1/(4*ab2*sa2));
        double mx2 = (-1/(std::sqrt(1-w2*w2)))*(nx2/(4*ab2*sb2));
        double mx3 = (-1/(std::sqrt(1-w3*w3)))*(nx3/(4*sa2*sb2));

        dx += (m*(1+e1)*mx1 + m*(acos(w1) - M_PI/3)*e1*k_s*mx1)/((1+e1)*(1+e1));
        dx += (m*(1+e2)*mx2 + m*(acos(w2) - M_PI/3)*e2*k_s*mx2)/((1+e2)*(1+e2));
        dx += (m*(1+e3)*mx3 + m*(acos(w3) - M_PI/3)*e3*k_s*mx3)/((1+e3)*(1+e3));

        // Calculate d/dx
        double ny1 = 4*CGAL::to_double(b.y() - a.y())*std::sqrt(sa2)*std::sqrt(ab2) - (sa2 + ab2 - sb2)*std::sqrt(ab2)*2*CGAL::to_double(s.y() - a.y())/std::sqrt(sa2);
        double ny2 = 4*CGAL::to_double(a.y() - b.y())*std::sqrt(sb2)*std::sqrt(ab2) - (sb2 + ab2 - sa2)*std::sqrt(ab2)*2*CGAL::to_double(s.y() - b.y())/std::sqrt(sb2);
        double ny3 = 4*CGAL::to_double(2*s.y() - a.y() - b.y())*std::sqrt(sa2)*std::sqrt(sb2) - (sa2 + sb2 - ab2)*(std::sqrt(sb2)*2*CGAL::to_double(s.y() - a.y())/std::sqrt(sa2) + std::sqrt(sa2)*2*CGAL::to_double(s.y() - b.y())/std::sqrt(sb2));

        double my1 = (-1/(std::sqrt(1-w1*w1)))*(ny1/(4*ab2*sa2));
        double my2 = (-1/(std::sqrt(1-w2*w2)))*(ny2/(4*ab2*sb2));
        double my3 = (-1/(std::sqrt(1-w3*w3)))*(ny3/(4*sa2*sb2));

        dy += (m*(1+e1)*my1 + m*(acos(w1) - M_PI/3)*e1*k_s*my1)/((1+e1)*(1+e1));
        dy += (m*(1+e2)*my2 + m*(acos(w2) - M_PI/3)*e2*k_s*my2)/((1+e2)*(1+e2));
        dy += (m*(1+e3)*my3 + m*(acos(w3) - M_PI/3)*e3*k_s*my3)/((1+e3)*(1+e3));
    }

    return Point(dx,dy);
}

double get_step_size(Problem* problem){
    Polygon polygon = problem->get_boundary();
    
    auto bbox = polygon.bbox();
    double max_size = std::sqrt(pow(bbox.x_span(),2) + pow(bbox.y_span(),2)) * 0.1;

    int count = 0;

    while (max_size >= 10.0) {
        max_size /= 10.0;
        count++;
    }
    while (max_size < 1.0) {
        max_size *= 10.0;
        count--;
    }

    double step_size = std::pow(10, count);
    return step_size;
}

Point locally_optimize_position(Point steiner, std::vector<Polygon>& triangles, Problem *problem){

    double step_size = 1e6;
    //double step_size = get_step_size(problem);
    //std::cout << "Step size: " << step_size << std::endl;
    //double step_size = compute_step_size(problem);
    //std::cout << "Computed step size: " << step_size << std::endl;

    double TOL = 1e-7;
    const int MAX_ITER = 1000;

    Point s = steiner;

    Segment constraint;

    int iteration = 0;

    std::vector<Polygon> neighborhood = find_neighborhood(s, triangles, false);

    try{
        // for now only optimize positions of interior points
        /*if(!is_interior_vertex(s, problem->get_boundary())){
            Segment boundaryConstraint = find_boundary_segment(s, problem->get_boundary());
            Point tmp = move_point_on_segment(s, boundaryConstraint, 10000, neighborhood);
            update_neighborhood(triangles, steiner, tmp);

            std::cout << "Steiner: " << s << " Moved steiner: " << tmp << std::endl;

        } else if (is_on_constraint(s, problem, &constraint)){
            Point tmp = move_point_on_segment(s, constraint, 10000, neighborhood);
            update_neighborhood(triangles, steiner, tmp);

        } else {*/
        if((is_interior_vertex_with_tolerance(s, problem->get_boundary())) && !is_on_constraint_with_tolerance(s, problem, &constraint)){
            Point gradient = calculate_gradient_refined_sigmoid(s, neighborhood);

            //std::cout << "gradient norm: " << norm(gradient) << std::endl;
            while(norm(gradient) > TOL && iteration < MAX_ITER){

                // * Armijo version
                /*Point tmp = line_search(triangles, {s}, {gradient}, step_size)[0];
                update_neighborhood(neighborhood, s, tmp);
                s = tmp;*/

                // * Use this code if not using Armijo
                Point tmp(s.x() - step_size*gradient.x(), s.y() - step_size*gradient.y());
                // Only proceed if the new point position is not outside the initial neighborhood
                if(is_in_the_neighborhood(tmp, neighborhood)){              
                    update_neighborhood(neighborhood, s, tmp);
                    s = tmp;
                }
                else{
                    std::cout << "I'm outside!!!\n";
                    step_size /= 10;
                    //break;
                }

                gradient = calculate_gradient_refined_sigmoid(s, neighborhood);

                iteration++;
            }

            //std::cout << "start point: " << steiner << " end point: " << s << std::endl;

            // update triangles with a new steiner point position
            update_neighborhood(triangles, steiner, s);

            //std::cout << "optimized gradient norm: " << norm(gradient) << std::endl;
        } 
    }catch(const CGAL::Assertion_exception& e){
        std::cout << "CGAL Exception 3 \n" << std::endl;
    }

    return s;
}