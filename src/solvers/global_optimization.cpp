#include "global_optimization.hpp"

double get_signed_area(const Point& p1, const Point& p2, const Point& p3){
    auto area = (p1.x() * (p2.y() - p3.y()) + p2.x() * (p3.y() - p1.y()) + p3.x() * (p1.y() - p2.y()))/2.0;
    return CGAL::to_double(area);
}

double norm_double(const Point& gradient) {
    double x = CGAL::to_double(gradient.x());
    double y = CGAL::to_double(gradient.y());

    return std::sqrt(x*x + y*y);
}

bool norm_is_smaller(std::vector<Point>& gradients, const double TOL){
    bool result = true;
    for(const Point& g : gradients){
        if(norm_double(g) > TOL){
            result = false;
        }
    }
    return result;
}

// Laplace function from "Harmonic triangulation"
Point calculate_gradient_point_laplace(const Point& s, std::vector<Polygon>& neighborhood, Problem* problem){
    double dx = 0.0;
    double dy = 0.0;

    for(Polygon triangle : neighborhood){
        int index1 = find_point_index(s, triangle, false);

        Point a = triangle.vertex((index1 + 1) % 3);
        Point b = triangle.vertex((index1 + 2) % 3);

        double ab2 = CGAL::to_double(CGAL::squared_distance(a, b));

        // Order of points s a b does not influence the derivative
        double area = get_signed_area(s, a, b);
        int coeff = (area < 0)? -1 : 1;

        // Calculate dx part 1
        dx += ab2*0.5*(-1*coeff*CGAL::to_double(a.y() - b.y()))/(area * area);

        // Calculate dy part 1
        dy += 0.5*ab2*(-1*coeff*CGAL::to_double(b.x() - a.x()))/(area * area);

        std::vector<Point> points = problem->get_points();

        // Check whether a is also a steiner point
        auto it1 = std::find(points.begin(), points.end(), a);
        if(it1 == points.end()){
            double sb2 = CGAL::to_double(CGAL::squared_distance(s, b));
            
            dx += (2*CGAL::to_double(s.x() - b.x())*std::abs(area) - sb2*0.5*(-1*coeff*CGAL::to_double(a.y() - b.y())))/(area * area);
            dy += (2*CGAL::to_double(s.y() - b.y())*std::abs(area) - sb2*0.5*(-1*coeff*CGAL::to_double(b.x() - a.x())))/(area * area);
        }

        // Check whether b is also a steiner point
        auto it2 = std::find(points.begin(), points.end(), b);
        if(it2 == points.end()){
            double sa2 = CGAL::to_double(CGAL::squared_distance(s, a));
            
            dx += (2*CGAL::to_double(s.x() - a.x())*std::abs(area) - sa2*0.5*(-1*coeff*CGAL::to_double(a.y() - b.y())))/(area * area);
            dy += (2*CGAL::to_double(s.y() - a.y())*std::abs(area) - sa2*0.5*(-1*coeff*CGAL::to_double(b.x() - a.x())))/(area * area);
        }
    }
    return Point(dx,dy);
}

std::vector<Point> calculate_gradient(std::vector<Point>& steiner_points, std::vector<Polygon>& triangles, Problem *problem, bool debug){
    std::vector<Point> gradients;
    Segment constraint;
    int count = 0;
    for(Point steiner : steiner_points){
        if(is_interior_vertex_with_tolerance(steiner, problem->get_boundary()) && !is_on_constraint_with_tolerance(steiner, problem, &constraint)){
            Point gradient(0,0);
            try{
                std::vector<Polygon> neighborhood = find_neighborhood(steiner, triangles, (debug && count==80));
                //gradient = calculate_gradient_point_laplace(steiner, neighborhood, problem);
                gradient = calculate_gradient_refined_sigmoid(steiner, neighborhood);
            }catch(const CGAL::Assertion_exception& e){
                std::cout << "CGAL Exception for point " << steiner << std::endl;
            }
            gradients.push_back(gradient);
        } else{
            gradients.push_back(Point(0,0));
        }
        count++;
    }
    return gradients;
}

std::vector<Point> globally_optimize_position(std::vector<Point>& steiner_points, std::vector<Polygon>& triangles, Problem *problem, bool debug){
    double TOL = 1e-7;
    const int MAX_ITER = 100;

    std::vector<double> step_sizes(steiner_points.size(), 1e-1);
    double step_size = 1e-1;

    int iteration = 0;

    std::vector<Point> s = steiner_points;

    try{
        std::vector<Point> gradients = calculate_gradient(s, triangles, problem, debug);

        while(!norm_is_smaller(gradients, TOL) && iteration < MAX_ITER){

            //* Armijo
            /*std::vector<Point> tmp = line_search(triangles, s, gradients, step_size);
            for(int i = 0; i<steiner_points.size(); i++){
                update_neighborhood(triangles, s[i], tmp[i]);
            }
            s = tmp;*/

            //* Gradient descent
            for(int i = 0; i<steiner_points.size(); i++){
                
                Point tmp(s[i].x() - step_sizes[i]*gradients[i].x(), s[i].y() - step_sizes[i]*gradients[i].y());  

                //check whether the new point position is valid
                std::vector<Polygon> neighborhood = find_neighborhood(s[i], triangles, false);
                if(is_in_the_neighborhood(tmp, neighborhood)){
                    update_neighborhood(triangles, s[i], tmp);
                    s[i] = tmp; 
                }  
                else{
                    step_sizes[i] = step_sizes[i]/10;
                    //std::cout << "Outside the neighborhood: index - " << i << std::endl; 
                }  
            }  

            gradients = calculate_gradient(s, triangles, problem, false);
            iteration++;
        }
    }catch(const CGAL::Assertion_exception& e){
        std::cout << "CGAL Exception 3 \n" << std::endl;
    }

    return s;
}