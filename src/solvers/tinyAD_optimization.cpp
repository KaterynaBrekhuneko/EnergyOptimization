#include "../../libs/TinyAD/include/TinyAD/ScalarFunction.hh"
#include "../../libs/TinyAD/include/TinyAD/Utils/LineSearch.hh"
#include "../../libs/TinyAD/include/TinyAD/Utils/NewtonDirection.hh"
#include "../../libs/TinyAD/include/TinyAD/Utils/NewtonDecrement.hh"

#include <iostream>
#include <vector>
#include <unordered_map>

#include <Eigen/Dense>

#include "tinyAD_optimization.hpp"

// Energy function: ln(1 + e^(x - π/2))
template <typename T>
T angle_cost_ln( Problem* problem,
        const Eigen::Vector2<T>& a,
        const Eigen::Vector2<T>& b,
        const Eigen::Vector2<T>& c)
{
    // Compute all three normalized triangle edge vectors
    Eigen::Vector2<T> ab = (b - a).normalized();
    Eigen::Vector2<T> bc = (c - b).normalized();
    Eigen::Vector2<T> ca = (a - c).normalized();

    // Compute cosine of all 3 angles (unsigned)
    T cos_angle_a = (-ca).dot(ab);
    T cos_angle_b = (-ab).dot(bc);
    T cos_angle_c = (-bc).dot(ca);
   
    if(abs(cos_angle_a) >= 1){
        std::cout << "Cos out of range for Point a: " << a[0] << " " << a[1] << std::endl;
    }
    if(abs(cos_angle_b) >= 1){
        std::cout << "Cos out of range for Point b: " << b[0] << " " << b[1] << std::endl;
    }
    if(abs(cos_angle_c) >= 1){
        std::cout << "Cos out of range for Point c: " << c[0] << " " << c[1] << std::endl;
    }
        
    // Calculate and return the sigmoid term
    T angle_a = acos(cos_angle_a);
    T angle_b = acos(cos_angle_b);
    T angle_c = acos(cos_angle_c);

    T ln_a = log(1+exp(angle_a - M_PI/2));
    T ln_b = log(1+exp(angle_b - M_PI/2));
    T ln_c = log(1+exp(angle_c - M_PI/2));

    return ln_a + ln_b + ln_c;
}

// Energy function: 1 / (1 + exp(-k(x - π/2)))
template <typename T>
T angle_cost_sigmoid( Problem* problem,
        const Eigen::Vector2<T>& a,
        const Eigen::Vector2<T>& b,
        const Eigen::Vector2<T>& c)
{
    double k = 1;

    // Compute all three normalized triangle edge vectors
    Eigen::Vector2<T> ab = (b - a).normalized();
    Eigen::Vector2<T> bc = (c - b).normalized();
    Eigen::Vector2<T> ca = (a - c).normalized();

    // Compute cosine of all 3 angles (unsigned)
    T cos_angle_a = (-ca).dot(ab);
    T cos_angle_b = (-ab).dot(bc);
    T cos_angle_c = (-bc).dot(ca);

    // Check if all values are in the interval [-1, 1]
    /*if(abs(cos_angle_a) >= 1 || abs(cos_angle_b) >= 1 || abs(cos_angle_c) >= 1){
        auto bbox = problem->get_boundary().bbox();
        auto scale = std::max(box.xmax() - box.xmin(), box.ymax() - box.ymin());
        std::cout << "Point: " << Point((p.x() - box.xmin()) * 560 / scale + 16, p.y() * 560 / scale + 272) << std::endl;

        TINYAD_WARNING("cosine out of range!");
    }*/
   
    if(abs(cos_angle_a) >= 1){
        std::cout << "Cos out of range for Point a: " << a[0] << " " << a[1] << std::endl;
    }
    if(abs(cos_angle_b) >= 1){
        std::cout << "Cos out of range for Point b: " << b[0] << " " << b[1] << std::endl;
    }
    if(abs(cos_angle_c) >= 1){
        std::cout << "Cos out of range for Point c: " << c[0] << " " << c[1] << std::endl;
    }
        
    // Calculate and return the sigmoid term
    T angle_a = acos(cos_angle_a);
    T angle_b = acos(cos_angle_b);
    T angle_c = acos(cos_angle_c);

    T sigmoid_a = 1/(1+exp(-k*(angle_a - M_PI/2)));
    T sigmoid_b = 1/(1+exp(-k*(angle_b - M_PI/2)));
    T sigmoid_c = 1/(1+exp(-k*(angle_c - M_PI/2)));

    return sigmoid_a + sigmoid_b + sigmoid_c;
}

// Energy function: s / (1 + exp(-k(x - π/2))) + m(x - π/2) / (1 + exp(-k_s(x - π/2)))
template <typename T>
T angle_cost_refined_sigmoid( Problem* problem,
        const Eigen::Vector2<T>& a,
        const Eigen::Vector2<T>& b,
        const Eigen::Vector2<T>& c)
{
    double k = 1;
    double s = 1;
    double m = 0.1;
    double k_s = 1;

    // Compute all three normalized triangle edge vectors
    Eigen::Vector2<T> ab = (b - a).normalized();
    Eigen::Vector2<T> bc = (c - b).normalized();
    Eigen::Vector2<T> ca = (a - c).normalized();

    // Compute cosine of all 3 angles (unsigned)
    T cos_angle_a = (-ca).dot(ab);
    T cos_angle_b = (-ab).dot(bc);
    T cos_angle_c = (-bc).dot(ca);

    if(abs(cos_angle_a) >= 1){
        std::cout << "Cos out of range for Point a: " << a[0] << " " << a[1] << std::endl;
    }
    if(abs(cos_angle_b) >= 1){
        std::cout << "Cos out of range for Point b: " << b[0] << " " << b[1] << std::endl;
    }
    if(abs(cos_angle_c) >= 1){
        std::cout << "Cos out of range for Point c: " << c[0] << " " << c[1] << std::endl;
    }
        
    // Calculate and return the sigmoid term
    T angle_a = acos(cos_angle_a);
    T angle_b = acos(cos_angle_b);
    T angle_c = acos(cos_angle_c);

    T sigmoid_a = s/(1+exp(-k*(angle_a - M_PI/2))) + (m*(angle_a - M_PI/2))/(1+exp(-k_s*(angle_a - M_PI/2)));
    T sigmoid_b = s/(1+exp(-k*(angle_b - M_PI/2))) + (m*(angle_b - M_PI/2))/(1+exp(-k_s*(angle_b - M_PI/2)));
    T sigmoid_c = s/(1+exp(-k*(angle_c - M_PI/2))) + (m*(angle_c - M_PI/2))/(1+exp(-k_s*(angle_c - M_PI/2)));

    return sigmoid_a + sigmoid_b + sigmoid_c;
}

template <typename T>
T cotangent(const Eigen::Vector2<T>& u, const Eigen::Vector2<T>& v) {
    T dot_product = u.dot(v);
    T determinant = u.x() * v.y() - u.y() * v.x();  // 2D cross product
    
    if (determinant == 0) {
        std::cerr << "Warning: Zero determinant, undefined cotangent." << std::endl;
        return std::numeric_limits<T>::infinity();
    }

    return dot_product / determinant;
}

// Energy function: 
template <typename T>
T angle_cost_laplace( Problem* problem,
        const Eigen::Vector2<T>& a,
        const Eigen::Vector2<T>& b,
        const Eigen::Vector2<T>& c)
{
    // Compute all three normalized triangle edge vectors
    Eigen::Vector2<T> ab = (b - a);
    Eigen::Vector2<T> bc = (c - b);
    Eigen::Vector2<T> ca = (a - c);

    // Compute cotangents of all three angles
    T cot_a = cotangent<T>(-ca, ab);
    T cot_b = cotangent<T>(-ab, bc);
    T cot_c = cotangent<T>(-bc, ca);

    return abs(cot_a) + abs(cot_b) + abs(cot_c);
}

int find_obtuse_angle(Eigen::Vector2d& a, Eigen::Vector2d& b, Eigen::Vector2d& c){
    std::array<Eigen::Vector2d, 3> points = {a, b, c};
    for (int i = 0; i < 3; i++) {
        auto v1 = points[(i + 1) % 3] - points[i];
        auto v2 = points[(i + 2) % 3] - points[i];
        if (v1.dot(v2) < 0) {
            return i;
        }
    }
    return -1;
}

bool is_on_constraint(Point& steiner, std::vector<Segment>& constraints, Segment* constraint){
    for(Segment c : constraints){
        double distance = std::sqrt(CGAL::to_double(CGAL::squared_distance(c, steiner)));
        if(c.has_on(steiner) /*|| distance < 1e-6*/){
            *constraint = c;
            return true;
        }
    }
    return false;
}

bool is_on_boundary(Point& steiner, Polygon& boundary, Segment* constraint){
    int n = boundary.size();
    for (int i = 0; i < n; ++i) {
        Segment edge(boundary.vertex(i), boundary.vertex((i + 1) % n)); 
        double distance = std::sqrt(CGAL::to_double(CGAL::squared_distance(edge, steiner)));
        if (edge.has_on(steiner) /*|| distance < 1e-6*/) {
            *constraint = edge;
            return true; 
        }
    }
    return false;
}

int find_point_idx(Point& p, std::vector<Point>& points, std::vector<Point>& steiner){
    for(int i = 0; i<points.size(); i++){
        if(p == points[i]){
            return i;
        }
    }

    for(int j = 0; j<steiner.size(); j++){
        if(p == steiner[j]){
            return j + points.size();
        }
    }

    throw std::runtime_error("Point not found in points or steiner!");
}

const bool is_constrained_point(Point& p, Polygon& boundary, std::vector<Segment>& constraints, Segment* constraint) {
    return is_on_constraint(p, constraints, constraint) || is_on_boundary(p, boundary, constraint);
}

std::tuple<std::vector<Point>, std::vector<int>, std::vector<Segment>> get_constrained_points(std::vector<Point>& steiner, Polygon& boundary, std::vector<Segment>& constraints){
    std::vector<Point> constrained_points;
    std::vector<Segment> constraints_of_points;
    std::vector<int> constrained_indices;
    for(int i = 0; i< steiner.size(); i++){
        Point s = steiner[i];
        Segment constraint;
        if(is_constrained_point(s, boundary, constraints, &constraint)){
            constrained_points.push_back(s);
            constrained_indices.push_back(i);
            constraints_of_points.push_back(constraint);
        }
    }
    return {constrained_points, constrained_indices, constraints_of_points};
}

double signed_area(const Eigen::Vector2d& a, const Eigen::Vector2d& b, const Eigen::Vector2d& c) {
    Eigen::Matrix2d M = TinyAD::col_mat(b - a, c - a);
    return M.determinant();
}

void flip_edge(Problem* problem, Eigen::MatrixXd& V, Eigen::MatrixXi& F, int tri1, int tri2, int v1, int v2, bool consider_energy) {
    // Get the vertices of the triangles
    int a = F(tri1, v1);
    int b = F(tri1, (v1 + 1) % 3);
    int c = F(tri1, (v1 + 2) % 3);
    int d = F(tri2, (v2 + 2) % 3);

    // Only flip if it benefits energy
    double energy_before = angle_cost_sigmoid<double>(problem, V.row(a), V.row(b), V.row(c)) + angle_cost_sigmoid<double>(problem, V.row(d), V.row(b), V.row(c));
    double energy_after = angle_cost_sigmoid<double>(problem, V.row(a), V.row(b), V.row(d)) + angle_cost_sigmoid<double>(problem, V.row(a), V.row(c), V.row(d));

    if(!consider_energy || energy_after <= energy_before){
        // Replace the triangles with the flipped configuration 
        //(find the correct point order so the triangle will not be considered flipped during optimization)
        if(signed_area(V.row(a), V.row(b), V.row(d)) > 0){
            F.row(tri1) = Eigen::Vector3i(a, b, d);
        } else if(signed_area(V.row(a), V.row(d), V.row(b)) > 0){
            F.row(tri1) = Eigen::Vector3i(a, d, b);
        } else {
            std::cout << "1: " << signed_area(V.row(a), V.row(b), V.row(d)) << " 2: " << signed_area(V.row(a), V.row(d), V.row(b)) << std::endl;
        }

        if(signed_area(V.row(a), V.row(c), V.row(d)) > 0){
            F.row(tri2) = Eigen::Vector3i(a, c, d);
        } else if(signed_area(V.row(a), V.row(d), V.row(c)) > 0){
            F.row(tri2) = Eigen::Vector3i(a, d, c);
        } else {
            std::cout << "1: " << signed_area(V.row(a), V.row(c), V.row(d)) << " 2: " << signed_area(V.row(a), V.row(d), V.row(c)) << std::endl;
        }
    }
}

// Main function to find and flip all obtuse triangles
void flip_obtuse_triangles(Problem* problem, Eigen::MatrixXd& V, Eigen::MatrixXi& F, bool consider_energy) {
    for (int i = 0; i < F.rows(); ++i) {
        // Triangle vertices
        Eigen::Vector2d a = V.row(F(i, 0));
        Eigen::Vector2d b = V.row(F(i, 1));
        Eigen::Vector2d c = V.row(F(i, 2));
        int obtuseIdx = find_obtuse_angle(a, b, c);

        if (obtuseIdx != -1) {
            // Find the neighboring triangle that shares this edge
            int v1 = F(i, (obtuseIdx + 1) % 3);
            int v2 = F(i, (obtuseIdx + 2) % 3);
            
            // Look for a neighboring triangle
            for (int j = 0; j < F.rows(); ++j) {
                if (j == i) continue;
                for (int k = 0; k < 3; ++k) {
                    if ((F(j, k) == v1 && F(j, (k + 1) % 3) == v2) ||
                        (F(j, k) == v2 && F(j, (k + 1) % 3) == v1)) {
                        // Perform the edge flip
                        flip_edge(problem, V, F, i, j, obtuseIdx, k, consider_energy);
                        break;
                    }
                }
            }
        }
    }   
}

void alter_direction(Eigen::VectorXd& d, Eigen::VectorXi& BS, Eigen::MatrixXd& BS_VAR, Eigen::MatrixXd& BS_TANGENT){
    for(int i = 0;  i < BS.size(); i++){
        Eigen::VectorXd dir_1 = d.row(2*BS[i]);
        Eigen::VectorXd dir_2 = d.row(2*BS[i]+1);
        Eigen::Vector2d dir;
        dir[0] = dir_1[0];
        dir[1] = dir_2[0];

        Eigen::Vector4d constraint = BS_TANGENT.row(i);
        Eigen::Vector2d start = constraint.segment<2>(0);
        Eigen::Vector2d end = constraint.segment<2>(2);
        Eigen::Vector2d seg = end - start;

        double projectionFactor = dir.dot(seg) / seg.squaredNorm();
        Eigen::Vector2d new_dir = projectionFactor * seg;
        d[2*BS[i]] = new_dir[0];
        d[2*BS[i]+1] = new_dir[1];
    }
}

void optimizeTinyAD(Problem* problem){
    auto points = problem->get_points(); 
    auto steiner = problem->get_steiner(); 
    auto triangles = problem->get_triangulation();
    auto boundary = problem->get_boundary();
    auto constraints = problem->get_constraints();

    // Convert all points to a vector
    Eigen::MatrixXd V(points.size() + steiner.size(), 2);
    for (size_t i = 0; i < points.size(); ++i) {
        V(i, 0) = CGAL::to_double(points[i].x()); 
        V(i, 1) = CGAL::to_double(points[i].y()); 
    }
    for (size_t i = points.size(); i < points.size() + steiner.size(); ++i) {
        V(i, 0) = CGAL::to_double(steiner[i - points.size()].x()); 
        V(i, 1) = CGAL::to_double(steiner[i - points.size()].y()); 
    }

    // Convert triangles to a matrix
    Eigen::MatrixXi F(triangles.size(), 3);
    for (size_t i = 0; i < triangles.size(); ++i) {
        const Polygon& polygon = triangles[i];

        Point a = polygon.vertex(0);
        Point b = polygon.vertex(1);
        Point c = polygon.vertex(2);

        F(i, 0) = find_point_idx(a, points, steiner);
        F(i, 1) = find_point_idx(b, points, steiner);
        F(i, 2) = find_point_idx(c, points, steiner);
    }

    // Define barrier terms for the constraints
    auto [constrained_points, constrained_indices, constraints_of_points] = get_constrained_points(steiner, boundary, constraints);
    Eigen::VectorXi B(points.size());
    Eigen::MatrixXd B_VAR(points.size(), 2);
    for (size_t i = 0; i < points.size(); ++i) {
        B_VAR(i, 0) = CGAL::to_double(points[i].x()); 
        B_VAR(i, 1) = CGAL::to_double(points[i].y()); 
        B(i) = i;
    }

    Eigen::VectorXi BS(constrained_points.size());
    Eigen::MatrixXd BS_VAR(constrained_points.size(), 2);
    Eigen::MatrixXd BS_TANGENT(constrained_points.size(), 4);
    for (size_t i = 0; i < constrained_points.size(); ++i) {
        BS_VAR(i, 0) = CGAL::to_double(constrained_points[i].x()); 
        BS_VAR(i, 1) = CGAL::to_double(constrained_points[i].y());
        BS(i) = constrained_indices[i] + points.size(); 

        Segment c = constraints_of_points[i];
        
        BS_TANGENT(i, 0) = CGAL::to_double(c.target().x());
        BS_TANGENT(i, 1) = CGAL::to_double(c.target().y());
        BS_TANGENT(i, 2) = CGAL::to_double(c.source().x());
        BS_TANGENT(i, 3) = CGAL::to_double(c.source().y());
    }

    find_minimum(problem, V, F, B, B_VAR, BS, BS_VAR, BS_TANGENT);
    update_problem_tinyAD(problem, V, F);

    //std::cout << "Num rows in V: " << V.rows() << " and in F: " << F.rows() << std::endl;
    /*flip_obtuse_triangles(problem, V, F, true);
    update_problem_tinyAD(problem, V, F);*/

    /*find_minimum(problem, V, F, B, B_VAR);
    update_problem_tinyAD(problem, V, F);*/
}

void find_minimum(Problem* problem, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& B, Eigen::MatrixXd& B_VAR, Eigen::VectorXi& BS, Eigen::MatrixXd& BS_VAR, Eigen::MatrixXd& BS_TANGENT){
    auto points = problem->get_points(); 
    auto steiner = problem->get_steiner(); 
    auto triangles = problem->get_triangulation();
    auto boundary = problem->get_boundary();
    auto constraints = problem->get_constraints();
    
    // Objective settings
    const double w_penalty = 1e5;
    const double epsilon = 1e-6;
    
    // Objective
    auto func = TinyAD::scalar_function<2>(TinyAD::range(V.rows()));

    // Add objective term per triangle
    func.add_elements<3>(TinyAD::range(F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element) {
        // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);

        Eigen::Index f_idx = element.handle;
        Eigen::Vector2<T> a = element.variables(F(f_idx, 0));
        Eigen::Vector2<T> b = element.variables(F(f_idx, 1));
        Eigen::Vector2<T> c = element.variables(F(f_idx, 2));

        // Triangle flipped?
        Eigen::Matrix2<T> M = TinyAD::col_mat(b - a, c - a);
        if (M.determinant() <= 0.0){ 
            //std::cout << F.row(f_idx) << std::endl;
            //std::cout << a << " | " << b << " | " << c << "\n" << std::endl;
            return (T)INFINITY;
        }

        return angle_cost_sigmoid<T>(problem, a, b, c);
    });

    // Add penalty term per constrained vertex
    func.add_elements<1>(TinyAD::range(B.size()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element) {
        using T = TINYAD_SCALAR_TYPE(element);
        Eigen::Vector2<T> p = element.variables(B[element.handle]);
        Eigen::Vector2d p_target = B_VAR.row(element.handle);
        return w_penalty * (p_target - p).squaredNorm();
    });

    func.add_elements<1>(TinyAD::range(BS.size()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element) {
        using T = TINYAD_SCALAR_TYPE(element);
        Eigen::Vector2<T> p = element.variables(BS[element.handle]);
        Eigen::Vector2d p_target = BS_VAR.row(element.handle);

        /*Eigen::Vector4d constraint = BS_TANGENT.row(element.handle);
        Eigen::Vector2d start = constraint.segment<2>(0);
        Eigen::Vector2d end = constraint.segment<2>(2);
        Eigen::Vector2d seg = end - start;
        Eigen::Vector2<T> pt = p - start;

        auto cross = seg.x() * pt.y() - seg.y() * pt.x();
        if (abs(cross) > epsilon){
            return w_penalty;
        } 

        if(!( p.x() >= std::min(start.x(), end.x())  &&  p.x() <= std::max(start.x(), end.x())  &&
              p.y() >= std::min(start.y(), end.y())  &&  p.y() <= std::max(start.y(), end.y()) )){

            return w_penalty;
        }

        return 0;*/
        
        return w_penalty * (p_target - p).squaredNorm();
        //return 0;
    });

    // Initialize x with the 2D vertex positions
    Eigen::VectorXd x = func.x_from_data([&] (int v_idx) {
        return V.row(v_idx);
    });

    /*std::cout << "BS:\n";
    std::cout << BS;
    std::cout << "\npoints size:\n";
    std::cout << points.size();*/

    TINYAD_DEBUG_OUT("Start energy: " << func.eval(x));

    // Projected Newton
    TinyAD::LinearSolver solver;
    const int max_iters = 1000;
    const double convergence_eps = 1e-6;
    for (int iter = 0; iter < max_iters; ++iter)
    {
        auto [f, g, H_proj] = func.eval_with_hessian_proj(x);

        int s = x.cols();
        int s1 = x.rows();
        int s2 = V.cols();
        int s3 = V.rows();

        Eigen::VectorXd d = TinyAD::newton_direction(g, H_proj, solver);

        double newton_decrement = TinyAD::newton_decrement<double>(d, g);
        //alter_direction(d, BS, BS_VAR, BS_TANGENT);
        //TINYAD_DEBUG_OUT("Energy | Newton decrement in iteration " << iter << ": " << f << " | " << newton_decrement);
        if(newton_decrement < convergence_eps)
            break;
        /*if(iter == 9){
            func.x_to_data(x, [&] (int v_idx, const Eigen::Vector2d& p) {
                if(v_idx >= points.size()){
                    Point start_p = steiner[v_idx - points.size()];
                    if(!is_constrained_point(start_p, boundary, constraints)){
                        V.row(v_idx) = p;
                    }
                }
            });
            update_problem_tinyAD(problem, V, F);
        }*/
        //x = TinyAD::line_search(x, d, f, g, func, 1.0, 0.8, 256);
        x = TinyAD::line_search(x, d, f, g, func);
    }
    TINYAD_DEBUG_OUT("Final energy: " << func.eval(x));

    // Write final vertex positions to mesh.
    func.x_to_data(x, [&] (int v_idx, const Eigen::Vector2d& p) {
        if(v_idx >= points.size()){
            Point start_p = steiner[v_idx - points.size()];
            Segment c;
            if(!is_constrained_point(start_p, boundary, constraints, &c)){
                V.row(v_idx) = p;
            }
        }
    });
}

void update_problem_tinyAD(Problem* problem, Eigen::MatrixXd& V, Eigen::MatrixXi& F){
    auto points = problem->get_points(); 
    auto steiner = problem->get_steiner(); 

    problem->clear_solution();
    // Update points and steiner
    for (size_t i = 0; i < points.size(); ++i) {
        points[i] = Point(V(i, 0), V(i, 1)); // Convert Eigen row to CGAL::Point_2
    }
    for (size_t i = 0; i < steiner.size(); ++i) {
        Point new_s = Point(V(points.size() + i, 0), V(points.size() + i, 1));
        steiner[i] = new_s;
        problem->add_steiner(new_s);
    }

    for (int i = 0; i < F.rows(); ++i) {
        Polygon poly;
        for (int j = 0; j < 3; ++j) {
            int idx = F(i, j);
            Point p(V(idx, 0), V(idx, 1)); // Get updated point from V
            poly.push_back(p);
        }
        problem->add_triangle(poly);
    }
}
