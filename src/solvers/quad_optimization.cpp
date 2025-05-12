#include "quad_optimization.hpp"

#include "../../libs/TinyAD/include/TinyAD/ScalarFunction.hh"
#include "../../libs/TinyAD/include/TinyAD/Utils/LineSearch.hh"
#include "../../libs/TinyAD/include/TinyAD/Utils/NewtonDirection.hh"
#include "../../libs/TinyAD/include/TinyAD/Utils/NewtonDecrement.hh"

#include <Eigen/Dense>

// Energy function: 1 / (1 + exp(-k(x - π/2)))
template <typename T>
T angle_cost_sigmoid_quad( Problem* problem,
        const Eigen::Vector2<T>& a,
        const Eigen::Vector2<T>& b,
        const Eigen::Vector2<T>& c,
        const Eigen::Vector2<T>& d)
{
    double k = 1;

    // Compute all three normalized triangle edge vectors
    Eigen::Vector2<T> ab = (b - a).normalized();
    Eigen::Vector2<T> bc = (c - b).normalized();
    Eigen::Vector2<T> cd = (d - c).normalized();
    Eigen::Vector2<T> da = (a - d).normalized();

    // Compute cosine of all 3 angles (unsigned)
    T cos_angle_a = (-da).dot(ab);
    T cos_angle_b = (-ab).dot(bc);
    T cos_angle_c = (-bc).dot(cd);
    T cos_angle_d = (-cd).dot(da);
   
    if(abs(cos_angle_a) >= 1){
        std::cout << "Cos out of range for Point a: " << a[0] << " " << a[1] << std::endl;
    }
    if(abs(cos_angle_b) >= 1){
        std::cout << "Cos out of range for Point b: " << b[0] << " " << b[1] << std::endl;
    }
    if(abs(cos_angle_c) >= 1){
        std::cout << "Cos out of range for Point c: " << c[0] << " " << c[1] << std::endl;
    }
    if(abs(cos_angle_d) >= 1){
        std::cout << "Cos out of range for Point c: " << c[0] << " " << c[1] << std::endl;
    }
        
    // Calculate and return the sigmoid term
    T angle_a = acos(cos_angle_a);
    T angle_b = acos(cos_angle_b);
    T angle_c = acos(cos_angle_c);
    T angle_d = acos(cos_angle_d);

    T sigmoid_a = 1/(1+exp(-k*(angle_a - M_PI/2)));
    T sigmoid_b = 1/(1+exp(-k*(angle_b - M_PI/2)));
    T sigmoid_c = 1/(1+exp(-k*(angle_c - M_PI/2)));
    T sigmoid_d = 1/(1+exp(-k*(angle_d - M_PI/2)));

    return sigmoid_a + sigmoid_b + sigmoid_c + sigmoid_d;
}

// Energy function: 1 / (1 + exp(-k(x - π/2)))
template <typename T>
T angle_cost_symmetric_quad( Problem* problem,
        const Eigen::Vector2<T>& a,
        const Eigen::Vector2<T>& b,
        const Eigen::Vector2<T>& c,
        const Eigen::Vector2<T>& d)
{
    double k = 1;

    // Compute all three normalized triangle edge vectors
    Eigen::Vector2<T> ab = (b - a).normalized();
    Eigen::Vector2<T> bc = (c - b).normalized();
    Eigen::Vector2<T> cd = (d - c).normalized();
    Eigen::Vector2<T> da = (a - d).normalized();

    // Compute cosine of all 3 angles (unsigned)
    T cos_angle_a = (-da).dot(ab);
    T cos_angle_b = (-ab).dot(bc);
    T cos_angle_c = (-bc).dot(cd);
    T cos_angle_d = (-cd).dot(da);
   
    if(abs(cos_angle_a) >= 1){
        std::cout << "Cos out of range for Point a: " << a[0] << " " << a[1] << std::endl;
    }
    if(abs(cos_angle_b) >= 1){
        std::cout << "Cos out of range for Point b: " << b[0] << " " << b[1] << std::endl;
    }
    if(abs(cos_angle_c) >= 1){
        std::cout << "Cos out of range for Point c: " << c[0] << " " << c[1] << std::endl;
    }
    if(abs(cos_angle_d) >= 1){
        std::cout << "Cos out of range for Point c: " << c[0] << " " << c[1] << std::endl;
    }
        
    // Calculate and return the sigmoid term
    T angle_a = acos(cos_angle_a);
    T angle_b = acos(cos_angle_b);
    T angle_c = acos(cos_angle_c);
    T angle_d = acos(cos_angle_d);

    T value_a = (angle_a - M_PI/2)*(angle_a - M_PI/2);
    T value_b = (angle_b - M_PI/2)*(angle_b - M_PI/2);
    T value_c = (angle_c - M_PI/2)*(angle_c - M_PI/2);
    T value_d = (angle_d - M_PI/2)*(angle_d - M_PI/2);

    T angle_energy = value_a + value_b + value_c + value_d;

    /*T diag_ac = (c - a).norm();
    T diag_bd = (d - b).norm();
    T length_energy = pow(diag_ac - diag_bd, 2);*/

    T length_energy = (ab.norm() - cd.norm())*(ab.norm() - cd.norm()) + (bc.norm() - da.norm())*(bc.norm() - da.norm()) + (ab.norm() - cd.norm())*(bc.norm() - da.norm());

    return length_energy + angle_energy;
}

bool is_on_constraint_quad(Point& steiner, std::vector<Segment>& constraints, Segment* constraint){
    for(Segment c : constraints){
        double distance = std::sqrt(CGAL::to_double(CGAL::squared_distance(c, steiner)));
        if(c.has_on(steiner) || distance < 1e-6){
            *constraint = c;
            return true;
        }
    }
    return false;
}

bool is_on_boundary_quad(Point& steiner, Polygon& boundary, Segment* constraint){
    int n = boundary.size();
    for (int i = 0; i < n; ++i) {
        Segment edge(boundary.vertex(i), boundary.vertex((i + 1) % n)); 
        double distance = std::sqrt(CGAL::to_double(CGAL::squared_distance(edge, steiner)));
        if (edge.has_on(steiner) || distance < 1e-6) {
            *constraint = edge;
            return true; 
        }
    }
    return false;
}

const bool is_constrained_point_quad(Point& p, Polygon& boundary, std::vector<Segment>& constraints, Segment* constraint) {
    return is_on_constraint_quad(p, constraints, constraint) || is_on_boundary_quad(p, boundary, constraint);
}

std::tuple<std::vector<Point>, std::vector<int>, std::vector<Segment>> get_constrained_points_quad(std::vector<Point>& steiner, Polygon& boundary, std::vector<Segment>& constraints){
    std::vector<Point> constrained_points;
    std::vector<Segment> constraints_of_points;
    std::vector<int> constrained_indices;
    for(int i = 0; i< steiner.size(); i++){
        Point s = steiner[i];
        Segment constraint;
        if(is_constrained_point_quad(s, boundary, constraints, &constraint)){
            constrained_points.push_back(s);
            constrained_indices.push_back(i);
            constraints_of_points.push_back(constraint);
        }
    }
    return {constrained_points, constrained_indices, constraints_of_points};
}

std::vector<Polygon> optimize_TinyAD_quad(Problem* problem, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::VectorXi B, Eigen::MatrixXd B_VAR){
    auto points = problem->get_points(); 
    auto boundary = problem->get_boundary();
    auto constraints = problem->get_constraints();

    // Get Steiner points from V:
    int num_steiner = V.rows() - points.size() - 1;
    Eigen::MatrixXd S = V.bottomRows(num_steiner);
    std::vector<Point> steiner;
    for(int i = 0; i < num_steiner; i++){
        Point s(S(i, 0), S(i, 1));
        steiner.push_back(s);
    }

    auto [constrained_points, constrained_indices, constraints_of_points] = get_constrained_points_quad(steiner, boundary, constraints);
    Eigen::VectorXi BS(constrained_points.size());
    Eigen::MatrixXd BS_VAR(constrained_points.size(), 2);
    for (size_t i = 0; i < constrained_points.size(); ++i) {
        BS_VAR(i, 0) = CGAL::to_double(constrained_points[i].x()); 
        BS_VAR(i, 1) = CGAL::to_double(constrained_points[i].y());
        BS(i) = constrained_indices[i] + points.size() + 1; 
    }

    find_minimum_quad(problem, V, F, B, B_VAR, BS, BS_VAR);
    std::vector<Polygon> quads = update_problem_quad(problem, V, F);

    return quads;
}

void find_minimum_quad(Problem* problem, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi& B, Eigen::MatrixXd& B_VAR, Eigen::VectorXi& BS, Eigen::MatrixXd& BS_VAR){
    auto points = problem->get_points(); 
    auto boundary = problem->get_boundary();
    auto constraints = problem->get_constraints();
    
    // Objective settings
    const double w_penalty = 1e5;
    const double epsilon = 1e-6;
    
    // Objective
    auto func = TinyAD::scalar_function<2>(TinyAD::range(V.rows()));

    // Add objective term per triangle
    func.add_elements<4>(TinyAD::range(F.rows()), [&] (auto& element) -> TINYAD_SCALAR_TYPE(element) {
        // Evaluate element using either double or TinyAD::Double
        using T = TINYAD_SCALAR_TYPE(element);

        Eigen::Index f_idx = element.handle;
        auto i0 = F(f_idx, 0);
        auto i1 = F(f_idx, 1);
        auto i2 = F(f_idx, 2);
        auto i3 = F(f_idx, 3);
        Eigen::Vector2<T> a = element.variables(F(f_idx, 0));
        Eigen::Vector2<T> b = element.variables(F(f_idx, 1));
        Eigen::Vector2<T> c = element.variables(F(f_idx, 2));
        Eigen::Vector2<T> d = element.variables(F(f_idx, 3));

        // Triangle flipped?
        /*Eigen::Matrix2<T> M = TinyAD::col_mat(b - a, c - a);
        if (M.determinant() <= 0.0){ 
            //std::cout << F.row(f_idx) << std::endl;
            //std::cout << a << " | " << b << " | " << c << "\n" << std::endl;
            return (T)INFINITY;
        }*/

        Eigen::Matrix2<T> M1 = TinyAD::col_mat(b - a, c - a);
        Eigen::Matrix2<T> M2 = TinyAD::col_mat(c - a, d - a);
        if (M1.determinant() <= 0.0 || M2.determinant() <= 0.0){
            return (T)INFINITY;
        }

        return angle_cost_symmetric_quad<T>(problem, a, b, c, d);
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
        return w_penalty * (p_target - p).squaredNorm();
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
    const int max_iters = 10;
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
        //x = TinyAD::line_search(x, d, f, g, func, 1.0, 0.8, 256);
        x = TinyAD::line_search(x, d, f, g, func);
    }
    TINYAD_DEBUG_OUT("Final energy: " << func.eval(x));

    // Write final vertex positions to mesh.
    func.x_to_data(x, [&] (int v_idx, const Eigen::Vector2d& p) {
        if(v_idx > points.size()){
            Point start_p(V(v_idx, 0), V(v_idx, 1));
            Segment c;
            if(!is_constrained_point_quad(start_p, boundary, constraints, &c)){
                V.row(v_idx) = p;
            }
        }
    });
}

std::vector<Polygon> update_problem_quad(Problem* problem, Eigen::MatrixXd& V, Eigen::MatrixXi& F){
    std::vector<Polygon> quads;

    auto points = problem->get_points(); 
    auto steiner = problem->get_steiner(); 

    // Update points and steiner
    for (size_t i = 0; i < steiner.size(); ++i) {
        Point new_s = Point(V(points.size() + i, 0), V(points.size() + i, 1));
        steiner[i] = new_s;
        problem->add_steiner(new_s);
    }

    for (int i = 0; i < F.rows(); ++i) {
        Polygon poly;
        for (int j = 0; j < 4; ++j) {
            int idx = F(i, j);
            Point p(V(idx, 0), V(idx, 1)); // Get updated point from V
            poly.push_back(p);
        }
        quads.push_back(poly);
    }

    return quads;
}