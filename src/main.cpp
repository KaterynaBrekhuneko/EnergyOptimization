#include "main.hpp"
#include "local_optimization.hpp"
#include "global_optimization.hpp"
#include "tinyAD_optimization.hpp"
#include "delaunay_refinement.hpp"
#include "quad_mesher.hpp"

#include <gmsh.h>

#include <chrono>
#include <ctime>

namespace fs = std::filesystem;

int num_new_steiner = 0;

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

std::pair<std::vector<Polygon>, std::vector<Polygon>> flip(Polygon& triangle, std::vector<Polygon>& triangulation, bool ignore){
    Segment longestEdge = find_longest_edge(triangle);

    std::vector<Polygon> removed_triangles;
    std::vector<Polygon> new_triangles;

    for (const auto& t : triangulation) {
        if (contains_edge(t, longestEdge.source(), longestEdge.target())) {
            removed_triangles.push_back(t);
        }
    }

    if(removed_triangles.size() == 2){
        /*triangulation.erase(std::remove_if(triangulation.begin(), triangulation.end(),
        [&](const Polygon& tri) { return contains_edge(tri, longestEdge.source(), longestEdge.target()); }),
        triangulation.end());*/

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
            new_triangles.push_back(t1);
            new_triangles.push_back(t2);
            return {removed_triangles, new_triangles};
        } else {
            //new_triangles.push_back(removed_triangles[0]);
            //new_triangles.push_back(removed_triangles[1]);
            removed_triangles.clear();
            new_triangles.clear();
            return {removed_triangles, new_triangles};
        }
    }

    //removed_triangles.clear();
    return {removed_triangles, new_triangles};
}

void perform_edge_flips(Problem *problem, bool ignore){
    std::vector<Polygon> triangulation = problem->get_triangulation();

    std::vector<Polygon> removed_triangles;
    std::vector<Polygon> new_triangles;

    std::vector<Polygon> to_flip = triangulation;
    for(Polygon& triangle : to_flip){
        if(is_obtuse_triangle(triangle)){
            auto [current_remove, current_new] = flip(triangle, triangulation, ignore);

            //add and remove if has not been done before
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
    
    /*problem->clear_solution();
    for(int i = 0; i < new_steiner.size(); i++){
        problem->add_steiner(new_steiner[i]);
    }*/
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
                    num_new_steiner++;
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

void write_to_csv_obtuse(std::string filepath, std::vector<Mesh_Statistics> all_stats) {
    std::ofstream out(filepath);
    if (!out) throw std::runtime_error("Could not open file for writing: " + filepath);

    out << "instance,steiner_meshing,steiner_optimized,obtuse_meshing,obtuse_optimized\n";
    for (Mesh_Statistics stats : all_stats) {
        out << stats.get_name() << ","
            << stats.get_steiner_after_meshing() << ","
            << stats.get_steiner_after_optimization() << ","
            << stats.get_obtuse_after_meshing() << ","
            << stats.get_obtuse_after_optimization() << "\n";
    }
}

void write_to_csv_equilateral(std::string filepath, std::vector<Mesh_Statistics> all_stats) {
    std::ofstream out(filepath);
    if (!out) throw std::runtime_error("Could not open file for writing: " + filepath);

    out << "name,deviation,aspect_ratio,min_angle,max_angle,num_steiner,max_length\n";
    for (Mesh_Statistics stats : all_stats) {
        out << stats.get_name() << ","
            << stats.get_deviation() << ","
            << stats.get_aspect_ratio() << ","
            << stats.get_min_angle() << ","
            << stats.get_max_angle() << ","
            << stats.get_steiner_after_meshing() <<","
            << stats.get_max_length() <<"\n";
    }
}

void print_current_time(){
    auto now = std::chrono::system_clock::now();              // get current time_point
    std::time_t now_c = std::chrono::system_clock::to_time_t(now); // convert to time_t

    std::tm* local_time = std::localtime(&now_c);             // convert to local time

    std::cout << "Current time: "
              << std::put_time(local_time, "%Y-%m-%d %H:%M:%S")
              << std::endl;
}

int main(int argc, char **argv){
    //* Calculate mean angle deviation for equilateral meshes
    print_current_time();
    std::vector<Mesh_Statistics> all_stats;

    //Problem *problem = new Problem("../instances/ortho/ortho_10_d2723dcc.instance.json");
    //uniform_mesh(problem);

    std::vector<fs::directory_entry> entries;
    for (const auto& entry : fs::directory_iterator("../instances/point-set-current")) {
        entries.push_back(entry);
    }
    std::sort(entries.begin(), entries.end(),
              [](const fs::directory_entry& a, const fs::directory_entry& b) {
                  return a.path().filename() < b.path().filename(); // alphabetical by name
              });
          

    /*for (const auto& entry : entries) {
        if (entry.is_regular_file() && entry.path().extension() == ".json") {
            Problem *problem = new Problem("../instances/point-set-current/" + entry.path().filename().string());*/

            Problem *problem = new Problem("../instances/simple-exterior/simple-polygon-exterior_10_40642b31.instance.json");
            std::string name = problem->get_name();
            //if(!(name == "simple-polygon_100_4b4ba391" || name == "simple-polygon_10_297edd18" || name == "simple-polygon_20_4bd3c2e5")){
            //if(!(name == "point-set_10_c04b0024" || name == "point-set_20_54ab0b47" || name == "point-set_40_1b92b629" || name == "point-set_60_ac318d72" || name == "point-set_80_1675b331")){
                std::cout << BLUE << "\ncurrent instance: " << problem->get_name() << RESET << std::endl;
                //Mesh_Statistics stats  = uniform_mesh(problem);
                auto start = std::chrono::high_resolution_clock::now();
                classic_delaunay_refinement(problem);
                auto end = std::chrono::high_resolution_clock::now();

                std::chrono::duration<double> duration = end - start;
                auto total_seconds = duration.count();
                int minutes = total_seconds/60;
                int seconds = total_seconds - minutes*60;
                std::cout << GREEN << "Execution time: " << minutes << " minutes " << seconds << " seconds" << RESET << std::endl;

                /*all_stats.push_back(stats);
                write_to_csv_obtuse("../results/modified3_initial_ln_simple_exterior.csv", all_stats);*/ 
            //}
        //}
    //}

    //write_to_csv_obtuse("../results/modified3_initial_ln_simple_exterior.csv", all_stats);

    //print_current_time();

    /*all_stats.clear();
    entries.clear();
    for (const auto& entry : fs::directory_iterator("../instances/ortho")) {
        entries.push_back(entry);
    }
    std::sort(entries.begin(), entries.end(),
              [](const fs::directory_entry& a, const fs::directory_entry& b) {
                  return a.path().filename() < b.path().filename(); // alphabetical by name
              });
          

    for (const auto& entry : entries) {
        if (entry.is_regular_file() && entry.path().extension() == ".json") {
            Problem *problem = new Problem("../instances/ortho/" + entry.path().filename().string());

            //Problem *problem = new Problem("../instances/point-set/point-set_10_7451a2a9.instance.json");
            std::string name = problem->get_name();
            std::cout << "current instance: " << problem->get_name() << std::endl;
            Mesh_Statistics stats  = uniform_mesh(problem);

            all_stats.push_back(stats);
            write_to_csv_obtuse("../results/modified3_initial_ln_ortho.csv", all_stats); 
        }
    }

    write_to_csv_obtuse("../results/modified3_initial_ln_ortho.csv", all_stats);

    print_current_time();

    all_stats.clear();
    entries.clear();
    for (const auto& entry : fs::directory_iterator("../instances/point-set")) {
        entries.push_back(entry);
    }
    std::sort(entries.begin(), entries.end(),
              [](const fs::directory_entry& a, const fs::directory_entry& b) {
                  return a.path().filename() < b.path().filename(); // alphabetical by name
              });
          

    for (const auto& entry : entries) {
        if (entry.is_regular_file() && entry.path().extension() == ".json") {
            Problem *problem = new Problem("../instances/point-set/" + entry.path().filename().string());

            //Problem *problem = new Problem("../instances/point-set/point-set_10_7451a2a9.instance.json");
            std::string name = problem->get_name();
            std::cout << "current instance: " << problem->get_name() << std::endl;
            Mesh_Statistics stats  = uniform_mesh(problem);

            all_stats.push_back(stats);
            write_to_csv_obtuse("../results/modified3_initial_ln_point_set.csv", all_stats); 
        }
    }

    write_to_csv_obtuse("../results/modified3_initial_ln_point_set.csv", all_stats);

    print_current_time();

    all_stats.clear();
    entries.clear();
    for (const auto& entry : fs::directory_iterator("../instances/simple")) {
        entries.push_back(entry);
    }
    std::sort(entries.begin(), entries.end(),
              [](const fs::directory_entry& a, const fs::directory_entry& b) {
                  return a.path().filename() < b.path().filename(); // alphabetical by name
              });
          

    for (const auto& entry : entries) {
        if (entry.is_regular_file() && entry.path().extension() == ".json") {
            Problem *problem = new Problem("../instances/simple/" + entry.path().filename().string());

            //Problem *problem = new Problem("../instances/point-set/point-set_10_7451a2a9.instance.json");
            std::string name = problem->get_name();
            std::cout << "current instance: " << problem->get_name() << std::endl;
            Mesh_Statistics stats  = uniform_mesh(problem);

            all_stats.push_back(stats);
            write_to_csv_obtuse("../results/modified3_initial_ln_simple.csv", all_stats); 
        }
    }

    write_to_csv_obtuse("../results/modified3_initial_ln_simple.csv", all_stats);

    print_current_time();*/

    /*locally_optimize_solution(problem);
    //globally_optimize_solution(problem);
    std::cout << "Num obtuse after optimization 1: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution({});

    locally_optimize_obtuse(problem);
    //globally_optimize_obtuse(problem);
    std::cout << "Num obtuse after optimization 2: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution({});

    locally_optimize_obtuse(problem);
    //globally_optimize_obtuse(problem);
    std::cout << "Num obtuse after optimization 3: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution({});

    perform_edge_flips(problem, false);
    std::cout << "Num obtuse after flips 1: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution({});

    fix_boundary(problem);
    std::cout << "Num obtuse after fixing boundary 1: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution({});

    locally_optimize_obtuse(problem);
    //globally_optimize_obtuse(problem);
    std::cout << "Num obtuse after optimization 4: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution({});

    fix_boundary(problem);
    std::cout << "Num obtuse after fixing boundary 2: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution({});

    perform_edge_flips(problem, false);
    std::cout << "Num obtuse after flips 2: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution({});

    locally_optimize_obtuse(problem);
    //globally_optimize_obtuse(problem);
    std::cout << "Num obtuse after optimization 5: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution({});

    fix_boundary(problem);
    std::cout << "Num obtuse after fixing boundary 3: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution({});
    //std::cout << "new_steiner : " << problem->get_steiner().size() << std::endl;

    locally_optimize_obtuse(problem);
    //globally_optimize_obtuse(problem);
    std::cout << "Num obtuse after optimization 6: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution({});

    fix_boundary(problem);
    std::cout << "Num obtuse after fixing boundary 4: " << countObtuse(problem->get_triangulation()) << std::endl;
    problem->visualize_solution({});

    std::cout << "new_steiner : " << problem->get_steiner().size() << "\n\n" << std::endl;

    // Sort instance files in the directory in alphabetical order
    /*std::vector<fs::directory_entry> entries;
    for (const auto& entry : fs::directory_iterator("../instances/simple")) {
        entries.push_back(entry);
    }
    std::sort(entries.begin(), entries.end(),
              [](const fs::directory_entry& a, const fs::directory_entry& b) {
                  return a.path().filename() < b.path().filename(); // alphabetical by name
              });
          

    for (const auto& entry : entries) {
        if (entry.is_regular_file() && entry.path().extension() == ".json") {
            Problem *problem = new Problem("../instances/simple/" + entry.path().filename().string());

            //Problem *problem = new Problem("../instances/simple/simple-polygon_20_4bd3c2e5.instance.json");
            std::cout << "current instance: " << problem->get_name() << std::endl;

            Mesh_Statistics stats  = refine(problem);

            //Mesh_Statistics stats  = mesh_equilateral(problem);
            /*all_stats.push_back(stats);
            write_to_csv_equilateral("../results/mean_deviations_initial_ortho.csv", all_stats);
        }
    }

    write_to_csv_equilateral("../results/mean_deviations_initial_ortho.csv", all_stats);

    print_current_time();*/

    /*all_stats.clear();
    // Sort instance files in the directory in alphabetical order
    entries.clear();
    for (const auto& entry : fs::directory_iterator("../instances/simple")) {
        entries.push_back(entry);
    }
    std::sort(entries.begin(), entries.end(),
              [](const fs::directory_entry& a, const fs::directory_entry& b) {
                  return a.path().filename() < b.path().filename(); // alphabetical by name
              });
          

    for (const auto& entry : entries) {
        if (entry.is_regular_file() && entry.path().extension() == ".json") {
            Problem *problem = new Problem("../instances/simple/" + entry.path().filename().string());

            std::cout << "current instance: " << problem->get_name() << std::endl;

            Mesh_Statistics stats  = refine(problem);

            //Mesh_Statistics stats  = mesh_equilateral(problem);
            all_stats.push_back(stats);
            write_to_csv_equilateral("../results/mean_deviations_initial_simple.csv", all_stats);
        }
    }

    write_to_csv_equilateral("../results/mean_deviations_initial_simple.csv", all_stats);

    print_current_time();

    all_stats.clear();
    // Sort instance files in the directory in alphabetical order
    entries.clear();
    for (const auto& entry : fs::directory_iterator("../instances/point-set")) {
        entries.push_back(entry);
    }
    std::sort(entries.begin(), entries.end(),
              [](const fs::directory_entry& a, const fs::directory_entry& b) {
                  return a.path().filename() < b.path().filename(); // alphabetical by name
              });
          

    for (const auto& entry : entries) {
        if (entry.is_regular_file() && entry.path().extension() == ".json") {
            Problem *problem = new Problem("../instances/point-set/" + entry.path().filename().string());

            std::cout << "current instance: " << problem->get_name() << std::endl;

            Mesh_Statistics stats  = refine(problem);

            //Mesh_Statistics stats  = mesh_equilateral(problem);
            all_stats.push_back(stats);
            write_to_csv_equilateral("../results/mean_deviations_initial_point_set.csv", all_stats);
        }
    }

    write_to_csv_equilateral("../results/mean_deviations_initial_point_set.csv", all_stats);

    print_current_time();

    all_stats.clear();
    // Sort instance files in the directory in alphabetical order
    entries.clear();
    for (const auto& entry : fs::directory_iterator("../instances/simple-exterior")) {
        entries.push_back(entry);
    }
    std::sort(entries.begin(), entries.end(),
              [](const fs::directory_entry& a, const fs::directory_entry& b) {
                  return a.path().filename() < b.path().filename(); // alphabetical by name
              });
          

    for (const auto& entry : entries) {
        if (entry.is_regular_file() && entry.path().extension() == ".json") {
            Problem *problem = new Problem("../instances/simple-exterior/" + entry.path().filename().string());

            std::cout << "current instance: " << problem->get_name() << std::endl;

            Mesh_Statistics stats  = refine(problem);

            //Mesh_Statistics stats  = mesh_equilateral(problem);
            all_stats.push_back(stats);
            write_to_csv_equilateral("../results/mean_deviations_initial_simple_exterior.csv", all_stats);
        }
    }

    write_to_csv_equilateral("../results/mean_deviations_initial_simple_exterior.csv", all_stats);

    print_current_time();*/


    //* Calculate angle stats for individual equilateral meshes
    /*Problem *problem = new Problem("../instances/point-set/point-set_80_837b0f11.instance.json");

    mesh_equilateral_single(problem);*/

    //Problem *problem = new Problem("../instances/challenge_instances_cgshop25/ortho_60_f744490d.instance.json");
    //problem->visualize_solution();
    //Problem *problem = new Problem(argv[1]);

    
    //* Test all instances
    /*print_current_time();
    std::vector<Mesh_Statistics> all_stats;

    // Sort instance files in the directory in alphabetical order
    std::vector<fs::directory_entry> entries;
    for (const auto& entry : fs::directory_iterator("../instances/ortho")) {
        entries.push_back(entry);
    }
    std::sort(entries.begin(), entries.end(),
              [](const fs::directory_entry& a, const fs::directory_entry& b) {
                  return a.path().filename() < b.path().filename(); // alphabetical by name
              });

    for (const auto& entry : entries) {
        if (entry.is_regular_file() && entry.path().extension() == ".json") {
            Problem *problem = new Problem("../instances/ortho/" + entry.path().filename().string());

            //Problem *problem = new Problem("../instances/ortho/ortho_60_f744490d.instance.json");
            std::cout << "current instance: " << problem->get_name() << std::endl;

            Mesh_Statistics stats = uniform_mesh(problem);
            all_stats.push_back(stats);
            write_to_csv_obtuse("../results/ortho_uniform_tinyAD.csv", all_stats);
        }
    }

    write_to_csv_obtuse("../results/ortho_uniform_tinyAD.csv", all_stats);
    print_current_time();*/

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

   /*std::vector<Mesh_Statistics> all_stats;

    for (const auto& entry : fs::directory_iterator("../instances/ortho")) {
        if (entry.is_regular_file() && entry.path().extension() == ".json") {
            Problem *problem = new Problem("../instances/ortho/" + entry.path().filename().string());*/
            //Problem *problem = new Problem("../instances/ortho/ortho_60_c423f527.instance.json");
            
            /*Solver* solver = new Mesh();
            solver->solve(problem);
            std::cout << "Instance name: " << problem->get_name() << std::endl;*/

            //step_by_step_mesh(problem);

            /*Mesh_Statistics stats;
            stats.set_name(problem->get_name());
            stats.set_steiner_after_meshing(problem->get_steiner().size());
            stats.set_obtuse_after_meshing(countObtuse(problem->get_triangulation()));*/
            
            /*std::cout << "Num obtuse before: " << countObtuse(problem->get_triangulation()) << std::endl;
            problem->visualize_solution();

            std::cout << "num_steiner: " << problem->get_steiner().size() << std::endl;

            locally_optimize_solution(problem);
            //globally_optimize_solution(problem);
            std::cout << "Num obtuse after optimization 1: " << countObtuse(problem->get_triangulation()) << std::endl;
            problem->visualize_solution();

            locally_optimize_obtuse(problem);
            //globally_optimize_obtuse(problem);
            std::cout << "Num obtuse after optimization 2: " << countObtuse(problem->get_triangulation()) << std::endl;
            problem->visualize_solution();

            locally_optimize_obtuse(problem);
            //globally_optimize_obtuse(problem);
            std::cout << "Num obtuse after optimization 3: " << countObtuse(problem->get_triangulation()) << std::endl;
            problem->visualize_solution();

            perform_edge_flips(problem, false);
            std::cout << "Num obtuse after flips 1: " << countObtuse(problem->get_triangulation()) << std::endl;
            problem->visualize_solution();

            fix_boundary(problem);
            std::cout << "Num obtuse after fixing boundary 1: " << countObtuse(problem->get_triangulation()) << std::endl;
            problem->visualize_solution();

            locally_optimize_obtuse(problem);
            //globally_optimize_obtuse(problem);
            std::cout << "Num obtuse after optimization 4: " << countObtuse(problem->get_triangulation()) << std::endl;
            problem->visualize_solution();

            fix_boundary(problem);
            std::cout << "Num obtuse after fixing boundary 2: " << countObtuse(problem->get_triangulation()) << std::endl;
            problem->visualize_solution();

            perform_edge_flips(problem, false);
            std::cout << "Num obtuse after flips 2: " << countObtuse(problem->get_triangulation()) << std::endl;
            problem->visualize_solution();

            locally_optimize_obtuse(problem);
            //globally_optimize_obtuse(problem);
            std::cout << "Num obtuse after optimization 5: " << countObtuse(problem->get_triangulation()) << std::endl;
            problem->visualize_solution();

            fix_boundary(problem);
            std::cout << "Num obtuse after fixing boundary 3: " << countObtuse(problem->get_triangulation()) << std::endl;
            problem->visualize_solution();
            //std::cout << "new_steiner : " << problem->get_steiner().size() << std::endl;

            locally_optimize_obtuse(problem);
            //globally_optimize_obtuse(problem);
            std::cout << "Num obtuse after optimization 6: " << countObtuse(problem->get_triangulation()) << std::endl;
            problem->visualize_solution();

            fix_boundary(problem);
            std::cout << "Num obtuse after fixing boundary 4: " << countObtuse(problem->get_triangulation()) << std::endl;
            problem->visualize_solution();

            std::cout << "new_steiner : " << problem->get_steiner().size() << "\n\n" << std::endl;*/

            /*stats.set_steiner_after_optimization(problem->get_steiner().size());
            stats.set_obtuse_after_optimization(countObtuse(problem->get_triangulation()));

            all_stats.push_back(stats);
            write_to_csv_obtuse("../results/ortho_uniform_step.csv", all_stats);
        }
    }

    write_to_csv_obtuse("../results/ortho_uniform_step.csv", all_stats);*/


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
