#include "main.hpp"

namespace fs = std::filesystem;

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

void write_to_csv_fully_non_obtuse(std::string filepath, std::vector<Mesh_Statistics> all_stats) {
    std::ofstream out(filepath);
    if (!out) throw std::runtime_error("Could not open file for writing: " + filepath);

    out << "name,#obtuse,#steiner,#steiner_opt,time\n";
    for (Mesh_Statistics stats : all_stats) {
        out << stats.get_name() << ","
            << stats.get_obtuse_after_meshing() << ","
            << stats.get_steiner_after_meshing() << ","
            << 0 << ","
            << stats.get_time() << "\n";
    }
}

void print_current_time(){
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);

    std::tm *local_time = std::localtime(&now_c);

    std::cout << "Current time: "
              << std::put_time(local_time, "%Y-%m-%d %H:%M:%S")
              << std::endl;
}

int main(int argc, char **argv){
    print_current_time();

    Problem *problem = new Problem("../instances/point-set/point-set_10_d009159f.instance.json");
    std::cout << "current instance: " << problem->get_name() << std::endl;
    problem->visualize_solution({});

    //* Quad Mesh
    /*Problem *problem = new Problem("../instances/simple-exterior/simple-polygon-exterior-20_60_53ad6d23.instance.json");
    uniform_mesh(problem);
    build_quad_mesh_medians(problem);*/

    //* Fully non-obtuse mesh
    /*print_current_time();
    std::vector<Mesh_Statistics> all_stats;

    //Problem *problem = new Problem("../instances/ortho/ortho_10_d2723dcc.instance.json");
    //uniform_mesh(problem);

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

            std::string name = problem->get_name();
            
                std::cout << BLUE << "\ncurrent instance: " << problem->get_name() << RESET << std::endl;
                //Mesh_Statistics stats  = uniform_mesh(problem);
                auto start = std::chrono::high_resolution_clock::now();
                Mesh_Statistics stats = offcenter_delaunay_refinement(problem);
                //Mesh_Statistics stats;
                //classic_delaunay_refinement(problem);
                auto end = std::chrono::high_resolution_clock::now();

                std::chrono::duration<double> duration = end - start;
                auto total_seconds = duration.count();
                int minutes = total_seconds/60;
                int seconds = total_seconds - minutes*60;
                std::cout << GREEN << "Execution time: " << minutes << " minutes " << seconds << " seconds" << RESET << std::endl;

                stats.set_time(total_seconds);

                //all_stats.push_back(stats);
                //write_to_csv_fully_non_obtuse("../results/fully_nonobtuse_ortho_offcenter.csv", all_stats);
        }
    }

    write_to_csv_fully_non_obtuse("../results/fully_nonobtuse_ortho_offcenter.csv", all_stats);*/


    //* Equilateral mesh
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

    
    //* Uniform non-obtuse
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
            //write_to_csv_obtuse("../results/ortho_uniform_tinyAD.csv", all_stats);
        }
    }

    write_to_csv_obtuse("../results/ortho_uniform_tinyAD.csv", all_stats);
    print_current_time();*/

    //* Local and global optimization with gradient descent
    /*std::vector<Mesh_Statistics> all_stats;

    for (const auto& entry : fs::directory_iterator("../instances/ortho")) {
        if (entry.is_regular_file() && entry.path().extension() == ".json") {
            Problem *problem = new Problem("../instances/ortho/" + entry.path().filename().string());
            
            Polygon boundary = problem->get_boundary();
            auto bbox = boundary.bbox();
            double max_size = std::sqrt(pow(bbox.x_span(),2) + pow(bbox.y_span(),2)) * 0.1;

            std::vector<Point> points = problem->get_points();
            std::set<Point> point_set(points.begin(), points.end());
            CDT cdt = problem->generate_CDT();
            problem->update_problem(cdt, point_set);
            problem->generate_uniform_mesh(cdt, max_size);

            Mesh_Statistics stats;
            stats.set_name(problem->get_name());
            stats.set_steiner_after_meshing(problem->get_steiner().size());
            stats.set_obtuse_after_meshing(count_obtuse_triangles(problem));
            
            std::cout << "Num obtuse before: " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});

            std::cout << "num_steiner: " << problem->get_steiner().size() << std::endl;

            locally_optimize_solution(problem);
            //globally_optimize_solution(problem);
            std::cout << "Num obtuse after optimization 1: " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});

            locally_optimize_obtuse(problem);
            //globally_optimize_obtuse(problem);
            std::cout << "Num obtuse after optimization 2: " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});

            locally_optimize_obtuse(problem);
            //globally_optimize_obtuse(problem);
            std::cout << "Num obtuse after optimization 3: " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});

            perform_edge_flips(problem, false);
            std::cout << "Num obtuse after flips 1: " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});

            fix_boundary(problem);
            std::cout << "Num obtuse after fixing boundary 1: " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});

            locally_optimize_obtuse(problem);
            //globally_optimize_obtuse(problem);
            std::cout << "Num obtuse after optimization 4: " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});

            fix_boundary(problem);
            std::cout << "Num obtuse after fixing boundary 2: " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});

            perform_edge_flips(problem, false);
            std::cout << "Num obtuse after flips 2: " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});

            locally_optimize_obtuse(problem);
            //globally_optimize_obtuse(problem);
            std::cout << "Num obtuse after optimization 5: " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});

            fix_boundary(problem);
            std::cout << "Num obtuse after fixing boundary 3: " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});
            //std::cout << "new_steiner : " << problem->get_steiner().size() << std::endl;

            locally_optimize_obtuse(problem);
            //globally_optimize_obtuse(problem);
            std::cout << "Num obtuse after optimization 6: " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});

            fix_boundary(problem);
            std::cout << "Num obtuse after fixing boundary 4: " << count_obtuse_triangles(problem) << std::endl;
            problem->visualize_solution({});

            std::cout << "new_steiner : " << problem->get_steiner().size() << "\n\n" << std::endl;

            stats.set_steiner_after_optimization(problem->get_steiner().size());
            stats.set_obtuse_after_optimization(count_obtuse_triangles(problem));

            all_stats.push_back(stats);
            write_to_csv_obtuse("../results/ortho_uniform_step.csv", all_stats);
        }
    }

    write_to_csv_obtuse("../results/ortho_uniform_step.csv", all_stats);*/

    return 0;
}
