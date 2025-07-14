#include "problem.hpp"
#include <filesystem>

Problem::Problem(std::string file_name)
//Problem::Problem(char *file_name)
{
    std::ifstream f(file_name);

    if (f.fail()){
        std::cerr << "file not found!\n";
        return;
    }

    json data = json::parse(f);

    name = data["instance_uid"];

    int num_points = data["num_points"];
    points.reserve(num_points);

    for (int i = 0; i < num_points; i++){
        int x = data["points_x"][i];
        int y = data["points_y"][i];
        points.push_back(Point(x, y));
    }

    for (int index : data["region_boundary"]){
        boundary_indices.push_back(index);
        boundary.push_back(points[index]);
    }

    int num_constraints = data["num_constraints"];
    constraints.reserve(num_constraints);

    for (int i = 0; i < num_constraints; i++){
        auto constraint = data["additional_constraints"][i];
        constraints_indices.push_back(constraint);
        constraints.push_back(Segment(points[constraint[0]], points[constraint[1]]));
    }

    for (int i = 0; i < data["steiner_x"].size(); i++){
        int x = data["steiner_x"][i];
        int y = data["steiner_y"][i];
        steiner.push_back(Point(x, y));
    }

    num_obtuse = -1;
}

void Problem::write_problem_to_json(const std::string& output_file) {
    json data;

    // Basic fields
    data["instance_uid"] = name;

    // Points
    data["num_points"] = points.size();
    data["points_x"] = json::array();
    data["points_y"] = json::array();
    for (const auto& p : points) {
        data["points_x"].emplace_back(CGAL::to_double(p.x()));
        data["points_y"].emplace_back(CGAL::to_double(p.y()));
    }

    // Boundary indices
    data["region_boundary"] = boundary_indices;

    // Constraints
    data["num_constraints"] = constraints.size();
    data["additional_constraints"] = json::array();
    for (const auto& constraint : constraints_indices) {
        data["additional_constraints"].push_back(constraint);
    }

    // Steiner points
    data["steiner_x"] = json::array();
    data["steiner_y"] = json::array();
    for (const auto& s : steiner) {
        data["steiner_x"].emplace_back(CGAL::to_double(s.x()));
        data["steiner_y"].emplace_back(CGAL::to_double(s.y()));
    }

    // Write to file
    std::ofstream out(output_file);
    if (!out) {
        std::cerr << "Could not open output file for writing.\n";
        return;
    }
    out << std::setw(2) << data << std::endl;
}

void Problem::remove_triangle(Polygon t){
    triangulation.erase(
        std::remove_if(triangulation.begin(), triangulation.end(),
            [&t](const Polygon& poly) { return poly == t; }),
        triangulation.end()
    );
}

void Problem::update_triangulation(Point s, Point new_s){
    for(auto& triangle : triangulation){
        for (size_t i = 0; i < triangle.size(); ++i) {
            if (triangle.vertex(i) == s) {  
                triangle[i] = new_s;  
            }
        }
    }

    for(int i = 0; i< steiner.size(); i++){
        Point steiner_point = steiner[i];
        if(steiner_point == s){
            steiner[i] = new_s;
        }
    }
}

void Problem::update_steiner(Point s, Point new_s){
    // If the new_s is already present just remove s otherwise replace s with new_s
    bool is_present = false;
    for(int i = 0; i< steiner.size(); i++){
        Point steiner_point = steiner[i];
        if(steiner_point == new_s){
            is_present = true;
            break;
        }
    }
    for(int i = 0; i< points.size(); i++){
        Point point = points[i];
        if(point == new_s){
            is_present = true;
            break;
        }
    }

    if(is_present){
        steiner.erase(std::remove(steiner.begin(), steiner.end(), s), steiner.end());
    } else {
        for(int i = 0; i< steiner.size(); i++){
            Point steiner_point = steiner[i];
            if(steiner_point == s){
                steiner[i] = new_s;
                break;
            }
        }

        //std::cout << "sss" << std::endl;
    }
}

// TODO: Test & maybe delete
Point parse_to_point(const std::string& x_str, const std::string& y_str) {
    auto parse_rational = [](const std::string& str) -> K::FT {
        std::size_t pos = str.find('/');
        if (pos == std::string::npos) {
            throw std::invalid_argument("String does not represent a rational number");
        }
        std::string numerator_str = str.substr(0, pos);
        std::string denominator_str = str.substr(pos + 1);
        double num_d = std::stod(numerator_str);
        double denom_d = std::stod(denominator_str);
        K::FT num(num_d);
        K::FT den(denom_d);

        return num / den;
    };

    K::FT x = parse_rational(x_str);
    K::FT y = parse_rational(y_str);

    return Point(x, y);
}

// TODO do i need this?
void Problem::load_solution()
{
    std::string file_name = "../solutions/" + name + ".solution.json";
    std::ifstream f(file_name);

    if (f.fail()){
        std::cerr << "file not found!\n";
        return;
    }

    json data = json::parse(f);

    int num_steiner = data["steiner_points_x"].size();

    for (int i = 0; i < num_steiner; i++) {
        std::string x = data["steiner_points_x"][i];
        std::string y = data["steiner_points_y"][i];

        steiner.push_back(parse_to_point(x, y));
    }
}

int find_obtuse_angle(Polygon& triangle){
    auto a2 = CGAL::squared_distance(triangle.vertex(0), triangle.vertex(1)); // side opposite to vertex (x1, y1)
    auto b2 = CGAL::squared_distance(triangle.vertex(1), triangle.vertex(2)); // side opposite to vertex (x2, y2)
    auto c2 = CGAL::squared_distance(triangle.vertex(0), triangle.vertex(2)); // side opposite to vertex (x3, y3)

    if (a2 > b2 + c2)
        return 2; // Angle at (x1, y1) is obtuse
    if (b2 > a2 + c2)
        return 0; // Angle at (x2, y2) is obtuse
    if (c2 > a2 + b2)
        return 1; // Angle at (x3, y3) is obtuse

    return -1; // No obtuse angles
}

bool is_obtuse_triangle(Polygon& triangle){
    return find_obtuse_angle(triangle) != -1;
}

int count_obtuse_triangles(Problem* problem){
    int count = 0;
    for(Polygon& triangle : problem->get_triangulation()){
        if(is_obtuse_triangle(triangle)){
            count++;
        }
    }
    return count;
}

double squared_distance(Point& a, Point& b){
    return CGAL::to_double(CGAL::squared_distance(a, b));
}
double distance(Point& a, Point& b){
    return std::sqrt(CGAL::to_double(CGAL::squared_distance(a, b)));
}

void Problem::visualize_solution(std::vector<Polygon> problematic_triangles){
    // assuming run from build folder
    std::string path = "../offcenter_delaunay_refinement/ipe/" + name + ".ipe";

    std::vector<Segment> triangle_segments;
    std::vector<Polygon> obtuse_triangles;

    for (Polygon triangle : triangulation){
        for (Segment seg : triangle.edges()){
            triangle_segments.push_back(seg);
        }
        if (is_obtuse_triangle(triangle)){
            obtuse_triangles.push_back(triangle);
        }
    }

    //steiner = {};

    to_IPE(path, points, constraints, {}, boundary.vertices(), {}, triangle_segments, obtuse_triangles, problematic_triangles);
}

void Problem::visualize_solution_voronoi(std::vector<std::pair<Point, Point>> voronoi_edges){
    // assuming run from build folder
    std::string path = "../offcenter_delaunay_refinement/ipe/" + name + ".ipe";

    std::vector<Segment> triangle_segments;
    std::vector<Polygon> obtuse_triangles;

    for (Polygon triangle : triangulation){
        for (Segment seg : triangle.edges()){
            triangle_segments.push_back(seg);
        }
        if (is_obtuse_triangle(triangle)){
            obtuse_triangles.push_back(triangle);
        }
    }

    //steiner = {};

    to_IPE_voronoi(path, points, constraints, {}, boundary.vertices(), {}, triangle_segments, obtuse_triangles, {}, voronoi_edges);
}

void Problem::save_intermidiate_result(std::string path){
    std::vector<Segment> triangle_segments;
    std::vector<Polygon> obtuse_triangles;

    for (Polygon triangle : triangulation){
        for (Segment seg : triangle.edges()){
            triangle_segments.push_back(seg);
        }
        if (is_obtuse_triangle(triangle)){
            obtuse_triangles.push_back(triangle);
        }
    }

    to_SVG(path, points, constraints, boundary, steiner, triangulation, obtuse_triangles);
}

void to_IPE(std::string path, std::vector<Point> points, std::vector<Segment> constraints, std::vector<Segment> new_constraints, std::vector<Point> boundary, std::vector<Point> steiner, std::vector<Segment> triangulation, std::vector<Polygon> obtuse_triangles, std::vector<Polygon> problematic_triangles){
    std::ofstream o(path);

    auto xmin = boundary[0].x();
    auto xmax = boundary[0].x();
    auto ymin = boundary[0].y();
    auto ymax = boundary[0].y();

    for (Point p : boundary){
        xmin = std::min(xmin, p.x());
        xmax = std::max(xmax, p.x());
        ymin = std::min(ymin, p.y());
        ymax = std::max(ymax, p.y());
    }
    auto scale = std::max(xmax - xmin, ymax - ymin);
    
    double scale_d = CGAL::to_double(scale);
    double xmin_d = CGAL::to_double(xmin);
    double xmax_d = CGAL::to_double(xmax);
    double ymin_d = CGAL::to_double(ymin);
    double ymax_d = CGAL::to_double(ymax);

    // Header of the IPE File
    o << "<?xml version=\"1.0\"?>\n";
    o << "<!DOCTYPE ipe SYSTEM \"ipe.dtd\">\n";
    o << "<ipe version=\"70218\" creator=\"Ipe 7.2.24\">\n";
    o << "<info created=\"D:20221020151441\" modified=\"D:20221020151441\"/>\n";
    o << "<ipestyle name=\"basic\">\n";
    o << "<symbol name=\"mark/disk(sx)\" transformations=\"translations\">\n";
    o << "<path fill=\"sym-stroke\">\n";
    o << "0.6 0 0 0.6 0 0 e\n";
    o << "</path>\n";
    o << "</symbol>\n";
    o << "<anglesize name=\"22.5 deg\" value=\"22.5\"/>\n";
    o << "<anglesize name=\"30 deg\" value=\"30\"/>\n";
    o << "<anglesize name=\"45 deg\" value=\"45\"/>\n";
    o << "<anglesize name=\"60 deg\" value=\"60\"/>\n";
    o << "<anglesize name=\"90 deg\" value=\"90\"/>\n";
    o << "<arrowsize name=\"large\" value=\"10\"/>\n";
    o << "<arrowsize name=\"small\" value=\"5\"/>\n";
    o << "<arrowsize name=\"tiny\" value=\"3\"/>\n";

    o << "<color name=\"blue\" value=\"0 0 1\"/>\n";
    o << "<color name=\"gray\" value=\"0.745\"/>\n";
    o << "<color name=\"green\" value=\"0 1 0\"/>\n";
    o << "<color name=\"red\" value=\"1 0 0\"/>\n";
    o << "<color name=\"pink\" value=\"1 0.753 0.796\"/>\n";
    o << "<pen name=\"heavier\" value=\"0.8\"/>\n";
    o << "<pen name=\"fat\" value=\"1.4\"/>\n";
    o << "<pen name=\"ultrafat\" value=\"2\"/>\n";

    o << "<gridsize name=\"16 pts (~6 mm)\" value=\"16\"/>\n";
    o << "<gridsize name=\"32 pts (~12 mm)\" value=\"32\"/>\n";
    o << "<gridsize name=\"4 pts\" value=\"4\"/>\n";
    o << "<gridsize name=\"8 pts (~3 mm)\" value=\"8\"/>\n";
    o << "<opacity name=\"10%\" value=\"0.1\"/>\n";
    o << "<opacity name=\"25%\" value=\"0.25\"/>\n";
    o << "<opacity name=\"50%\" value=\"0.5\"/>\n";
    o << "<opacity name=\"75%\" value=\"0.75\"/>\n";
    o << "<symbolsize name=\"large\" value=\"5\"/>\n";
    o << "<symbolsize name=\"small\" value=\"2\"/>\n";
    o << "<symbolsize name=\"tiny\" value=\"1.1\"/>\n";
    o << "<textsize name=\"huge\" value=\"\\huge\"/>\n";
    o << "<textsize name=\"large\" value=\"\\large\"/>\n";
    o << "<textsize name=\"small\" value=\"\\small\"/>\n";
    o << "<textsize name=\"tiny\" value=\"\tiny\"/>\n";
    o << "<tiling name=\"falling\" angle=\"-60\" step=\"4\" width=\"1\"/>\n";
    o << "<tiling name=\"rising\" angle=\"30\" step=\"4\" width=\"1\"/>\n";
    o << "</ipestyle>\n";
    o << "<page>\n";
    o << "<layer name=\"hull\"/>\n";
    o << "<layer name=\"constraints\"/>\n";
    o << "<layer name=\"triangulation\"/>\n";
    o << "<layer name=\"obtuse\"/>\n";
    o << "<layer name=\"blub\"/>\n";
    o << "<layer name=\"problems\"/>\n";
    o << "<view layers=\"hull constraints triangulation obtuse blub\" active=\"triangulation\"/>\n";

    for (Polygon p : problematic_triangles){
        Point a = p[0];
        Point b = p[1];
        Point c = p[2];
        o << "<path layer=\"problems\" fill=\"red\" stroke-opacity=\"opaque\">\n";
        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272) << " m \n";
        o << ((b.x() - xmin) * 560 / scale + 16) << " " << (b.y() * 560 / scale + 272) << " l \n";
        o << ((c.x() - xmin) * 560 / scale + 16) << " " << (c.y() * 560 / scale + 272) << " l \n";
        o << "</path>\n";
    }

    for (Polygon p : obtuse_triangles){
        Point a = p[0];
        Point b = p[1];
        Point c = p[2];
        o << "<path layer=\"obtuse\" fill=\"pink\" stroke-opacity=\"opaque\">\n";
        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272) << " m \n";
        o << ((b.x() - xmin) * 560 / scale + 16) << " " << (b.y() * 560 / scale + 272) << " l \n";
        o << ((c.x() - xmin) * 560 / scale + 16) << " " << (c.y() * 560 / scale + 272) << " l \n";
        o << "</path>\n";
    }

    for (Segment v : triangulation){
        o << "<path layer=\"triangulation\" stroke=\"black\">\n";

        Point a = v.start();
        Point b = v.end();

        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272) << " m \n";
        o << ((b.x() - xmin) * 560 / scale + 16) << " " << (b.y() * 560 / scale + 272) << " l \n";
        o << "</path>\n";
    }

    for (int i = 0; i < boundary.size(); i++){
        Point a = boundary[i];
        Point b = boundary[(i + 1) % boundary.size()];

        o << "<path layer=\"hull\" stroke=\"blue\" pen=\"fat\">\n";

        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272) << " m \n";
        o << ((b.x() - xmin) * 560 / scale + 16) << " " << (b.y() * 560 / scale + 272) << " l \n";
        o << "</path>\n";
    }

    for (Segment v : constraints){
        Point a = v.start();
        Point b = v.end();

        o << "<path layer=\"constraints\" stroke=\"red\" pen=\"fat\">\n";
        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272) << " m \n";
        o << ((b.x() - xmin) * 560 / scale + 16) << " " << (b.y() * 560 / scale + 272) << " l \n";
        o << "</path>\n";
    }

    for (Segment v : new_constraints){
        Point a = v.start();
        Point b = v.end();

        o << "<path layer=\"constraints\" stroke=\"green\" pen=\"fat\">\n";
        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272) << " m \n";
        o << ((b.x() - xmin) * 560 / scale + 16) << " " << (b.y() * 560 / scale + 272) << " l \n";
        o << "</path>\n";
    }

    for (Point a : points){
        o << "<use layer=\"points\" name=\"mark/disk(sx)\" pos=\"";
        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272);
        o << "\" size=\"normal\" stroke=\"black\"/>\n";
    }

    for (Point a : steiner){
        o << "<use layer=\"blub\" name=\"mark/disk(sx)\" pos=\"";
        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272);
        o << "\" size=\"normal\" stroke=\"red\"/>\n";
    }

    o << "</page>\n";
    o << "</ipe>\n";

    o.close();

    std::string systemCom = "ipe " + path + " > /dev/null 2>&1";
    int systemRet = system(systemCom.c_str());

    if (systemRet == -1){
        printf("Could not open IPE");
    }
}

void to_SVG(std::string path, std::vector<Point> points, std::vector<Segment> constraints, Polygon boundary, std::vector<Point> steiner, std::vector<Polygon> triangles, std::vector<Polygon> obtuse_triangles) {
    std::ofstream o(path);
    
    // Compute bounding box
    auto xmin = boundary[0].x();
    auto xmax = boundary[0].x();
    auto ymin = boundary[0].y();
    auto ymax = boundary[0].y();

    for (Point p : boundary) {
        xmin = std::min(xmin, p.x());
        xmax = std::max(xmax, p.x());
        ymin = std::min(ymin, p.y());
        ymax = std::max(ymax, p.y());
    }
    auto scale = std::max(xmax - xmin, ymax - ymin);

    double width = 600;
    double height = 600;

    o << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
    o << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"" << width << "\" height=\"" << height << "\">\n";

    // Add a white background rectangle
    o << "<rect width=\"100%\" height=\"100%\" fill=\"white\" />\n";

    // Helper lambda to scale and shift points
    auto transform_x = [&](double x) { return ((x - xmin) * 560.0 / scale) + 16.0; };
    auto transform_y = [&](double y) { return height - (((y - ymin) * 560.0 / scale) + 16.0); }; // invert Y axis

    // Fill obtuse triangles with pink
    for (const auto& poly : obtuse_triangles) {
        o << "<polygon points=\"";
        for (const auto& p : poly) {
            o << transform_x(CGAL::to_double(p.x())) << "," << transform_y(CGAL::to_double(p.y())) << " ";
        }
        o << "\" fill=\"pink\" stroke=\"none\" opacity=\"1\" />\n";
    }

    // Draw triangles (black edges)
    for (auto poly : triangles) {
        for (std::size_t i = 0; i < poly.size(); ++i) {
            const Point& a = poly[i];
            const Point& b = poly[(i + 1) % poly.size()];

            o << "<line x1=\"" << transform_x(CGAL::to_double(a.x())) << "\" y1=\"" << transform_y(CGAL::to_double(a.y())) 
              << "\" x2=\"" << transform_x(CGAL::to_double(b.x())) << "\" y2=\"" << transform_y(CGAL::to_double(b.y())) 
              << "\" stroke=\"black\" stroke-width=\"1\" />\n";
        }
    }

    // Draw boundary (blue, thicker)
    for (std::size_t i = 0; i < boundary.size(); ++i) {
        Point a = boundary[i];
        Point b = boundary[(i + 1) % boundary.size()];

        o << "<line x1=\"" << transform_x(CGAL::to_double(a.x())) << "\" y1=\"" << transform_y(CGAL::to_double(a.y())) 
          << "\" x2=\"" << transform_x(CGAL::to_double(b.x())) << "\" y2=\"" << transform_y(CGAL::to_double(b.y())) 
          << "\" stroke=\"blue\" stroke-width=\"3\" />\n";
    }

    // Draw constraint edges (red, thicker)
    for (Segment v : constraints) {
        Point a = v.start();
        Point b = v.end();

        o << "<line x1=\"" << transform_x(CGAL::to_double(a.x())) << "\" y1=\"" << transform_y(CGAL::to_double(a.y())) 
          << "\" x2=\"" << transform_x(CGAL::to_double(b.x())) << "\" y2=\"" << transform_y(CGAL::to_double(b.y())) 
          << "\" stroke=\"red\" stroke-width=\"2\" />\n";
    }

    // Draw points (black disks)
    for (Point a : points) {
        o << "<circle cx=\"" << transform_x(CGAL::to_double(a.x())) << "\" cy=\"" << transform_y(CGAL::to_double(a.y())) 
          << "\" r=\"2\" fill=\"black\" />\n";
    }

    // Draw Steiner points (red disks)
    /*for (Point a : steiner) {
        o << "<circle cx=\"" << transform_x(CGAL::to_double(a.x())) << "\" cy=\"" << transform_y(CGAL::to_double(a.y())) 
          << "\" r=\"2\" fill=\"red\" />\n";
    }*/

    o << "</svg>\n";
    o.close();
}

void to_IPE_voronoi(std::string path, std::vector<Point> points, std::vector<Segment> constraints, std::vector<Segment> new_constraints, std::vector<Point> boundary, std::vector<Point> steiner, std::vector<Segment> triangulation, std::vector<Polygon> obtuse_triangles, std::vector<Polygon> problematic_triangles, std::vector<std::pair<Point, Point>> voronoi_edges){
    std::ofstream o(path);

    auto xmin = boundary[0].x();
    auto xmax = boundary[0].x();
    auto ymin = boundary[0].y();
    auto ymax = boundary[0].y();

    for (Point p : boundary){
        xmin = std::min(xmin, p.x());
        xmax = std::max(xmax, p.x());
        ymin = std::min(ymin, p.y());
        ymax = std::max(ymax, p.y());
    }
    auto scale = std::max(xmax - xmin, ymax - ymin);

    // Header of the IPE File
    o << "<?xml version=\"1.0\"?>\n";
    o << "<!DOCTYPE ipe SYSTEM \"ipe.dtd\">\n";
    o << "<ipe version=\"70218\" creator=\"Ipe 7.2.24\">\n";
    o << "<info created=\"D:20221020151441\" modified=\"D:20221020151441\"/>\n";
    o << "<ipestyle name=\"basic\">\n";
    o << "<symbol name=\"mark/disk(sx)\" transformations=\"translations\">\n";
    o << "<path fill=\"sym-stroke\">\n";
    o << "0.6 0 0 0.6 0 0 e\n";
    o << "</path>\n";
    o << "</symbol>\n";
    o << "<anglesize name=\"22.5 deg\" value=\"22.5\"/>\n";
    o << "<anglesize name=\"30 deg\" value=\"30\"/>\n";
    o << "<anglesize name=\"45 deg\" value=\"45\"/>\n";
    o << "<anglesize name=\"60 deg\" value=\"60\"/>\n";
    o << "<anglesize name=\"90 deg\" value=\"90\"/>\n";
    o << "<arrowsize name=\"large\" value=\"10\"/>\n";
    o << "<arrowsize name=\"small\" value=\"5\"/>\n";
    o << "<arrowsize name=\"tiny\" value=\"3\"/>\n";

    o << "<color name=\"blue\" value=\"0 0 1\"/>\n";
    o << "<color name=\"gray\" value=\"0.745\"/>\n";
    o << "<color name=\"green\" value=\"0 1 0\"/>\n";
    o << "<color name=\"red\" value=\"1 0 0\"/>\n";
    o << "<color name=\"pink\" value=\"1 0.753 0.796\"/>\n";
    o << "<pen name=\"heavier\" value=\"0.8\"/>\n";
    o << "<pen name=\"fat\" value=\"1.4\"/>\n";
    o << "<pen name=\"ultrafat\" value=\"2\"/>\n";

    o << "<gridsize name=\"16 pts (~6 mm)\" value=\"16\"/>\n";
    o << "<gridsize name=\"32 pts (~12 mm)\" value=\"32\"/>\n";
    o << "<gridsize name=\"4 pts\" value=\"4\"/>\n";
    o << "<gridsize name=\"8 pts (~3 mm)\" value=\"8\"/>\n";
    o << "<opacity name=\"10%\" value=\"0.1\"/>\n";
    o << "<opacity name=\"25%\" value=\"0.25\"/>\n";
    o << "<opacity name=\"50%\" value=\"0.5\"/>\n";
    o << "<opacity name=\"75%\" value=\"0.75\"/>\n";
    o << "<symbolsize name=\"large\" value=\"5\"/>\n";
    o << "<symbolsize name=\"small\" value=\"2\"/>\n";
    o << "<symbolsize name=\"tiny\" value=\"1.1\"/>\n";
    o << "<textsize name=\"huge\" value=\"\\huge\"/>\n";
    o << "<textsize name=\"large\" value=\"\\large\"/>\n";
    o << "<textsize name=\"small\" value=\"\\small\"/>\n";
    o << "<textsize name=\"tiny\" value=\"\tiny\"/>\n";
    o << "<tiling name=\"falling\" angle=\"-60\" step=\"4\" width=\"1\"/>\n";
    o << "<tiling name=\"rising\" angle=\"30\" step=\"4\" width=\"1\"/>\n";
    o << "</ipestyle>\n";
    o << "<page>\n";
    o << "<layer name=\"hull\"/>\n";
    o << "<layer name=\"constraints\"/>\n";
    o << "<layer name=\"triangulation\"/>\n";
    o << "<layer name=\"obtuse\"/>\n";
    o << "<layer name=\"blub\"/>\n";
    o << "<layer name=\"problems\"/>\n";
    o << "<view layers=\"hull constraints triangulation obtuse blub\" active=\"triangulation\"/>\n";

    for (Polygon p : problematic_triangles){
        Point a = p[0];
        Point b = p[1];
        Point c = p[2];
        o << "<path layer=\"problems\" fill=\"red\" stroke-opacity=\"opaque\">\n";
        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272) << " m \n";
        o << ((b.x() - xmin) * 560 / scale + 16) << " " << (b.y() * 560 / scale + 272) << " l \n";
        o << ((c.x() - xmin) * 560 / scale + 16) << " " << (c.y() * 560 / scale + 272) << " l \n";
        o << "</path>\n";
    }

    for (Polygon p : obtuse_triangles){
        Point a = p[0];
        Point b = p[1];
        Point c = p[2];
        o << "<path layer=\"obtuse\" fill=\"pink\" stroke-opacity=\"opaque\">\n";
        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272) << " m \n";
        o << ((b.x() - xmin) * 560 / scale + 16) << " " << (b.y() * 560 / scale + 272) << " l \n";
        o << ((c.x() - xmin) * 560 / scale + 16) << " " << (c.y() * 560 / scale + 272) << " l \n";
        o << "</path>\n";
    }

    for (Segment v : triangulation){
        o << "<path layer=\"triangulation\" stroke=\"black\">\n";

        Point a = v.start();
        Point b = v.end();

        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272) << " m \n";
        o << ((b.x() - xmin) * 560 / scale + 16) << " " << (b.y() * 560 / scale + 272) << " l \n";
        o << "</path>\n";
    }

    for (int i = 0; i < boundary.size(); i++){
        Point a = boundary[i];
        Point b = boundary[(i + 1) % boundary.size()];

        o << "<path layer=\"hull\" stroke=\"blue\" pen=\"fat\">\n";

        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272) << " m \n";
        o << ((b.x() - xmin) * 560 / scale + 16) << " " << (b.y() * 560 / scale + 272) << " l \n";
        o << "</path>\n";
    }

    for (Segment v : constraints){
        Point a = v.start();
        Point b = v.end();

        o << "<path layer=\"constraints\" stroke=\"red\" pen=\"fat\">\n";
        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272) << " m \n";
        o << ((b.x() - xmin) * 560 / scale + 16) << " " << (b.y() * 560 / scale + 272) << " l \n";
        o << "</path>\n";
    }

    for (Segment v : new_constraints){
        Point a = v.start();
        Point b = v.end();

        o << "<path layer=\"constraints\" stroke=\"green\" pen=\"fat\">\n";
        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272) << " m \n";
        o << ((b.x() - xmin) * 560 / scale + 16) << " " << (b.y() * 560 / scale + 272) << " l \n";
        o << "</path>\n";
    }

    for (Point a : points){
        o << "<use layer=\"points\" name=\"mark/disk(sx)\" pos=\"";
        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272);
        o << "\" size=\"normal\" stroke=\"black\"/>\n";
    }

    for (Point a : steiner){
        o << "<use layer=\"blub\" name=\"mark/disk(sx)\" pos=\"";
        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272);
        o << "\" size=\"normal\" stroke=\"red\"/>\n";
    }

    for(auto edge : voronoi_edges){
        Point a = edge.first;
        Point b = edge.second;

        o << "<path layer=\"voronoi\" stroke=\"green\" pen=\"heavier\">\n";
        o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272) << " m\n";
        o << ((b.x() - xmin) * 560 / scale + 16) << " " << (b.y() * 560 / scale + 272) << " l\n";
        o << "</path>\n";
    }

    o << "</page>\n";
    o << "</ipe>\n";

    o.close();

    std::string systemCom = "ipe " + path + " > /dev/null 2>&1";
    int systemRet = system(systemCom.c_str());

    if (systemRet == -1){
        printf("Could not open IPE");
    }
}