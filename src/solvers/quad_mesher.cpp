#include "quad_mesher.hpp"
#include "quad_optimization.hpp"

#include <CGAL/intersections.h>
#include <vector>
#include <gmsh.h>

#include <fstream>
#include <algorithm>


typedef CGAL::Bbox_2 Bbox;

void real_coord_to_ipe(Point p, Bbox box)
{
    auto scale = std::max(box.xmax() - box.xmin(), box.ymax() - box.ymin());

    auto x = (p.x() - box.xmin())*560/scale + 16;
    auto y = p.y()*560/scale + 272;

    std::cout << GREEN << "x: " << x << " y: " << y << RESET << std::endl;
}

double shortest_edge_length(const std::vector<Polygon>& quads, Problem* problem) {
    double min1 = std::numeric_limits<double>::infinity();
    double min2 = std::numeric_limits<double>::infinity();
    double min3 = std::numeric_limits<double>::infinity();

    Segment edge1, edge2, edge3;

    for (const auto& poly : quads) {
        int n = poly.size();
        for (int i = 0; i < n; ++i) {
            const Point& p1 = poly[i];
            const Point& p2 = poly[(i + 1) % n]; // wrap around
            double len = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p2)));

            if (len < min1) {
                min3 = min2; edge3 = edge2;
                min2 = min1; edge2 = edge1;
                min1 = len;  edge1 = Segment(p1, p2);
            } else if (len < min2) {
                min3 = min2; edge3 = edge2;
                min2 = len;  edge2 = Segment(p1, p2);
            } else if (len < min3) {
                min3 = len;  edge3 = Segment(p1, p2);
            }
        }
    }

    Bbox bbox = problem->get_boundary().bbox();

    // Print shortest edge
    std::cout << GREEN << "Shortest edge:\n";
    std::cout << "x: " << edge1[0].x() << " y: " << edge1[0].y() << "\n";
    std::cout << "x: " << edge1[1].x() << " y: " << edge1[1].y() << RESET << std::endl;
    real_coord_to_ipe(edge1[0], bbox);
    real_coord_to_ipe(edge1[1], bbox);


    // Print third shortest edge
    std::cout << RED << "Third shortest edge:\n";
    std::cout << "x: " << edge3[0].x() << " y: " << edge3[0].y() << "\n";
    std::cout << "x: " << edge3[1].x() << " y: " << edge3[1].y() << RESET << std::endl;
    real_coord_to_ipe(edge3[0], bbox);
    real_coord_to_ipe(edge3[1], bbox);

    return min3;
}


double get_sizing_quad(Problem* problem) {
    Polygon boundary = problem->get_boundary();
    auto bbox = boundary.bbox();
    double max_size = std::sqrt(std::pow(bbox.x_span(), 2) + std::pow(bbox.y_span(), 2)) * 0.1;
    
    return max_size / 2;
}

int find_point_idx_medians(Point& p, std::vector<Point>& points, std::vector<Point>& steiner){
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

double angle(const Point& prev, const Point& curr, const Point& next) {
    Vector v1 = prev - curr;
    Vector v2 = next - curr;
    double dot = CGAL::to_double(v1 * v2);
    double det = CGAL::to_double(v1.x() * v2.y() - v1.y() * v2.x());
    return std::atan2(std::abs(det), dot) * 180.0 / CGAL_PI;
}

double average_angle_deviation(const std::vector<Polygon>& quads) {
    double total_deviation = 0.0;
    int total_angles = 0;

    for (const auto& quad : quads) {
        if (quad.size() != 4) continue; // skip non-quads

        for (int i = 0; i < 4; ++i) {
            const Point& prev = quad[(i + 3) % 4];
            const Point& curr = quad[i];
            const Point& next = quad[(i + 1) % 4];
            double theta = angle(prev, curr, next);
            total_deviation += std::abs(theta - 90.0);
            ++total_angles;
        }
    }

    return total_angles > 0 ? total_deviation / total_angles : 0.0;
}

double edge_length_quad(const Point& a, const Point& b) {
    return std::sqrt(CGAL::to_double(CGAL::squared_distance(a, b)));
}

void compute_angles(const std::vector<Polygon>& quads, std::string path){
    std::vector<double> all_angles;

    for (const Polygon& quad : quads) {
        if (quad.size() != 4) continue; // skip non-triangles

        const Point& A = quad[0];
        const Point& B = quad[1];
        const Point& C = quad[2];
        const Point& D = quad[3];

        all_angles.push_back(angle(D, A, B));
        all_angles.push_back(angle(A, B, C));
        all_angles.push_back(angle(B, C, D));
        all_angles.push_back(angle(C, D, A));
    }

    std::ofstream out(path);
    out << "angle\n"; 
    for (double angle : all_angles) {
        out << angle << "\n";
    }
    out.close();
}

std::vector<double> compute_aspect_ratios(const std::vector<Polygon>& quads, std::string path) {
    std::vector<double> aspect_ratios;

    for (const auto& quad : quads) {
        if (quad.size() != 4) continue; // skip invalid quads

        std::vector<double> lengths;
        for (int i = 0; i < 4; ++i) {
            lengths.push_back(edge_length_quad(quad[i], quad[(i + 1) % 4]));
        }

        double min_len = *std::min_element(lengths.begin(), lengths.end());
        double max_len = *std::max_element(lengths.begin(), lengths.end());

        if (min_len > 0.0) {
            aspect_ratios.push_back(max_len / min_len);
        } else {
            aspect_ratios.push_back(std::numeric_limits<double>::infinity()); // degenerate case
        }
    }

    std::ofstream out(path);
    out << "ratio\n"; 
    for (double ratio : aspect_ratios) {
        out << ratio << "\n";
    }
    out.close();

    return aspect_ratios;
}

double average_ar(const std::vector<double>& ars) {
    double sum = 0.0;
    for (double ar : ars) sum += ar;
    return ars.empty() ? 0.0 : sum / ars.size();
}

double max_ar(const std::vector<double>& ars) {
    return ars.empty() ? 0.0 : *std::max_element(ars.begin(), ars.end());
}

std::tuple<Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi, Eigen::MatrixXd>
convert_quad_mesh_medians(Problem* problem, std::vector<Polygon>& quads){

    auto points = problem->get_points(); 
    auto steiner = problem->get_steiner(); 
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

    // Convert triangles into a matrix
    Eigen::MatrixXi F(quads.size(), 4);
    for (size_t i = 0; i < quads.size(); ++i) {
        const Polygon& polygon = quads[i];

        Point a = polygon.vertex(0);
        Point b = polygon.vertex(1);
        Point c = polygon.vertex(2);
        Point d = polygon.vertex(3);

        F(i, 0) = find_point_idx_medians(a, points, steiner);
        F(i, 1) = find_point_idx_medians(b, points, steiner);
        F(i, 2) = find_point_idx_medians(c, points, steiner);
        F(i, 3) = find_point_idx_medians(d, points, steiner);
    }

    // Define barrier terms for the constraints
    Eigen::VectorXi B(points.size());
    Eigen::MatrixXd B_VAR(points.size(), 2);
    for (size_t i = 0; i < points.size(); ++i) {
        B_VAR(i, 0) = CGAL::to_double(points[i].x()); 
        B_VAR(i, 1) = CGAL::to_double(points[i].y()); 
        B(i) = i;
    }
    return {V, F, B, B_VAR};
}

std::vector<Polygon> generate_quads_medians(Problem* problem, CDT& cdt){
    std::vector<Polygon> quads;
    std::set<Point> steiner_point_set;

    for ( Face_handle f : cdt.finite_face_handles()) {
        if (!problem->triangle_is_inside(cdt, f)) continue;

        Point a = f->vertex(0)->point();
        Point b = f->vertex(1)->point();
        Point c = f->vertex(2)->point();
        Point p1 = CGAL::midpoint(a, b);
        Point p2 = CGAL::midpoint(b, c);
        Point p3 = CGAL::midpoint(a, c);

        Point o;

        Line line1(a, p2);
        Line line2(b, p3);

        // Compute the intersection
        auto result = CGAL::intersection(line1, line2);

        if (const Point* ip = boost::get<Point>(&*result)) {
            o = *ip;
        } else {
            throw std::runtime_error("Lines do not intersect in a single point");
        }

        Polygon q1, q2, q3;

        q1.push_back(p3);
        q1.push_back(a);
        q1.push_back(p1);
        q1.push_back(o);

        q2.push_back(p1);
        q2.push_back(b);
        q2.push_back(p2);
        q2.push_back(o);

        q3.push_back(p2);
        q3.push_back(c);
        q3.push_back(p3);
        q3.push_back(o);

        quads.push_back(q1);
        quads.push_back(q2);
        quads.push_back(q3);

        steiner_point_set.emplace(o);
        steiner_point_set.emplace(p1);
        steiner_point_set.emplace(p2);
        steiner_point_set.emplace(p3);
    }
    for(auto steiner_point : steiner_point_set){
        problem->add_steiner(steiner_point);
    }
    return quads;
}

void build_quad_mesh_medians(Problem* problem){
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

    Mesher mesher(cdt);
    double max_size = get_sizing_quad(problem);
    mesher.set_criteria(Criteria(0, max_size*2));
    mesher.refine_mesh();

    problem->update_problem(cdt, point_set);
    std::cout << GREEN << "num steiner: " << problem->get_steiner().size() << RESET << std::endl;
    auto quads = generate_quads_medians(problem, cdt);

    to_IPE_quad("../solutions/ipe/QUAD-" + problem->get_name() + ".ipe", problem->get_points(), problem->get_constraints(), problem->get_boundary(), {}, quads);

    //* Stats before
    std::cout << RED << "Angle deviation before: " << average_angle_deviation(quads) << RESET << std::endl;
    auto ars = compute_aspect_ratios(quads, "../results/aspect_ratios_quad_medians_initial_angles.csv");
    std::cout << RED << "Average AR before = " << average_ar(ars) << RESET << std::endl;
    std::cout << RED << "Max AR before = " << max_ar(ars) << RESET << std::endl;
    compute_angles(quads, "../results/angles_quad_medians_initial_angles.csv");

    auto [V, F, B, B_VAR] = convert_quad_mesh_medians(problem, quads);

    auto optimized_quads = optimize_TinyAD_quad(problem, V, F, B, B_VAR);

    std::cout << GREEN << "shortest edge length: " << shortest_edge_length(quads, problem) << std::endl;

    //* Stats after
    std::cout << RED << "Angle deviation after: " << average_angle_deviation(optimized_quads) << RESET << std::endl;
    ars = compute_aspect_ratios(optimized_quads, "../results/aspect_ratios_quad_medians_optimized_angles.csv");
    std::cout << RED << "Average AR after = " << average_ar(ars) << RESET << std::endl;
    std::cout << RED << "Max AR after = " << max_ar(ars) << RESET << std::endl;
    compute_angles(optimized_quads, "../results/angles_quad_medians_optimized_angles.csv");

    to_IPE_quad("../solutions/ipe/QUAD-" + problem->get_name() + ".ipe", problem->get_points(), problem->get_constraints(), problem->get_boundary(), {}, optimized_quads);
}

// GMSH Quad mesh---------------------------------------------------------------------------------------------

double get_sq_length_of_shortest_edge(Polygon& boundary){
    double min_length = DBL_MAX;
    for(int i = 0; i<boundary.size(); i++){
        Point a = boundary[i];
        Point b = boundary[(i+1)%boundary.size()];

        double length = CGAL::to_double(CGAL::squared_distance(a, b));
        if(length < min_length){
            min_length = length;
        }
    }
    return min_length;
}

double get_index_of_longest_edge(Polygon& boundary){
    double max_length = 0;
    int index;
    for(int i = 0; i<boundary.size(); i++){
        Point a = boundary[i];
        Point b = boundary[(i+1)%boundary.size()];

        double length = CGAL::to_double(CGAL::squared_distance(a, b));
        if(length > max_length){
            max_length = length;
            index = i;
        }
    }
    return index;
}

std::tuple<Point, int> calculate_additional_point(Problem* problem){
    Polygon boundary = problem->get_boundary();

    double min_length = get_sq_length_of_shortest_edge(boundary);

    int i = get_index_of_longest_edge(boundary);

    Point a = boundary[i];
    Point b = boundary[(i+1)%boundary.size()];

    // Check whether the shortest edge is the last boundary edge
    /*Point a = boundary[boundary.size() - 1];
    Point b = boundary[0];
    if((std::abs(CGAL::to_double(CGAL::squared_distance(a, b)) - min_length)) < 1e-10){
        a = boundary[0];
        b = boundary[1];
        std::cout << RED << "Added point is on the first boundary Segment!" << RESET << std::endl;
    }*/

    min_length = std::sqrt(min_length);

    Vector ab = a - b;
    double ab_len = std::sqrt(CGAL::to_double(ab.squared_length()));
    double scale = min_length / ab_len;
    FT exact_scale(scale);
    Vector direction = ab * exact_scale;
    Point c = a - direction;

    Line line(a, b);
    Point projected_point = line.projection(c);
    return {projected_point, i};
}

std::vector<int> get_ccw_sorted_values(const std::map<Point, int>& point_map) {
    std::vector<std::pair<Point, int>> entries(point_map.begin(), point_map.end());

    // Step 1: Compute centroid
    double cx = 0.0, cy = 0.0;
    for (const auto& [p, _] : entries) {
        cx += CGAL::to_double(p.x());
        cy += CGAL::to_double(p.y());
    }
    cx /= entries.size();
    cy /= entries.size();
    Point centroid(cx, cy);

    // Step 2: Sort by angle relative to centroid
    std::sort(entries.begin(), entries.end(),
        [&centroid](const std::pair<Point, int>& a, const std::pair<Point, int>& b) {
            double angle_a = std::atan2(CGAL::to_double(a.first.y() - centroid.y()), CGAL::to_double(a.first.x() - centroid.x()));
            double angle_b = std::atan2(CGAL::to_double(b.first.y() - centroid.y()), CGAL::to_double(b.first.x() - centroid.x()));
            return angle_a < angle_b;
        });

    // Step 3: Extract values in sorted order
    std::vector<int> sorted_values;
    for (const auto& [_, val] : entries) {
        sorted_values.push_back(val);
    }

    return sorted_values;
}

std::tuple<std::vector<Polygon>, Eigen::MatrixXd, Eigen::MatrixXi, Eigen::VectorXi, Eigen::MatrixXd>
convert_quad_mesh(Problem* problem){

    std::vector<Polygon> quads;
    // Init structures for optimization
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::VectorXi B(problem->get_points().size());
    Eigen::MatrixXd B_VAR(problem->get_points().size(), 2);

    try {
        std::vector<std::size_t> nodeTags;
        std::vector<double> coords;
        std::vector<double> paramCoords;
        gmsh::model::mesh::getNodes(nodeTags, coords, paramCoords);

        V.resize(nodeTags.size(), 2);

        std::unordered_map<std::size_t, std::array<double, 3>> nodeIdToCoord;
        for (std::size_t i = 0; i < nodeTags.size(); ++i) {
            nodeIdToCoord[nodeTags[i]] = {
                coords[3 * i + 0],
                coords[3 * i + 1],
                coords[3 * i + 2]
            };
            V(i, 0) = coords[3 * i + 0];
            V(i, 1) = coords[3 * i + 1];
            if(i < problem->get_points().size()){
                B(i) = i;
                B_VAR(i, 0) = coords[3 * i + 0];
                B_VAR(i, 1) = coords[3 * i + 1];
            }
        }

        // STEP 2: Loop over all 2D entities
        std::vector<std::pair<int, int>> entities;
        gmsh::model::getEntities(entities, 2); // dimension 2 = surface elements

        for (const auto &[dim, tag] : entities) {
            std::vector<int> elementTypes;
            std::vector<std::vector<std::size_t>> elementTags;
            std::vector<std::vector<std::size_t>> nodeTagsPerElement;

            gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTagsPerElement, dim, tag);

            for (std::size_t i = 0; i < elementTypes.size(); ++i) {
                int elementType = elementTypes[i];
                std::string name;
                int dimDummy;
                int order;
                int numNodesPerElem;
                std::vector<double> dummyCoords;
                int numPrimaryNodes;

                gmsh::model::mesh::getElementProperties(
                    elementType, name, dimDummy, order, numNodesPerElem, dummyCoords, numPrimaryNodes
                );

                //std::cout << "Element type: " << name << " (nodes per face: " << numNodesPerElem << ")\n";

                const auto &nodeListFlat = nodeTagsPerElement[i];

                // Init faces for optimization
                int totalFaceCount = nodeListFlat.size();
                F.resize(totalFaceCount/4, 4);

                for (std::size_t j = 0; j < nodeListFlat.size(); j += numNodesPerElem) {
                    //std::cout << "Face:\n";

                    Polygon q;
                    int nodeId0 = nodeListFlat[j + 0];
                    int nodeId1 = nodeListFlat[j + 1];
                    int nodeId2 = nodeListFlat[j + 2];
                    int nodeId3 = nodeListFlat[j + 3];
                    const auto &xyz0 = nodeIdToCoord[nodeId0];
                    const auto &xyz1 = nodeIdToCoord[nodeId1];
                    const auto &xyz2 = nodeIdToCoord[nodeId2];
                    const auto &xyz3 = nodeIdToCoord[nodeId3];

                    std::map<Point, int> current_quad;
                    current_quad[Point(xyz0[0], xyz0[1])] = nodeId0;
                    current_quad[Point(xyz1[0], xyz1[1])] = nodeId1;
                    current_quad[Point(xyz2[0], xyz2[1])] = nodeId2;
                    current_quad[Point(xyz3[0], xyz3[1])] = nodeId3;

                    std::vector<int> ccw = get_ccw_sorted_values(current_quad);
                    //std::cout << RED << numNodesPerElem << RESET << std::endl;
                    F(j/4, 0) = ccw[0] - 1;
                    F(j/4, 1) = ccw[1] - 1;
                    F(j/4, 2) = ccw[2] - 1;
                    F(j/4, 3) = ccw[3] - 1;


                    for (int k = 0; k < numNodesPerElem; ++k) {
                        std::size_t nodeId = nodeListFlat[j + k];
                        const auto &xyz = nodeIdToCoord[nodeId];
                        //std::cout << "  Vertex " << nodeId << ": ("<< xyz[0] << ", " << xyz[1] << ", " << xyz[2] << ")\n";

                        Point p(xyz[0], xyz[1]);
                        q.push_back(p);

                        //std::cout << "  j: " << j << " k: " << k << std::endl;
                        //F(j/4,k) = nodeId - 1;
                    }

                    Point q0 = q[0];
                    
                    quads.push_back(q);
                }
            }
        }

    } catch (const std::exception &e) {
        std::cerr << "Gmsh error: " << e.what() << "\n";
    }
    gmsh::finalize();

    return {quads, V, F, B, B_VAR};
}

void build_quad_mesh_gmsh(Problem* problem){
    try{
        gmsh::initialize();
        gmsh::model::add(problem->get_name());

        double lc = get_sizing_quad(problem);  // uniform mesh size

        auto points = problem->get_points();
        auto boundary = problem->get_boundary();
        auto constraints = problem->get_constraints();
        auto boundary_indices = problem->get_boundary_indices();
        auto constraints_indices = problem->get_constraints_indices();

        Eigen::MatrixXd VARS(points.size(), 2);
        for (int i = 0; i < points.size(); ++i) {
            VARS(i, 0) = CGAL::to_double(points[i].x()); 
            VARS(i, 1) = CGAL::to_double(points[i].y()); 
        }

        for(int i = 0; i < points.size(); ++i){
            int p = gmsh::model::geo::addPoint(VARS(i, 0), VARS(i, 1), 0, lc);
        }

        // Add additional point to make the num of triangles in the triangle mesh even
        //auto [c, index] = calculate_additional_point(problem);
        //std::cout << RED << "final point: " << c.x() << " " << c.y() << RESET << std::endl;
        //int p = gmsh::model::geo::addPoint(CGAL::to_double(c.x()), CGAL::to_double(c.y()), 0, lc);
        //boundary_indices.push_back(p-1);
        //boundary_indices.insert(boundary_indices.begin() + index + 1, p-1);
        //problem->add_steiner(c);

        //to_IPE_quad("../solutions/ipe/OPTIMIZED_QUAD-" + problem->get_name() + ".ipe", problem->get_points(), problem->get_constraints(), problem->get_boundary(), problem->get_steiner(), {});

        std::vector<int> boundary_segments;
        for(int i = 0; i < boundary_indices.size(); i++){
            int p1 = boundary_indices[i] + 1;
            int p2 = boundary_indices[(i+1)%boundary_indices.size()] + 1;
            int l = gmsh::model::geo::addLine(p1, p2);
            boundary_segments.push_back(l);
            //std::cout << GREEN << p1 << " " << p2 << " | " << RESET;
        }
        //std::cout << std::endl;
        int cl = gmsh::model::geo::addCurveLoop(boundary_segments);
        int pl = gmsh::model::geo::addPlaneSurface({cl});

        gmsh::model::geo::synchronize();

        // embed input points into the mesh
        for(int i = 0; i < points.size(); ++i){
            gmsh::model::mesh::embed(0, {i+1}, 2, pl);
        }

        // embed constraints into the mesh
        std::vector<int> constraints_lines;
        for(int i = 0; i<constraints.size(); i++){
            auto [i1, i2] = constraints_indices[i];
            int constraint_line = gmsh::model::geo::addLine(i1+1, i2+1);
            constraints_lines.push_back(constraint_line);
        }
        gmsh::model::geo::synchronize();
        for(int i = 0; i<constraints.size(); i++){
            gmsh::model::mesh::embed(1, {constraints_lines[i]}, 2, pl);
        }

        gmsh::option::setNumber("Mesh.Optimize", 0);
        gmsh::option::setNumber("Mesh.OptimizeNetgen", 0);
        gmsh::option::setNumber("Mesh.Smoothing", 0);

        gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 2);
        gmsh::model::mesh::setRecombine(2, pl);

        gmsh::model::mesh::generate(2);

        // Print result to result.txt
        add_text_and_overwrite(problem, "../result.txt");
        
        auto [quads, V, F, B, B_VAR] = convert_quad_mesh(problem);

        //* Stats before
        std::cout << RED << "Angle deviation before: " << average_angle_deviation(quads) << RESET << std::endl;
        auto ars = compute_aspect_ratios(quads, "../results/aspect_ratios_quad_blossom_initial_2.csv");
        std::cout << RED << "Average AR before = " << average_ar(ars) << RESET << std::endl;
        std::cout << RED << "Max AR before = " << max_ar(ars) << RESET << std::endl;
        compute_angles(quads, "../results/angles_quad_blossom_initial_2.csv");

        to_IPE_quad("../solutions/ipe/QUAD-" + problem->get_name() + ".ipe", problem->get_points(), problem->get_constraints(), problem->get_boundary(), {}, quads);
        //to_SVG_quad("../quad_meshes/QUAD-" + problem->get_name() + ".svg", problem->get_points(), problem->get_constraints(), problem->get_boundary(), {}, quads);

        auto optimized_quads = optimize_TinyAD_quad(problem, V, F, B, B_VAR);

        std::cout << GREEN << shortest_edge_length(optimized_quads, problem) << RESET << std::endl;

        //* Stats after
        std::cout << RED << "Angle deviation after: " << average_angle_deviation(optimized_quads) << RESET << std::endl;
        ars = compute_aspect_ratios(optimized_quads, "../results/aspect_ratios_quad_blossom_optimized_angles_2.csv");
        std::cout << RED << "Average AR after = " << average_ar(ars) << RESET << std::endl;
        std::cout << RED << "Max AR after = " << max_ar(ars) << RESET << std::endl;
        compute_angles(optimized_quads, "../results/angles_quad_blossom_optimized_angles_2.csv");

        to_IPE_quad("../solutions/ipe/OPTIMIZED_QUAD-" + problem->get_name() + ".ipe", problem->get_points(), problem->get_constraints(), problem->get_boundary(), {}, optimized_quads);
    } catch (const std::exception &e) {
        std::cerr << GREEN << "Could not mesh instance " << problem->get_name() << RESET << std::endl;
    }
}

void to_IPE_quad(std::string path, std::vector<Point> points, std::vector<Segment> constraints, Polygon boundary, std::vector<Point> steiner, std::vector<Polygon> quads){
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
    o << "<view layers=\"hull constraints triangulation obtuse blub\" active=\"triangulation\"/>\n";


    for (auto poly : quads) {
        for (std::size_t i = 0; i < poly.size(); ++i) {
            const Point& a = poly[i];
            const Point& b = poly[(i + 1) % poly.size()];

            o << "<path layer=\"triangulation\" stroke=\"black\">\n";
            o << ((a.x() - xmin) * 560 / scale + 16) << " " << (a.y() * 560 / scale + 272) << " m \n";
            o << ((b.x() - xmin) * 560 / scale + 16) << " " << (b.y() * 560 / scale + 272) << " l \n";
            o << "</path>\n";
        }
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
#if __APPLE__
    systemCom = "open -W /Applications/Ipe.app " + path;
#endif
    int systemRet = system(systemCom.c_str());

    if (systemRet == -1){
        printf("Could not open IPE");
    }
}

void to_SVG_quad(std::string path, std::vector<Point> points, std::vector<Segment> constraints, Polygon boundary, std::vector<Point> steiner, std::vector<Polygon> quads) {
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

    // Draw quads (black edges)
    for (auto poly : quads) {
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
    for (Point a : steiner) {
        o << "<circle cx=\"" << transform_x(CGAL::to_double(a.x())) << "\" cy=\"" << transform_y(CGAL::to_double(a.y())) 
          << "\" r=\"2\" fill=\"red\" />\n";
    }

    o << "</svg>\n";
    o.close();
}

void add_text_and_overwrite(Problem* problem, const std::string& path) {
    // Step 1: Read the existing content
    std::ifstream input_file(path);
    std::stringstream buffer;

    if (input_file.is_open()) {
        buffer << input_file.rdbuf(); // read entire file into buffer
        input_file.close();
    } // if file does not exist, we just continue with empty buffer

    // Step 2: Modify the content
    buffer << "\n" + problem->get_name() + "\n"; // append new text

    // Step 3: Loop over all 2D entities
    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities, 2); // dimension 2 = surface elements

    for (const auto &[dim, tag] : entities) {
        std::vector<int> elementTypes;
        std::vector<std::vector<std::size_t>> elementTags;
        std::vector<std::vector<std::size_t>> nodeTagsPerElement;

        gmsh::model::mesh::getElements(elementTypes, elementTags, nodeTagsPerElement, dim, tag);

        for (std::size_t i = 0; i < elementTypes.size(); ++i) {
            int elementType = elementTypes[i];
            std::string name;
            int dimDummy;
            int order;
            int numNodesPerElem;
            std::vector<double> dummyCoords;
            int numPrimaryNodes;

            gmsh::model::mesh::getElementProperties(
                elementType, name, dimDummy, order, numNodesPerElem, dummyCoords, numPrimaryNodes
            );

            //std::cout << "Element type: " << name << " (nodes per face: " << numNodesPerElem << ")\n";
            const auto &nodeListFlat = nodeTagsPerElement[i];

            // Init faces for optimization
            int totalFaceCount = nodeListFlat.size()/4;

            buffer << name << ": " << totalFaceCount << "\n"; 
        }
    }

    // Step 4: Write back (overwrite)
    std::ofstream output_file(path, std::ios::out | std::ios::trunc); // trunc = overwrite
    if (!output_file.is_open()) {
        throw std::runtime_error("Failed to open file for writing: " + path);
    }

    output_file << buffer.str(); // write the updated content
    output_file.close();
}