#include "quad_mesher.hpp"

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Boolean_set_operations_2/oriented_side.h>
#include <CGAL/intersections.h>
#include <vector>
#include <gmsh.h>

typedef CGAL::Delaunay_mesh_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;

typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;

typedef CDT::Face_handle Face_handle;
typedef CDT::Edge Edge;

double get_sizing_quad(Problem* problem) {
    Polygon boundary = problem->get_boundary();
    auto bbox = boundary.bbox();
    double max_size = std::sqrt(std::pow(bbox.x_span(), 2) + std::pow(bbox.y_span(), 2)) * 0.1;
    
    return max_size / 2;
}

std::vector<Polygon> generate_quads_medians(Problem* problem, CDT& cdt){
    std::vector<Polygon> quads;

    for ( Face_handle f : cdt.finite_face_handles()) {
        if (!problem->triangle_is_inside<CDT, Face_handle>(cdt, f)) continue;

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

    CDT cdt = problem->generate_CDT<CDT>();
    problem->update_problem<CDT, Face_handle>(cdt, point_set);

    Mesher mesher(cdt);
    double max_size = get_sizing_quad(problem);
    mesher.set_criteria(Criteria(0, max_size*2));
    mesher.refine_mesh();

    auto quads = generate_quads_medians(problem, cdt);

    to_IPE_quad("../solutions/ipe/QUAD-" + problem->get_name() + ".ipe", problem->get_points(), problem->get_constraints(), problem->get_boundary(), {}, quads);
}

std::vector<Polygon> convert_quad_mesh(){
    std::vector<Polygon> quads;
    try {
        std::vector<std::size_t> nodeTags;
        std::vector<double> coords;
        std::vector<double> paramCoords;
        gmsh::model::mesh::getNodes(nodeTags, coords, paramCoords);

        std::unordered_map<std::size_t, std::array<double, 3>> nodeIdToCoord;
        for (std::size_t i = 0; i < nodeTags.size(); ++i) {
            nodeIdToCoord[nodeTags[i]] = {
                coords[3 * i + 0],
                coords[3 * i + 1],
                coords[3 * i + 2]
            };
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

                std::cout << "Element type: " << name << " (nodes per face: " << numNodesPerElem << ")\n";

                const auto &nodeListFlat = nodeTagsPerElement[i];
                for (std::size_t j = 0; j < nodeListFlat.size(); j += numNodesPerElem) {
                    std::cout << "Face:\n";

                    Polygon q;
                    for (int k = 0; k < numNodesPerElem; ++k) {
                        std::size_t nodeId = nodeListFlat[j + k];
                        const auto &xyz = nodeIdToCoord[nodeId];
                        std::cout << "  Vertex " << nodeId << ": ("<< xyz[0] << ", " << xyz[1] << ", " << xyz[2] << ")\n";

                        Point p(xyz[0], xyz[1]);
                        q.push_back(p);
                    }
                    quads.push_back(q);
                }
            }
        }

    } catch (const std::exception &e) {
        std::cerr << "Gmsh error: " << e.what() << "\n";
    }
    gmsh::finalize();

    return quads;
}

void build_quad_mesh_gmsh(Problem* problem){
    gmsh::initialize();
    gmsh::model::add(problem->get_name());

    double lc = get_sizing_quad(problem)*2;  // uniform mesh size

    auto points = problem->get_points();
    auto boundary_indices = problem->get_boundary_indices();

    Eigen::MatrixXd V(points.size(), 2);
    for (int i = 0; i < points.size(); ++i) {
        V(i, 0) = CGAL::to_double(points[i].x()); 
        V(i, 1) = CGAL::to_double(points[i].y()); 
    }

    for(int i = 0; i < points.size(); ++i){
        int p = gmsh::model::geo::addPoint(V(i, 0), V(i, 1), 0, lc);
    }

    // Add additional point to make the num of triangles in the triangle mesch even
    double additional_x = 0;
    double additional_y = 982072;
    int p = gmsh::model::geo::addPoint(additional_x, additional_y, 0, lc);
    boundary_indices.push_back(p-1);

    std::vector<int> boundary;
    for(int i = 0; i < boundary_indices.size(); ++i){
        int p1 = boundary_indices[i] + 1;
        int p2 = boundary_indices[(i+1)%boundary_indices.size()] + 1;
        int l = gmsh::model::geo::addLine(p1, p2);
        boundary.push_back(l);
    }
    int cl = gmsh::model::geo::addCurveLoop(boundary);
    int pl = gmsh::model::geo::addPlaneSurface({cl});

    gmsh::model::geo::synchronize();

    gmsh::option::setNumber("Mesh.Optimize", 0);
    gmsh::option::setNumber("Mesh.OptimizeNetgen", 0);
    gmsh::option::setNumber("Mesh.Smoothing", 0);

    gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 1);
    gmsh::model::mesh::setRecombine(2, pl);

    gmsh::model::mesh::generate(2);
    gmsh::write("aaa.msh");
    
    std::vector<Polygon> quads = convert_quad_mesh();
    to_IPE_quad("../solutions/ipe/QUAD-" + problem->get_name() + ".ipe", problem->get_points(), problem->get_constraints(), problem->get_boundary(), {}, quads);

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