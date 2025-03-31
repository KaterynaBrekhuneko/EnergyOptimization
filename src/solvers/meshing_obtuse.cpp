#include "meshing_obtuse.hpp"

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_vertex_base_2.h>
#include <CGAL/Boolean_set_operations_2/oriented_side.h>
#include <CGAL/enum.h>
#include <CGAL/Double_map.h>

typedef CGAL::Delaunay_mesh_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;

typedef CDT::Geom_traits Geom_traits;
typedef CDT::Face_handle Face_handle;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Edge Edge;

typedef CGAL::Double_map<Face_handle, std::pair<double, double>, ScoreComparator> Bad_Faces;

Problem* problem;
CDT cdt;
Bad_Faces bad_faces;
std::deque<Edge> bad_constraints;

// What CGAL uses
struct Face_compare {
    bool operator()(const Face_handle& fh1, const Face_handle& fh2) const {
      if(fh1->vertex(0)->point() < fh2->vertex(0)->point())
        return true;
      else if(fh1->vertex(0)->point() == fh2->vertex(0)->point()) {
        if(fh1->vertex(1)->point() < fh2->vertex(1)->point())
          return true;
        else if(fh1->vertex(1)->point() == fh2->vertex(1)->point() &&
                fh1->vertex(2)->point() < fh2->vertex(2)->point())
          return true;
      }
      return false;
    }
};

struct ScoreComparator {
    bool operator()(const std::pair<double, double>& a, const std::pair<double, double>& b) const {
        // Sort descending by angle score (higher obtuseness first)
        if (a.first != b.first)
            return a.first > b.first;

        // If angle scores are equal, sort by side length descending
        return a.second > b.second;
    }
};

bool triangle_is_inside(const Polygon& boundary, const Face_handle& triangle){
    Vector a = Vector(0, 0);
    for (int i = 0; i < 3; i++) {
        CDT::Point p = cdt.point(triangle->vertex(i));
        a += Vector(p.x(), p.y());
    }
    Point c = Point(a.x() / 3, a.y() / 3);
    return CGAL::oriented_side(c, boundary) == CGAL::POSITIVE;
}

bool is_encroached(const Edge& edge){
    const Face_handle fh = edge.first;
    int i = edge.second;

    const Vertex_handle& va = fh->vertex(cdt. cw(i));
    const Vertex_handle& vb = fh->vertex(cdt.ccw(i));

    const Vertex_handle& vi = fh->vertex(i);
    const Vertex_handle& mvi = cdt.tds().mirror_vertex(fh, i);

    // Chech whether vi or mvi are encroached by edge
    typedef typename Geom_traits::Angle_2 Angle_2;

    const Angle_2 angle = cdt.geom_traits().angle_2_object();

    const Point& a = va->point();
    const Point& b = vb->point();

    return  (cdt.is_infinite(vi) || (angle(a, vi->point(), b) == CGAL::ACUTE))
            &&
            (cdt.is_infinite(mvi) || (angle(a, mvi->point(), b) == CGAL::ACUTE));
}

// Adds encroached constrained edges to the queue
void find_bad_constraints(){
    for(const Edge& edge : cdt.constrained_edges()){
        if(is_encroached(edge)){
            bad_constraints.push_back(edge);
        }
    }
}

// Adds obtuse triangles and triangles with a too large side length to the queue
void find_bad_faces(double max_length){
    for(const auto& fh : cdt.finite_face_handles()){
        if(triangle_is_inside(problem->get_boundary(), fh)){
            const Point& pa = fh->vertex(0)->point();
            const Point& pb = fh->vertex(1)->point();
            const Point& pc = fh->vertex(2)->point();

            double
                a = CGAL::to_double(CGAL::squared_distance(pb, pc)),
                b = CGAL::to_double(CGAL::squared_distance(pc, pa)),
                c = CGAL::to_double(CGAL::squared_distance(pa, pb));
            
            double max_sq = std::max({a, b, c});
            double obtuseness = max_sq - (a + b + c - max_sq);
            double side_length_ratio = max_sq / max_sq;

            std::pair<double, double> quality = {obtuseness, side_length_ratio};

            if(side_length_ratio >= 1 || obtuseness > 1){
                bad_faces.insert(fh, quality);
            }
        }
    }
}

//! max_length has to be squared
void init(Problem* prob, double max_length){
    problem = prob;
    cdt = problem->generate_CDT<CDT>();
    find_bad_constraints();
    find_bad_faces(max_length);
}

bool step_by_step_refine(){
    // Break if there are no elements to refine
    if(!is_refinement_possible()) return false;

    if(!bad_constraints.empty()){
        Edge& constraint = bad_constraints.back();
        bad_constraints.pop_back();
        refine_constraint(constraint);
    } else {
        auto it = bad_faces.front();
        auto entry = *it.base();
        std::pair<double, double> quality = entry.first;
        Face_handle fh = entry.second;

        bad_faces.pop_front();

        refine_face(fh, quality);
    }
}

void refine_constraint(Edge& constraint){

}

void refine_face(Face_handle fh, std::pair<double, double> quality){

}

bool is_refinement_possible(){
    return (bad_constraints.empty() || bad_faces.empty());
}
