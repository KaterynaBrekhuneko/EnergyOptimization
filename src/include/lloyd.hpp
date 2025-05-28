#include "problem.hpp"
#include "solver.hpp"

class Lloyd : public Solver
{
public:
    void solve(Problem* prob);
};

bool is_on_constraint_lloyd(const Point& steiner, Problem* problem);
bool is_on_boundary_lloyd(Point& p, Polygon& boundary);

template <typename CDT, typename Vertex_handle, typename Vertex_circulator>
Point compute_voronoi_centroid(const Vertex_handle& v, CDT& cdt) {
    std::vector<Point> neighbors;
    Vertex_circulator vc = cdt.incident_vertices(v);
    if (vc != 0) {
        do {
            if (!cdt.is_infinite(vc)) {
                neighbors.push_back(vc->point());
            }
        } while (++vc != cdt.incident_vertices(v));
    }

    if (neighbors.empty()) return v->point();

    double sum_x = 0, sum_y = 0;
    for (const auto& neighbor : neighbors) {
        sum_x += CGAL::to_double(neighbor.x());
        sum_y += CGAL::to_double(neighbor.y());
    }

    return Point(sum_x / neighbors.size(), sum_y / neighbors.size());
}

template <typename CDT, typename Vertex_handle, typename Vertex_circulator>
void iterate_lloyd(CDT& cdt, Problem *problem, std::vector<Segment>& constraints, int iterations) {
    std::vector<Point> points = problem->get_points();
    Polygon boundary = problem->get_boundary();

    std::vector<Segment> all_constraints = constraints;
    for (int i = 0; i < boundary.size(); i++) all_constraints.push_back(Segment(boundary[i], boundary[(i+1)%boundary.size()]));

    for (int iter = 0; iter < iterations; ++iter) {
        std::vector<Point> new_positions;

        for (auto v = cdt.finite_vertices_begin(); v != cdt.finite_vertices_end(); ++v) {
            Point current = v->point();

            if(std::find(points.begin(), points.end(), current) == points.end() && !is_on_boundary_lloyd(current, boundary) && !is_on_constraint_lloyd(current, problem)){
                Point new_point = compute_voronoi_centroid<CDT, Vertex_handle, Vertex_circulator>(v, cdt);
                if(boundary.bounded_side(new_point) != CGAL::ON_UNBOUNDED_SIDE){
                    new_positions.push_back(new_point);
                } else {
                    new_positions.push_back(v->point());
                }
            } else {
                new_positions.push_back(v->point());
            }

        }

        // Clear and reinsert points while preserving constraints
        cdt.clear();
        cdt.insert_constraints(all_constraints.begin(), all_constraints.end());
        cdt.insert(new_positions.begin(), new_positions.end());
    }
}