// -----------------------------------------------------------------------------
//
//  Gmsh C++ tutorial 11
//
//  Unstructured quadrangular meshes
//
// -----------------------------------------------------------------------------

#include <set>
#include <gmsh.h>

#include <iostream>

int main(int argc, char **argv)
{
  gmsh::initialize();

  gmsh::model::add("t11");

  // We have seen in tutorials `t3.cpp' and `t6.cpp' that extruded and
  // transfinite meshes can be "recombined" into quads, prisms or
  // hexahedra. Unstructured meshes can be recombined in the same way. Let's
  // define a simple geometry with an analytical mesh size field:

  double lc = 0.2;  // uniform mesh size

  int p1 = gmsh::model::geo::addPoint(-1.25, -.5, 0, lc);
  int p2 = gmsh::model::geo::addPoint(1.25, -.5, 0, lc);
  int p3 = gmsh::model::geo::addPoint(1.25, 1.25, 0, lc);
  int p4 = gmsh::model::geo::addPoint(-1.25, 1.25, 0, lc);

  std::cout << "points: " << p1 << " " << p2 << " " << p3 << " " << p4 << " " << std::endl;

  int l1 = gmsh::model::geo::addLine(p1, p2);
  int l2 = gmsh::model::geo::addLine(p2, p3);
  int l3 = gmsh::model::geo::addLine(p3, p4);
  int l4 = gmsh::model::geo::addLine(p4, p1);

  int cl = gmsh::model::geo::addCurveLoop({l1, l2, l3, l4});
  int pl = gmsh::model::geo::addPlaneSurface({cl});

  gmsh::model::geo::synchronize();

  /*gmsh::model::mesh::field::add("MathEval", 1);
  gmsh::model::mesh::field::setString(
    1, "F", "0.01*(1.0+30.*(y-x*x)*(y-x*x) + (1-x)*(1-x))");
  gmsh::model::mesh::field::setAsBackgroundMesh(1);*/

  //* To generate quadrangles instead of triangles, we can simply add
  //*gmsh::model::mesh::setRecombine(2, pl);

  // If we'd had several surfaces, we could have used the global option
  // "Mesh.RecombineAll":
  //
  //gmsh::option::setNumber("Mesh.RecombineAll", 1);

  // The default recombination algorithm is called "Blossom": it uses a minimum
  // cost perfect matching algorithm to generate fully quadrilateral meshes from
  // triangulations. More details about the algorithm can be found in the
  // following paper: J.-F. Remacle, J. Lambrechts, B. Seny, E. Marchandise,
  // A. Johnen and C. Geuzaine, "Blossom-Quad: a non-uniform quadrilateral mesh
  // generator using a minimum cost perfect matching algorithm", International
  // Journal for Numerical Methods in Engineering 89, pp. 1102-1119, 2012.

  // For even better 2D (planar) quadrilateral meshes, you can try the
  // experimental "Frontal-Delaunay for quads" meshing algorithm, which is a
  // triangulation algorithm that enables to create right triangles almost
  // everywhere: J.-F. Remacle, F. Henrotte, T. Carrier-Baudouin, E. Bechet,
  // E. Marchandise, C. Geuzaine and T. Mouton. A frontal Delaunay quad mesh
  // generator using the L^inf norm. International Journal for Numerical Methods
  // in Engineering, 94, pp. 494-512, 2013. Uncomment the following line to try
  // the Frontal-Delaunay algorithms for quads:
  //
  // gmsh::option::setNumber("Mesh.Algorithm", 8);

  // The default recombination algorithm might leave some triangles in the mesh,
  // if recombining all the triangles leads to badly shaped quads. In such
  // cases, to generate full-quad meshes, you can either subdivide the resulting
  // hybrid mesh (with `Mesh.SubdivisionAlgorithm' set to 1), or use the
  // full-quad recombination algorithm, which will automatically perform a
  // coarser mesh followed by recombination, smoothing and
  // subdivision. Uncomment the following line to try the full-quad algorithm:
  //
  // gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 2); // or 3

  // You can also set the subdivision step alone, with
  //
  //gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 1);

  //gmsh::model::mesh::refine();

  //gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 1);
  //gmsh::model::mesh::setRecombine(2, pl);

  //gmsh::model::mesh::setAlgorithm(2, 8);
  gmsh::option::setNumber("Mesh.Optimize", 0);
  gmsh::option::setNumber("Mesh.OptimizeNetgen", 0);
  gmsh::option::setNumber("Mesh.Smoothing", 0);

  gmsh::option::setNumber("Mesh.RecombinationAlgorithm", 3);
  gmsh::model::mesh::setRecombine(2, pl);

  gmsh::model::mesh::generate(2);
  gmsh::write("t11.msh");

  // Note that you could also apply the recombination algorithm and/or the
  // subdivision step explicitly after meshing, as follows:
  //
  /*gmsh::model::mesh::generate(2);
  gmsh::model::mesh::recombine();
  gmsh::option::setNumber("Mesh.SubdivisionAlgorithm", 1);
  gmsh::model::mesh::refine();
  gmsh::write("t11.msh");*/

  // Launch the GUI to see the results:
  std::set<std::string> args(argv, argv + argc);
  if(!args.count("-nopopup")) gmsh::fltk::run();

  gmsh::finalize();
  return 0;
}
