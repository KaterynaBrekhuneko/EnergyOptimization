# Mesh Optimization for PSLGs

This project provides a framework for generating and optimizing high-quality triangular and quadrilateral meshes for planar straight-line graphs (PSLGs), with a particular focus on eliminating obtuse angles using angle-based energy functions. It was developed in the context of research related to the [2025 CG:SHOP Challenge](https://cgshop.ibr.cs.tu-bs.de/competition/cg-shop-2025/#problem-description).

## Build Instructions

### Dependencies

- C++17
- [CGAL](https://www.cgal.org/) (version ≥ 3.28)
- [Boost](https://www.boost.org/) (version ≥ 1.83)
- [Eigen3](https://eigen.tuxfamily.org/) (version ≥ 3.4.0)
- [nlohmann/json](https://github.com/nlohmann/json)
- [TinyAD](https://github.com/patr-schm/TinyAD) (for newton-based optimization)

### Compile and run (with CMake)

```bash
mkdir build && cd build
cmake ../src
cmake --build .

./energy

