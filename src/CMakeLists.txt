cmake_minimum_required(VERSION 3.27)
set (CMAKE_CXX_STANDARD 23)
project(MNOT)
file (GLOB_RECURSE sources *.c *.cpp *.tpp)

# CGAL and its components
find_package( CGAL QUIET COMPONENTS )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )
  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# Eigen
find_package( Eigen3 REQUIRED )

# TinyAD
set(TINYAD_DIR "/home/kateryna/MNOT/libs/TinyAD")
include_directories(${TINYAD_DIR}/src)

include_directories(include tools)

add_executable(mnot ${sources})

set(GMSH_INCLUDE_DIR /home/kateryna/ENERGYOPTIMIZATION/libs/GMSH/api)
find_library(GMSH_LIBRARY gmsh
    PATHS /home/kateryna/ENERGYOPTIMIZATION/libs/GMSH/build
)

# Include directories
target_include_directories(mnot PRIVATE
    ${GMSH_INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIRS}
    ${TINYAD_DIR}/src
)

# Link libraries
target_link_libraries(mnot PRIVATE
    ${GMSH_LIBRARY}
    CGAL::CGAL
)

# target_include_directories(mnot PRIVATE ${EIGEN3_INCLUDE_DIRS})
# target_include_directories(mnot PRIVATE ${TINYAD_DIR}/src)
# target_link_libraries(mnot CGAL::CGAL)
