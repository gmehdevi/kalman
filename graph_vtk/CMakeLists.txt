cmake_minimum_required(VERSION 3.10)
project(graph_vtk)

# Find the VTK package
find_package(VTK 9.3 REQUIRED)

# Add Boost if required
find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# Create your executable
add_executable(graph_vtk main.cpp)

# Link VTK and Boost libraries
target_link_libraries(graph_vtk ${VTK_LIBRARIES} ${Boost_LIBRARIES})

# Auto-initialize VTK modules
if (VTK_VERSION VERSION_GREATER_EQUAL "8.90")
  vtk_module_autoinit(
    TARGETS graph_vtk
    MODULES ${VTK_LIBRARIES}
  )
endif()
