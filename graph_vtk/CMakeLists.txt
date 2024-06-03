cmake_minimum_required(VERSION 3.10)
project(graph_vtk)


set(VTK_DIR "/home/gmehdevi/sgoinfre/VTK-9.3.0-build")
find_package(VTK REQUIRED)

set(Boost_INCLUDE_DIR "/home/gmehdevi/code/kalman/include/boost")

include(${VTK_USE_FILE})  
# Create your executable
add_executable(graph_vtk main.cpp)

# Link VTK and Boost libraries
target_link_libraries(graph_vtk ${VTK_LIBRARIES} ${Boost_LIBRARIES})

# Auto-initialize VTK modules
if (VTK_VERSION VERSION_GREATER_EQUAL "9.0")
  vtk_module_autoinit(
    TARGETS graph_vtk
    MODULES ${VTK_LIBRARIES}
  )
endif()
