cmake_minimum_required(VERSION 3.18.0)

project(desal CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wreturn-type -Wempty-body -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wreturn-type -Wempty-body -Wno-deprecated -Wl,-no-as-needed")


find_package(Gperftools REQUIRED)

#set(CMAKE_BULD_TYPE "RelWithDebInfo")

#ADD_DEFINITIONS(-DNDEBUG)

include_directories(include)

list(APPEND src
    "src/point.cpp"
    "src/face.cpp"
    "src/cell.cpp"
    "src/mesh.cpp"
    "src/patch.cpp"
    "src/triangle.cpp"
    "src/gauss_seidel.cpp"
    "src/gradient.cpp"
    "src/read_mesh_gmsh.cpp")

add_library(desal SHARED ${src})
target_link_libraries(desal ${GPERFTOOLS_PROFILER})

enable_testing()
add_subdirectory(unittest)
