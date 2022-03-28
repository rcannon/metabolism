cmake_minimum_required(VERSION 3.0.0)
project(metab VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
list(APPEND CMAKE_INCLUDE_PATH "include")

include(CTest)
enable_testing()

add_subdirectory(eigen)
add_subdirectory(ifopt)

#add_library(metab metab.cpp)
add_executable(Feasible_point_relaxed feasible_point_relaxed.cc)
add_executable(Maximum_entropy_relaxed maximum_entropy_relaxed.cc)


target_link_libraries(Feasible_point_relaxed PUBLIC eigen ifopt)
target_link_libraries(Maximum_entropy_relaxed PUBLIC eigen ifopt)

target_include_directories(Feasible_point_relaxed PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/eigen"
                          "${PROJECT_SOURCE_DIR}/ifopt"
                          "/share/apps/ipopt/3.14"
                          )

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
