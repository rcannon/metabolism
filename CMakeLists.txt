cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
#set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

#set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/PREREQS_INSTALL")

project(metabolism)

find_package(Eigen3 REQUIRED)

set(ifopt_DIR "${PROJECT_SOURCE_DIR}/PREREQS_INSTALL/share/ifopt/cmake")
find_package(ifopt REQUIRED)

#include_directories(${Eigen3_INCLUDE_DIRS})
#include_directories(${ifopt_INCLUDE_DIRS})

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(metabolism src/main.cc)

target_link_libraries(metabolism PUBLIC Eigen3::Eigen)
target_link_libraries(metabolism PUBLIC ifopt::ifopt_ipopt)

# when expanding the project with more files,
# only add the .cc files to this;
# the new header files should automatically
# be accounted for.
target_sources(metabolism PRIVATE 
        ${PROJECT_SOURCE_DIR}/src/ifopt_constraint_class.cc 
        ${PROJECT_SOURCE_DIR}/src/ifopt_cost_class.cc
        ${PROJECT_SOURCE_DIR}/src/ifopt_variable_class.cc
        ${PROJECT_SOURCE_DIR}/src/maximum_entropy_relaxed.cc
        ${PROJECT_SOURCE_DIR}/src/read_data_files.cc
)