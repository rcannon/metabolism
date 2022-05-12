cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
#set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

#set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/INSTALL")

project(metabolism)

find_package(Eigen3 REQUIRED)
#include_directories(${EIGEN_INCLUDE_DIRS})

set(ifopt_DIR "${PROJECT_SOURCE_DIR}/INSTALL/share/ifopt/cmake")
find_package(ifopt REQUIRED)
#include_directories(${ifopt_INCLUDE_DIRS})

include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
)

add_executable(metabolism main.cc)
target_link_libraries(metabolism PUBLIC Eigen3::Eigen)
target_link_libraries(metabolism PUBLIC ifopt::ifopt_ipopt)