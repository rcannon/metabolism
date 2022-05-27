cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# comment line below for production
set(CMAKE_BUILD_TYPE Debug)

project(metabolism)

find_package(Eigen3 REQUIRED)

set(ifopt_DIR "${PROJECT_SOURCE_DIR}/PREREQS_INSTALL/share/ifopt/cmake")
find_package(ifopt REQUIRED)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(run src/main.cc)
#add_executable(test src/metabolism_test.cc)

target_link_libraries(run PUBLIC Eigen3::Eigen)
target_link_libraries(run PUBLIC ifopt::ifopt_ipopt)

#target_link_libraries(test PUBLIC Eigen3::Eigen)
#target_link_libraries(test PUBLIC ifopt::ifopt_ipopt)

# when expanding the project with more files,
# only add the .cc files to this;
# the new header files should automatically
# be accounted for.
target_sources(run PRIVATE 
        ${PROJECT_SOURCE_DIR}/src/ifopt_constraint_classes.cc 
        ${PROJECT_SOURCE_DIR}/src/ifopt_cost_class.cc
        ${PROJECT_SOURCE_DIR}/src/ifopt_variable_class.cc
        ${PROJECT_SOURCE_DIR}/src/maximum_entropy_relaxed.cc
        ${PROJECT_SOURCE_DIR}/src/read_data_files.cc
        ${PROJECT_SOURCE_DIR}/src/metabolism.cc
        ${PROJECT_SOURCE_DIR}/src/find_initial_values.cc
        ${PROJECT_SOURCE_DIR}/src/write_vector_to_csv.cc
)

#target_sources(test PRIVATE 
#        ${PROJECT_SOURCE_DIR}/src/ifopt_constraint_classes.cc 
#        ${PROJECT_SOURCE_DIR}/src/ifopt_cost_class.cc
#        ${PROJECT_SOURCE_DIR}/src/ifopt_variable_class.cc
#        ${PROJECT_SOURCE_DIR}/src/maximum_entropy_relaxed.cc
#        ${PROJECT_SOURCE_DIR}/src/read_data_files.cc
#        ${PROJECT_SOURCE_DIR}/src/metabolism.cc
#        ${PROJECT_SOURCE_DIR}/src/find_initial_values.cc
#)