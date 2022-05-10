cmake_minimum_required(VERSION 3.0.0)
project(metab VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
list(APPEND CMAKE_INCLUDE_PATH "include")

include(CTest)
enable_testing()

add_subdirectory(eigen)
add_subdirectory(ifopt)

add_executable(Run_maximum_entropy run_maximum_entropy.cc)

target_link_libraries(Run_maximum_entropy PUBLIC eigen ifopt)

target_include_directories(Run_maximum_entropy PUBLIC
                          "${PROJECT_BINARY_DIR}"
                          "${PROJECT_SOURCE_DIR}/eigen"
                          "${PROJECT_SOURCE_DIR}/ifopt"
                          "/share/apps/ipopt/3.14"
                          )

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
