cmake_minimum_required(VERSION 3.10)
project(rule_based_optimal_control)

# Use C++17 (change if needed)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add the include/ directory to the include path
include_directories(include)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

set(IPOPT_ROOT "/usr/local/")
include_directories(${IPOPT_ROOT}/include/coin-or)
link_directories(${IPOPT_ROOT}/lib)

add_executable(${PROJECT_NAME}
    src/rule_based_optimal_control.cpp
    include/rule_based_optimal_control.h
)
target_link_libraries(${PROJECT_NAME} ipopt)