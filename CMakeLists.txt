cmake_minimum_required(VERSION 3.14)
project(sgtl_test)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp)
add_executable(test ${SOURCE_FILES})

# Add libraries
include_directories(eigen-3.3.9) # Eigen - header-only
include_directories(sgtl_lib)
add_subdirectory(sgtl_lib)

target_link_libraries(test sgtl)