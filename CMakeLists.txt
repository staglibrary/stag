cmake_minimum_required(VERSION 3.14)
project(stag_test)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp)
add_executable(stagt ${SOURCE_FILES})

# Add libraries
include_directories(eigen-3.3.9) # Eigen - header-only
include_directories(stag_lib)
add_subdirectory(stag_lib)

target_link_libraries(stagt stag)