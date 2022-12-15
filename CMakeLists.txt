cmake_minimum_required(VERSION 3.14)
project(stag_test)

# Require at least C++20
set(CMAKE_CXX_STANDARD 20)

# Add the cmake folder
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Download googletest
include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Find and include the Eigen library
find_package(Eigen3 3.1 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Google test
enable_testing()

# Add the test files to the test object
add_executable(
        stagt
        test/graph_test.cpp
        test/stag_test.cpp
        test/utility_test.cpp
        test/graphio_test.cpp
        test/random_test.cpp
        test/cluster_test.cpp test/spectrum_test.cpp)

# Add the test data files - these need to be copied to the binary directory for the testing application to work
# correctly
configure_file(test/data/test1.edgelist test/data/test1.edgelist COPYONLY)
configure_file(test/data/test2.edgelist test/data/test2.edgelist COPYONLY)
configure_file(test/data/test3.edgelist test/data/test3.edgelist COPYONLY)
configure_file(test/data/test4.edgelist test/data/test4.edgelist COPYONLY)
configure_file(test/data/test5.edgelist test/data/test5.edgelist COPYONLY)

# Add libraries
include_directories(stag_lib)
add_subdirectory(stag_lib stag_lib)

target_link_libraries(
        stagt
        GTest::gtest_main
        stag
)

# Do the googletest magic
include(GoogleTest)
gtest_discover_tests(stagt)

# Build the docs
add_subdirectory("docs")

