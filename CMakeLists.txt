cmake_minimum_required(VERSION 3.14)
project(stag_test)

# Require at least C++20
set(CMAKE_CXX_STANDARD 20)

# Add the cmake folder
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Default to release build if CMAKE_BUILD_TYPE is not configured
if (CMAKE_BUILD_TYPE STREQUAL "")
    message(STATUS "[stag] Defaulting to release build")
    set(CMAKE_BUILD_TYPE "Release")
else()
    message(STATUS "[stag] Build type: ${CMAKE_BUILD_TYPE}")
endif()

# Download googletest
include(FetchContent)
FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.12.1
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Make google test available, and disable testing
FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
    FetchContent_Populate(googletest)
    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# Find and include the Eigen library
find_package(Eigen3 3.1 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# Find and include the Spectra library
find_package(Spectra 1.0.1 REQUIRED)
get_property(SPECTRA_INCLUDE_DIR TARGET Spectra::Spectra
        PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
include_directories(${SPECTRA_INCLUDE_DIR})

# Google test
enable_testing()

# Find and include the STAG library
include_directories(stag_lib)
add_subdirectory(stag_lib stag_lib)

# Add the test files to the test object
add_executable(
        stagt EXCLUDE_FROM_ALL
        test/graph_test.cpp
        test/stag_test.cpp
        test/utility_test.cpp
        test/graphio_test.cpp
        test/random_test.cpp
        test/cluster_test.cpp
        test/spectrum_test.cpp)

target_link_libraries(
        stagt
        GTest::gtest_main
        stag
)

add_executable(
        stagt_static EXCLUDE_FROM_ALL
        test/graph_test.cpp
        test/stag_test.cpp
        test/utility_test.cpp
        test/graphio_test.cpp
        test/random_test.cpp
        test/cluster_test.cpp
        test/spectrum_test.cpp)

target_link_libraries(
        stagt_static
        GTest::gtest_main
        stag_static
)

# Add the test data files - these need to be copied to the binary directory for the testing application to work
# correctly
configure_file(test/data/test1.edgelist test/data/test1.edgelist COPYONLY)
configure_file(test/data/test2.edgelist test/data/test2.edgelist COPYONLY)
configure_file(test/data/test3.edgelist test/data/test3.edgelist COPYONLY)
configure_file(test/data/test4.edgelist test/data/test4.edgelist COPYONLY)
configure_file(test/data/test5.edgelist test/data/test5.edgelist COPYONLY)
configure_file(test/data/test6.edgelist test/data/test6.edgelist COPYONLY)
configure_file(test/data/test1.adjacencylist test/data/test1.adjacencylist COPYONLY)
configure_file(test/data/test2.adjacencylist test/data/test2.adjacencylist COPYONLY)
configure_file(test/data/test3.adjacencylist test/data/test3.adjacencylist COPYONLY)
configure_file(test/data/test6.adjacencylist test/data/test6.adjacencylist COPYONLY)
configure_file(test/data/test7.adjacencylist test/data/test7.adjacencylist COPYONLY)
configure_file(test/data/test8.adjacencylist test/data/test8.adjacencylist COPYONLY)
configure_file(test/data/badgraph.adjacencylist test/data/badgraph.adjacencylist COPYONLY)
configure_file(test/data/hugegraph.adjacencylist test/data/hugegraph.adjacencylist COPYONLY)

# Create the STAG tools
add_executable(
        stag_edge2adj
        stagtools/edge2adj.cpp
)
target_link_libraries(
        stag_edge2adj
        stag_static
)

add_executable(
        stag_adj2edge
        stagtools/adj2edge.cpp
)

target_link_libraries(
        stag_adj2edge
        stag_static
)

add_executable(
        stag_sbm
        stagtools/sbm.cpp
)

target_link_libraries(
        stag_sbm
        stag_static
)

# Do the googletest magic
include(GoogleTest)
gtest_discover_tests(stagt)

# Build the docs
add_subdirectory("docs")

# Configure the install target for the stag tools.
install(TARGETS stag_edge2adj
        RUNTIME DESTINATION bin)

install(TARGETS stag_adj2edge
        RUNTIME DESTINATION bin)

install(TARGETS stag_sbm
        RUNTIME DESTINATION bin)
