cmake_minimum_required(VERSION 3.14)
project(stag_test)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)

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

# Google test
enable_testing()

# Add the test files to the test object
add_executable(
        stagt
        test/graph_test.cpp
        test/stag_test.cpp)

# Add libraries
include_directories(eigen-3.3.9) # Eigen - header-only
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
