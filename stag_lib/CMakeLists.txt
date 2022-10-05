project(stag)

set(HEADER_FILES
        stag.h
        utility.h)

set(SOURCE_FILES
        stag.cpp
        utility.cpp)

# Add libraries
find_package(Eigen3 NO_MODULE REQUIRED)

# Define the STAG library
add_library(stag SHARED ${SOURCE_FILES} ${HEADER_FILES})

# Create an object file named 'stag.so'.
set_target_properties(stag PROPERTIES PREFIX "")
