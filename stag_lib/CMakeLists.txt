project(stag)

set(HEADER_FILES
        stag.h
        )

set(SOURCE_FILES
        stag.cpp
        )

# Add libraries
find_package(Eigen3 NO_MODULE REQUIRED)

# Define the STAG library
add_library(stag STATIC ${SOURCE_FILES} ${HEADER_FILES})
