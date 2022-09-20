project(sgtl)

set(HEADER_FILES
        sgtl.h
        )

set(SOURCE_FILES
        sgtl.cpp
        )

# Add libraries
find_package(Eigen3 NO_MODULE REQUIRED)

# Define the SGTL library
add_library(sgtl STATIC ${SOURCE_FILES} ${HEADER_FILES})
