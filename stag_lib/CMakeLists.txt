project(stag)

set(HEADER_FILES
        stag.h
        graph.h
        utility.h
        graphio.h
        random.h
        cluster.h)

set(SOURCE_FILES
        graph.cpp
        utility.cpp
        graphio.cpp
        random.cpp
        cluster.cpp)

# Add libraries
include_directories(../eigen-3.3.9) # Eigen - header-only

# Define the STAG library
add_library(stag SHARED ${SOURCE_FILES} ${HEADER_FILES})

# Create an object file named 'stag.so'.
set_target_properties(stag PROPERTIES PREFIX "")
