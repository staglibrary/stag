project(stag)

set(HEADER_FILES
        stag.h
        graph.h
        utility.h
        graphio.h
        random.h
        cluster.h
        spectrum.h
        KMeansRex/KMeansRexCore.h
        KMeansRex/KMeansRexCoreInterface.h
        )

set(SOURCE_FILES
        graph.cpp
        utility.cpp
        graphio.cpp
        random.cpp
        cluster.cpp
        spectrum.cpp
        KMeansRex/KMeansRexCore.cpp
        KMeansRex/mersenneTwister2002.c
        )

# Find and include the Eigen library
find_package(Eigen3 3.1 REQUIRED)
message(STATUS "[stag] Found Eigen!")
include_directories(${EIGEN3_INCLUDE_DIR})

# Find and include the Spectra library
find_package(Spectra 1.0.1 REQUIRED)
message(STATUS "[stag] Found Spectra!")
include_directories(${SPECTRA_INCLUDE_DIR})

# Define the STAG library
add_library(stag SHARED ${SOURCE_FILES} ${HEADER_FILES})
set_target_properties(stag PROPERTIES PREFIX "")
