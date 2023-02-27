project(stag)

set(STAG_PUBLIC_HEADERS
        stag.h
        graph.h
        utility.h
        graphio.h
        random.h
        cluster.h
        spectrum.h
        )

set(HEADER_FILES
        ${STAG_PUBLIC_HEADERS}
        KMeansRex/KMeansRexCore.h
        KMeansRex/KMeansRexCoreInterface.h
        KMeansRex/mersenneTwister2002.h
        )

set(SOURCE_FILES
        graph.cpp
        utility.cpp
        graphio.cpp
        random.cpp
        cluster.cpp
        spectrum.cpp
        KMeansRex/KMeansRexCore.cpp
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
set_target_properties(stag PROPERTIES PUBLIC_HEADER "${STAG_PUBLIC_HEADERS}")

# Configure the install target
include(CMakePackageConfigHelpers)

configure_package_config_file("${CMAKE_CURRENT_LIST_DIR}/stagConfig.cmake.in"
                              "${CMAKE_CURRENT_BINARY_DIR}/stagConfig.cmake"
                              INSTALL_DESTINATION lib/stag)

install(TARGETS stag
        LIBRARY DESTINATION lib/stag
        PUBLIC_HEADER DESTINATION include/stag)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/stagConfig.cmake"
        DESTINATION lib/stag)
