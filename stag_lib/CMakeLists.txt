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

#####################################
# Define the STAG library
#####################################
# We would like to define two version of the library: shared and static.
# In order to do this, we first build an 'object' library which compiles the
# sources into object files.
add_library(stag_object OBJECT ${SOURCE_FILES} ${HEADER_FILES})
target_compile_options(stag_object PRIVATE -Werror -Wall -Wextra)
set_property(TARGET stag_object PROPERTY POSITION_INDEPENDENT_CODE 1)

# Build the shared library which will be installed
add_library(stag SHARED $<TARGET_OBJECTS:stag_object>)
target_compile_options(stag PRIVATE -Werror -Wall -Wextra)
set_target_properties(stag PROPERTIES PREFIX "")
set_target_properties(stag PROPERTIES PUBLIC_HEADER "${STAG_PUBLIC_HEADERS}")

# Define a static version of the STAG library for the stagtools
add_library(stag_static STATIC $<TARGET_OBJECTS:stag_object>)
target_compile_options(stag_static PRIVATE -Werror -Wall -Wextra)

######################################
# Install
######################################
include(CMakePackageConfigHelpers)

configure_package_config_file("${CMAKE_CURRENT_LIST_DIR}/stagConfig.cmake.in"
                              "${CMAKE_CURRENT_BINARY_DIR}/stagConfig.cmake"
                              INSTALL_DESTINATION lib/stag)

install(TARGETS stag
        LIBRARY DESTINATION lib/stag
        PUBLIC_HEADER DESTINATION include/stag)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/stagConfig.cmake"
        DESTINATION lib/stag)
