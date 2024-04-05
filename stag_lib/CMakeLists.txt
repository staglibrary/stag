project(stag)

set(STAG_PUBLIC_HEADERS
        stag.h
        Graph/graph.h
        utility.h
        Graph/graphio.h
        Graph/random.h
        Cluster/cluster.h
        Graph/spectrum.h
        KDE/lsh.h
        KDE/kde.h
        definitions.h
        KDE/data.h
        KDE/LSHTable.h
        )

set(HEADER_FILES
        ${STAG_PUBLIC_HEADERS}
        KMeansRex/KMeansRexCore.h
        KMeansRex/KMeansRexCoreInterface.h
        KMeansRex/mersenneTwister2002.h
        multithreading/ctpl_stl.h
        )

set(SOURCE_FILES
        utility.cpp
        Graph/graph.cpp
        Graph/graphio.cpp
        Graph/random.cpp
        Graph/spectrum.cpp
        Cluster/cluster.cpp
        KDE/lsh.cpp
        KDE/LSHTable.cpp
        KDE/kde.cpp
        KDE/data.cpp
        KMeansRex/KMeansRexCore.cpp
        )

# Generate correct files on windows
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)

# Find and include the Eigen library
find_package(Eigen3 3.1 REQUIRED)
message(STATUS "[stag] Found Eigen: ${EIGEN3_INCLUDE_DIR}")
include_directories(${EIGEN3_INCLUDE_DIR})

# Find and include the Spectra library
find_package(Spectra 1.0.1 REQUIRED)
get_property(SPECTRA_INCLUDE_DIR TARGET Spectra::Spectra
             PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "[stag] Found Spectra: ${SPECTRA_INCLUDE_DIR}")
include_directories(${SPECTRA_INCLUDE_DIR})

#####################################
# Define the STAG library
#####################################
# We would like to define two version of the library: shared and static.
# In order to do this, we first build an 'object' library which compiles the
# sources into object files.
add_library(stag_object OBJECT ${SOURCE_FILES} ${HEADER_FILES})
set_property(TARGET stag_object PROPERTY POSITION_INDEPENDENT_CODE 1)

# Build the shared library which will be installed
add_library(stag SHARED $<TARGET_OBJECTS:stag_object>)
set_target_properties(stag PROPERTIES PREFIX "")
set_target_properties(stag PROPERTIES PUBLIC_HEADER "${STAG_PUBLIC_HEADERS}")

# Define a static version of the STAG library for the stagtools
add_library(stag_static STATIC $<TARGET_OBJECTS:stag_object>)

# Configure the compiler options to error on warnings
if(MSVC)
    message(STATUS "[stag] Compiling for Windows.")
    set(_COMPILE_OPTS /W4)
else()
    if (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RELEASE")
        message(STATUS "[stag] Enabling compiler optimisations.")
        set(_COMPILE_OPTS -Werror -Wall -Wextra -O3)
    else()
        message(STATUS "[stag] Disabling compiler optimisations.")
        set(_COMPILE_OPTS -Werror -Wall -Wextra)
    endif()
endif()
message(STATUS "[stag] Set compiler options: ${_COMPILE_OPTS}")
target_compile_options(stag_object PRIVATE "${_COMPILE_OPTS}")
target_compile_options(stag PRIVATE "${_COMPILE_OPTS}")
target_compile_options(stag_static PRIVATE "${_COMPILE_OPTS}")

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
