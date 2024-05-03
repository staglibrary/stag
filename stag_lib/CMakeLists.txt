project(stag)

set(STAG_PUBLIC_HEADERS
        stag.h
        graph.h
        utility.h
        graphio.h
        random.h
        cluster.h
        spectrum.h
        lsh.h
        kde.h
        definitions.h
        data.h
        )

set(HEADER_FILES
        ${STAG_PUBLIC_HEADERS}
        KMeansRex/KMeansRexCore.h
        KMeansRex/KMeansRexCoreInterface.h
        KMeansRex/mersenneTwister2002.h
        multithreading/ctpl_stl.h
        indicators/indicators.hpp
        )

set(SOURCE_FILES
        utility.cpp
        graph.cpp
        graphio.cpp
        random.cpp
        spectrum.cpp
        cluster.cpp
        lsh.cpp
        kde.cpp
        data.cpp
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

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
message(STATUS "[stag] Found Threads")

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
target_link_libraries(stag PRIVATE Threads::Threads)

# Define a static version of the STAG library for the stagtools
add_library(stag_static STATIC $<TARGET_OBJECTS:stag_object>)
target_link_libraries(stag_static PRIVATE Threads::Threads)

# Configure the compiler options to error on warnings
if(MSVC)
    message(STATUS "[stag] Compiling for Windows.")
    set(_COMPILE_OPTS /W4)
else()
    if (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RELEASE")
        message(STATUS "[stag] Enabling compiler optimisations.")
        set(_COMPILE_OPTS -Wall -Wextra -O3)
    else()
        message(STATUS "[stag] Disabling compiler optimisations.")
        set(_COMPILE_OPTS -Wall -Wextra)
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

# Install the STAG shared library.
install(TARGETS stag
        LIBRARY DESTINATION lib/stag)

# Install the STAG headers, maintaining the directory structure.
install(DIRECTORY "${CMAKE_SOURCE_DIR}/stag_lib/"
        DESTINATION "include/stag"
        FILES_MATCHING PATTERN "*.h")

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/stagConfig.cmake"
        DESTINATION lib/stag)
