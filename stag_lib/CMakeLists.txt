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

# Define the STAG library
add_library(stag SHARED ${SOURCE_FILES} ${HEADER_FILES})
set_target_properties(stag PROPERTIES PREFIX "")

#################
# Add libraries #
#################
# Eigen - header-only
include_directories(../eigen-3.3.9)

# libneo4j-client
set(libneo4j-build-dir "${CMAKE_BINARY_DIR}/libneo4j-client/")
set(libneo4j-lib-file "${libneo4j-build-dir}/lib/libneo4j-client.so")
add_custom_command(
        OUTPUT ${libneo4j-lib-file}
        COMMAND ./autogen.sh
        COMMAND ./configure --disable-tools --prefix=${libneo4j-build-dir}
        COMMAND make install
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/libneo4j-client/
        )
add_custom_target(build-libneo4j DEPENDS ${libneo4j-lib-file})
add_library(libneo4j-client SHARED IMPORTED)
set_target_properties(libneo4j-client PROPERTIES
        IMPORTED_LOCATION ${libneo4j-lib-file}
        INTERFACE_INCLUDE_DIRECTORIES "${libneo4j-build-dir}/include/")
add_dependencies(libneo4j-client build-libneo4j)

# The stag library depends on the libneo4j build
add_dependencies(stag libneo4j-client)
target_link_libraries(stag libneo4j-client)
