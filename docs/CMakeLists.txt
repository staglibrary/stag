find_program(DOXYGEN_EXECUTABLE doxygen)

#############################################################
# Doxygen build
#############################################################

# Find all the public headers
file(GLOB_RECURSE STAG_PUBLIC_HEADERS ${CMAKE_SOURCE_DIR}/stag_lib/*.h)

set(DOXYGEN_DOCS_DIR ${PROJECT_SOURCE_DIR}/docs)
set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/stag_lib)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/xml/index.xml)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

#Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us
add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                   DEPENDS ${STAG_PUBLIC_HEADERS} ${DOXYGEN_DOCS_DIR}
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                   COMMENT "Generating docs with doxygen")

add_custom_target(Doxygen DEPENDS ${DOXYGEN_INDEX_FILE})

#############################################################
# Documentation target
#############################################################
add_custom_target(Documentation DEPENDS ${DOXYGEN_INDEX_FILE})
