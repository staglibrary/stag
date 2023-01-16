find_program(SPHINX_EXECUTABLE sphinx-build)
find_program(DOXYGEN_EXECUTABLE doxygen)

#############################################################
# Doxygen build
#############################################################

# Find all the public headers
file(GLOB_RECURSE STAG_PUBLIC_HEADERS ${CMAKE_SOURCE_DIR}/stag_lib/*.h)

set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/stag_lib)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/docs/doxygen)
set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/xml/index.xml)
set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

#Replace variables inside @@ with the current values
configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR}) #Doxygen won't create this for us
add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
                   DEPENDS ${STAG_PUBLIC_HEADERS}
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                   COMMENT "Generating docs with doxygen")

add_custom_target(Doxygen DEPENDS ${DOXYGEN_INDEX_FILE})

#############################################################
# Sphinx build
#############################################################

set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs/sphinx)
set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

# All documentation files should go here
set(DOCUMENTATION_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/index.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/getting-started.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/reference.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/graphs.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/graphio.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/cluster.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/random.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/types.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/utility.rst
    ${CMAKE_CURRENT_SOURCE_DIR}/_static/css/custom.css
    ${CMAKE_CURRENT_SOURCE_DIR}/_static/css/style.css
    ${CMAKE_CURRENT_SOURCE_DIR}/_static/js/custom.js
    ${DOXYGEN_INDEX_FILE}
   )

# Only regenerate Sphinx when:
# - Doxygen has rerun
# - Our doc files have been updated
# - The Sphinx config has been updated
add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                   COMMAND 
                     ${SPHINX_EXECUTABLE} -a -b html
                     # Tell Breathe where to find the Doxygen output
                     -Dbreathe_projects.STAG=${DOXYGEN_OUTPUT_DIR}/xml
                   ${SPHINX_SOURCE} ${SPHINX_BUILD}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   DEPENDS ${DOCUMENTATION_FILES}
                   MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
                   COMMENT "Generating documentation with Sphinx")

#############################################################
# Documentation target
#############################################################
add_custom_target(Documentation ALL DEPENDS ${SPHINX_INDEX_FILE})
