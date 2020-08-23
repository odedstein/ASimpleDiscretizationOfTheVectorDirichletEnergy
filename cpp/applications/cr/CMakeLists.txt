project(cr)


# Setup MEX compilation, as in gptoolbox
set( CMAKE_SHARED_LIBRARY_SUFFIX ".${Matlab_MEX_EXTENSION}" )  # set suffix to .mexa64
set( CMAKE_SHARED_LIBRARY_PREFIX )      # remove the "lib" prefix
set( CMAKE_SHARED_LIBRARY_CXX_FLAGS )       # remove the -fPIC option. mex does not accept the "-fPIC" option
set( CMAKE_POSITION_INDEPENDENT_CODE ON)
set( CMAKE_INSTALL_RPATH "\$ORIGIN" )
add_definitions(-DPROJECT_NAME="${PROJECT_NAME}")

# Compile cr_vector_laplacian_and_mass
set(SRCFILES "cr-vector-laplacian-and-mass.cpp")
list(APPEND SRCFILES ${GLOBAL_CPPS})
add_library("cr_vector_laplacian_and_mass" SHARED ${SRCFILES})
target_include_directories("cr_vector_laplacian_and_mass" PRIVATE ${PR_INCLUDE_DIRS})
target_link_libraries("cr_vector_laplacian_and_mass" PUBLIC ${LINKED_LIBS})
configure_file("cr_vector_laplacian_and_mass.m" "${CMAKE_CURRENT_BINARY_DIR}/cr_vector_laplacian_and_mass.m" COPYONLY)

# Compile cr_vector_laplacian
set(SRCFILES "cr-vector-laplacian.cpp")
list(APPEND SRCFILES ${GLOBAL_CPPS})
add_library("cr_vector_laplacian" SHARED ${SRCFILES})
target_include_directories("cr_vector_laplacian" PRIVATE ${PR_INCLUDE_DIRS})
target_link_libraries("cr_vector_laplacian" PUBLIC ${LINKED_LIBS})
configure_file("cr_vector_laplacian.m" "${CMAKE_CURRENT_BINARY_DIR}/cr_vector_laplacian.m" COPYONLY)

# Compile cr_twisted_vector_laplacian
set(SRCFILES "cr-twisted-vector-laplacian.cpp")
list(APPEND SRCFILES ${GLOBAL_CPPS})
add_library("cr_twisted_vector_laplacian" SHARED ${SRCFILES})
target_include_directories("cr_twisted_vector_laplacian" PRIVATE ${PR_INCLUDE_DIRS})
target_link_libraries("cr_twisted_vector_laplacian" PUBLIC ${LINKED_LIBS})
configure_file("cr_twisted_vector_laplacian.m" "${CMAKE_CURRENT_BINARY_DIR}/cr_twisted_vector_laplacian.m" COPYONLY)

# Compile cr_vector_mass
set(SRCFILES "cr-vector-mass.cpp")
list(APPEND SRCFILES ${GLOBAL_CPPS})
add_library("cr_vector_mass" SHARED ${SRCFILES})
target_include_directories("cr_vector_mass" PRIVATE ${PR_INCLUDE_DIRS})
target_link_libraries("cr_vector_mass" PUBLIC ${LINKED_LIBS})
configure_file("cr_vector_mass.m" "${CMAKE_CURRENT_BINARY_DIR}/cr_vector_mass.m" COPYONLY)

# Compile cr_shape_operator
set(SRCFILES "cr-shape-operator.cpp")
list(APPEND SRCFILES ${GLOBAL_CPPS})
add_library("cr_shape_operator" SHARED ${SRCFILES})
target_include_directories("cr_shape_operator" PRIVATE ${PR_INCLUDE_DIRS})
target_link_libraries("cr_shape_operator" PUBLIC ${LINKED_LIBS})
configure_file("cr_shape_operator.m" "${CMAKE_CURRENT_BINARY_DIR}/cr_shape_operator.m" COPYONLY)

# Compile cr_edges
set(SRCFILES "cr-edges.cpp")
list(APPEND SRCFILES ${GLOBAL_CPPS})
add_library("cr_edges" SHARED ${SRCFILES})
target_include_directories("cr_edges" PRIVATE ${PR_INCLUDE_DIRS})
target_link_libraries("cr_edges" PUBLIC ${LINKED_LIBS})
configure_file("cr_edges.m" "${CMAKE_CURRENT_BINARY_DIR}/cr_edges.m" COPYONLY)

