cmake_minimum_required(VERSION 3.7)
project(curvature-hessian)
set(CMAKE_CXX_STANDARD 14)

# Max OS X settings
set(CMAKE_OSX_DEPLOYMENT_TARGET 10.13)

# Find libraries
list(APPEND CMAKE_MODULE_PATH
    "${PROJECT_SOURCE_DIR}/external/cmake")
set(MATLAB_ADDITIONAL_VERSIONS "R2019b=9.7" "R2018b=9.5" "R2020a=9.8")
find_package(MATLAB COMPONENTS MEX_COMPILER MX_LIBRARY ENG_LIBRARY)

# Eigen
find_package(Eigen3 3.3.4 REQUIRED NO_MODULE HINTS "external/eigen/cmake")

# Copy over libigl sources to libigl
file(GLOB files "${PROJECT_SOURCE_DIR}/include/libigl/*")
foreach(file ${files})
    get_filename_component(filename ${file} NAME)
    configure_file("${file}" "${PROJECT_SOURCE_DIR}/external/libigl/include/igl/${filename}" COPYONLY)
endforeach()

# libigl
list(APPEND CMAKE_MODULE_PATH
    "${PROJECT_SOURCE_DIR}/external/libigl/cmake")
option(LIBIGL_USE_STATIC_LIBRARY "Use libigl as static library" OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"         ON)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)
include(${PROJECT_SOURCE_DIR}/external/libigl/cmake/libigl.cmake)
set(LINKED_LIBS igl::core)
if(LIBIGL_WITH_OPENGL_GLFW)
    set(LINKED_LIBS ${LINKED_LIBS} igl::opengl_glfw)
endif()
if(LIBIGL_WITH_MATLAB)
    add_definitions("-DMEX")
    link_directories("${Matlab_MEX_COMPILER}/../")
    set(LINKED_LIBS ${LINKED_LIBS} igl::matlab)
endif()

# include code
set(PR_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include/tools" "${LIBIGL_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/include/other_papers" "${PROJECT_SOURCE_DIR}/external/geometry-central/include" "/usr/local/include")

# Static library
if(LIBIGL_USE_STATIC_LIBRARY)
    file(GLOB SOURCES_TOOLS "${PROJECT_SOURCE_DIR}/include/tools/*.cpp")
    add_library("covariant_tools" ${SOURCES_TOOLS})
    target_include_directories(covariant_tools PUBLIC ${PR_INCLUDE_DIRS} ${LIBIGL_SOURCE_DIR})
    target_link_libraries(covariant_tools ${LINKED_LIBS})
    set(LINKED_LIBS ${LINKED_LIBS} covariant_tools)
endif()

# Store location of data directory
set(DATA_PATH "${PROJECT_SOURCE_DIR}/data" CACHE PATH "location of data")
add_definitions(-DDATA_PATH="${DATA_PATH}/")

# Add your projects
add_subdirectory("applications")

