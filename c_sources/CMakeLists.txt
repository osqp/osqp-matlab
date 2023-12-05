cmake_minimum_required( VERSION 3.18 )

project( osqp-matlab )


# Find the MATLAB installation directory and various settings for it
find_package( Matlab REQUIRED COMPONENTS MEX_COMPILER)

message( STATUS "Matlab root is " ${Matlab_ROOT_DIR} )

include_directories( ${Matlab_INCLUDE_DIRS} )

# The mex interface uses C++11
set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

# Some options
option( MATLAB_CUDA "Build CUDA backend" OFF )

if( CMAKE_COMPILER_IS_GNUCXX )
    # Add debug symbols and optimizations to the build
    set( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O2" )
    set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2" )
endif()

# Insist on the pre-2018 complex data API so that mxGetPr will work correctly
add_compile_definitions( MATLAB_MEXSRC_RELEASE=R2017b )
message( STATUS "Using Matlab pre-2018a API for mxGetPr compatibility" )

# Some parts of the main libraries need to know we are building for MATLAB
# So pass a definition to let them know
add_compile_definitions( MATLAB )

# Configure the settings of the OSQP library
if( MATLAB_CUDA )
    message( STATUS "Building CUDA-based mex file" )
    set( OSQP_ALGEBRA_BACKEND    "cuda" CACHE STRING "Build CUDA mex file" )
    set( OSQP_CODEGEN            OFF CACHE BOOL "Disable codegen" )
    set( OSQP_ENABLE_DERIVATIVES OFF CACHE BOOL "Disable derivatives" )

    # Give CMake hints for where to find libraries
    set( CUDAToolkit_ROOT "${MATLAB_CUDA_ROOT}" CACHE STRING "CUDA Toolkit root directory" )
    set( CUDA_TOOLKIT_ROOT_DIR "${MATLAB_CUDA_ROOT}" CACHE STRING "CUDA Toolkit root directory" )
    set( CUDA_HOME "${MATLAB_CUDA_ROOT}" CACHE STRING "CUDA Toolkit root directory" )
else()
    message( STATUS "Building builtin algebramex file")
    set( OSQP_ALGEBRA_BACKEND    "builtin" CACHE STRING "Build builtin algebra")
    set( OSQP_CODEGEN            ON CACHE BOOL "Enable codegen" )
    set( OSQP_ENABLE_DERIVATIVES ON CACHE BOOL "Enable derivatives" )
endif()

set( OSQP_USE_FLOAT          OFF CACHE BOOL "Use doubles" )
set( OSQP_ENABLE_PRINTING    ON CACHE BOOL "Enable printing" )
set( OSQP_ENABLE_PROFILING   ON CACHE BOOL "Enable profiling" )
set( OSQP_ENABLE_INTERRUPT   ON CACHE BOOL "Enable interrupt" )

# Configure the build outputs for the OSQP library
set( OSQP_BUILD_STATIC_LIB ON CACHE BOOL "Build static library" )
set( OSQP_BUILD_SHARED_LIB OFF CACHE BOOL "Disable shared library" )
set( OSQP_BUILD_UNITTESTS  OFF CACHE BOOL "Disable unit tests" )
set( OSQP_BUILD_DEMO_EXE   OFF CACHE BOOL "Disable demo executable" )

# Add the custom MATLAB memory handling
set( OSQP_CUSTOM_MEMORY "memory_matlab.h" )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

# Actually pull in the OSQP targets
add_subdirectory( osqp_sources )

# Add OSQP include directories that are needed
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/osqp_sources/include/public
                     ${CMAKE_CURRENT_SOURCE_DIR}/osqp_sources/include/private )

# Interrupt support requires the libut library
get_filename_component( Matlab_LIB_DIR ${Matlab_MEX_LIBRARY} DIRECTORY )

if( WIN32 )
    set( UT_LIBRARY "${Matlab_LIB_DIR}/libut.lib" )
else()
    set( UT_LIBRARY "-L${Matlab_LIB_DIR} -lut" )
endif()

if( MATLAB_CUDA )
    set( ADDITIONAL_LIBS "-L${MATLAB_CUDA_LIBS}" )
else()
    set( ADDITIONAL_LIBS "")
endif()

matlab_add_mex( NAME osqp_mex
                SRC ${CMAKE_CURRENT_SOURCE_DIR}/osqp_mex.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/interrupt_matlab.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/memory_matlab.c
                    ${CMAKE_CURRENT_SOURCE_DIR}/osqp_struct_info.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/osqp_struct_settings.cpp
                LINK_TO osqpstatic
                        ${UT_LIBRARY}
                        ${ADDITIONAL_LIBS}
                # Force compilation in the traditional C API (equivalent to the -R2017b flag)
                R2017b )