###########################################################################
#                 CMake Build File for ICM-BIDEM
#
#    Written By: Zach Cobell
#
###########################################################################
#
# The CMake build system enable ICM-HYDRO to be deployed and built
# in a cross platform environment. 
#
###########################################################################
#...Set the default build type
IF(DEFINED CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
        build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
        Release RelWithDebInfo MinSizeRel.")
ELSE()
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
        options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
        RelWithDebInfo MinSizeRel.")
ENDIF()
###########################################################################

###########################################################################
#  GENERAL OPTIONS
###########################################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)
PROJECT(ICM_BIDEM Fortran)
ENABLE_LANGUAGE( C )

###########################################################################


###########################################################################
#...Put the static libraries in the CMakeFiles folder so they don't
#   contaminate the build directory
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles)
###########################################################################


###########################################################################
#...Set the version string (used for libraries)
SET(ICMBIDEM_VERSION_MAJOR 23)
SET(ICMBIDEM_VERSION_MINOR  0)
SET(ICMBIDEM_VERSION_PATCH  0)
SET(ICMBIDEM_VERSION_STRING ${ICMBIDEM_VERSION_MAJOR}.${ICMBIDEM_VERSION_MINOR}.${ICMBIDEM_VERSION_PATCH})
###########################################################################


###########################################################################
#...Determine architecture specific parameters
INCLUDE(${CMAKE_SOURCE_DIR}/cmake/architecture.cmake)
###########################################################################


###########################################################################
# ICM_BIDEM
###########################################################################
SET( BIDEM hydro_v${ICMBIDEM_VERSION_STRING} )
ADD_EXECUTABLE( ${BIDEM} CMP_BI.f90
                         mod_global.f90
                         mod_util.f90 )
SET_TARGET_PROPERTIES( ${BIDEM} PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/CMakeFiles/mod)
SET_TARGET_PROPERTIES( ${BIDEM} PROPERTIES COMPILE_FLAGS "${Fortran_LINELENGTH_FLAG} ${Fortran_COMPILER_SPECIFIC_FLAG}" )
