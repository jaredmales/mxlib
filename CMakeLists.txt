cmake_minimum_required(VERSION 3.17)
set(CMAKE_VERBOSE_MAKEFILE ON)

include(FindPkgConfig)
include(CheckLanguage)
find_package(PkgConfig REQUIRED)

project(mxlib LANGUAGES CXX C)
include_directories(include)

#########################################
# Build the list of source files
#
# Each subdirectory populates OBJLIB_SOURCES with its contents
#########################################
add_subdirectory(source)

#########################################
# Generate GIT status headers
#########################################

# UNCOMP is git status of the headers
add_custom_command(
  OUTPUT include/mxlib_uncomp_version.h
  COMMAND bash ./gengithead.sh ./ ./include/mxlib_uncomp_version.h MXLIB_UNCOMP
  DEPENDS ${OBJLIB_SOURCES} include/mxlib_uncomp_version.h ./gengithead.sh
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

# COMP is git status of the compiled source
add_custom_command(
  OUTPUT source/mxlib_comp_version.h
  COMMAND bash ./gengithead.sh ./ ./source/mxlib_comp_version.h MXLIB_COMP
  DEPENDS ${OBJLIB_SOURCES} source/mxlib_comp_version.h ./gengithead.sh
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

##########################################################################
# The Object Library
#
# In addition to the sources, this depends on both the COMP and UNCOMP GIT status headers
##########################################################################
add_library(objlib OBJECT ${OBJLIB_SOURCES} source/mxlib_comp_version.h include/mxlib_uncomp_version.h)

# shared libraries need PIC
set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)

###############################################################
# The actual library
###############################################################
# shared and static libraries built from the same object files
add_library(mxlib-shared SHARED $<TARGET_OBJECTS:objlib> $<TARGET_OBJECTS:levmar> $<TARGET_OBJECTS:sofa>)
add_library(mxlib-static STATIC $<TARGET_OBJECTS:objlib> $<TARGET_OBJECTS:levmar> $<TARGET_OBJECTS:sofa>)

#This makes it so we have libmxlib.so and libmxlib.a, without the -shared/-static
SET_TARGET_PROPERTIES(mxlib-shared PROPERTIES OUTPUT_NAME mxlib CLEAN_DIRECT_OUTPUT 1)
SET_TARGET_PROPERTIES(mxlib-static PROPERTIES OUTPUT_NAME mxlib CLEAN_DIRECT_OUTPUT 1)


### TODO
#  cxxversion
#
#  OPTIMIZATIONS      X
#  CXXFLAGS           X
#  CXXVERSION         X
#  CFLAGS             X
#  EIGEN              X
#  CFITSIO            X
#  GSL                X
#  levmar             X
#  SOFA               X
#
# Options to implement:
#                 linux   darwin
#  NEED_CUDA        -       -
#  NEED_ISIO        X       ?
#  USE_FFT_FROM     X
#  USE_BLAS_FROM
#       MKL         X       -
#       openblas    X       -
#       openblas-serial -   -
#       openblas-openmp -   -
#
#
#
#

set(MXLIB_CVERSION "--std=c99" CACHE STRING "c version")
set(MXLIB_CXXVERSION "--std=c++20" CACHE STRING "c++ version")

set(MXLIB_OPTIMIZE "-O3 --fast-math" CACHE STRING "optimization settings")

set(MXLIB_DEFINES "-D_XOPEN_SOURCE=700" CACHE STRING "compiler command line defines")

set(MXLIB_CFLAGS "-march=native -mfma -DNDEBUG" CACHE STRING "c compiler flags")
set(MXLIB_CXXFLAGS "-march=native -mfma -DNDEBUG" CACHE STRING "c++ compiler flags")

option(MXLIB_USE_OPENMP "Use OpenMP" ON)

option(MXLIB_USE_CUDA  "Use CUDA library" ON)
option(MXLIB_USE_ISIO  "Use ImageStreamIO library" ON)

set(MXLIB_USE_FFT_FROM "fftw" CACHE STRING "Which library to use for the FFT interface")

set(MXLIB_USE_BLAS_FROM "mkl" CACHE STRING "Which library to use for BLAS and Lapack (mkl, openblas)")

#These options are for MKL, which has a cmake config but it's hard to work with
set(MXLIB_MKL_INCLUDE_DIRECTORY "$ENV{MKLROOT}/include" CACHE STRING "MKL include directory")
set(MXLIB_MKL_DEFINITIONS "-DMXLIB_MKL" CACHE STRING "MKL defines")
set(MXLIB_MKL_CFLAGS "-m64" CACHE STRING "MKL cflags")
set(MXLIB_MKL_LINK_DIRECTORIES "$ENV{MKLROOT}/lib/intel64" CACHE STRING "MKL link directories")
set(MXLIB_MKL_LDFLAGS "-Wl,--no-as-needed" CACHE STRING "MKL linker flags")
set(MXLIB_MKL_LIBS "mkl_intel_lp64;mkl_sequential;mkl_core;pthread;m;dl" CACHE STRING "MKL libraries")

############################################
# Compiler Setup
############################################

add_compile_definitions(${MXLIB_DEFINES})

set(CMAKE_C_FLAGS "${MXLIB_CVERSION} ${MXLIB_OPTIMIZE} ${MXLIB_CFLAGS}")
set(CMAKE_CXX_FLAGS "${MXLIB_CXXVERSION} ${MXLIB_OPTIMIZE} ${MXLIB_CXXFLAGS}")

set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} ${CMAKE_CXX_FLAGS} ${MXLIB_DEFINES}")

############################################
# OpenMP setup
############################################

if(MXLIB_USE_OPENMP)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} -fopenmp")

    message("mxlib: OpenMP support enabled")

    #This doesn't work on 24.04 as of 07/2025:
    #find_package(OpenMP)

    #if(OpenMP_C_FOUND)

    #    target_compile_features(objlib PUBLIC OpenMP)
    #    target_compile_features(mxlib-shared PUBLIC OpenMP)
    #    target_compile_features(mxlib-static PUBLIC OpenMP)

    #    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} -I${OpenMP_C_INCLUDE_DIRS} ${OpenMP_C_FLAGS}")

    #    message("mxlib: OpenMP support enabled")
    #else()
    #    message("mxlib: OpenMP support not found")
    #endif()

else()
    message("mxlib: OpenMP support not enabled")
endif()

############################################
# boost setup
############################################

find_package(Boost REQUIRED COMPONENTS system filesystem)

include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(objlib PUBLIC ${Boost_LIBRARIES})
target_link_libraries (mxlib-shared PUBLIC ${Boost_LIBRARIES})
target_link_libraries (mxlib-static PUBLIC ${Boost_LIBRARIES})

# don't add anything to MXLIB_PC_CFLAGS b/c we don't link boost for apps

############################################
# GSL setup
############################################

pkg_check_modules(GSL REQUIRED IMPORTED_TARGET gsl)
target_link_libraries(objlib PUBLIC PkgConfig::GSL)
target_link_libraries (mxlib-shared PUBLIC PkgConfig::GSL)
target_link_libraries (mxlib-static PUBLIC PkgConfig::GSL)

set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} ${GSL_CFLAGS}")

############################################
# Eigen setup
############################################

pkg_check_modules(EIGEN3 REQUIRED IMPORTED_TARGET eigen3)
target_link_libraries(objlib PUBLIC PkgConfig::EIGEN3)
target_link_libraries (mxlib-shared PUBLIC PkgConfig::EIGEN3)
target_link_libraries (mxlib-static PUBLIC PkgConfig::EIGEN3)

#Don't let Eigen use openmp
add_compile_definitions(-DEIGEN_DONT_PARALLELIZE)

set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} ${EIGEN3_CFLAGS} -DEIGEN_DONT_PARALLELIZE")
set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${EIGEN3_LDFLAGS}")

############################################
# BLAS & Lapack setup
############################################

if(${MXLIB_USE_BLAS_FROM} STREQUAL "mkl")
    #MKL has a cmake config but it is hard to customize, so we do it manually
    include_directories(${MXLIB_MKL_INCLUDE_DIRECTORY})
    add_compile_definitions(${MXLIB_MKL_DEFINITIONS})
    add_compile_options(${MXLIB_MKL_CFLAGS})

    add_link_options(${MXLIB_MKL_LDFLAGS})
    target_link_directories(objlib PUBLIC ${MXLIB_MKL_LINK_DIRECTORIES})
    target_link_directories(mxlib-shared PUBLIC ${MXLIB_MKL_LINK_DIRECTORIES})
    target_link_directories(mxlib-static PUBLIC ${MXLIB_MKL_LINK_DIRECTORIES})
    target_link_libraries(objlib PUBLIC ${MXLIB_MKL_LIBS})
    target_link_libraries(mxlib-shared PUBLIC ${MXLIB_MKL_LIBS})
    target_link_libraries(mxlib-static PUBLIC ${MXLIB_MKL_LIBS})


    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} -I${MXLIB_MKL_INCLUDE_DIRECTORY} ${MXLIB_MKL_DEFINITIONS} ${MXLIB_MKL_CFLAGS}" )
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${MXLIB_MKL_LDFLAGS}")
elseif(${MXLIB_USE_BLAS_FROM} STREQUAL "openblas")

    pkg_check_modules(OPENBLAS REQUIRED IMPORTED_TARGET openblas)
    target_link_libraries(objlib PUBLIC PkgConfig::OPENBLAS)
    target_link_libraries (mxlib-shared PUBLIC PkgConfig::OPENBLAS)
    target_link_libraries (mxlib-static PUBLIC PkgConfig::OPENBLAS)

    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} ${OPENBLAS_CFLAGS}")

    pkg_check_modules(LAPACK REQUIRED IMPORTED_TARGET lapack)
    target_link_libraries(objlib PUBLIC PkgConfig::LAPACK)
    target_link_libraries (mxlib-shared PUBLIC PkgConfig::LAPACK)
    target_link_libraries (mxlib-static PUBLIC PkgConfig::LAPACK)

    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} ${LAPACK_CFLAGS}")

else(NOT ${MXLIB_USE_BLAS_FROM} STREQUAL "")
    message(FATAL_ERROR "MXLIB_USE_BLAS_FROM = ${MXLIB_USE_BLAS_FROM} is not supported.")
endif()

if(NOT ${MXLIB_USE_BLAS_FROM} STREQUAL "")
  message("mxlib: BLAS & Lapack support enabled with ${MXLIB_USE_BLAS_FROM}")
else()
  message("mxlib: BLAS & Lapack support not enabled")
endif()

############################################
# FFTW setup
############################################
if(${MXLIB_USE_FFT_FROM} STREQUAL "fftw")

    pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)
    include_directories(${FFTW3_INCLUDE_DIRS})
    target_link_libraries(objlib PUBLIC PkgConfig::FFTW3)
    target_link_libraries(mxlib-shared PUBLIC PkgConfig::FFTW3)
    target_link_libraries(mxlib-static PUBLIC PkgConfig::FFTW3)
    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} -I${FFTW3_INCLUDE_DIRS}")

    pkg_check_modules(FFTW3F REQUIRED IMPORTED_TARGET fftw3f)
    include_directories(${FFTW3F_INCLUDE_DIRS})
    target_link_libraries(objlib PUBLIC PkgConfig::FFTW3F)
    target_link_libraries(mxlib-shared PUBLIC PkgConfig::FFTW3F)
    target_link_libraries(mxlib-static PUBLIC PkgConfig::FFTW3F)
    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} -I${FFTW3F_INCLUDE_DIRS}")

    pkg_check_modules(FFTW3L REQUIRED IMPORTED_TARGET fftw3l)
    include_directories(${FFTW3L_INCLUDE_DIRS})
    target_link_libraries(objlib PUBLIC PkgConfig::FFTW3L)
    target_link_libraries(mxlib-shared PUBLIC PkgConfig::FFTW3L)
    target_link_libraries(mxlib-static PUBLIC PkgConfig::FFTW3L)
    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} -I${FFTW3F_INCLUDE_DIRS}")

    pkg_check_modules(FFTW3Q REQUIRED IMPORTED_TARGET fftw3q)
    include_directories(${FFTW3Q_INCLUDE_DIRS})
    target_link_libraries(objlib PUBLIC PkgConfig::FFTW3Q)
    target_link_libraries(mxlib-shared PUBLIC PkgConfig::FFTW3Q)
    target_link_libraries(mxlib-static PUBLIC PkgConfig::FFTW3Q)

    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} -I${FFTW3Q_INCLUDE_DIRS}")
else()
    message(FATAL_ERROR "mxlib: MXLIB_USE_FFT_FROM = ${MXLIB_USE_FFT_FROM} is not supported.")
endif()

if(NOT ${MXLIB_USE_FFT_FROM} STREQUAL "")
  message("mxlib: FFT support enabled with ${MXLIB_USE_FFT_FROM}")
else()
  message("mxlib: FFT support not enabled (must have fftw3, fftw3f, fftw3l, and fftw3q)")
endif()

############################################
# cfitsio setup
############################################

pkg_check_modules(CFITSIO REQUIRED IMPORTED_TARGET cfitsio)

include_directories(${CFITSIO_INCLUDE_DIRS})
target_link_libraries(objlib PUBLIC PkgConfig::CFITSIO)
target_link_libraries(mxlib-shared PUBLIC PkgConfig::CFITSIO)
target_link_libraries(mxlib-static PUBLIC PkgConfig::CFITSIO)

set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} ${CFITSIO_CFLAGS}")



############################################
# ImageStreamIO (ISIO) setup
############################################
if(MXLIB_USE_ISIO)

    pkg_check_modules(ISIO REQUIRED IMPORTED_TARGET ImageStreamIO)
    target_link_libraries(objlib PUBLIC PkgConfig::ISIO)
    target_link_libraries(mxlib-shared PUBLIC PkgConfig::ISIO)
    target_link_libraries(mxlib-static PUBLIC PkgConfig::ISIO)

    include_directories(${ISIO_INCLUDE_DIRS})
    add_compile_definitions("-DMXLIB_MILK")

    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} ${ISIO_CFLAGS} -DMXLIB_MILK")

    message("mxlib: ImageStreamIO support enabled")
else()
    message("mxlib: ImageStreamIO support not enabled")
endif()

############################################
# CUDA setup
############################################
if(MXLIB_USE_CUDA)
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
  else()
    message(WARNING "No CUDA support found")
    set(MXLIB_USE_CUDA OFF)
  endif()
endif(MXLIB_USE_CUDA)

if(MXLIB_USE_CUDA)
  message("mxlib: CUDA support enabled but not yet supported")
else()
  message("mxlib: CUDA support not enabled")
endif(MXLIB_USE_CUDA)


###############################################################
# Generate PC file CFLAGS and LDFLAGS
###############################################################

#list-ify to remove duplicates then de-list-ify
string(REPLACE " " ";" MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS}")
list(REMOVE_DUPLICATES MXLIB_PC_CFLAGS)
string(REPLACE ";" " " MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS}")

string(REPLACE " " ";" MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS}")
list(REMOVE_DUPLICATES MXLIB_PC_LDFLAGS)
string(REPLACE ";" " " MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS}")

message("cflags= ${MXLIB_PC_CFLAGS}")
message("ldflags= ${MXLIB_PC_LDFLAGS}")



