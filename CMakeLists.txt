#################################################################
#                      cmake configuration for mxlib
#
#################################################################
## TODO
#    - make FFT optional
#    - make BLAS & Lapack optional
#    - make CFITSIO optional
#    - deal with HASQUAD (should be MXLIB_HASQUAD)
#    - enable build without GSL
#    - enable build without CFITSIO
#    - test build without fftw
#    - test build without BLAS
#    - test build without ISIO
#    - test build without CUDA
# (i.o.w. minimum depends should be boost and Eigen.  maybe we don't even need Eigen)
#
#
## Compliance with Makefile build;
#                    linux   darwin
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
#                     linux   darwin
#  NEED_CUDA            X
#  NEED_ISIO            X       ?
#  USE_FFT_FROM         X
#  USE_BLAS_FROM:      ------------
#       MKL             X       -
#       openblas        X       -
#       openblas-serial -       -
#       openblas-openmp -       -
#
# New Todos:
#

#########################################
# basic setup
#
#########################################

cmake_minimum_required(VERSION 3.24)

#set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

include(FindPkgConfig)
include(CheckLanguage)
find_package(PkgConfig REQUIRED)

project(mxlib LANGUAGES CXX C
              VERSION "0.2.0"
              DESCRIPTION "c++ library for astronomy signal processing"
              HOMEPAGE_URL "https://jaredmales.github.io/mxlib-doc/")

## macro linklibdirs_from_list
 # Turn list of library directories into a -L delimited string for linking
 #
 # First: LINK_LIBDIRS: output string
 # Second: LIB_DIRS_LIST: input list to process
 #
macro(linklibdirs_from_list LINK_LIBDIRS LIBDIR_LIST)

    string(REPLACE ";" " -L" ${LINK_LIBDIRS} "${${LIBDIR_LIST}}")
    set(${LINK_LIBDIRS} "-L${${LINK_LIBDIRS}}")

endmacro()

## Extract the library name from a path
 # e.g. /path/to/libname.so produces name
 #
 # LIBNAME: output contains the library name
 # LIBPATH: the path to the library
 #
macro(libname_from_path LIBNAME LIBPATH)

    get_filename_component(${LIBNAME} ${${LIBPATH}} NAME)
    get_filename_component(LIBEXT ${${LIBPATH}} LAST_EXT)

    string(REPLACE "lib" "" ${LIBNAME} ${${LIBNAME}})
    string(REPLACE ${LIBEXT} "" ${LIBNAME} ${${LIBNAME}})

endmacro()

## macro liblist_strip_namespace
 # Strip namespace from library names in a list and add a prefix
 # assumes that the same prefix is added to all libs in the list
 # E.g. this can be used to strip Boost:: and replace it with boost_
 #
 # LIBRARY_LIST: input/output list of libraries
 # prefix: input the prefix to add
 #
macro(liblist_strip_namespace LIBRARY_LIST prefix)

    #strip namespace from beginning
    string(REGEX REPLACE "^[^;]*::" "${prefix}" ${LIBRARY_LIST} "${${LIBRARY_LIST}}")

    #strip remaining namespaces
    string(REGEX REPLACE "[ ;].*::" ";${prefix}" ${LIBRARY_LIST} "${${LIBRARY_LIST}}")

endmacro()

## Function linklibs_from_list
 # Turn list of library names into a -l delimited string for linking
 # If the library is specified with a full path we don't mess with it
 #
 # First: LINK_LIBRARIES: output string
 # Second: LIBRARY_LIST: input list to process
 #
macro(linklibs_from_list LINK_LIBRARIES LIBRARY_LIST)

    foreach(lib ${${LIBRARY_LIST}})

         cmake_path(HAS_PARENT_PATH lib LIB_IS_PATH)

         #message("lib: ${lib} ${LIB_IS_PATH}")

         if(NOT LIB_IS_PATH)
            string(PREPEND lib "-l")
         endif()

         set(LIBRARY_LIST_OUT "${LIBRARY_LIST_OUT} ${lib}")

    endforeach()

    #message("LIBRARY_LIST_OUT: ${LIBRARY_LIST_OUT}")

    set(${LINK_LIBRARIES} "${LIBRARY_LIST_OUT}")

    #message("LINK_LIBRARIES: ${${LINK_LIBRARIES}}")

endmacro()

## Function dump_cmake_variables
 # Dump all the current variables for debugging
 #
function(dump_cmake_variables)
    get_cmake_property(_variableNames VARIABLES)
    list (SORT _variableNames)
    foreach (_variableName ${_variableNames})
        if ((NOT DEFINED ARGV0) OR _variableName MATCHES ${ARGV0})
            message(STATUS "${_variableName}=${${_variableName}}")
        endif()
    endforeach()
endfunction()

#########################################
# mxlib settings
#
# each setting is prefixed with MXLIB_
#########################################

set(MXLIB_CVERSION "--std=c99" CACHE STRING "mxlib c version")
set(MXLIB_CXXVERSION "--std=c++20" CACHE STRING "mxlib c++ version (minimum is c++20)")

#Jared's setup: "-march=native -mfma -DNDEBUG"
set(MXLIB_CFLAGS "-march=native" CACHE STRING "mxlib c compiler flags")
set(MXLIB_CXXFLAGS "-march=native" CACHE STRING "mxlib c++ compiler flags")

set(MXLIB_OPTIMIZE "-O3 --fast-math" CACHE STRING "mxlib optimization settings")

set(MXLIB_DEFINES "-D_XOPEN_SOURCE=700" CACHE STRING "mxlib compiler command line defines")

option(MXLIB_USE_OPENMP "Whether or not to use OpenMP in mxlib and other libraries" ON )

option(MXLIB_USE_CUDA  "Whether or not to use CUDA library for mxlib" ON)
option(MXLIB_USE_ISIO  "Whether or not to use the ImageStreamIO library for mxlib" ON)
option(MXLIB_BUILD_TESTS "Whether or not to build mxlib tests" OFF)
option(MXLIB_BUILD_COVERAGE "Whether or not to enable coverage instrumentation for C/C++ builds" OFF)
set(MXLIB_COVERAGE_TEST_TIMEOUT 300 CACHE STRING "CTest timeout (seconds) for coverage target")
option(MXLIB_BUILD_DOCS "Whether or not to build Doxygen docs as part of the default build" OFF)

set(MXLIB_USE_FFT_FROM "fftw" CACHE STRING "Which library to use for the FFT interface in mxlib")

set(MXLIB_USE_BLAS_FROM "openblas" CACHE STRING "Which library to use for BLAS and Lapack in mxlib.  Supported choices: mkl, openblas.")

#These options are for MKL, which has a cmake config but it's hard to work with
set(MXLIB_MKL_INCLUDE_DIRECTORY "$ENV{MKLROOT}/include" CACHE STRING "MKL include directory")
set(MXLIB_MKL_DEFINITIONS "-DMXLIB_MKL" CACHE STRING "MKL compiler defines")
set(MXLIB_MKL_CFLAGS "-m64" CACHE STRING "MKL compiiler cflags (used for cxx too)")
set(MXLIB_MKL_LIBDIRS "$ENV{MKLROOT}/lib/intel64" CACHE STRING "list of MKL link directories")
set(MXLIB_MKL_LDFLAGS "-Wl,--no-as-needed" CACHE STRING "MKL linker flags")
set(MXLIB_MKL_LIBS "mkl_intel_lp64;mkl_sequential;mkl_core;pthread;m;dl" CACHE STRING "list of MKL libraries")

set(MXLIB_CUDA_ARCH "native" CACHE STRING "the cuda architecture")

#########################################
# Validate config settings
#
# Only things needed for conditional includes
#########################################

# CUDA
if(MXLIB_USE_CUDA)

    check_language(CUDA)

    find_package(CUDAToolkit)

    if(CUDAToolkit_FOUND)
	    message("found cudatoolkit")
    endif()

    if(CMAKE_CUDA_COMPILER)
	    message("found cuda compiler")
    endif()

    if(CMAKE_CUDA_COMPILER AND CUDAToolkit_FOUND)

    else()
        message(WARNING "No CUDA support found, setting MXLIB_USE_CUDA to off")
        set(MXLIB_USE_CUDA OFF)
    endif()
endif()


#########################################
# header files
#
#
#########################################

#first generate the errno header
add_executable(generate_errno_info ${CMAKE_SOURCE_DIR}/include/error/generators/generate_errno_info.cpp)

add_custom_command(
    OUTPUT  ${CMAKE_SOURCE_DIR}/include/error/generators/errno_info.hpp
    COMMAND generate_errno_info
    DEPENDS generate_errno_info include/error/generators/c++11_errnos.txt
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/include/error/generators/
    #DEPENDS_EXPLICIT_ONLY
)

#then generate the FITS header
add_executable(generate_fits_status_info ${CMAKE_SOURCE_DIR}/include/error/generators/generate_fits_status_info.cpp)

add_custom_command(
    OUTPUT  ${CMAKE_SOURCE_DIR}/include/error/generators/fits_status_info.hpp
    COMMAND generate_fits_status_info
    DEPENDS generate_fits_status_info include/error/generators/fits_error_codes.txt
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/include/error/generators/
    #DEPENDS_EXPLICIT_ONLY
)

#then generate the CUDA header
if(MXLIB_USE_CUDA)
    add_executable(generate_cudaError_info ${CMAKE_SOURCE_DIR}/include/error/generators/generate_cudaError_info.cpp)
    target_link_libraries(generate_cudaError_info PUBLIC CUDA::cudart)

    add_custom_command(
        OUTPUT  ${CMAKE_SOURCE_DIR}/include/error/generators/cudaError_info.hpp
        COMMAND generate_cudaError_info
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/include/error/generators/
        #DEPENDS_EXPLICIT_ONLY
        )
endif()


set(ERROR_T_SRCS ${CMAKE_SOURCE_DIR}/include/error/generators/generate_error_t.cpp
                                ${CMAKE_SOURCE_DIR}/include/error/generators/errno_info.hpp
                                ${CMAKE_SOURCE_DIR}/include/error/generators/fits_status_info.hpp
                                )
if(MXLIB_USE_CUDA)
set(ERROR_T_SRCS ${ERROR_T_SRCS} ${CMAKE_SOURCE_DIR}/include/error/generators/cudaError_info.hpp)
endif()

#now generate errno_t.hpp
add_executable(generate_error_t ${ERROR_T_SRCS})

add_custom_command(
    OUTPUT  ${CMAKE_SOURCE_DIR}/include/error/error_t.hpp
    COMMAND generate_error_t
    DEPENDS generate_error_t include/error/generators/errno_info.hpp include/error/generators/error_t_list.txt
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/include/error/generators/
    #DEPENDS_EXPLICIT_ONLY
)

include_directories(include)

add_subdirectory(include)


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
  COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/gengithead.sh ${CMAKE_CURRENT_SOURCE_DIR} include/mxlib_uncomp_version.h MXLIB_UNCOMP
  DEPENDS ${OBJLIB_INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR}/gengithead.sh
  #WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

# COMP is git status of the compiled source
add_custom_command(
  OUTPUT source/mxlib_comp_version.h
  COMMAND bash ./gengithead.sh ./ ./source/mxlib_comp_version.h MXLIB_COMP
  DEPENDS ${OBJLIB_SOURCES} ./gengithead.sh
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

##########################################################################
# The Object Library
#
# In addition to the sources, this depends on both the COMP and UNCOMP GIT status headers
##########################################################################
add_library(objlib OBJECT ${OBJLIB_SOURCES} ${OBJLIB_INCLUDES} source/mxlib_comp_version.h include/mxlib_uncomp_version.h include/error/error_t.hpp)

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

install (DIRECTORY include/ DESTINATION include/mx FILES_MATCHING PATTERN "CMakeLists.txt" EXCLUDE PATTERN "*.hpp" PATTERN "*.h" PATTERN "*.inc" PATTERN "generators" EXCLUDE)
install (FILES include/mxlib.hpp include/mxError.hpp include/mxException.hpp DESTINATION include/mx)
install (FILES include/error/error_t.hpp DESTINATION include/mx/error)

#We install the generated version behind cmake's back so it does get treated as a changed depend
install (FILES ${CMAKE_CURRENT_BINARY_DIR}/include/mxlib_uncomp_version.h DESTINATION include/mx)

install (TARGETS mxlib-shared EXPORT mxlib-sharedTargets DESTINATION lib)
install (TARGETS mxlib-static EXPORT mxlib-staticTargets DESTINATION lib)

install(
    EXPORT mxlib-sharedTargets
    FILE mxlib-sharedTargets.cmake
    DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/mxlib"
)
############################################
# Compiler Setup
############################################

add_compile_definitions(${MXLIB_DEFINES})

set(CMAKE_C_FLAGS "${MXLIB_CVERSION} ${MXLIB_OPTIMIZE} ${MXLIB_CFLAGS}")
set(CMAKE_CXX_FLAGS "${MXLIB_CXXVERSION} ${MXLIB_OPTIMIZE} ${MXLIB_CXXFLAGS}")

set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} ${CMAKE_CXX_FLAGS} ${MXLIB_DEFINES}")

if(MXLIB_BUILD_COVERAGE)
    # Match the legacy Make coverage behavior: --coverage and -O0 for instrumented builds.
    add_compile_definitions(MXLIB_BUILD_COVERAGE=1)
    add_compile_options(
        $<$<COMPILE_LANGUAGE:C>:--coverage>
        $<$<COMPILE_LANGUAGE:C>:-O0>
        $<$<COMPILE_LANGUAGE:CXX>:--coverage>
        $<$<COMPILE_LANGUAGE:CXX>:-O0>
    )
    add_link_options(--coverage)
endif()

############################################
# OpenMP setup
############################################

if(MXLIB_USE_OPENMP)

    find_package(OpenMP)

    if(NOT (OpenMP_C_FOUND AND OpenMP_CXX_FOUND))

        message("mxlib: OpenMP support not found")

        set(MXLIB_USE_OPENMP OFF)
    else()

        target_link_libraries(objlib PUBLIC OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
        target_link_libraries(mxlib-shared PUBLIC OpenMP::OpenMP_C OpenMP::OpenMP_CXX)
        target_link_libraries(mxlib-static PUBLIC OpenMP::OpenMP_C OpenMP::OpenMP_CXX)

        set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} ${OpenMP_CXX_FLAGS}")
        set(MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS} ${OpenMP_CXX_FLAGS}")

        #this is so dumb, the following does not produce a linkable system:
        # linklibdirs_from_list(OpenMP_LIBDIRS_FLAGS OpenMP_CXX_LIBRARY)
        # set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${OpenMP_LIBDIRS_FLAGS}")

        # linklibs_from_list(OpenMP_LIBRARIES_FLAGS OpenMP_CXX_LIB_NAMES)
        # set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${OpenMP_LIBRARIES_FLAGS}")
        #but this does
        set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} -fopenmp")
        #which does not appear to be supported by any variable produced in FindOpenMP

        message("mxlib: OpenMP support enabled")

    endif()

else()
    message("mxlib: OpenMP support not enabled")
endif()

############################################
# boost setup
############################################

find_package(Boost REQUIRED COMPONENTS system) # filesystem)

include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(objlib PUBLIC ${Boost_LIBRARIES})
target_link_libraries (mxlib-shared PUBLIC ${Boost_LIBRARIES})
target_link_libraries (mxlib-static PUBLIC ${Boost_LIBRARIES})

# don't add anything to MXLIB_PC_CFLAGS b/c we don't link boost for apps

linklibdirs_from_list(Boost_LIBDIRS_FLAGS Boost_LIBRARY_DIRS)
set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${Boost_LIBDIRS_FLAGS}")

liblist_strip_namespace(Boost_LIBRARIES "boost_")
linklibs_from_list(Boost_LIBRARIES_FLAGS Boost_LIBRARIES)
set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${Boost_LIBRARIES_FLAGS}")


############################################
# Eigen setup
############################################

pkg_check_modules(EIGEN3 REQUIRED IMPORTED_TARGET eigen3)
target_link_libraries(objlib PUBLIC PkgConfig::EIGEN3)
target_link_libraries (mxlib-shared PUBLIC PkgConfig::EIGEN3)
target_link_libraries (mxlib-static PUBLIC PkgConfig::EIGEN3)

#Don't let Eigen use openmp
add_compile_definitions(-DEIGEN_DONT_PARALLELIZE)

set(MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS} -I${EIGEN3_INCLUDE_DIRS} ${EIGEN3_CFLAGS} ${EIGEN3_CFLAGS_OTHER}")

set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} -I${EIGEN3_INCLUDE_DIRS} ${EIGEN3_CFLAGS} ${EIGEN3_CFLAGS_OTHER} -DEIGEN_DONT_PARALLELIZE")
set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${EIGEN3_LDFLAGS}")

############################################
# FFTW setup
############################################
if(${MXLIB_USE_FFT_FROM} STREQUAL "fftw")

    #double
    pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)
    include_directories(${FFTW3_INCLUDE_DIRS})
    target_link_libraries(objlib PUBLIC PkgConfig::FFTW3)
    target_link_libraries(mxlib-shared PUBLIC PkgConfig::FFTW3)
    target_link_libraries(mxlib-static PUBLIC PkgConfig::FFTW3)

    set(MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS} ${FFTW3_CFLAGS} ${FFTW3_CFLAGS_OTHER}")
    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} -I${FFTW3_INCLUDE_DIRS}")

    linklibdirs_from_list(FFTW3_LIBDIRS_FLAGS FFTW3_LIBRARY_DIRS)
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${FFTW3_LIBDIRS_FLAGS}")

    linklibs_from_list(FFTW3_LIBRARIES_FLAGS FFTW3_LIBRARIES)
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${FFTW3_LIBRARIES_FLAGS}")

    #float
    pkg_check_modules(FFTW3F REQUIRED IMPORTED_TARGET fftw3f)
    include_directories(${FFTW3F_INCLUDE_DIRS})
    target_link_libraries(objlib PUBLIC PkgConfig::FFTW3F)
    target_link_libraries(mxlib-shared PUBLIC PkgConfig::FFTW3F)
    target_link_libraries(mxlib-static PUBLIC PkgConfig::FFTW3F)

    set(MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS} ${FFTW3F_CFLAGS} ${FFTW3F_CFLAGS_OTHER}")
    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} -I${FFTW3F_INCLUDE_DIRS}")

    linklibdirs_from_list(FFTW3F_LIBDIRS_FLAGS FFTW3F_LIBRARY_DIRS)
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${FFTW3F_LIBDIRS_FLAGS}")

    linklibs_from_list(FFTW3F_LIBRARIES_FLAGS FFTW3F_LIBRARIES)
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${FFTW3F_LIBRARIES_FLAGS}")

    #long double
    pkg_check_modules(FFTW3L REQUIRED IMPORTED_TARGET fftw3l)
    include_directories(${FFTW3L_INCLUDE_DIRS})
    target_link_libraries(objlib PUBLIC PkgConfig::FFTW3L)
    target_link_libraries(mxlib-shared PUBLIC PkgConfig::FFTW3L)
    target_link_libraries(mxlib-static PUBLIC PkgConfig::FFTW3L)

    set(MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS} ${FFTW3L_CFLAGS} ${FFTW3L_CFLAGS_OTHER}")
    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} -I${FFTW3L_INCLUDE_DIRS}")

    linklibdirs_from_list(FFTW3L_LIBDIRS_FLAGS FFTW3L_LIBRARY_DIRS)
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${FFTWL_LIBDIRS_FLAGS}")

    linklibs_from_list(FFTW3L_LIBRARIES_FLAGS FFTW3L_LIBRARIES)
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${FFTW3L_LIBRARIES_FLAGS}")

    #quad/f128
    pkg_check_modules(FFTW3Q IMPORTED_TARGET fftw3q)

    if(FFTW3Q_FOUND)
        add_compile_definitions("-DHAS_QUAD")
        include_directories(${FFTW3Q_INCLUDE_DIRS})
        target_link_libraries(objlib PUBLIC PkgConfig::FFTW3Q)
        target_link_libraries(mxlib-shared PUBLIC PkgConfig::FFTW3Q)
        target_link_libraries(mxlib-static PUBLIC PkgConfig::FFTW3Q)

        set(MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS} ${FFTW3Q_CFLAGS} ${FFTW3Q_CFLAGS_OTHER}")
        set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} -DHAS_QUAD -I${FFTW3Q_INCLUDE_DIRS}")

        linklibdirs_from_list(FFTW3Q_LIBDIRS_FLAGS FFTW3Q_LIBRARY_DIRS)
        set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${FFTW3Q_LIBDIRS_FLAGS}")

        linklibs_from_list(FFTW3Q_LIBRARIES_FLAGS FFTW3Q_LIBRARIES)
        set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${FFTW3Q_LIBRARIES_FLAGS}")
    else()
        message("mxlib: quad support in FFTW not found")
    endif()


else()
    message(FATAL_ERROR "mxlib: MXLIB_USE_FFT_FROM = ${MXLIB_USE_FFT_FROM} is not supported.")
endif()

if(NOT ${MXLIB_USE_FFT_FROM} STREQUAL "")
  message("mxlib: FFT support enabled with ${MXLIB_USE_FFT_FROM}")
else()
  message("mxlib: FFT support not enabled (must have fftw3, fftw3f, fftw3l, and fftw3q)")
endif()


############################################
# BLAS & Lapack setup
############################################

if(${MXLIB_USE_BLAS_FROM} STREQUAL "mkl")
    #MKL has a cmake config but it is hard to customize, so we do it manually
    include_directories(${MXLIB_MKL_INCLUDE_DIRECTORY})
    add_compile_definitions(${MXLIB_MKL_DEFINITIONS})
    add_compile_options(${MXLIB_MKL_CFLAGS})

    add_link_options(${MXLIB_MKL_LDFLAGS})
    target_link_directories(objlib PUBLIC ${MXLIB_MKL_LIBDIRS})
    target_link_directories(mxlib-shared PUBLIC ${MXLIB_MKL_LIBDIRS})
    target_link_directories(mxlib-static PUBLIC ${MXLIB_MKL_LIBDIRS})
    target_link_libraries(objlib PUBLIC ${MXLIB_MKL_LIBS})
    target_link_libraries(mxlib-shared PUBLIC ${MXLIB_MKL_LIBS})
    target_link_libraries(mxlib-static PUBLIC ${MXLIB_MKL_LIBS})


    set(MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS} ${MXLIB_MKL_CFLAGS}")

    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} -I${MXLIB_MKL_INCLUDE_DIRECTORY} ${MXLIB_MKL_DEFINITIONS} ${MXLIB_MKL_CFLAGS}" )

    linklibdirs_from_list(MXLIB_MKL_LIBDIRS_FLAGS MXLIB_MKL_LIBDIRS)
    linklibs_from_list(MXLIB_MKL_LIBS_FLAGS MXLIB_MKL_LIBS)
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${MXLIB_MKL_LDFLAGS} ${MXLIB_MKL_LIBDIRS_FLAGS} ${MXLIB_MKL_LIBS_FLAGS}")

elseif(${MXLIB_USE_BLAS_FROM} STREQUAL "openblas")

    pkg_check_modules(OPENBLAS REQUIRED IMPORTED_TARGET openblas)
    target_link_libraries(objlib PUBLIC PkgConfig::OPENBLAS)
    target_link_libraries (mxlib-shared PUBLIC PkgConfig::OPENBLAS)
    target_link_libraries (mxlib-static PUBLIC PkgConfig::OPENBLAS)

    set(MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS} ${OPENBLAS_CFLAGS} ${OPENBLAS_CFLAGS_OTHER}")
    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} ${OPENBLAS_CFLAGS} ${OPENBLAS_CFLAGS_OTHER}")

    linklibdirs_from_list(OPENBLAS_LIBDIRS_FLAGS OPENBLAS_LIBRARY_DIRS)
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${OPENBLAS_LIBDIRS_FLAGS}")

    linklibs_from_list(OPENBLAS_LIBRARIES_FLAGS OPENBLAS_LIBRARIES)
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${OPENBLAS_LIBRARIES_FLAGS}")


    pkg_check_modules(LAPACK REQUIRED IMPORTED_TARGET lapack)
    target_link_libraries(objlib PUBLIC PkgConfig::LAPACK)
    target_link_libraries (mxlib-shared PUBLIC PkgConfig::LAPACK)
    target_link_libraries (mxlib-static PUBLIC PkgConfig::LAPACK)

    set(MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS} ${LAPACK_CFLAGS} ${LAPACK_CFLAGS_OTHER}")
    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} ${LAPACK_CFLAGS} ${LAPACK_CFLAGS_OTHER}")

    linklibdirs_from_list(LAPACK_LIBDIRS_FLAGS LAPACK_LIBRARY_DIRS)
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${LAPACK_LIBDIRS_FLAGS}")

    linklibs_from_list(LAPACK_LIBRARIES_FLAGS LAPACK_LIBRARIES)
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${LAPACK_LIBRARIES_FLAGS}")

else(NOT ${MXLIB_USE_BLAS_FROM} STREQUAL "")

    message(FATAL_ERROR "MXLIB_USE_BLAS_FROM = ${MXLIB_USE_BLAS_FROM} is not supported.")

endif()

if(NOT ${MXLIB_USE_BLAS_FROM} STREQUAL "")
  message("mxlib: BLAS & Lapack support enabled with ${MXLIB_USE_BLAS_FROM}")
else()
  message("mxlib: BLAS & Lapack support not enabled")
endif()

############################################
# GSL setup
#
# After Blas to ensure gslcblas is used as last restort
############################################

pkg_check_modules(GSL REQUIRED IMPORTED_TARGET gsl)
target_link_libraries(objlib PUBLIC PkgConfig::GSL)
target_link_libraries (mxlib-shared PUBLIC PkgConfig::GSL)
target_link_libraries (mxlib-static PUBLIC PkgConfig::GSL)

set(MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS} ${GSL_CFLAGS} ${GSL_CFLAGS_OTHER}")

set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} ${GSL_CFLAGS}")

linklibdirs_from_list(GSL_LIBDIRS_FLAGS GSL_LIBRARY_DIRS)
set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${GSL_LIBDIRS_FLAGS}")

linklibs_from_list(GSL_LIBRARIES_FLAGS GSL_LIBRARIES)
set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${GSL_LIBRARIES_FLAGS}")


############################################
# cfitsio setup
############################################

pkg_check_modules(CFITSIO REQUIRED IMPORTED_TARGET cfitsio)

include_directories(${CFITSIO_INCLUDE_DIRS})
target_link_libraries(objlib PUBLIC PkgConfig::CFITSIO)
target_link_libraries(mxlib-shared PUBLIC PkgConfig::CFITSIO)
target_link_libraries(mxlib-static PUBLIC PkgConfig::CFITSIO)

set(MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS} ${CFITSIO_CFLAGS} ${CFITSIO_CFLAGS_OTHER}")
set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} ${CFITSIO_CFLAGS}")

linklibdirs_from_list(CFITSIO_LIBDIRS_FLAGS CFITSIO_LIBRARY_DIRS)
set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${CFITSIO_LIBDIRS_FLAGS}")

linklibs_from_list(CFITSIO_LIBRARIES_FLAGS CFITSIO_LIBRARIES)
set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${CFITSIO_LIBRARIES_FLAGS}")

############################################
# ImageStreamIO (ISIO) setup
############################################
if(MXLIB_USE_ISIO)

    pkg_check_modules(ISIO IMPORTED_TARGET ImageStreamIO)

    if(ISIO_FOUND)
        target_link_libraries(objlib PUBLIC PkgConfig::ISIO)
        target_link_libraries(mxlib-shared PUBLIC PkgConfig::ISIO)
        target_link_libraries(mxlib-static PUBLIC PkgConfig::ISIO)

        include_directories(${ISIO_INCLUDE_DIRS})
        add_compile_definitions("-DMXLIB_MILK")

        set(MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS} ${ISIO_CFLAGS} ${ISIO_CFLAGS_OTHER} -DMXLIB_MILK")
        set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} ${ISIO_CFLAGS} -DMXLIB_MILK")

        linklibdirs_from_list(ISIO_LIBDIRS_FLAGS ISIO_LIBRARY_DIRS)
        set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${ISIO_LIBDIRS_FLAGS}")

        linklibs_from_list(ISIO_LIBRARIES_FLAGS ISIO_LIBRARIES)
        set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} ${ISIO_LIBRARIES_FLAGS}")
    else()
        message(WARNING "No ImageStreamIO support found, setting MXLIB_USE_ISIO to off")
        set(MXLIB_USE_ISIO OFF)
    endif()
endif()

if(MXLIB_USE_ISIO)
    message("mxlib: ImageStreamIO support enabled")
else()
    message("mxlib: ImageStreamIO support not enabled")
endif()

############################################
# CUDA setup
############################################
if(MXLIB_USE_CUDA)

    set_property( TARGET objlib PROPERTY CUDA_ARCHITECTURES ${MXLIB_CUDA_ARCH})
    set_property( TARGET mxlib-shared PROPERTY CUDA_ARCHITECTURES ${MXLIB_CUDA_ARCH})
    set_property( TARGET mxlib-static PROPERTY CUDA_ARCHITECTURES ${MXLIB_CUDA_ARCH})

    enable_language(CUDA)

    include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

    target_link_libraries(mxlib-shared PUBLIC CUDA::cudart)
    target_link_libraries(mxlib-static PUBLIC CUDA::cudart_static)

    libname_from_path(CUDA_cudart_LIBNAME CUDA_cudart_LIBRARY)

    target_link_libraries(mxlib-shared PUBLIC CUDA::cublas)
    target_link_libraries(mxlib-static PUBLIC CUDA::cublas_static)

    libname_from_path(CUDA_cublas_LIBNAME CUDA_cublas_LIBRARY)

    target_link_libraries(mxlib-shared PUBLIC CUDA::cufft)
    target_link_libraries(mxlib-static PUBLIC CUDA::cufft_static)

    libname_from_path(CUDA_cufft_LIBNAME CUDA_cufft_LIBRARY)

    target_link_libraries(mxlib-shared PUBLIC CUDA::curand)
    target_link_libraries(mxlib-static PUBLIC CUDA::curand_static)

    libname_from_path(CUDA_curand_LIBNAME CUDA_curand_LIBRARY)

    target_link_libraries(mxlib-shared PUBLIC CUDA::cusolver)
    target_link_libraries(mxlib-static PUBLIC CUDA::cusolver_static)

    libname_from_path(CUDA_cusolver_LIBNAME CUDA_cusolver_LIBRARY)

    set(MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS} ${MXLIB_CXXVERSION} ${MXLIB_CXXFLAGS} ${MXLIB_OPTIMIZE} ${MXLIB_DEFINES}")

    set(MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS} -DEIGEN_NO_CUDA -DMXLIB_CUDA -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")

    add_compile_definitions(-DEIGEN_NO_CUDA)
    add_compile_definitions(-DMXLIB_CUDA)
    add_compile_definitions(-DMXLIB_MILK)

    #Strip the includes, extra spaces, and prepend -Xcompiler
    string(REGEX REPLACE "-I[^ ]*" "" MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS}")
    string(REGEX REPLACE "[ ]+" " " MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS}")
    string(REPLACE " " " -Xcompiler=" MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS}")

    #and add to the CUDA flags

    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${MXLIB_NVCCXX_FLAGS}")
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --std c++20")

    set(CMAKE_CUDA_STANDARD 20)
    #set(CMAKE_CXX_STANDARD_REQUIRED ON)

    set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} -DEIGEN_NO_CUDA -DMXLIB_CUDA -I${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")

    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} -L${CUDAToolkit_LIBRARY_DIR}")
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} -l${CUDA_cudart_LIBNAME}")
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} -l${CUDA_cublas_LIBNAME}")
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} -l${CUDA_cufft_LIBNAME}")
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} -l${CUDA_curand_LIBNAME}")
    set(MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS} -l${CUDA_cusolver_LIBNAME}")

    if(MXLIB_USE_ISIO)
        # add flag so that MILK gets CUDA support
        set(MXLIB_NVCCXX_FLAGS "${MXLIB_NVCCXX_FLAGS} -DHAVE_CUDA")
        set(MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS} -DHAVE_CUDA")
    endif()

    message("mxlib: CUDA support enabled")


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
#remove stray flags caused by empty strings
string(REGEX REPLACE "-I |-I$" " " MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS}")
string(REGEX REPLACE "-D |-D$ " " " MXLIB_PC_CFLAGS "${MXLIB_PC_CFLAGS}")

string(REPLACE " " ";" MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS}")
list(REMOVE_DUPLICATES MXLIB_PC_LDFLAGS)
string(REPLACE ";" " " MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS}")
#remove stray flags caused by empty strings
string(REGEX REPLACE "-l |-l$ " " " MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS}")
string(REGEX REPLACE "-L |-L$ " " " MXLIB_PC_LDFLAGS "${MXLIB_PC_LDFLAGS}")

message("cflags= ${MXLIB_PC_CFLAGS}")
message("ldflags= ${MXLIB_PC_LDFLAGS}")

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/pkgconfig/${PROJECT_NAME}.pc.in
                ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc
                @ONLY )

install( FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc
         DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig )

install( PROGRAMS gengithead.sh DESTINATION bin)

############################################
# Tests
############################################

if(MXLIB_BUILD_TESTS)
    include(CTest)
    add_subdirectory(tests)
else()
    add_subdirectory(tests EXCLUDE_FROM_ALL)
endif()

add_custom_target(tests
    DEPENDS mxlibTests
)

set(MXLIB_COVERAGE_BUILD_DIR ${CMAKE_BINARY_DIR}/coverage-build)
set(MXLIB_COVERAGE_INFO ${CMAKE_BINARY_DIR}/coverage.info)
set(MXLIB_COVERAGE_FILTERED_INFO ${CMAKE_BINARY_DIR}/coverage_filtered.info)
set(MXLIB_DOC_OUTPUT_DIR ${CMAKE_BINARY_DIR}/doc)
set(MXLIB_DOC_HTML_DIR ${MXLIB_DOC_OUTPUT_DIR}/html)
set(MXLIB_COVERAGE_REPORT_DIR ${MXLIB_DOC_HTML_DIR}/coverage)

find_package(Doxygen QUIET)

if(Doxygen_FOUND)
    set(MXLIB_DOXYGEN_CONFIG_IN ${CMAKE_SOURCE_DIR}/doc/mxlib_doxygen.in)
    set(MXLIB_DOXYGEN_CONFIG ${CMAKE_BINARY_DIR}/mxlib_doxygen)
    set(MXLIB_DOXYGEN_OUTPUT_DIRECTORY ${MXLIB_DOC_OUTPUT_DIR})

    configure_file(${MXLIB_DOXYGEN_CONFIG_IN} ${MXLIB_DOXYGEN_CONFIG} @ONLY)

    if(MXLIB_BUILD_DOCS)
        set(_mxlib_docs_all ALL)
    else()
        set(_mxlib_docs_all)
    endif()

    add_custom_target(docs ${_mxlib_docs_all}
        COMMAND ${CMAKE_COMMAND} -E make_directory ${MXLIB_DOC_OUTPUT_DIR}
        COMMAND ${DOXYGEN_EXECUTABLE} ${MXLIB_DOXYGEN_CONFIG}
        COMMAND ${CMAKE_COMMAND} -DCOVERAGE_REPORT_DIR=${MXLIB_COVERAGE_REPORT_DIR} -P ${CMAKE_SOURCE_DIR}/cmake/ensure_coverage_placeholder.cmake
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/doc
        DEPENDS ${MXLIB_DOXYGEN_CONFIG}
                ${CMAKE_SOURCE_DIR}/doc/mxlib_header.html
                ${CMAKE_SOURCE_DIR}/doc/mxlibDoxygenLayout.xml
                ${CMAKE_SOURCE_DIR}/doc/main.dox
                ${CMAKE_SOURCE_DIR}/doc/coverage_report.dox
                ${CMAKE_SOURCE_DIR}/cmake/ensure_coverage_placeholder.cmake
        USES_TERMINAL
    )
else()
    message(WARNING "doxygen not found. The docs target will not be available.")
endif()

find_program(MXLIB_LCOV_EXECUTABLE lcov)
find_program(MXLIB_GENHTML_EXECUTABLE genhtml)

if(MXLIB_LCOV_EXECUTABLE AND MXLIB_GENHTML_EXECUTABLE)
    add_custom_target(coverage
        COMMAND ${CMAKE_COMMAND} -S ${CMAKE_SOURCE_DIR} -B ${MXLIB_COVERAGE_BUILD_DIR} -DMXLIB_BUILD_TESTS=ON -DMXLIB_BUILD_COVERAGE=ON -DMXLIB_USE_CUDA=OFF -DMXLIB_USE_OPENMP=${MXLIB_USE_OPENMP} -DMXLIB_USE_ISIO=${MXLIB_USE_ISIO} -DMXLIB_USE_FFT_FROM=${MXLIB_USE_FFT_FROM} -DMXLIB_USE_BLAS_FROM=${MXLIB_USE_BLAS_FROM}
        COMMAND ${CMAKE_COMMAND} --build ${MXLIB_COVERAGE_BUILD_DIR} --target mxlibTests -j
        COMMAND ${MXLIB_LCOV_EXECUTABLE} --directory ${MXLIB_COVERAGE_BUILD_DIR} --zerocounters
        COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${MXLIB_COVERAGE_BUILD_DIR} --output-on-failure --timeout ${MXLIB_COVERAGE_TEST_TIMEOUT}
        COMMAND ${MXLIB_LCOV_EXECUTABLE} --directory ${MXLIB_COVERAGE_BUILD_DIR} --capture --output-file ${MXLIB_COVERAGE_INFO}
        COMMAND ${MXLIB_LCOV_EXECUTABLE} --remove ${MXLIB_COVERAGE_INFO} ${CMAKE_SOURCE_DIR}/tests/* /usr/* /sys/* /tty/* --ignore-errors unused,unused --ignore-errors inconsistent,inconsistent --output-file ${MXLIB_COVERAGE_FILTERED_INFO}
        COMMAND ${CMAKE_COMMAND} -E make_directory ${MXLIB_COVERAGE_REPORT_DIR}
        COMMAND ${MXLIB_GENHTML_EXECUTABLE} ${MXLIB_COVERAGE_FILTERED_INFO} --output-directory ${MXLIB_COVERAGE_REPORT_DIR} --title mxlib --hierarchical --merge-aliases --suppress-aliases --filter function --demangle-cpp --css-file ${CMAKE_SOURCE_DIR}/tests/coverage/gcov.css
        USES_TERMINAL
    )
    if(TARGET docs)
        add_dependencies(coverage docs)
    endif()
else()
    message(WARNING "lcov/genhtml not found. The coverage target will not be available.")
endif()

add_custom_target(coverage_clean
    COMMAND ${CMAKE_COMMAND} -E rm -f ${MXLIB_COVERAGE_INFO}
    COMMAND ${CMAKE_COMMAND} -E rm -f ${MXLIB_COVERAGE_FILTERED_INFO}
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${MXLIB_COVERAGE_REPORT_DIR}
    COMMAND ${CMAKE_COMMAND} -E rm -rf ${MXLIB_COVERAGE_BUILD_DIR}
)

if(TARGET docs)
    add_custom_target(docs_clean
        COMMAND ${CMAKE_COMMAND} -E rm -rf ${MXLIB_DOC_OUTPUT_DIR}
    )
endif()

#dump_cmake_variables(".")
