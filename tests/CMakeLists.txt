#################################################################
#                      cmake configuration for mxlib tests
#################################################################

set(MXLIB_TEST_MAIN_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/testMain.cpp)

# Catch2's main translation unit is identical for every per-source test
# executable. Compile it once and reuse the object to keep the complete header
# inventory practical to build.
add_library(mxlibTestMain OBJECT ${MXLIB_TEST_MAIN_SOURCE})
target_link_libraries(mxlibTestMain PRIVATE mxlib-shared)
if(NOT MXLIB_BUILD_TESTS)
    set_target_properties(mxlibTestMain PROPERTIES EXCLUDE_FROM_ALL TRUE EXCLUDE_FROM_DEFAULT_BUILD TRUE)
endif()

# Build the production inventory from the source tree so adding a public header
# or implementation also requires adding its corresponding test source.
file(GLOB_RECURSE MXLIB_PRODUCTION_SOURCES CONFIGURE_DEPENDS
    ${PROJECT_SOURCE_DIR}/include/*.hpp
    ${PROJECT_SOURCE_DIR}/include/*.h
    ${PROJECT_SOURCE_DIR}/include/*.inc
    ${PROJECT_SOURCE_DIR}/include/*.cpp
    ${PROJECT_SOURCE_DIR}/source/*.hpp
    ${PROJECT_SOURCE_DIR}/source/*.h
    ${PROJECT_SOURCE_DIR}/source/*.inc
    ${PROJECT_SOURCE_DIR}/source/*.cpp
)

# These files are vendored or generated outputs. Authored generator-tool
# implementations remain in the strict project-owned source inventory.
list(FILTER MXLIB_PRODUCTION_SOURCES EXCLUDE REGEX "/include/app/inih/")
list(FILTER MXLIB_PRODUCTION_SOURCES EXCLUDE REGEX "/include/app/optionparser/")
list(FILTER MXLIB_PRODUCTION_SOURCES EXCLUDE REGEX "/include/astro/sofa/")
list(FILTER MXLIB_PRODUCTION_SOURCES EXCLUDE REGEX "/include/error/generators/[^/]+_info\\.hpp$")
list(FILTER MXLIB_PRODUCTION_SOURCES EXCLUDE REGEX "/include/error/error_t\\.hpp$")
list(FILTER MXLIB_PRODUCTION_SOURCES EXCLUDE REGEX "/include/mxlib_uncomp_version\\.h$")
list(FILTER MXLIB_PRODUCTION_SOURCES EXCLUDE REGEX "/source/app/optionparser\\.cpp$")
list(FILTER MXLIB_PRODUCTION_SOURCES EXCLUDE REGEX "/source/mxlib_(comp|uncomp)_version\\.h$")
list(FILTER MXLIB_PRODUCTION_SOURCES EXCLUDE REGEX "/source/vendor/")
list(SORT MXLIB_PRODUCTION_SOURCES)

function(mxlib_test_source_for_production PRODUCTION_SOURCE OUTPUT_VARIABLE)
    file(RELATIVE_PATH _mxlib_production_rel "${PROJECT_SOURCE_DIR}" "${PRODUCTION_SOURCE}")

    # mamajek.inc is a non-standalone data initializer used by astroSpectrum.
    if(_mxlib_production_rel STREQUAL "include/astro/mamajek.inc")
        set(_mxlib_expected_test "${CMAKE_CURRENT_SOURCE_DIR}/include/astro/astroSpectrum_test.cpp")
    elseif(_mxlib_production_rel MATCHES "^(include|source)/(.+)\\.(hpp|h|inc|cpp)$")
        set(_mxlib_expected_test "${CMAKE_CURRENT_SOURCE_DIR}/include/${CMAKE_MATCH_2}_test.cpp")
    else()
        message(FATAL_ERROR "No mxlib test-path convention for production source: ${PRODUCTION_SOURCE}")
    endif()

    set(${OUTPUT_VARIABLE} "${_mxlib_expected_test}" PARENT_SCOPE)
endfunction()

set(MXLIB_EXPECTED_TEST_SOURCES)
set(MXLIB_MISSING_TEST_MAPPINGS)
foreach(_mxlib_production_source IN LISTS MXLIB_PRODUCTION_SOURCES)
    mxlib_test_source_for_production("${_mxlib_production_source}" _mxlib_expected_test)
    list(APPEND MXLIB_EXPECTED_TEST_SOURCES "${_mxlib_expected_test}")

    if(NOT EXISTS "${_mxlib_expected_test}")
        file(RELATIVE_PATH _mxlib_production_rel "${PROJECT_SOURCE_DIR}" "${_mxlib_production_source}")
        file(RELATIVE_PATH _mxlib_test_rel "${CMAKE_CURRENT_SOURCE_DIR}" "${_mxlib_expected_test}")
        list(APPEND MXLIB_MISSING_TEST_MAPPINGS "${_mxlib_production_rel} -> tests/${_mxlib_test_rel}")
    endif()
endforeach()
list(REMOVE_DUPLICATES MXLIB_EXPECTED_TEST_SOURCES)
list(SORT MXLIB_EXPECTED_TEST_SOURCES)

if(MXLIB_MISSING_TEST_MAPPINGS)
    list(JOIN MXLIB_MISSING_TEST_MAPPINGS "\n  " _mxlib_missing_test_message)
    message(FATAL_ERROR "Missing mxlib production test sources:\n  ${_mxlib_missing_test_message}")
endif()

file(GLOB_RECURSE MXLIB_DISCOVERED_TEST_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/include/*_test.cpp
)
list(SORT MXLIB_DISCOVERED_TEST_SOURCES)

# Generated error_t.hpp has a test but is intentionally absent from the
# project-owned production-source inventory above.
set(MXLIB_EXTRA_TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/error/error_t_test.cpp
)

set(MXLIB_CUDA_TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/sim/deformableMirror_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/sim/directPhaseReconstructor_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/sim/simulatedAOSystem_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/sim/wooferTweeterDM_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/sim/wooferTweeterFilter_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/sim/wooferTweeterReconstructor_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/error/generators/generate_cudaError_info_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/cuda/cublasHandle_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/cuda/cudaPtr_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/cuda/cusolverDnHandle_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/cuda/cusolverDnParams_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/cuda/syevdxT_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/cuda/templateCublas_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/cuda/templateCuda_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/cuda/templateCufft_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/cuda/templateCurand_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/ft/fftTcuda_test.cpp
)

set(MXLIB_ISIO_TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/sim/turbAtmosphere_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/sim/turbSubHarmonic_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/improc/milkImage_test.cpp
)

set(MXLIB_CUDA_ISIO_TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/sim/pyramidSensor_test.cpp
)

# These ownership placeholders remain part of the strict production inventory
# and generated test documentation, but are intentionally not compiled. The AO
# headers are obsolete duplicate implementations whose dependencies and APIs no
# longer exist. psdFilterCuda additionally depends on the obsolete
# templateCudaPtr and CUDA-samples helper headers, removed pointwiseMul and
# scalarMul APIs, and the pre-rank psdFilter specialization.
set(MXLIB_LEGACY_HEADER_TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/sim/directPhaseReconstructorD_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/sim/directPhaseReconstructorOrtho_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/sim/pyramidSensorSepQuad_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/sigproc/psdFilterCuda_test.cpp
)

set(MXLIB_ALL_GATED_TEST_SOURCES
    ${MXLIB_CUDA_TEST_SOURCES}
    ${MXLIB_ISIO_TEST_SOURCES}
    ${MXLIB_CUDA_ISIO_TEST_SOURCES}
)
list(LENGTH MXLIB_ALL_GATED_TEST_SOURCES _mxlib_gated_count)
list(REMOVE_DUPLICATES MXLIB_ALL_GATED_TEST_SOURCES)
list(LENGTH MXLIB_ALL_GATED_TEST_SOURCES _mxlib_unique_gated_count)
if(NOT _mxlib_gated_count EQUAL _mxlib_unique_gated_count)
    message(FATAL_ERROR "An mxlib test source appears in more than one feature-gated test set")
endif()

set(MXLIB_ACCOUNTED_TEST_SOURCES
    ${MXLIB_EXPECTED_TEST_SOURCES}
    ${MXLIB_EXTRA_TEST_SOURCES}
)
list(REMOVE_DUPLICATES MXLIB_ACCOUNTED_TEST_SOURCES)

foreach(_mxlib_expected_test IN LISTS MXLIB_ACCOUNTED_TEST_SOURCES)
    if(NOT _mxlib_expected_test IN_LIST MXLIB_DISCOVERED_TEST_SOURCES)
        message(FATAL_ERROR "Expected mxlib test source was not discovered: ${_mxlib_expected_test}")
    endif()
endforeach()

foreach(_mxlib_discovered_test IN LISTS MXLIB_DISCOVERED_TEST_SOURCES)
    if(NOT _mxlib_discovered_test IN_LIST MXLIB_ACCOUNTED_TEST_SOURCES)
        message(FATAL_ERROR "Test source has no production-source mapping or approved exception: ${_mxlib_discovered_test}")
    endif()

    file(STRINGS "${_mxlib_discovered_test}" _mxlib_test_case_lines REGEX "^[ \t]*TEST_CASE[ \t]*\\(")
    if(NOT _mxlib_test_case_lines)
        message(FATAL_ERROR "mxlib Catch2 test source has no TEST_CASE declaration: ${_mxlib_discovered_test}")
    endif()
endforeach()

foreach(_mxlib_gated_test IN LISTS MXLIB_ALL_GATED_TEST_SOURCES)
    if(NOT _mxlib_gated_test IN_LIST MXLIB_DISCOVERED_TEST_SOURCES)
        message(FATAL_ERROR "Feature-gated mxlib test source was not discovered: ${_mxlib_gated_test}")
    endif()
endforeach()

foreach(_mxlib_legacy_test IN LISTS MXLIB_LEGACY_HEADER_TEST_SOURCES)
    if(NOT _mxlib_legacy_test IN_LIST MXLIB_DISCOVERED_TEST_SOURCES)
        message(FATAL_ERROR "Legacy-header mxlib test source was not discovered: ${_mxlib_legacy_test}")
    endif()
    if(_mxlib_legacy_test IN_LIST MXLIB_ALL_GATED_TEST_SOURCES)
        message(FATAL_ERROR "Legacy-header mxlib test source is also feature-gated: ${_mxlib_legacy_test}")
    endif()
endforeach()

# Standalone developer programs are inventoried separately and are not linked
# into Catch2 executables.
set(MXLIB_MANUAL_TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/manual/directPhaseSensor.cpp
)
file(GLOB_RECURSE MXLIB_DISCOVERED_MANUAL_TEST_SOURCES CONFIGURE_DEPENDS
    ${CMAKE_CURRENT_SOURCE_DIR}/manual/*.cpp
)
foreach(_mxlib_manual_test IN LISTS MXLIB_MANUAL_TEST_SOURCES)
    if(NOT _mxlib_manual_test IN_LIST MXLIB_DISCOVERED_MANUAL_TEST_SOURCES)
        message(FATAL_ERROR "Inventoried mxlib manual test source was not discovered: ${_mxlib_manual_test}")
    endif()
endforeach()
foreach(_mxlib_manual_test IN LISTS MXLIB_DISCOVERED_MANUAL_TEST_SOURCES)
    if(NOT _mxlib_manual_test IN_LIST MXLIB_MANUAL_TEST_SOURCES)
        message(FATAL_ERROR "Uninventoried mxlib manual test source: ${_mxlib_manual_test}")
    endif()
endforeach()

# Every non-gated discovered test is built by default; optional tests are added
# only when their required dependencies are configured.
set(MXLIB_TEST_SOURCES ${MXLIB_DISCOVERED_TEST_SOURCES})
list(REMOVE_ITEM MXLIB_TEST_SOURCES ${MXLIB_ALL_GATED_TEST_SOURCES})
list(REMOVE_ITEM MXLIB_TEST_SOURCES ${MXLIB_LEGACY_HEADER_TEST_SOURCES})
if(MXLIB_USE_CUDA)
    list(APPEND MXLIB_TEST_SOURCES ${MXLIB_CUDA_TEST_SOURCES})
endif()
if(MXLIB_USE_ISIO)
    list(APPEND MXLIB_TEST_SOURCES ${MXLIB_ISIO_TEST_SOURCES})
endif()
if(MXLIB_USE_CUDA AND MXLIB_USE_ISIO)
    list(APPEND MXLIB_TEST_SOURCES ${MXLIB_CUDA_ISIO_TEST_SOURCES})
endif()
list(REMOVE_DUPLICATES MXLIB_TEST_SOURCES)
list(SORT MXLIB_TEST_SOURCES)

set(MXLIB_TEST_TARGETS)

function(mxlib_add_test_executable TEST_SOURCE)
    file(RELATIVE_PATH _mxlib_test_rel "${CMAKE_CURRENT_SOURCE_DIR}" "${TEST_SOURCE}")
    string(REPLACE "/" "_" _mxlib_test_name "${_mxlib_test_rel}")
    string(REPLACE "." "_" _mxlib_test_name "${_mxlib_test_name}")
    set(_mxlib_test_target "mxlibTest_${_mxlib_test_name}")

    add_executable(${_mxlib_test_target} $<TARGET_OBJECTS:mxlibTestMain> ${TEST_SOURCE})
    target_link_libraries(${_mxlib_test_target} PRIVATE mxlib-shared)

    if(NOT MXLIB_BUILD_TESTS)
        set_target_properties(${_mxlib_test_target} PROPERTIES EXCLUDE_FROM_ALL TRUE EXCLUDE_FROM_DEFAULT_BUILD TRUE)
    endif()

    if(MXLIB_BUILD_TESTS AND BUILD_TESTING)
        add_test(NAME ${_mxlib_test_target} COMMAND ${_mxlib_test_target})
    endif()

    set(MXLIB_TEST_TARGETS ${MXLIB_TEST_TARGETS} ${_mxlib_test_target} PARENT_SCOPE)
endfunction()

foreach(_mxlib_test_source IN LISTS MXLIB_TEST_SOURCES)
    mxlib_add_test_executable("${_mxlib_test_source}")
endforeach()

add_custom_target(mxlibTests DEPENDS ${MXLIB_TEST_TARGETS})

if(MXLIB_BUILD_TESTS AND BUILD_TESTING)
    add_custom_target(mxlibTestRun
        COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_BINARY_DIR} --output-on-failure
        DEPENDS mxlibTests
        USES_TERMINAL
    )
endif()

set(MXLIB_ONE_TEST
    ""
    CACHE STRING
          "Single test source to build, relative to tests/ (e.g. include/math/geo_test.cpp), similar to tests/Makefile.one"
)

if(MXLIB_ONE_TEST)
    if(IS_ABSOLUTE "${MXLIB_ONE_TEST}")
        set(MXLIB_ONE_TEST_SOURCE "${MXLIB_ONE_TEST}")
    else()
        set(MXLIB_ONE_TEST_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/${MXLIB_ONE_TEST}")
    endif()

    if(MXLIB_ONE_TEST_SOURCE MATCHES "\\.o$")
        string(REGEX REPLACE "\\.o$" ".cpp" MXLIB_ONE_TEST_SOURCE "${MXLIB_ONE_TEST_SOURCE}")
    elseif(MXLIB_ONE_TEST_SOURCE MATCHES "\\.hpp$")
        string(REGEX REPLACE "\\.hpp$" ".cpp" MXLIB_ONE_TEST_SOURCE "${MXLIB_ONE_TEST_SOURCE}")
    elseif(NOT MXLIB_ONE_TEST_SOURCE MATCHES "\\.cpp$")
        string(APPEND MXLIB_ONE_TEST_SOURCE ".cpp")
    endif()

    if(NOT EXISTS "${MXLIB_ONE_TEST_SOURCE}")
        message(FATAL_ERROR "MXLIB_ONE_TEST source does not exist: ${MXLIB_ONE_TEST_SOURCE}")
    endif()

    add_executable(mxlibTestOne $<TARGET_OBJECTS:mxlibTestMain> ${MXLIB_ONE_TEST_SOURCE})
    target_link_libraries(mxlibTestOne PRIVATE mxlib-shared)
    if(NOT MXLIB_BUILD_TESTS)
        set_target_properties(mxlibTestOne PROPERTIES EXCLUDE_FROM_ALL TRUE EXCLUDE_FROM_DEFAULT_BUILD TRUE)
    endif()

    add_custom_target(mxlibTestOneRun
        COMMAND $<TARGET_FILE:mxlibTestOne>
        DEPENDS mxlibTestOne
        USES_TERMINAL
    )

    if(MXLIB_BUILD_TESTS AND BUILD_TESTING)
        add_test(NAME mxlibTestOne COMMAND mxlibTestOne)
    endif()
endif()
