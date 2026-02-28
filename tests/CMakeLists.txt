#################################################################
#                      cmake configuration for mxlib tests
#################################################################

set(MXLIB_TEST_MAIN_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/testMain.cpp)

set(MXLIB_TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/analysis/aoAtmosphere_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/analysis/aoSystem_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/analysis/aoPSDs_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/app/appConfigurator_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/astro/astroDynamics_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/error/error_t_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/error/error_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/error/exception_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/improc/imageMasks_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ioutils/fileUtils_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ioutils/readColumns_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ioutils/stringUtils_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ioutils/fits/fitsHeaderCard_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ioutils/fits/fitsFile_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ioutils/fits/fitsHeader_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/geo_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/func/moffat_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/ft/fftT_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/ft/fftwEnvironment_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/ft/ftTypes_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/templateBLAS_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/templateLapack_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/math/randomT_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/sigproc/circularBuffer_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/sigproc/psdUtils_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/sigproc/psdFilter_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/sigproc/zernike_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/sigproc/signalWindows_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/improc/imageTransforms_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/improc/imageUtils_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/improc/imageFilters_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/sys/timeUtils_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/improc/imageXCorrFFT_test.cpp
)

if(MXLIB_USE_CUDA AND MXLIB_ENABLE_GPU_RUNTIME_TESTS)
    list(APPEND MXLIB_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/include/ao/sim/pyramidSensor_test.cpp)
    list(APPEND MXLIB_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/include/wfp/fraunhoferPropagator_test.cpp)
    list(APPEND MXLIB_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/include/wfp/imagingUtils_test.cpp)
    list(APPEND MXLIB_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/include/math/cuda/templateCublas_test.cpp)
    list(APPEND MXLIB_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/include/math/ft/fftTcuda_test.cpp)
endif()

set(MXLIB_TEST_TARGETS)

function(mxlib_add_test_executable TEST_SOURCE)
    file(RELATIVE_PATH _mxlib_test_rel "${CMAKE_CURRENT_SOURCE_DIR}" "${TEST_SOURCE}")
    string(REPLACE "/" "_" _mxlib_test_name "${_mxlib_test_rel}")
    string(REPLACE "." "_" _mxlib_test_name "${_mxlib_test_name}")
    set(_mxlib_test_target "mxlibTest_${_mxlib_test_name}")

    add_executable(${_mxlib_test_target} ${MXLIB_TEST_MAIN_SOURCE} ${TEST_SOURCE})
    target_link_libraries(${_mxlib_test_target} PRIVATE mxlib-shared)

    if(NOT MXLIB_BUILD_TESTS)
        set_target_properties(${_mxlib_test_target} PROPERTIES EXCLUDE_FROM_ALL TRUE EXCLUDE_FROM_DEFAULT_BUILD TRUE)
    endif()

    if(MXLIB_BUILD_TESTS AND BUILD_TESTING)
        add_test(NAME ${_mxlib_test_target} COMMAND ${_mxlib_test_target})
        if(_mxlib_test_target STREQUAL "mxlibTest_include_sigproc_psdFilter_test_cpp" AND MXLIB_BUILD_COVERAGE)
            # psdFilter coverage runtime can vary with system load; allow extra time to avoid transient timeouts.
            math(EXPR _mxlib_psd_filter_timeout "${MXLIB_COVERAGE_TEST_TIMEOUT} * 3")
            set_tests_properties(${_mxlib_test_target} PROPERTIES TIMEOUT ${_mxlib_psd_filter_timeout})
        endif()
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

    add_executable(mxlibTestOne ${MXLIB_TEST_MAIN_SOURCE} ${MXLIB_ONE_TEST_SOURCE})
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
