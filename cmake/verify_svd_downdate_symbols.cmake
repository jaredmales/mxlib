# Verify that libmxlib exposes only the ABI-v2 SVD-deletion shared-library boundary.

if(NOT DEFINED MXLIB_NM_EXECUTABLE OR NOT DEFINED MXLIB_SHARED_LIBRARY)
    message(FATAL_ERROR "The SVD-deletion symbol check requires nm and the shared-library path")
endif()

execute_process(
    COMMAND "${MXLIB_NM_EXECUTABLE}" -D -C --defined-only "${MXLIB_SHARED_LIBRARY}"
    RESULT_VARIABLE MXLIB_NM_RESULT
    OUTPUT_VARIABLE MXLIB_DYNAMIC_SYMBOLS
    ERROR_VARIABLE MXLIB_NM_ERROR
)
if(NOT MXLIB_NM_RESULT EQUAL 0)
    message(FATAL_ERROR "nm failed while checking ${MXLIB_SHARED_LIBRARY}: ${MXLIB_NM_ERROR}")
endif()

set(MXLIB_REQUIRED_SVD_DELETION_SYMBOLS
    "mx::math::svdDeletionResult<double, mx::math::svdDeletionAbiV2Tag>::svdDeletionResult()"
    "mx::math::svdDeletionWorkspace<double, mx::math::svdDeletionAbiV2Tag>::svdDeletionWorkspace()"
    "mx::math::detail::svdRemoveRowsAbiV2<double>"
)
foreach(MXLIB_REQUIRED_SYMBOL IN LISTS MXLIB_REQUIRED_SVD_DELETION_SYMBOLS)
    string(FIND "${MXLIB_DYNAMIC_SYMBOLS}" "${MXLIB_REQUIRED_SYMBOL}" MXLIB_REQUIRED_SYMBOL_POSITION)
    if(MXLIB_REQUIRED_SYMBOL_POSITION EQUAL -1)
        message(FATAL_ERROR "Required ABI-v2 symbol is missing: ${MXLIB_REQUIRED_SYMBOL}")
    endif()
endforeach()

foreach(MXLIB_SCALAR_TYPE IN ITEMS float double)
    foreach(MXLIB_HANDLE_TYPE IN ITEMS svdDeletionResult svdDeletionWorkspace)
        set(MXLIB_UNTAGGED_HANDLE "mx::math::${MXLIB_HANDLE_TYPE}<${MXLIB_SCALAR_TYPE}>::")
        string(FIND "${MXLIB_DYNAMIC_SYMBOLS}" "${MXLIB_UNTAGGED_HANDLE}" MXLIB_UNTAGGED_HANDLE_POSITION)
        if(NOT MXLIB_UNTAGGED_HANDLE_POSITION EQUAL -1)
            message(FATAL_ERROR "Obsolete untagged SVD-deletion handle symbol remains: ${MXLIB_UNTAGGED_HANDLE}")
        endif()
    endforeach()

    set(MXLIB_TAGGED_RESULT
        "mx::math::svdDeletionResult<${MXLIB_SCALAR_TYPE}, mx::math::svdDeletionAbiV2Tag>")
    set(MXLIB_TAGGED_WORKSPACE
        "mx::math::svdDeletionWorkspace<${MXLIB_SCALAR_TYPE}, mx::math::svdDeletionAbiV2Tag>")
    set(MXLIB_FORBIDDEN_HANDLE_ADAPTERS
        "${MXLIB_TAGGED_RESULT}::prepare("
        "${MXLIB_TAGGED_RESULT}::singularValues()"
        "${MXLIB_TAGGED_RESULT}::squaredSingularValues()"
        "${MXLIB_TAGGED_RESULT}::rotation()"
        "${MXLIB_TAGGED_WORKSPACE}::prepare("
    )
    foreach(MXLIB_FORBIDDEN_ADAPTER IN LISTS MXLIB_FORBIDDEN_HANDLE_ADAPTERS)
        string(FIND "${MXLIB_DYNAMIC_SYMBOLS}" "${MXLIB_FORBIDDEN_ADAPTER}" MXLIB_FORBIDDEN_ADAPTER_POSITION)
        if(NOT MXLIB_FORBIDDEN_ADAPTER_POSITION EQUAL -1)
            message(FATAL_ERROR "Consumer-side SVD-deletion adapter was exported: ${MXLIB_FORBIDDEN_ADAPTER}")
        endif()
    endforeach()
endforeach()

string(REGEX MATCH "[^\n]*(svdDeletion|SvdDeletion|svdRemove)[^\n]*Eigen::Ref[^\n]*"
             MXLIB_OBSOLETE_EIGEN_SYMBOL "${MXLIB_DYNAMIC_SYMBOLS}")
if(MXLIB_OBSOLETE_EIGEN_SYMBOL)
    message(FATAL_ERROR "Obsolete Eigen-facing SVD-deletion symbol remains: ${MXLIB_OBSOLETE_EIGEN_SYMBOL}")
endif()
