
##########################################################################
# The SOFA Library
#
##########################################################################

set(SOFA_DIR "20231011")

file(GLOB LIBSRCS CONFIGURE_DEPENDS "${SOFA_DIR}/c/src/*.h" "${SOFA_DIR}/c/src/*.c")

add_library(sofa OBJECT ${LIBSRCS})
include_directories(${SOFA_DIR}/c/src/)
set_property(TARGET sofa PROPERTY POSITION_INDEPENDENT_CODE 1)

target_compile_options(sofa PUBLIC -pedantic -Wall -O -fPIC)

set_target_properties(sofa PROPERTIES LINKER_LANGUAGE C)
