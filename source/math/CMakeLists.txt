add_subdirectory(cuda)
add_subdirectory(fit)
add_subdirectory(ft)
add_subdirectory(func)
add_subdirectory(plot)

set(OBJLIB_SOURCES ${OBJLIB_SOURCES} source/math/gslInterpolator.cpp
                                     source/math/templateBLAS.cpp
                                     source/math/templateLapack.cpp PARENT_SCOPE)

