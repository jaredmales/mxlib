
if(MXLIB_USE_CUDA)

    set(OBJLIB_SOURCES  ${OBJLIB_SOURCES}
                        source/math/cuda/templateCublas.cu
                        source/math/cuda/templateCufft.cpp
                        source/math/cuda/templateCurand.cpp PARENT_SCOPE)

endif()
