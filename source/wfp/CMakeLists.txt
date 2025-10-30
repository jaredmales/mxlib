
set(OBJLIB_SOURCES ${OBJLIB_SOURCES} source/wfp/imagingUtils.cpp PARENT_SCOPE)

if(MXLIB_USE_CUDA)

set(OBJLIB_SOURCES ${OBJLIB_SOURCES} source/wfp/imagingUtils.cu source/wfp/imagingUtils.cpp PARENT_SCOPE)

endif()
