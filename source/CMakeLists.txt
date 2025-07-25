add_subdirectory(ao)
add_subdirectory(app)
add_subdirectory(improc)
add_subdirectory(ioutils)
add_subdirectory(ipc)
add_subdirectory(math)
add_subdirectory(sigproc)
add_subdirectory(sys)

set(OBJLIB_SOURCES ${OBJLIB_SOURCES} source/mxError.cpp source/mxlib.cpp PARENT_SCOPE)

add_subdirectory(vendor/levmar-2.6)
add_subdirectory(vendor/sofa)
