add_subdirectory(ao)
add_subdirectory(app)
add_subdirectory(astro)
add_subdirectory(error)
add_subdirectory(base)
add_subdirectory(improc)
add_subdirectory(ioutils)
add_subdirectory(ipc)
add_subdirectory(math)
add_subdirectory(meta)
add_subdirectory(sigproc)
add_subdirectory(sys)
add_subdirectory(wfp)

set(OBJLIB_INCLUDES ${OBJLIB_INCLUDES} include/mxlib.hpp
                                       include/mxError.hpp
                                       include/mxException.hpp
                                       include/mxlib_uncomp_version.h
                                       PARENT_SCOPE )
