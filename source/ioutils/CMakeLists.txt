add_subdirectory(fits)

set(OBJLIB_SOURCES ${OBJLIB_SOURCES}
                   source/ioutils/fileUtils.cpp
                   source/ioutils/stringUtils.cpp
                   source/ioutils/textTable.cpp PARENT_SCOPE)

