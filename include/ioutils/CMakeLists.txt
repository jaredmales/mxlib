add_subdirectory(fits)

set(OBJLIB_INCLUDES ${OBJLIB_INCLUDES} include/ioutils/binVector.hpp
                                       include/ioutils/fileUtils.hpp
                                       include/ioutils/ioutils.hpp
                                       include/ioutils/pout.hpp
                                       include/ioutils/rawBinary.hpp
                                       include/ioutils/readColumns.hpp
                                       include/ioutils/stringUtils.hpp
                                       include/ioutils/textTable.hpp
                                       PARENT_SCOPE )
