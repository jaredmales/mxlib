add_subdirectory(analysis)
add_subdirectory(sim)

set(OBJLIB_INCLUDES ${OBJLIB_INCLUDES} include/ao/aoPaths.hpp
                                     include/ao/basis.hpp
                                     include/ao/fourierBasis.hpp
                                     include/ao/influenceFunctions.hpp
                                     include/ao/pupil.hpp
                                     include/ao/zernikeBasis.hpp
                                     PARENT_SCOPE)
