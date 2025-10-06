add_subdirectory(cuda)
add_subdirectory(fit)
add_subdirectory(ft)
add_subdirectory(func)
add_subdirectory(plot)

set(OBJLIB_INCLUDES ${OBJLIB_INCLUDES} include/math/constants.hpp
                                       include/math/eigenLapack.hpp
                                       include/math/geo.hpp
                                       include/math/gslInterpolation.hpp
                                       include/math/gslInterpolator.hpp
                                       include/math/histogramUniform.hpp
                                       include/math/logInterpolator.hpp
                                       include/math/logRadProfIntegrator.hpp
                                       include/math/math.hpp
                                       include/math/point2D.hpp
                                       include/math/radprofIntegral.hpp
                                       include/math/randomSeed.hpp
                                       include/math/randomT.hpp
                                       include/math/roots.hpp
                                       include/math/templateBLAS.hpp
                                       include/math/templateLapack.hpp
                                       include/math/vectorUtils.hpp
                                       PARENT_SCOPE )
