add_subdirectory(sofa)

set(OBJLIB_INCLUDES ${OBJLIB_INCLUDES} include/astro/astroDynamics.hpp
                                       include/astro/astroSpectrum.hpp
                                       include/astro/cahoyAlbedos.hpp
                                       include/astro/constants.hpp
                                       include/astro/mamajek.inc
                                       include/astro/phoenixSpectrum.hpp
                                       include/astro/stars.hpp
                                       include/astro/astroSpectra.hpp
                                       include/astro/blackbody.hpp
                                       include/astro/kepler.hpp
                                       include/astro/orbitUtils.hpp
                                       include/astro/planets.hpp
                                       include/astro/sofa.hpp
                                       include/astro/units.hpp
                                       PARENT_SCOPE )
