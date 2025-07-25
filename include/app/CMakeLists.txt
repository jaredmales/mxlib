add_subdirectory(optionparser)

set(OBJLIB_INCLUDES ${OBJLIB_INCLUDES} include/app/appConfigurator.hpp
                                       include/app/application.hpp
                                       include/app/clOptions.hpp
                                       include/app/configTarget.hpp
                                       include/app/iniFile.hpp
                                       include/app/ini.hpp
                                       PARENT_SCOPE )
