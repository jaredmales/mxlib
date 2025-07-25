##########################################################################
# The levmar Library
#
##########################################################################

set(LIBSRCS lm.c Axb.c misc.c lmlec.c lmbc.c lmblec.c lmbleic.c)

add_library(levmar OBJECT ${LIBSRCS})
set_property(TARGET levmar PROPERTY POSITION_INDEPENDENT_CODE 1)

set(CONFIGFLAGS "-ULINSOLVERS_RETAIN_MEMORY")

#don't use target_compile_definitions b/c it adds -D, which makes -U not work
add_definitions(${CONFIGFLAGS} ${ARCHFLAGS})

target_compile_options(levmar PUBLIC -O3 -ffast-math -funroll-loops -Wall -fPIC)

set_target_properties(levmar PROPERTIES LINKER_LANGUAGE C)

