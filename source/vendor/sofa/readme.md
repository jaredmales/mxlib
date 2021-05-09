This is a copy of the IAU's Standards of Fundamental Astronomy (SOFA) library.

See: http://www.iausofa.org/

The source is unmodified, with the exception of the makefile.  The modifications to it are:

-Addition of the following include statements:
   ```
   -include ../../../../../../mk/MxLib.mk
   -include ../../../../../../mk/Common.mk
   ```
- edit `INSTALL_PATH`, `SOFA_LIB_DIR`, `SOFA_INC_DIR`
- Addition of `-fPIC` is in `CFLAGF` and `CFLAGX` for use in a shared library

`libsofa_c` is statically linked into the `libmxlib.so` binary.

Additionally the headers `sofa.h` and `sofam.h` are copied to `include/astro/sofa`.

You can use libsofa via mxlib by including `<mx/astro/sofa/sofa.h>` and linking `libmxlib.so`.  For c++ you can include `<mx/astro/sofa.hpp>`, which imposes the namespace `sofa`.

To install libsofa_c.a independently, go to the XXXXXX/c/src directory and type `make install`.


