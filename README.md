# mxlib

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/a3171d09105445bcb7b0ea29487f3256)](https://www.codacy.com/app/jaredmales/mxlib?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jaredmales/mxlib&amp;utm_campaign=Badge_Grade) [![Build (Linux, macOS)](https://github.com/jaredmales/mxlib/actions/workflows/build-ubuntu-latest.yml/badge.svg)](https://github.com/jaredmales/mxlib/actions/workflows/build-ubuntu-latest.yml)

This is the C/C++ library of Jared Males. It contains code I have developed for data analysis and other tasks, primarily related to astronomy.

The documentation is located here: https://jaredmales.github.io/mxlib-doc/

See the [User's Guide](https://jaredmales.github.io/mxlib-doc/modules.html) for [installation instructions](https://jaredmales.github.io/mxlib-doc/group__installation.html)

## CMake Tests

Build tests on demand (always available, even if `MXLIB_BUILD_TESTS=OFF`):

```bash
cmake -S . -B _build
cmake --build _build --target tests -j
```

Configure with tests enabled:

```bash
cmake -S . -B _build -DMXLIB_BUILD_TESTS=ON
```

With `MXLIB_BUILD_TESTS=ON`, test executables are part of the default build.
With `MXLIB_BUILD_TESTS=OFF`, they are skipped by default and built only via `tests`/`mxlibTest` targets.

Build the aggregate test executable:

```bash
cmake --build _build --target mxlibTest -j
```

Run tests:

```bash
ctest --test-dir _build --output-on-failure
```

Run the aggregate test executable directly:

```bash
cmake --build _build --target mxlibTestRun
```

Build and run a single test source (Makefile.one equivalent):

```bash
cmake -S . -B _build -DMXLIB_BUILD_TESTS=ON -DMXLIB_ONE_TEST=include/math/geo_test.cpp
cmake --build _build --target mxlibTestOne -j
cmake --build _build --target mxlibTestOneRun
```
