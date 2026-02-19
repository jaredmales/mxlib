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
With `MXLIB_BUILD_TESTS=OFF`, they are skipped by default and built only via `tests`/`mxlibTests` targets.

Build all test executables:

```bash
cmake --build _build --target mxlibTests -j
```

Run tests:

```bash
ctest --test-dir _build --output-on-failure
```

Run the CTest test suite directly:

```bash
cmake --build _build --target mxlibTestRun
```

Build and run a single test source (Makefile.one equivalent):

```bash
cmake -S . -B _build -DMXLIB_BUILD_TESTS=ON -DMXLIB_ONE_TEST=include/math/geo_test.cpp
cmake --build _build --target mxlibTestOne -j
cmake --build _build --target mxlibTestOneRun
```

## Coverage

Coverage generation is integrated into CMake and modeled after the MagAOX flow.

Prerequisites:

```bash
lcov --version
genhtml --version
```

Generate an HTML report:

```bash
cmake -S . -B _build
cmake --build _build --target coverage
```

Optional: tune coverage test timeout (default `300` seconds):

```bash
cmake -S . -B _build -DMXLIB_COVERAGE_TEST_TIMEOUT=600
```

Coverage artifacts are written under `_build/`:

- `_build/coverage.info`
- `_build/coverage_filtered.info`
- `_build/coverage_report/index.html`

Clean coverage artifacts:

```bash
cmake --build _build --target coverage_clean
```

Convenience scripts are also available:

- `tests/coverage/make_coverage`
- `tests/coverage/update_coverage`
