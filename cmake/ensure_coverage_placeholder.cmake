if(NOT DEFINED COVERAGE_REPORT_DIR OR COVERAGE_REPORT_DIR STREQUAL "")
    message(FATAL_ERROR "COVERAGE_REPORT_DIR is required")
endif()

set(_coverage_index "${COVERAGE_REPORT_DIR}/index.html")

# Never overwrite a real coverage report.
if(EXISTS "${_coverage_index}")
    message(STATUS "Coverage index already exists at ${_coverage_index}; leaving it unchanged.")
    return()
endif()

file(MAKE_DIRECTORY "${COVERAGE_REPORT_DIR}")

file(WRITE "${_coverage_index}" "<!doctype html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <title>mxlib coverage report</title>
  <link href=\"../doxygen.css\" rel=\"stylesheet\" type=\"text/css\">
  <link href=\"../doxygen-awesome.css\" rel=\"stylesheet\" type=\"text/css\">
  <link href=\"../doxygen-awesome-sidebar-only.css\" rel=\"stylesheet\" type=\"text/css\">
  <link href=\"../doxygen-awesome-sidebar-only-darkmode-toggle.css\" rel=\"stylesheet\" type=\"text/css\">
  <style>
    body {
      font-family: var(--font-family, sans-serif);
      color: var(--page-foreground-color, #2f4153);
      background: var(--page-background-color, #ffffff);
      margin: 2rem;
    }
    h1 { margin: 0 0 0.75rem 0; }
    p { margin: 0.5rem 0; }
    code {
      font-family: var(--font-family-monospace, monospace);
      color: var(--page-foreground-color, #2f4153);
    }
  </style>
</head>
<body>
  <h1>Coverage report has not been generated.</h1>
  <p>Run <code>cmake --build _build --target coverage</code> to generate it.</p>
</body>
</html>
")

message(STATUS "Wrote coverage placeholder at ${_coverage_index}")
