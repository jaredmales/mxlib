Follow these code style and documentation rules exactly.

1) File-Level Documentation
- Each header/source should have a top Doxygen file block:
  ```
  /** \file
   * \brief
   *
   * \ingroup
   */
  ```
- the group for a file is usually xxxxxx_files, which is a subgroup of xxxxxx.
- do not include the \author tag

2) Include Guards and Includes
- Match existing project include-guard naming convention.
- Keep include ordering consistent with project style.
- Do not introduce new include style unless project already uses it.

3) Function Declaration Documentation
- Every public or otherwise Doxygen-discoverable function declaration must have a brief `///` summary.
- Add a short `/** ... */` details block only when needed.
- Document every parameter inline at declaration site using:
  - `type name /**< [in] description */`
- Apply to normal methods, constructors, slots, and signals.
- Keep return-value docs where project uses them.
- Document file-local and other deliberately hidden implementation helpers in the same concise prose style, but use
  ordinary `//` and `/* ... */` comments so they do not enter the generated API documentation.

4) Member Variable Documentation
- Document non-trivial class members with `///`.
- Describe role/ownership/state, not just type.

5) Naming and Structure
- Use project member naming convention (e.g. `m_` prefix).
- Keep declaration ordering/grouping stable:
  - public/protected/private
  - slots/signals grouped consistently.
- Leave a blank line between declarations for readability.

6) Header vs Source Placement
- Keep non-trivial definitions out of headers.
- Move implementations to `.cpp` unless intentionally inline.

7) Editing Discipline
- Preserve existing behavior unless explicitly requested.
- When renaming members/APIs, update all dependent call sites.
- Keep changes minimal and scoped.

8) Formatting and Verification
- Run `clang-format` on touched files.
- Ensure docs and naming are consistent after formatting.
- Report any places where project style is ambiguous before making assumptions.

9) Doxygen Named Section Ordering
- For classes that expose configuration via member data + accessors, keep named sections split into:
  - `... - Data` for protected/private member state
  - `...` (without `- Data`) for public access functions
- Place the `... - Data` section before the corresponding public accessor section.

10) Header Declaration Parameter Docs
- In headers, prefer inline parameter documentation on declarations (`type name /**< ... */`) rather than separate `\param` lists, unless there is a specific reason to deviate.

11) PR Prompt Attribution
- At the top of PR descriptions, include an explicit attribution line when work was performed with Codex.
- Preferred format:
  - `This work was performed by GPT-5.3-Codex in response to the prompt: "...".`
- Include the primary user prompt verbatim (or a faithful condensed version if it is extremely long).

12) Unit Test Documentation
- Use Catch2 `TEST_CASE` for every top-level test. Do not add or retain `SCENARIO` or `SCENARIO_METHOD` declarations.
- Add a brief Doxygen block immediately before every Catch2 `TEST_CASE`, which matches the TEST_CASE name exactly. This is so that the brief matches the Catch2 test name to run it.
- Use \ingroup to place this TEST_CASE under the appropriate unit testing group
- Example:
  ```
  /// Test error_t boolean operators and functions
  /**
   * \ingroup error_error_unit_tests
   */
  TEST_CASE( "Test error_t boolean operators and functions", "[error::error]" )
  ```
- When additional explanation is needed, in the main documentation block State the behavior being verified and identify the real production API under test.
- Let Doxygen discover real calls in the test body so the test appears in each production API's `Referenced by` list.
- Do not use `\test` or prose-only `\ref` commands to manufacture test-to-API links.
- Keep Catch2 assertion macros lexically inside the `TEST_CASE` or `SECTION` body; do not place `REQUIRE`, `CHECK`,
  or related assertions in helper or wrapper functions.
- Helpers may prepare inputs and compute comparisons, but return named values to the test body for local assertion.
  Prefer numeric actual/error and expected/tolerance members over a single opaque boolean so Catch2 reports useful
  expanded values.
- Use `CAPTURE` or `INFO` at the local assertion site to identify loop parameters, backend choices, and other case
  context.

13) Preserve Doxygen Links Through Test Harnesses
- Preserve Doxygen links to the real production APIs when test fixtures, wrappers, namespaces, macros, or private-access techniques prevent automatic symbol linking.
- Add explicit Doxygen-only code references to the production symbols inside the relevant test body when direct calls are otherwise hidden.
- Guard reference-only code with `#ifdef __DOXY_ONLY__` so it need not compile, and use raw calls or member references that Doxygen can add to the production symbol's `Referenced by` list.
- Hide harness-only helpers from generated documentation with `\cond` and `\endcond` when they would dominate or obscure production API links.
- Disable `clang-format` around non-compiling Doxygen-only reference blocks when necessary.

14) Function Argument Ordering
- For functions that write caller-owned results, place output-only arguments before the semantic inputs.
- This keeps call sites analogous to `y = f(x)`, with the value being written on the left and the inputs on the right.
- Reusable workspaces and other stateful input/output support objects may follow the semantic inputs when that keeps the
  numerical call readable and matches the surrounding API family.
- Preserve required ordering for language operators, callbacks, overrides, and external or compatibility-constrained APIs.

When you finish:
- Summarize what changed.
- List affected files.
- Note any follow-up items or potential edge cases.
