Follow these code style and documentation rules exactly.

1) File-Level Documentation
- Each header/source should have a top Doxygen file block:
  - \file
  - \brief
  - \author (if project uses it)

2) Include Guards and Includes
- Match existing project include-guard naming convention.
- Keep include ordering consistent with project style.
- Do not introduce new include style unless project already uses it.

3) Function Declaration Documentation
- Every function declaration must have a brief `///` summary.
- Add a short `/** ... */` details block only when needed.
- Document every parameter inline at declaration site using:
  - `type name /**< [in] description */`
- Apply to normal methods, constructors, slots, and signals.
- Keep return-value docs where project uses them.

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

12) Unit Test Doxygen Structure
- For files under `tests/include/...`, attach each `TEST_CASE` to the matching `*_unit_tests` Doxygen group so test docs mirror the `include/...` folder hierarchy.
- If a needed subgroup is missing, add it to `doc/groupdefs.dox` under `\defgroup unit_tests` in the same hierarchy style as neighboring groups.
- When Doxygen misses a symbol used only in templated/overloaded contexts, add a guarded reference block using `MXLIBTEST_DOXYGEN_REF` near the relevant test case:
  - reference only the target API symbols
  - keep it side-effect free
  - match existing project pattern (`#if` or `#ifdef`) used in nearby tests.

13) Test Doc Cleanup Rules
- Do not use `\test` commands or `\anchor` tags in headers/sources for unit-test linkage.
- Remove legacy test-reference text that only existed to support `\test` entries (including dangling list items like `Tests:` or `[test docs]` notes).
- In test files, Doxygen `\brief` text for a `TEST_CASE` should exactly match the Catch2 test name string; do not prefix with `Scenario:` or `Test Case`.

14) Catch2 Structure for This Repo
- Use `TEST_CASE` + `SECTION`; do not introduce or keep `SCENARIO`, `GIVEN`, `WHEN`, `AND_GIVEN`, or `AND_WHEN`.
- When converting existing tests, preserve behavior and assertion intent; only reshape the structure/macros.

15) Long-Running Coverage Tests
- For known long stochastic tests (for example `psdFilter_test`), include lightweight periodic progress output so timeout failures show where runtime was spent.
- Prefer a per-test CTest timeout override in `tests/CMakeLists.txt` for coverage builds instead of globally increasing all test timeouts.

When you finish:
- Summarize what changed.
- List affected files.
- Note any follow-up items or potential edge cases.
