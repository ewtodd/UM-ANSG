# Mandatory Rules
<!---->
## In all languages
<!---->
- Prefer slightly verbose, self-explanatory code over terse code that needs
  comments to be understood.
- Keep comments to only what explains something non-obvious.
This governs
  inline commentary inside function bodies.
  Doxygen blocks on public
  declarations are exempt and have their own rules under "API documentation
  (Doxygen)" below.
- Never embed a literal `\n` inside a string or print argument.
A line break
  is always its own explicit statement.
  In C++/ROOT, use
  `std::cout << ... << std::endl;`.
  In Python, split output into separate
  `print()` calls, and use a bare `print()` for a blank line rather than
  appending `\n`.
<!---->
## Project layout
<!---->
- This is an umbrella repo: one flake, one dev shell per project.
The
  projects are `73Ge`, `78mBr`, `YAG-PreProcessing`, and `YAP-PSD`.
- The 73Ge shell starts in `73Ge/czt`, the CZT analysis; the PSI CdTe
  analysis of the same line is in `73Ge/cdte`, and header-only code both
  use is in `73Ge/include`.
  The 78mBr shell starts in `78mBr/macros`; the
  YAG and YAP shells start in the project root.
  Each shell wires that project's
  `include/` into `CPLUS_INCLUDE_PATH`/`ROOT_INCLUDE_PATH` and `lib/` into
  `LD_LIBRARY_PATH`.
- Backend per project: `73Ge` uses the CUDA Analysis-Utilities build
  (`au.cuda` + `au.rootCuda`, `AU_ROOFIT_BACKEND_CUDA=1`); the other three
  use the default CPU build (`au.default`). The CUDA and CPU builds are
  distinct derivations pinned to one Analysis-Utilities commit.
- Run `nix develop .#<Project>` from the repo root. Never `nix-shell`.
<!---->
## Secrets (agenix)
<!---->
- `73Ge/include/BEF.hpp`, `73Ge/src/BEF.cpp`, and
  `73Ge/czt/ConvertBEF.cpp` are agenix-encrypted sources that the 73Ge
  shell decrypts into the tree on entry. They are gitignored by name.
- Never print, stage, or commit them. Do not `git add .` in a shell where
  they have been decrypted.
<!---->
## Nix
<!---->
- Always use flakes and flake-based commands (`nix run`, `nix shell`, etc.).
  Never use the old `nix-shell` approach.
- If you are confused, stop and ask for help.
This is especially critical in
  Nix.
- Follow the existing style of the surrounding modules.
- `utils` (Analysis-Utilities) is a pinned flake input. Bump it with
  `nix flake lock --update-input utils`, never by hand-editing the lock
  file.
<!---->
## C++ / ROOT
<!---->
- Use ROOT data types, and pick the *correct* one for the actual need rather
  than defaulting blindly: `Int_t` for ordinary ints, `Long64_t` for entry
  counts and large/64-bit values, `Double_t` for floating point, `TString`
  for string convenience, and so on.
  Match the width and signedness the code
  actually requires.
- Do not use modern C++ features: no `auto`, no smart pointers, no
  range-based (`for (x : c)`) iteration.
  Use explicit types and classic
  indexed/iterator loops.
- Lambdas are permitted where they are short, local, and improve readability
  over the alternatives: sort comparators defined next to the std::sort call,
  thread workers whose explicit parameter lists would be longer and harder to
  scan than an inline capture.
  A lambda that spans more than about 5 lines or
  captures by reference outside an immediately obvious scope (e.g. stored in a
  std::function returned from the function) should still be a named function.
  When in doubt, write a named function.
- In performance-critical code, always gate logging behind a compile- or
  run-time toggle so it can be disabled.
  The `std::endl` flush is therefore
  never a concern on hot paths.
<!---->
## Python
<!---->
- In Python that uses ROOT, never use matplotlib.
Look at nearby files for the
  established plotting approach, or ask which is preferred.
<!---->
## Explanations
<!---->
- For non-trivial changes, explain thoroughly what changed and why.
Do not
  over-summarize or truncate the reasoning.
  Trivial edits can stay terse.
<!---->
## Commits
<!---->
- A Conventional Commits summary line with scope (`feat(73Ge): …`,
  `fix(78mBr): …`, `docs: …`), a blank line, then a body that explains
  what changed and why.
- Agent commits: `--author="son-of-anton-bot <307402699+son-of-anton-bot@users.noreply.github.com>"`.
- End an agent commit's body with a `Model:` trailer naming the model 
  that wrote it: son-of-anton instances use `$SOA_MODEL` (e.g. `Model:   radiance/qwen3.8-27b`); other agents give their own model ID
  (e.g. `Model: claude-opus-5-5`).
- A `pre-commit` hook in `.githooks/` clang-formats staged C/C++ with
  `LLVM, BreakStringLiterals: false` and re-stages the result.
  `BreakStringLiterals` is deliberately off, so do not reformat string
  literals by hand.
