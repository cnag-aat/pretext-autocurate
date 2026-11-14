# Notes for Claude Code in this repository

- Always run shell commands under **bash**, not zsh (the macOS default): bioinformatics
  scripts and one-liners assume bash. Wrap commands as `bash -c '...'` or use `#!/bin/bash`
  scripts. (zsh does not word-split unquoted variables, so `$ARGS` holding several options
  arrives as one argument.)
- Start from `docs/HANDOFF.md` (state of the code, conventions, next steps), then the tail
  of `docs/PROGRESS.md`.
- Tests: `python3 -m unittest discover -s tests` (needs the data in `seed/fixtures/`, not in
  the repository; about 12 minutes).
- External tools live in conda envs `pretext-autocurate` (minimap2, miniprot, FastGA) and
  `pretext-autocurate-x86` (Red, PretextGraph; Intel builds, run under Rosetta).
