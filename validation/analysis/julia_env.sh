# shellcheck shell=bash
# julia_env.sh — resolve $JULIA to Julia >= 1.12 (the repository's analysis/
# environment). Source it. A bare `julia` on a cluster's PATH can be much older,
# so resolve on the version: $JULIA, juliaup's launcher, then `julia` on PATH.

_julia_ok() {                       # $1 = candidate; true if >= 1.12
  local v major minor
  v=$("$1" --version 2>/dev/null) || return 1
  v=${v##* }; major=${v%%.*}; minor=${v#*.}; minor=${minor%%.*}
  [ "$major" -gt 1 ] || { [ "$major" -eq 1 ] && [ "$minor" -ge 12 ]; }
}

resolve_julia() {
  local c
  [ -n "${JULIA:-}" ] && _julia_ok "$JULIA" && return 0
  for c in "$HOME/.juliaup/bin/julia" "$(command -v julia 2>/dev/null || true)"; do
    [ -n "$c" ] && [ -x "$c" ] || continue
    if _julia_ok "$c"; then JULIA="$c"; export JULIA; return 0; fi
  done
  echo "error: no Julia >= 1.12 found. Set JULIA=/path/to/julia." >&2
  return 1
}
