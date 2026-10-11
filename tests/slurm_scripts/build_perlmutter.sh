#!/bin/bash
# Run manually in a writable shared checkout, with the intended Python active.
set -euo pipefail
if ! type module >/dev/null 2>&1; then
    echo 'This build helper needs the Perlmutter module environment.' >&2; exit 2
fi
module load PrgEnv-gnu
module load cray-mpich
export CRAYPE_LINK_TYPE=dynamic
export CTREEBALLS_ROOT="${CTREEBALLS_ROOT:-$PWD}"
exec "${CBALLS_PYTHON:-python3}" "$CTREEBALLS_ROOT/tests/python/build_perlmutter.py" "$@"
