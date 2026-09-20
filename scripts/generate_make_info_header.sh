#!/bin/sh

set -eu

# Avoid three subprocesses per exported setting.
exec "${PYTHON:-python3}" "$(dirname "$0")/build_support.py" make-info "$1"
