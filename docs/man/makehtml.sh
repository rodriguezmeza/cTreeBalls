#!/bin/sh
set -eu
cd "$(dirname "$0")"
command -v mandoc >/dev/null 2>&1 || {
    echo "Install mandoc to regenerate the manual HTML." >&2
    exit 1
}
cp cballs.m cballs.1
mandoc -Thtml cballs.1 > cballs.html.tmp
mv cballs.html.tmp cballs.html
