#!/bin/sh

set -eu

output=$1
temporary="${output}.tmp"

{
    printf '%s\n' '#ifndef CBALLS_MAKE_INFO_H'
    printf '%s\n' '#define CBALLS_MAKE_INFO_H'
    printf '\n'

    for name in ${CBALLS_MAKE_INFO_VARIABLES}; do
        value=$(printenv "${name}" 2>/dev/null || :)
        escaped=$(printf '%s' "${value}" | sed 's/\\/\\\\/g; s/"/\\"/g')
        printf '#define CBALLS_MAKE_%s "%s"\n' "${name}" "${escaped}"
    done

    printf '\n%s\n' '#endif'
} > "${temporary}"

if test -f "${output}" && cmp -s "${temporary}" "${output}"; then
    rm -f "${temporary}"
else
    mv "${temporary}" "${output}"
fi
