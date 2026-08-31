#!/usr/bin/env bash
# One small DR1 file (12 forests; about 154 KiB), not a recursive survey download.
set -euo pipefail
destination="${1:-data/desi-lya-example}"
mkdir -p "$destination"
wget --continue --directory-prefix="$destination" \
  'https://data.desi.lbl.gov/public/dr1/vac/dr1/lya-deltas/v1.0/delta-lya-0-0/Delta/delta-1019.fits.gz'
