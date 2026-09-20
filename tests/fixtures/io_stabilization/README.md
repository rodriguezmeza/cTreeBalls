# Gadget/FITS input-output regression fixtures

These are the exact small, synthetic fixtures recorded during the active-addon
review (September 2026); they contain no external scientific observations.

* `single.snap`: native little-endian Gadget format 1, eight float32 positions,
  BoxSize 10. Tests box initialization and scalar-field assignment.
* `multi.snap.0`, `multi.snap.1`: the same catalog split into two four-particle
  files, with the first z coordinate changed from 3 to 11. Its wrapped position
  must be (1, 2, 1), catching the former z-from-y error.
* `radec.fits`: eight binary-table rows, double columns KAPPA, RA, DEC, WEIGHT.
  KAPPA and WEIGHT are 1; RA is linspace(0.05, 0.8, 8), DEC linspace(0.3, 0.8, 8).
  Angles are radians. Selecting weight column 99 must fail; the same fixture
  also supplies valid input for failed and successful export tests.

Run `make test-io-stabilization` for CLI checks and
`make test-io-stabilization-cython` for CLI plus Python-interface checks.
The CLI tests use only the Python standard library. Cython checks additionally
need the locally built extension and its NumPy dependency. The test target
requires a 3D double build with GADGETIOON=1, CFITSIOON=1,
OCTREE2BALLSOMPON=1 and OCTREE3PCF3DOMPON=1 (for optional XYZ weight/LOS columns).

The test script creates malformed/boundary variants in temporary directories;
it never changes these recorded files. FITS output is independently decoded as
big-endian doubles and checked numerically without requiring Astropy.

SHA-256:

* `single.snap`: `bcde609fff0fa78869a8b5dcbe25da63d1e8c94571422dac4655d341e2df7d21`
* `multi.snap.0`: `6a07a83ae558949d1eb69009421b5b88b04ca2c1e843ddc5fd54a9ddf3f1e606`
* `multi.snap.1`: `c5a70a4e1d18dd623c4603a790dccde8e1c7d5d2dffbaf6b2dcda3f124fe6349`
* `radec.fits`: `0cf3f7d880199b72002b8bb683648fed29d0851182d8c13ac3d647703e067324`
