# Ly-alpha forest all-engines driver

`lya_corr_all_engines.py` reads a forest catalog once and runs the active
cTreeBalls forest and physical-3D multipole methods. `LYAFORESTOMPON=1`
provides nine OpenMP names; `LYAFORESTMPION=1` provides eight MPI counterparts.
`OCTREE3PCF3DOMPON=1` and `OCTREE3PCF3DMPION=1` add two physical-3D methods.

List the methods compiled into the current extension:

```sh
python3 tests/python/lya_corr_all_engines.py --list-engines
```

## Input

NPZ catalogs contain `positions`, `delta`, `weights`, and integer
`forest_ids`. Six-column ASCII input is `x y z delta weight forest_id`.
DESI/PICCA image-layout FITS files can be supplied with `--fits`; the driver
converts absorption redshift to comoving Mpc/h using the explicitly reported
fiducial cosmology.

## Examples

```sh
python3 tests/python/lya_corr_all_engines.py \
  --catalog pixels.npz --engine all-omp --statistics 2pcf \
  --threads 8 --output Output_lya_2pcf
```

```sh
python3 tests/python/lya_corr_all_engines.py \
  --catalog pixels.npz \
  --engine lya-1d-tree-3pcf-omp octree-3pcf-3d-omp \
  --statistics 3pcf --threads 8 --output Output_lya_3pcf
```

The radial and physical-3D estimators have different geometry and are timed
and reported together, not asserted to be numerically identical.

```sh
python3 tests/python/lya_corr_all_engines.py \
  --catalog pixels.npz \
  --engine lya-2pcf-mpi lya-1d-tree-2pcf-mpi \
  --statistics 2pcf --mpi-ranks 2 --threads 4 \
  --output Output_lya_mpi
```

Use repeated `--mpi-extra-arg` arguments. Each launcher option and its value
should be a separate argument.

## Scientific contracts

The anisotropic methods bin parallel and transverse separation and exclude
same-forest pairs. Radial methods use signed or absolute line-of-sight lags as
documented by each engine. `lya-1d-tree-same-los-2pcf-omp` is deliberately
different: it accepts pairs within one forest, normalizes each occupied
forest/bin, and then averages forests equally.

Three-point forest methods require three distinct forest IDs. The physical-3D
multipole methods use `exclude-all-same-los`; they are not interchangeable with
the anisotropic five-dimensional estimator.

`summary.json` records catalog provenance, timings, products, comparisons, and
plot paths. Empty denominator bins publish finite zero in native output and
appear as missing data in plots.
