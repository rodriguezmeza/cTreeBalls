# Active source profile

This testing branch contains the enabled Makefile addon set. The core octree method is always available. Shared implementation dependencies live in support/ and do not register standalone methods. The Cython binding sources stay in python/; Python analysis, benchmark and regression scripts live in tests/python/. Inactive addons and addons/python_env are excluded from the branch.

The development checkout includes the enabled bundled GSL and CFITSIO sources. A public source distribution excludes those libraries and sets `GSLINTERNAL=0, CFITSIOLIBON=0` in its staged Makefiles.

## Shipped addon directories

- `addons/addons_include`
- `addons/balltree_2balls_mpi`
- `addons/balltree_2balls_omp`
- `addons/balltree_shear_sphere_2balls_mpi`
- `addons/balltree_shear_sphere_2balls_omp`
- `addons/cfitsio`
- `addons/cfitsiolib`
- `addons/class_lib`
- `addons/cmdline_defs_settings`
- `addons/gadget_io`
- `addons/gsl`
- `addons/iolib`
- `addons/kdtree_2balls_mpi`
- `addons/kdtree_2balls_omp`
- `addons/kdtree_box_omp`
- `addons/kdtree_shear_sphere_2balls_mpi`
- `addons/kdtree_shear_sphere_2balls_omp`
- `addons/lya_forest_mpi`
- `addons/lya_forest_omp`
- `addons/neighbor_boxes_omp`
- `addons/octree_2balls_mpi`
- `addons/octree_2balls_omp`
- `addons/octree_3pcf_3d_mpi`
- `addons/octree_3pcf_3d_omp`
- `addons/octree_shear_sphere_2balls_mpi`
- `addons/octree_shear_sphere_2balls_omp`
- `addons/pxd`

## Shared implementation dependencies

- `support/kd_tree`
- `support/native_octree_pair`
- `support/pca_tree`
- `support/scalar_octree`
- `support/scalar_octree_mpi`
- `support/shear_binary`
- `support/shear_kernel`
- `support/spherical_shear`

## Enabled runtime methods

| Name | ID | Geometry and estimator |
| --- | ---: | --- |
| `octree-sincos-omp` | 24 | 2D/3D scalar; octree; OpenMP; standard 2PCF and sine/cosine 3PCF multipoles |
| `kdtree-2balls-omp` | 197 | scalar; median k-d tree; dual-node OpenMP scan; 2PCF and LogMultipole angular 3PCF, masks, and complex edge correction |
| `kdtree-2balls-mpi` | 198 | scalar; median k-d tree; deterministic MPI+OpenMP dual-node scan; distributed 2PCF and LogMultipole angular 3PCF, masks, and complex edge correction |
| `balltree-2balls-omp` | 174 | scalar; FCFC PCA ball tree; dual-node dual/triple-node OpenMP scan; 2PCF and angular-multipole 3PCF with auto- and cross-catalog support |
| `balltree-2balls-mpi` | 179 | scalar; FCFC PCA ball tree; deterministic MPI+OpenMP dual/triple-node scan; distributed 2PCF and angular-multipole 3PCF with auto- and cross-catalog support |
| `octree-2balls-omp` | 176 | scalar; native octree; dual-node 2PCF and LogMultipole 3PCF; 2PCF and angular-multipole 3PCF with auto- and cross-catalog support |
| `octree-2balls-mpi` | 177 | scalar; native octree; deterministic MPI+OpenMP frontier; distributed 2PCF and LogMultipole angular-multipole 3PCF |
| `octree-shear-sphere-2balls-omp` | 200 | full-sky spin-2 shear; dual-node dual-node octree; OpenMP; xi+/xi- and four Gamma-x shear 3PCF components |
| `octree-shear-sphere-2balls-mpi` | 203 | full-sky spin-2 shear; native octree; adaptive MPI+OpenMP frontier; distributed xi+/xi- and four Gamma-x shear 3PCF components |
| `kdtree-shear-sphere-2balls-omp` | 201 | full-sky spin-2 shear; dual-node median KD tree; OpenMP; xi+/xi- and four Gamma-x shear 3PCF components |
| `kdtree-shear-sphere-2balls-mpi` | 204 | full-sky spin-2 shear; median KD tree; adaptive MPI+OpenMP frontier; distributed xi+/xi- and four Gamma-x shear 3PCF components |
| `balltree-shear-sphere-2balls-omp` | 202 | full-sky spin-2 shear; dual-node FCFC PCA ball tree; OpenMP; xi+/xi- and four Gamma-x shear 3PCF components |
| `balltree-shear-sphere-2balls-mpi` | 205 | full-sky spin-2 shear; FCFC PCA ball tree; adaptive MPI+OpenMP frontier; distributed xi+/xi- and four Gamma-x shear 3PCF components |
| `kdtree-box-omp` | 71 | periodic Cartesian box; k-d tree; OpenMP; box 2PCF |
| `neighbor-boxes-omp` | 73 | periodic Cartesian box; linked boxes; OpenMP; periodic pair counts and unweighted density correlation function |
| `octree-3pcf-3d-omp` | 166 | 3D scalar; exact octree leaves; OpenMP; spherical-harmonic 2PCF/3PCF multipoles |
| `octree-3pcf-3d-mpi` | 192 | 3D scalar; exact octree leaves; MPI+OpenMP pivot blocks; spherical-harmonic 2PCF/3PCF; data/random survey estimator and edge correction |
| `lya-2pcf-omp` | 169 | 3D observer-centered forest pixels; exact default, opt-in forest cell pairs; OpenMP; weighted anisotropic 2PCF |
| `lya-3pcf-omp` | 170 | 3D observer-centered forest pixels; exact default, opt-in bounded cell-bin approximation; OpenMP; weighted five-dimensional 3PCF |
| `lya-2pcf-3pcf-omp` | 171 | 3D observer-centered forest pixels; exact default, opt-in forest cell pairs; OpenMP; weighted 2PCF and 3PCF; shared discovery or independently optimized passes |
| `lya-los-tree-2pcf-omp` | 206 | 3D observer-centered forest pixels; exact default, opt-in forest cell pairs; OpenMP; exact weighted anisotropic 2PCF |
| `lya-los-tree-3pcf-omp` | 207 | 3D observer-centered forest pixels; exact default, opt-in bounded cell-bin approximation; OpenMP; exact weighted five-dimensional 3PCF |
| `lya-los-tree-2pcf-3pcf-omp` | 208 | 3D observer-centered forest pixels; exact default, opt-in forest cell pairs; OpenMP; weighted 2PCF and 3PCF; shared discovery or independently optimized passes |
| `lya-1d-2pcf-omp` | 180 | radial Lyman-alpha pixels; sorted 1D range scan; OpenMP; weighted radial-only 2PCF |
| `lya-1d-3pcf-omp` | 181 | radial Lyman-alpha pixels; sorted 1D range scan; OpenMP; weighted radial-only 3PCF |
| `lya-1d-2pcf-3pcf-omp` | 182 | radial Lyman-alpha pixels; sorted 1D range scan; OpenMP; weighted radial-only 2PCF and 3PCF |
| `lya-1d-tree-2pcf-omp` | 183 | radial Lyman-alpha pixels; exact 1D interval tree; OpenMP; weighted radial-only 2PCF |
| `lya-1d-tree-3pcf-omp` | 193 | radial Lyman-alpha pixels; exact 1D interval tree; OpenMP; weighted radial-only 3PCF |
| `lya-1d-tree-same-los-2pcf-omp` | 195 | per-LOS radial Lyman-alpha pixels; exact 1D interval trees; OpenMP; equal-LOS mean of weighted within-forest radial 2PCFs |
| `lya-2pcf-mpi` | 185 | Lyman-alpha pixels; 3D octree; MPI+OpenMP; weighted anisotropic 2PCF |
| `lya-3pcf-mpi` | 186 | Lyman-alpha pixels; 3D octree; MPI+OpenMP; weighted five-dimensional 3PCF |
| `lya-2pcf-3pcf-mpi` | 187 | Lyman-alpha pixels; 3D octree; MPI+OpenMP; weighted 2PCF and 3PCF |
| `lya-1d-2pcf-mpi` | 188 | Lyman-alpha pixels; radial range scan; MPI+OpenMP; weighted radial-only 2PCF |
| `lya-1d-3pcf-mpi` | 189 | Lyman-alpha pixels; radial range scan; MPI+OpenMP; weighted radial-only 3PCF |
| `lya-1d-2pcf-3pcf-mpi` | 190 | Lyman-alpha pixels; radial range scan; MPI+OpenMP; weighted radial-only 2PCF and 3PCF |
| `lya-1d-tree-2pcf-mpi` | 191 | Lyman-alpha pixels; radial interval tree; MPI+OpenMP; weighted radial-only 2PCF with same-forest subtraction |
| `lya-1d-tree-3pcf-mpi` | 194 | Lyman-alpha pixels; radial interval tree; MPI+OpenMP; weighted radial-only 3PCF with exact three-forest exclusion |
| `lya-anisotropic-multipole-3pcf-omp` | 209 | 3D observer-centered forest pixels; exact radial and LOS polar bins; LOS-tree discovery; OpenMP; anisotropic Legendre raw moments; approximate mu reconstruction; optional exact five-dimensional bins |
| `lya-los-tree-2pcf-mpi` | 210 | 3D observer-centered forest pixels; exact default, opt-in forest cell pairs; replicated catalog MPI + OpenMP; exact weighted anisotropic 2PCF |
| `lya-los-tree-3pcf-mpi` | 211 | 3D observer-centered forest pixels; exact default, opt-in bounded cell-bin approximation; replicated catalog MPI + OpenMP; exact weighted five-dimensional 3PCF |
| `lya-los-tree-2pcf-3pcf-mpi` | 212 | 3D observer-centered forest pixels; exact default, opt-in forest cell pairs; replicated catalog MPI + OpenMP; weighted 2PCF and 3PCF; shared discovery or independently optimized passes |

Use `options=make-info`, `options=print-options` and `options=print-search-methods` to inspect a particular compiled executable. The latter is authoritative after any supported build overrides.
