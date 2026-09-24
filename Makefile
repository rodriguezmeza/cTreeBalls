# ----- MAKE FILE FOR cballs CODE -----
# Mario A. Rodriguez-Meza, Ciudad de Mexico, 27.04.2023
#
#
# Nothing to do in this file. Make your settings in Makefile_settings file only
#
MACHINES_DIR = ./
MPIEXEC ?= mpiexec
.DEFAULT_GOAL := all
# Machine definitions and code settings. Edit this file according to your needs
include $(MACHINES_DIR)Makefile_settings

#
# Nothing to do below
#

EXECPREFIX = cballs
ifndef SEARCHMETHODISDEFINED
EXEC = $(EXECPREFIX)
else
ifdef DEFAULTSETON
EXEC = $(EXECPREFIX)
else
EXEC = $(EXECPREFIX)_$(SEARCHMETHODIS)
endif


endif

$(info )
$(info =====================================================)
$(info SEARCHMETHOD = [${SEARCHMETHOD}]  EXEC = [${EXEC}])
$(info External options = [${OPT2}])
$(info =====================================================)
$(info )

MAIN = main.o
PROFILE_EXEC = $(WRKDIR)/$(EXEC)
PROFILE_LIB = $(WRKDIR)/lib$(EXEC).a

OBJS = main.o cballsio.o cballs.o startrun.o testdata.o treeload.o \
	cballsutils.o search.o abi_check.o run_metadata.o engine_registry.o \
	runtime_context.o memory_catalog.o common_histogram.o smooth_pivots.o mpi_runtime.o

PYTHON_FILES = python/cyballs.pyx setup.py python/ccyballs.pxd.in
CBALLS_OBJECTS = $(sort $(OBJS) $(TOOLS) $(SOURCE) $(EXTERNAL) $(EXTERNALCXX))
CBALLS_COMPILE_CONFIG = $(WRKDIR)/compile-settings.json
CBALLS_LINK_CONFIG = $(WRKDIR)/link-settings.json
export CBALLS_OBJECTS VENDORED_SOURCE_PATTERNS

all: $(EXEC) lib$(EXEC).a cyballs

# Reader headers may be dropped by later legacy HEADERFILES assignments.
cballsio.o: $(wildcard $(MDIR)/addons/iolib/*.h) \
            $(wildcard $(MDIR)/addons/cfitsio/*.h)

search_balltree_2balls_omp.o search_balltree_2balls_mpi.o \
search_balltree_2balls_omp_3pcf.o search_balltree_2balls_mpi_3pcf.o \
search_octree_2balls_omp.o search_octree_2balls_mpi.o: \
	$(MDIR)/addons/balltree_2balls_omp/dual_node_edge_correction.h

# One archive writer for both native and Python targets, even with make -j.
# Starting from an empty archive removes members belonging to disabled addons.
$(PROFILE_LIB): $(CBALLS_OBJECTS) $(CBALLS_LINK_CONFIG)
	$(RM) $@.tmp
	$(AR) $@.tmp $(addprefix $(WRKDIR)/,$(CBALLS_OBJECTS))
	mv $@.tmp $@

.PHONY: lib$(EXEC).a
lib$(EXEC).a: $(PROFILE_LIB)
	@cmp -s $< $@ || cp $< $@

.PHONY: cyballs-static-lib
cyballs-static-lib: lib$(EXEC).a

.PHONY: $(EXEC)
$(EXEC): $(PROFILE_EXEC)
	@cmp -s $< $@ || cp $< $@

$(PROFILE_EXEC): $(OBJS) $(EXTERNAL) $(MAIN) $(CBALLS_LINK_CONFIG)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $@.tmp $(addprefix $(WRKDIR)/,$(sort $(OBJS) $(EXTERNAL) $(MAIN))) $(MLIBS) $(FITSIOLIBS)
	mv $@.tmp $@

#B to fix NDIM problem
.PHONY: print-cyballs-build-env
print-cyballs-build-env:
	@printf '__CBALLS_CC__=%s\n' '$(CC)'
	@printf '__CBALLS_MPICC__=%s\n' '$(MPICC)'
	@printf '__CBALLS_LIB__=%s\n' '$(EXEC)'
	@printf '__CBALLS_CPPFLAGS__=%s\n' '$(OPT2) $(INCLUDES)'
	@printf '__CBALLS_ADDONSON__=%s\n' '$(ADDONSON)'
	@printf '__CBALLS_CLASSLIBON__=%s\n' '$(CLASSLIBON)'
	@printf '__CBALLS_PXDON__=%s\n' '$(PXDON)'
	@printf '__CBALLS_CFITSIOON__=%s\n' '$(CFITSIOON)'
	@printf '__CBALLS_CFITSIOLIBON__=%s\n' '$(CFITSIOLIBON)'
	@printf '__CBALLS_USEGSL__=%s\n' '$(USEGSL)'
	@printf '__CBALLS_OPENMPMACHINE__=%s\n' '$(OPENMPMACHINE)'
	@printf '__CBALLS_SLEEFON__=%s\n' '$(SLEEF_ENABLED)'
	@printf '__CBALLS_SINGLEPON__=%s\n' '$(SINGLEPON)'
	@printf '__CBALLS_KDTREEMPION__=%s\n' '$(KDTREEMPION)'
	@printf '__CBALLS_KDTREE2BALLSOMPON__=%s\n' '$(KDTREE2BALLSOMPON)'
	@printf '__CBALLS_KDTREE2BALLSMPION__=%s\n' '$(KDTREE2BALLSMPION)'
	@printf '__CBALLS_KDTREESHEARSPHERE2BALLSOMPON__=%s\n' '$(KDTREESHEARSPHERE2BALLSOMPON)'
	@printf '__CBALLS_KDTREESHEARSPHERE2BALLSMPION__=%s\n' '$(KDTREESHEARSPHERE2BALLSMPION)'
	@printf '__CBALLS_BALLTREESHEARSPHERE2BALLSOMPON__=%s\n' '$(BALLTREESHEARSPHERE2BALLSOMPON)'
	@printf '__CBALLS_BALLTREESHEARSPHERE2BALLSMPION__=%s\n' '$(BALLTREESHEARSPHERE2BALLSMPION)'
	@printf '__CBALLS_OCTREESHEARSPHERE2BALLSMPION__=%s\n' '$(OCTREESHEARSPHERE2BALLSMPION)'
	@printf '__CBALLS_BALLTREE2BALLSOMP3PCFON__=%s\n' '$(BALLTREE2BALLSOMP3PCFON)'
	@printf '__CBALLS_BALLTREE2BALLSMPI3PCFON__=%s\n' '$(BALLTREE2BALLSMPI3PCFON)'
	@printf '__CBALLS_BALLTREE2BALLSMPION__=%s\n' '$(BALLTREE2BALLSMPION)'
	@printf '__CBALLS_OCTREE2BALLSOMPON__=%s\n' '$(OCTREE2BALLSOMPON)'
	@printf '__CBALLS_OCTREE2BALLSMPION__=%s\n' '$(OCTREE2BALLSMPION)'
	@printf '__CBALLS_BALLTREEMPION__=%s\n' '$(BALLTREEMPION)'
	@printf '__CBALLS_OCTREEGGGMPION__=%s\n' '$(OCTREEGGGMPION)'
	@printf '__CBALLS_OCTREEBALLS4MPION__=%s\n' '$(OCTREEBALLS4MPION)'
	@printf '__CBALLS_LYAFORESTOMPON__=%s\n' '$(LYAFORESTOMPON)'
	@printf '__CBALLS_LYAFORESTMPION__=%s\n' '$(LYAFORESTMPION)'
	@printf '__CBALLS_OCTREE3PCF3DOMPON__=%s\n' '$(OCTREE3PCF3DOMPON)'
	@printf '__CBALLS_OCTREE3PCF3DMPION__=%s\n' '$(OCTREE3PCF3DMPION)'
	@printf '__CBALLS_GSLINTERNAL__=%s\n' '$(GSLINTERNAL)'
	@printf '__CBALLS_MACOSX_DEPLOYMENT_TARGET__=%s\n' '$(MACOSX_DEPLOYMENT_TARGET)'
#E


#
#B to test cBalls under different profiles
#
.PHONY: test-default test-cell-production test-kdtree-box-frontier \
	test-balltree-2balls-omp test-balltree-2balls-mpi \
	test-balltree-2balls-3pcf test-octree-2balls-omp \
	test-octree-2balls-mask test-octree-2balls-mpi \
	test-kdtree-2balls-omp test-kdtree-2balls-mpi \
	test-octree-3pcf-3d-omp test-octree-3pcf-3d-mpi \
	test-lya-forest-omp test-lya-forest-mpi test-lya-forest-1d-omp \
	test-lya-corr-all-engines test-lya2pcf-reference \
	test-p3-cython test-sanitizer-smoke test-two-ball-edge \
	test-healpix-ordering test-parameter-parser \
	test-parameter-file-parser test-standalone-parser test-option-cache \
	test-sleef-vector-log \
	test-mixed-precision-profile test-singlep-search test-singlep \
	test-make-info test-make-info-profiles test-search-methods

test-default: all test-search-methods test-make-info test-option-cache

test-parameter-parser: lib$(EXEC).a
	mkdir -p $(WRKDIR)/tests
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(CCFLAG) \
		$(PROJECT_WARNING_FLAGS) $(INCLUDES) tests/test_parameter_parser.c \
		lib$(EXEC).a -o $(WRKDIR)/tests/test_parameter_parser \
		$(MLIBS) $(FITSIOLIBS)
	$(WRKDIR)/tests/test_parameter_parser

test-option-cache: lib$(EXEC).a
	mkdir -p $(WRKDIR)/tests
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(CCFLAG) \
		$(PROJECT_WARNING_FLAGS) $(INCLUDES) tests/test_options_cache.c \
		lib$(EXEC).a -o $(WRKDIR)/tests/test_options_cache \
		$(MLIBS) $(FITSIOLIBS)
	$(WRKDIR)/tests/test_options_cache

test-sleef-vector-log:
ifeq ($(SLEEF_ENABLED),1)
	mkdir -p $(WRKDIR)/tests
	$(CC) $(OPTFLAG) $(OMPFLAG) $(CCFLAG) $(SLEEF_CFLAGS) \
		-I$(MDIR)/addons/balltree_2balls_omp \
		tests/test_sleef_vector_log.c -o $(WRKDIR)/tests/test_sleef_vector_log \
		$(SLEEF_LDFLAGS) $(SLEEF_RPATH) $(SLEEF_LIBS) -lm
	$(WRKDIR)/tests/test_sleef_vector_log
else
	@echo "SKIP: SLEEF vector log backend is not enabled"
endif

test-make-info: $(EXEC)
	CBALLS=$(CURDIR)/$(EXEC) EXPECTED_SINGLEPON=$(SINGLEPON) \
		bash ./tests/make_tests/run_test_make_info

test-make-info-profiles:
	$(MAKE) SINGLEPON=0 test-make-info
	$(MAKE) SINGLEPON=1 test-make-info

test-search-methods: $(EXEC)
	CBALLS=$(CURDIR)/$(EXEC) bash ./tests/make_tests/run_test_search_methods

test-parameter-file-parser: $(EXEC)
	CBALLS=$(CURDIR)/$(EXEC) \
		bash ./tests/make_tests/run_test_standalone_parameter_parser

test-standalone-parser:
	$(MAKE) CLASSLIBON=0 EXEC=cballs_standalone_parser \
		cballs_standalone_parser
	CBALLS=$(CURDIR)/cballs_standalone_parser \
		bash ./tests/make_tests/run_test_standalone_parameter_parser

ifeq ($(CFITSIOON),1)
test-healpix-ordering: lib$(EXEC).a
	mkdir -p $(WRKDIR)/tests
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(CCFLAG) \
		$(PROJECT_WARNING_FLAGS) $(INCLUDES) tests/test_healpix_ordering.c \
		lib$(EXEC).a -o $(WRKDIR)/tests/test_healpix_ordering \
		$(MLIBS) $(FITSIOLIBS)
	$(WRKDIR)/tests/test_healpix_ordering
else
test-healpix-ordering:
	@echo "SKIP: CFITSIO support is disabled"
endif

test-cell-production: .base
	mkdir -p $(WRKDIR)/tests
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(CCFLAG) $(PROJECT_WARNING_FLAGS) $(INCLUDES) \
		tests/test_search_cell_aggregate.c -o $(WRKDIR)/tests/test_search_cell_aggregate $(MLIBS)
	$(WRKDIR)/tests/test_search_cell_aggregate
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(CCFLAG) $(PROJECT_WARNING_FLAGS) $(INCLUDES) \
		tests/test_cell_kappa_average.c -o $(WRKDIR)/tests/test_cell_kappa_average $(MLIBS)
	$(WRKDIR)/tests/test_cell_kappa_average
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(CCFLAG) $(PROJECT_WARNING_FLAGS) $(INCLUDES) \
		-DNOWKAvg tests/test_cell_kappa_average.c \
		-o $(WRKDIR)/tests/test_cell_kappa_average_nowkavg $(MLIBS)
	$(WRKDIR)/tests/test_cell_kappa_average_nowkavg
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(CCFLAG) $(PROJECT_WARNING_FLAGS) $(INCLUDES) \
		tests/test_mask_cell_state.c -o $(WRKDIR)/tests/test_mask_cell_state $(MLIBS)
	$(WRKDIR)/tests/test_mask_cell_state
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(CCFLAG) $(PROJECT_WARNING_FLAGS) $(INCLUDES) \
		tests/test_balls4_scan_contract.c -o $(WRKDIR)/tests/test_balls4_scan_contract $(MLIBS)
	$(WRKDIR)/tests/test_balls4_scan_contract

test-mixed-precision-profile: cyballs-static-lib
	mkdir -p $(WRKDIR)/tests
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(CCFLAG) $(PROJECT_WARNING_FLAGS) $(INCLUDES) \
		tests/test_mixed_precision.c lib$(EXEC).a \
		-o $(WRKDIR)/tests/test_mixed_precision $(MLIBS) $(FITSIOLIBS)
	$(WRKDIR)/tests/test_mixed_precision

test-singlep:
	$(MAKE) SINGLEPON=0 test-mixed-precision-profile
	$(MAKE) SINGLEPON=1 test-mixed-precision-profile
	$(MAKE) test-singlep-search

test-singlep-search:
	@tmp=$$(mktemp -d "$${TMPDIR:-/tmp}/ctreeballs-singlep-build.XXXXXX"); \
	status=0; \
	$(MAKE) -B SINGLEPON=0 cballs && cp cballs "$$tmp/cballs-double" && \
	$(MAKE) -B SINGLEPON=1 cballs && cp cballs "$$tmp/cballs-mixed" && \
	DOUBLE_CBALLS="$$tmp/cballs-double" \
	MIXED_CBALLS="$$tmp/cballs-mixed" \
		bash tests/make_tests/run_test_singlep_search || status=$$?; \
	$(MAKE) -B SINGLEPON=0 cballs cyballs-static-lib || status=$$?; \
	rm -rf "$$tmp"; \
	exit $$status

.PHONY: test-scalar-numerical-contract
test-scalar-numerical-contract: $(EXEC) cyballs
	$(PYTHON) tests/make_tests/test_scalar_numerical_contract.py

test-kdtree-box-frontier: $(EXEC)
	cd tests && CBALLS=$(CURDIR)/$(EXEC) \
		bash ./make_tests/run_test_kdtree_box_frontier

test-balltree-2balls-omp: $(EXEC)
	cd tests && CBALLS=$(CURDIR)/$(EXEC) \
		bash ./make_tests/run_test_balltree_2balls_omp
	cd tests && CBALLS=$(CURDIR)/$(EXEC) \
		bash ./make_tests/run_test_balltree_2balls_3pcf

test-balltree-2balls-mpi:
	$(MAKE) -B BALLTREE2BALLSMPION=1 TWOPCFON=1 TPCFON=1 cballs
	cd tests && CBALLS=$(CURDIR)/$(EXEC) MPIEXEC='$(MPIEXEC)' \
		bash ./make_tests/run_test_balltree_2balls_mpi

test-balltree-2balls-3pcf: $(EXEC)
	cd tests && CBALLS=$(CURDIR)/$(EXEC) \
		bash ./make_tests/run_test_balltree_2balls_3pcf

test-octree-2balls-omp: $(EXEC)
	cd tests && CBALLS=$(CURDIR)/$(EXEC) \
		bash ./make_tests/run_test_octree_2balls_omp

test-octree-2balls-mask: $(EXEC)
	$(PYTHON) tests/make_tests/test_octree_2balls_mask.py --cballs $(CURDIR)/$(EXEC)

.PHONY: test-two-ball-edge
test-two-ball-edge: $(EXEC)
	$(PYTHON) tests/make_tests/test_two_ball_edge_corrections.py \
		--cballs $(CURDIR)/$(EXEC) --dimension $(DEFDIMENSION)

.PHONY: test-kdtree-2balls-omp test-kdtree-2balls-mpi
test-kdtree-2balls-omp: $(EXEC)
	cd tests && CBALLS=$(CURDIR)/$(EXEC) \
		bash ./make_tests/run_test_kdtree_2balls_omp

test-kdtree-2balls-mpi:
	$(MAKE) -B KDTREE2BALLSOMPON=1 KDTREE2BALLSMPION=1 \
		TWOPCFON=1 TPCFON=1 cballs
	cd tests && CBALLS=$(CURDIR)/$(EXEC) MPIEXEC='$(MPIEXEC)' \
		bash ./make_tests/run_test_kdtree_2balls_mpi

test-octree-2balls-mpi:
	$(MAKE) -B OCTREE2BALLSMPION=1 cballs
	cd tests && CBALLS=$(CURDIR)/$(EXEC) MPIEXEC='$(MPIEXEC)' \
		bash ./make_tests/run_test_octree_2balls_mpi

.PHONY: test-octree-3pcf-3d-mpi
test-octree-3pcf-3d-mpi: $(EXEC)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tests/make_tests/test_octree_3pcf_3d_mpi.py \
		--cballs $(CURDIR)/$(EXEC) \
		--mpi-command "$(MPIEXEC) -n 2" $(CB3D_MPI_TEST_ARGS)

test-octree-3pcf-3d-omp: $(EXEC)
	CBALLS=$(CURDIR)/$(EXEC) \
		bash ./tests/make_tests/run_test_octree_3pcf_3d_omp

test-lya-forest-omp: $(EXEC)
	cd tests && CBALLS=$(CURDIR)/$(EXEC) \
		$(PYTHON) ./make_tests/test_lya_forest_omp.py

.PHONY: test-lya-corr-all-engines
test-lya-corr-all-engines: cyballs-static-lib
	CBALLS_STATIC_LIBRARY_READY=1 $(PYTHON) setup.py build_ext --inplace
	PYTHONPATH=$(CURDIR):$(CURDIR)/python:$(PYTHONPATH) \
		$(PYTHON) -m pytest -q tests/make_tests/test_lya_corr_all_engines.py

.PHONY: test-lya-forest-mpi
test-lya-forest-mpi: $(EXEC)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tests/make_tests/test_lya_forest_mpi.py \
		--cballs $(CURDIR)/$(EXEC) \
		--mpi-command "$(if $(MPIEXEC),$(MPIEXEC),mpiexec) -n 2" $(LYA_MPI_TEST_ARGS)

test-lya-forest-1d-omp: $(EXEC)
	cd tests && CBALLS=$(CURDIR)/$(EXEC) \
		$(PYTHON) ./make_tests/test_lya_forest_1d_omp.py

test-lya2pcf-reference: $(EXEC)
	cd tests && CBALLS=$(CURDIR)/$(EXEC) \
		LYA2PCF_SOURCE='$(LYA2PCF_SOURCE)' PYTHONDONTWRITEBYTECODE=1 \
		$(PYTHON) ./make_tests/test_lya2pcf_reference.py

test-p3-cython: cyballs
	$(PYTHON) tests/make_tests/test_p3_cython_startup.py

.PHONY: test-shear-sphere
test-shear-sphere: cyballs
	CBALLS_SHEAR_SPHERE_ENGINE=octree-shear-sphere-2balls-omp \
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tests/make_tests/test_shear_sphere_octree_omp.py

.PHONY: test-shear-sphere-2balls
test-shear-sphere-2balls: cyballs
	CBALLS_SHEAR_SPHERE_ENGINE=octree-shear-sphere-2balls-omp \
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tests/make_tests/test_shear_sphere_octree_omp.py
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tests/make_tests/test_shear_pivot_reuse.py

.PHONY: test-shear-sphere-kdtree-2balls
test-shear-sphere-kdtree-2balls: cyballs
	CBALLS_SHEAR_SPHERE_ENGINE=kdtree-shear-sphere-2balls-omp \
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tests/make_tests/test_shear_sphere_octree_omp.py

.PHONY: test-shear-sphere-balltree-2balls
test-shear-sphere-balltree-2balls: cyballs
	CBALLS_SHEAR_SPHERE_ENGINE=balltree-shear-sphere-2balls-omp \
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tests/make_tests/test_shear_sphere_octree_omp.py
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tests/make_tests/test_balltree_shear_build.py
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tests/make_tests/test_balltree_shear_pivot_reuse.py

.PHONY: test-shear-all-engines
test-shear-all-engines: cyballs
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) -m pytest -q tests/make_tests/test_shear_corr_all_engines.py

.PHONY: test-lya-los-tree test-benchmark-drivers
test-lya-los-tree: $(EXEC) cyballs
	CBALLS=$(CURDIR)/$(EXEC) $(PYTHON) -m pytest -q tests/make_tests/test_lya_forest_los_tree.py

test-benchmark-drivers: $(EXEC) cyballs
	$(PYTHON) -m pytest -q tests/make_tests/test_benchmark_timing.py \
		tests/make_tests/test_kappa_corr_all_engines.py tests/make_tests/test_kappa_corr_patch.py \
		tests/make_tests/test_shear_corr_all_engines.py tests/make_tests/test_lya_corr_all_engines.py \
		tests/make_tests/test_lya_analysis.py -k 'not mpi'

test-sanitizer-smoke: $(EXEC) test-cell-production
	cd tests && bash ./make_tests/run_test_sanitizer_smoke
#
#E to test cBalls under different profiles
#


ifeq ($(CLASSLIBON),1)
.PHONY: cyballs reinstall-cyballs
cyballs: cyballs-static-lib $(PYTHON_FILES)
	@CC='$(CC)' CBALLS_LIB='$(EXEC)' CBALLS_CPPFLAGS='$(OPT2) $(INCLUDES)' \
		CBALLS_STATIC_LIBRARY_READY=1 \
		$(PYTHON) scripts/incremental_cyballs.py --root $(CURDIR) --profile $(WRKDIR)

reinstall-cyballs: cyballs-static-lib $(PYTHON_FILES)
	@CC='$(CC)' CBALLS_LIB='$(EXEC)' CBALLS_CPPFLAGS='$(OPT2) $(INCLUDES)' \
		CBALLS_STATIC_LIBRARY_READY=1 \
		$(PYTHON) scripts/incremental_cyballs.py --root $(CURDIR) --profile $(WRKDIR) --force
else
cyballs:
endif

.PHONY : clean
clean: .base
	rm -rf $(WRKDIR);
	rm -f $(EXEC)
	rm -f cballs_standalone_parser
	rm -f lib$(EXEC).a
	rm -f $(MDIR)/python/ccyballs.pxd
	rm -f $(MDIR)/python/cyballs.c
	rm -f $(MDIR)/python/cyballs*.so $(MDIR)/python/cyballs*.pyd $(MDIR)/python/cyballs*.dylib
	rm -f $(MDIR)/cyballs*.so $(MDIR)/cyballs*.pyd $(MDIR)/cyballs*.dylib
	rm -rf $(MDIR)/python/build
	rm -rf cyballs.egg-info
	rm -rf build;
	rm -rf dist;

# Small recorded Gadget/FITS regressions; no optional Python FITS dependency.
.PHONY: test-io-stabilization test-io-stabilization-cython
ifeq ($(GADGETIOON)$(CFITSIOON)$(OCTREE2BALLSOMPON)$(OCTREE3PCF3DOMPON),1111)
test-io-stabilization: $(EXEC)
	$(PYTHON) tests/make_tests/test_io_stabilization.py --cballs $(CURDIR)/$(EXEC)

test-io-stabilization-cython: $(EXEC) cyballs-static-lib $(PYTHON_FILES)
	@test "$(CLASSLIBON)" = "1" || { echo "ERROR: requires CLASSLIBON=1"; exit 1; }
	CBALLS_STATIC_LIBRARY_READY=1 $(PYTHON) setup.py build_ext --inplace --force
	$(PYTHON) tests/make_tests/test_io_stabilization.py --cballs $(CURDIR)/$(EXEC) --cython
else
test-io-stabilization test-io-stabilization-cython:
	@echo "ERROR: requires GADGETIOON=1 CFITSIOON=1 OCTREE2BALLSOMPON=1 OCTREE3PCF3DOMPON=1"
	@exit 1
endif

# Active-profile regressions for count normalization and resource contracts.
.PHONY: test-runtime-stabilization-native test-runtime-stabilization test-runtime-stabilization-cython
test-runtime-stabilization-native: cyballs-static-lib
	mkdir -p $(WRKDIR)/tests
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(CCFLAG) \
		$(PROJECT_WARNING_FLAGS) $(INCLUDES) tests/test_runtime_stabilization.c \
		lib$(EXEC).a -o $(WRKDIR)/tests/test_runtime_stabilization \
		$(MLIBS) $(FITSIOLIBS)
	$(WRKDIR)/tests/test_runtime_stabilization

test-runtime-stabilization: $(EXEC) test-runtime-stabilization-native
	$(PYTHON) tests/make_tests/test_runtime_stabilization.py --cballs $(CURDIR)/$(EXEC)

test-runtime-stabilization-cython: $(EXEC) cyballs-static-lib test-runtime-stabilization-native $(PYTHON_FILES)
	CBALLS_STATIC_LIBRARY_READY=1 $(PYTHON) setup.py build_ext --inplace --force
	$(PYTHON) tests/make_tests/test_runtime_stabilization.py --cballs $(CURDIR)/$(EXEC) --cython

# Provenance, scalar-window support, and observer-relative mask regressions.
.PHONY: test-provenance-window-native test-provenance-window

test-provenance-window-native: cyballs-static-lib
	mkdir -p $(WRKDIR)/tests
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(CCFLAG) \
		$(PROJECT_WARNING_FLAGS) $(INCLUDES) tests/test_scalar_window.c \
		lib$(EXEC).a -o $(WRKDIR)/tests/test_scalar_window \
		$(MLIBS) $(FITSIOLIBS)
	$(WRKDIR)/tests/test_scalar_window

test-provenance-window: $(EXEC) test-provenance-window-native $(PYTHON_FILES)
	CBALLS_STATIC_LIBRARY_READY=1 $(PYTHON) setup.py build_ext --inplace --force
	PYTHONPATH=$(CURDIR):$(PYTHONPATH) $(PYTHON) -m pytest -q \
		tests/make_tests/test_provenance_window.py \
		tests/make_tests/test_kappa_corr_all_engines.py \
		tests/make_tests/test_two_ball_edge_cython.py -k 'not mpi'
	PYTHONPATH=$(CURDIR):$(PYTHONPATH) $(PYTHON) tests/make_tests/test_octree_2balls_mask.py \
		--cballs $(CURDIR)/$(EXEC) --cython
	$(PYTHON) tests/make_tests/test_two_ball_edge_corrections.py --cballs $(CURDIR)/$(EXEC) \
		--engine kdtree-2balls-omp --engine balltree-2balls-omp --engine octree-2balls-omp

# Resolve settings before stamping both native and Cython compilation units.
CBALLS_BUILD_HEADER = $(WRKDIR)/cballs_build_fingerprint.h
CBALLS_FINGERPRINT_VARIABLES = $(CBALLS_MAKE_INFO_VARIABLES) OBJS TOOLS SOURCE EXTERNALCXX EXEC
export CBALLS_FINGERPRINT_VARIABLES OBJS TOOLS SOURCE EXTERNALCXX EXEC
.PHONY: cballs-fingerprint-force print-build-fingerprint active-release-gate
cballs-fingerprint-force:

$(CBALLS_BUILD_HEADER): cballs-fingerprint-force scripts/build_fingerprint.py $(MAKEFILE_LIST) | .base
	@$(PYTHON) scripts/build_fingerprint.py --root $(CURDIR) --header $@

# Source provenance belongs to the metadata unit, not every engine.
run_metadata.o: $(CBALLS_BUILD_HEADER)
startrun.o: $(CBALLS_MAKE_INFO_HEADER)

.PHONY: cballs-config-force cballs-dependency-force
cballs-config-force cballs-dependency-force:

$(CBALLS_COMPILE_CONFIG): cballs-config-force scripts/build_support.py | .base
	@$(PYTHON) scripts/build_support.py compile $@

$(CBALLS_LINK_CONFIG): cballs-config-force scripts/build_support.py | .base
	@$(PYTHON) scripts/build_support.py link $@

$(CBALLS_OBJECTS): $(CBALLS_COMPILE_CONFIG)

# Bootstrap dependencies for old builds, or recover an accidentally removed .d.
$(foreach obj,$(CBALLS_OBJECTS),$(if $(wildcard $(WRKDIR)/$(obj:.o=.d)),,$(eval $(obj): cballs-dependency-force)))

print-build-fingerprint:
	@$(PYTHON) scripts/build_fingerprint.py --root $(CURDIR)

active-release-gate:
	$(PYTHON) scripts/active_release_gate.py $(RELEASE_GATE_ARGS)

.PHONY: test-incremental-build
test-incremental-build:
	$(PYTHON) tests/make_tests/test_incremental_build.py

.PHONY: test-two-ball-pivot-progress
test-two-ball-pivot-progress: $(EXEC)
	$(PYTHON) tests/make_tests/test_octree_2balls_progress.py --cballs $(CURDIR)/$(EXEC) \
		--engine octree-2balls-omp --engine kdtree-2balls-omp --engine balltree-2balls-omp

# Do not treat .d files as profile Makefiles in the fingerprint prerequisites.
-include $(wildcard $(addprefix $(WRKDIR)/,$(CBALLS_OBJECTS:.o=.d)))

.PHONY: check-capabilities generate-capabilities
check-capabilities:
	@$(PYTHON) scripts/generate_capabilities.py --check
generate-capabilities:
	@$(PYTHON) scripts/generate_capabilities.py
$(CBALLS_OBJECTS): | check-capabilities

.PHONY: test-resource-contracts test-runtime-context test-mpi-runtime-build affected-regressions benchmark-contracts
test-resource-contracts: cyballs-static-lib
	mkdir -p $(WRKDIR)/tests
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(CCFLAG) $(PROJECT_WARNING_FLAGS) $(INCLUDES) tests/test_resource_contracts.c lib$(EXEC).a -o $(WRKDIR)/tests/test_resource_contracts $(MLIBS) $(FITSIOLIBS)
	$(WRKDIR)/tests/test_resource_contracts
test-runtime-context: cyballs-static-lib
	mkdir -p $(WRKDIR)/tests
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(CCFLAG) $(PROJECT_WARNING_FLAGS) $(INCLUDES) tests/test_runtime_context.c lib$(EXEC).a -o $(WRKDIR)/tests/test_runtime_context $(MLIBS) $(FITSIOLIBS)
	$(WRKDIR)/tests/test_runtime_context
test-mpi-runtime-build: cyballs-static-lib
	mkdir -p $(WRKDIR)/tests
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) $(CCFLAG) $(PROJECT_WARNING_FLAGS) $(INCLUDES) tests/test_mpi_runtime.c lib$(EXEC).a -o $(WRKDIR)/tests/test_mpi_runtime $(MLIBS) $(FITSIOLIBS)
affected-regressions:
	$(PYTHON) scripts/affected_regressions.py $(REGRESSION_ARGS) --output affected-regressions.json
benchmark-contracts:
	$(PYTHON) scripts/benchmark_contracts.py $(BENCHMARK_ARGS)
