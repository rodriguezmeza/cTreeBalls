# ----- MAKE FILE FOR cballs CODE -----
# Mario A. Rodriguez-Meza, Ciudad de Mexico, 27.04.2023
#
#
# Nothing to do in this file. Make your settings in Makefile_settings file only
#
MACHINES_DIR = ./
MPIEXEC ?= mpiexec
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

OBJS = main.o cballsio.o cballs.o startrun.o testdata.o treeload.o \
	cballsutils.o search.o

PYTHON_FILES = python/cyballs.pyx setup.py python/ccyballs.pxd.in

all: $(EXEC) lib$(EXEC).a cyballs

# Reader headers may be dropped by later legacy HEADERFILES assignments.
cballsio.o: $(wildcard $(MDIR)/addons/iolib/*.h) \
            $(wildcard $(MDIR)/addons/cfitsio/*.h)

search_balltree_2balls_omp.o search_balltree_2balls_mpi.o \
search_balltree_2balls_omp_3pcf.o search_balltree_2balls_mpi_3pcf.o \
search_octree_2balls_omp.o search_octree_2balls_mpi.o: \
	$(MDIR)/addons/balltree_2balls_omp/treecorr_edge_correction.h

lib$(EXEC).a: $(OBJS) $(EXTERNAL)
	$(AR)  $@ $(addprefix $(WRKDIR)/, $(OBJS) $(TOOLS) $(SOURCE) $(EXTERNAL) $(EXTERNALCXX))

.PHONY: cyballs-static-lib
cyballs-static-lib: $(OBJS) $(EXTERNAL)
	$(RM) lib$(EXEC).a
	$(AR) lib$(EXEC).a $(addprefix $(WRKDIR)/, $(OBJS) $(TOOLS) $(SOURCE) $(EXTERNAL) $(EXTERNALCXX))

.PHONY: $(EXEC)
$(EXEC): $(PROFILE_EXEC)
	cp $(PROFILE_EXEC) $(EXEC)

$(PROFILE_EXEC): $(OBJS) $(EXTERNAL) $(MAIN)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $@ $(addprefix $(WRKDIR)/,$(notdir $^)) $(MLIBS) $(FITSIOLIBS)

#B to fix NDIM problem
.PHONY: print-cyballs-build-env
print-cyballs-build-env:
	@printf '__CBALLS_CC__=%s\n' '$(CC)'
	@printf '__CBALLS_LIB__=%s\n' '$(EXEC)'
	@printf '__CBALLS_CPPFLAGS__=%s\n' '$(OPT2) $(INCLUDES)'
	@printf '__CBALLS_ADDONSON__=%s\n' '$(ADDONSON)'
	@printf '__CBALLS_CLASSLIBON__=%s\n' '$(CLASSLIBON)'
	@printf '__CBALLS_PXDON__=%s\n' '$(PXDON)'
	@printf '__CBALLS_CFITSIOON__=%s\n' '$(CFITSIOON)'
	@printf '__CBALLS_CFITSIOLIBON__=%s\n' '$(CFITSIOLIBON)'
	@printf '__CBALLS_USEGSL__=%s\n' '$(USEGSL)'
	@printf '__CBALLS_OPENMPMACHINE__=%s\n' '$(OPENMPMACHINE)'
	@printf '__CBALLS_SINGLEPON__=%s\n' '$(SINGLEPON)'
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
.PHONY: test-default test-balls test-p0-regressions test-p1-regressions \
	test-p2-regressions test-p2-cython test-cell-production \
	test-openmp-determinism test-octree-ggg-fast-path test-kdtree-no-one-ball \
	test-balltree-omp test-balltree-2balls-omp test-balltree-2balls-mpi \
	test-balltree-2balls-3pcf \
	test-balltree-2balls-omp-3pcf test-balltree-2balls-mpi-3pcf \
	test-octree-2balls-omp test-octree-2balls-mask \
	test-balltree-mpi test-octree-2balls-mpi test-octree-ggg-mpi \
	test-octree-3pcf-3d-omp test-octree-balls4-no-smoothing \
		test-lya-forest-omp test-lya-forest-1d-omp \
	test-lya2pcf-reference \
	test-p3-cython test-sanitizer-smoke \
	test-balls0357-recovery \
	test-healpix-ordering test-parameter-parser \
	test-parameter-file-parser test-standalone-parser test-option-cache \
	test-mixed-precision-profile test-singlep-search test-singlep \
	test-make-info test-make-info-profiles test-search-methods

test-default: all
	cd tests && ./scripts/run_all_tests

test-balls:
	$(MAKE) -B BALLSON=1 OCTREESMOOTHINGON=0 all
	cd tests && MAKE_PROFILE_ARGS="BALLSON=1 OCTREESMOOTHINGON=0" ./make_tests/run_test_balls-omp

test-octree-ggg-omp:
	$(MAKE) -B OCTREEGGGOMPON=1 all
	cd tests && MAKE_PROFILE_ARGS="OCTREEGGGOMPON=1" ./make_tests/run_test_octree-ggg-omp

test-cython:
	$(MAKE) -B OCTREEGGGOMPON=1 all
	cd tests && MAKE_PROFILE_ARGS="OCTREEGGGOMPON=1" ./make_tests/run_test_cython

test-in-stop-out-run:
	$(MAKE) -B OCTREEGGGOMPON=1 all
	cd tests && MAKE_PROFILE_ARGS="OCTREEGGGOMPON=1" ./make_tests/run_test_in-stop-out_run

test-p0-regressions: $(EXEC) cyballs
	cd tests && bash ./make_tests/run_test_p0_regressions
	cd /tmp && $(PYTHON) $(CURDIR)/tests/make_tests/test_p0_cython_instances.py

test-p1-regressions: $(EXEC) cyballs test-healpix-ordering \
	test-parameter-parser test-parameter-file-parser test-standalone-parser
	bash ./tests/make_tests/run_test_p1_regressions

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

test-p2-regressions: $(EXEC)
	cd tests && bash ./make_tests/run_test_p2_regressions

test-p2-cython: cyballs
	$(PYTHON) tests/make_tests/test_p2_cython.py

test-balls0357-recovery:
	$(MAKE) -B BALLS0357ON=1 OCTREESMOOTHINGON=0 SMOOTHPIVOTON=1 all
	$(PYTHON) tests/make_tests/test_balls0357_recovery.py

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

test-octree-ggg-fast-path: $(EXEC) cyballs
	$(PYTHON) tests/make_tests/test_octree_ggg_fast_path.py
	$(PYTHON) tests/make_tests/test_octree_ggg_only_2pcf.py
	$(PYTHON) tests/make_tests/test_octree_ggg_edge_corrections.py

test-openmp-determinism: $(EXEC) cyballs
	cd tests && bash ./make_tests/run_test_openmp_determinism
	cd tests && CBALLS=$(CURDIR)/$(EXEC) \
		bash ./make_tests/run_test_kdtree_no_one_ball

test-kdtree-no-one-ball: $(EXEC)
	cd tests && CBALLS=$(CURDIR)/$(EXEC) \
		bash ./make_tests/run_test_kdtree_no_one_ball

test-balltree-omp: $(EXEC)
	cd tests && CBALLS=$(CURDIR)/$(EXEC) \
		bash ./make_tests/run_test_balltree_omp

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

test-balltree-2balls-omp-3pcf: $(EXEC)
	cd tests && CBALLS=$(CURDIR)/$(EXEC) \
		ENGINE=balltree-2balls-omp_3pcf \
		ENGINE_MACRO=BALLTREE2BALLSOMP3PCF \
		bash ./make_tests/run_test_balltree_2balls_3pcf

test-balltree-2balls-mpi-3pcf:
	$(MAKE) -B BALLTREE2BALLSMPI3PCFON=1 TPCFON=1 cballs
	cd tests && CBALLS=$(CURDIR)/$(EXEC) MPIEXEC='$(MPIEXEC)' \
		bash ./make_tests/run_test_balltree_2balls_mpi_3pcf

test-octree-2balls-omp: $(EXEC)
	cd tests && CBALLS=$(CURDIR)/$(EXEC) \
		bash ./make_tests/run_test_octree_2balls_omp

test-octree-2balls-mask: $(EXEC)
	$(PYTHON) tests/make_tests/test_octree_2balls_mask.py --cballs $(CURDIR)/$(EXEC)

.PHONY: test-two-ball-edge
test-two-ball-edge: $(EXEC)
	$(PYTHON) tests/make_tests/test_two_ball_edge_corrections.py \
		--cballs $(CURDIR)/$(EXEC) --dimension $(DEFDIMENSION)

test-octree-2balls-mpi:
	$(MAKE) -B OCTREE2BALLSMPION=1 cballs
	cd tests && CBALLS=$(CURDIR)/$(EXEC) MPIEXEC='$(MPIEXEC)' \
		bash ./make_tests/run_test_octree_2balls_mpi

test-balltree-mpi:
	$(MAKE) -B BALLTREEMPION=1 cballs
	cd tests && CBALLS=$(CURDIR)/$(EXEC) MPIEXEC='$(MPIEXEC)' \
		bash ./make_tests/run_test_balltree_mpi

test-octree-ggg-mpi:
	$(MAKE) -B OCTREEGGGMPION=1 cballs
	cd tests && CBALLS=$(CURDIR)/$(EXEC) MPIEXEC='$(MPIEXEC)' \
		bash ./make_tests/run_test_octree_ggg_mpi

.PHONY: test-octree-3pcf-3d-mpi
test-octree-3pcf-3d-mpi: $(EXEC)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tests/make_tests/test_octree_3pcf_3d_mpi.py \
		--cballs $(CURDIR)/$(EXEC) \
		--mpi-command "$(MPIEXEC) -n 2" $(CB3D_MPI_TEST_ARGS)

test-octree-3pcf-3d-omp: $(EXEC)
	CBALLS=$(CURDIR)/$(EXEC) \
		bash ./tests/make_tests/run_test_octree_3pcf_3d_omp

.PHONY: test-octree-balls4-no-smoothing test-octree-balls4-profile
test-octree-balls4-no-smoothing:
	$(MAKE) -B OCTREESMOOTHINGON=0 BALLSON=0 BALLS4SCANLEVON=0 \
		OCTREEBALLS4OMPON=1 test-octree-balls4-profile

test-octree-balls4-profile: $(EXEC) cyballs-static-lib
	CBALLS_STATIC_LIBRARY_READY=1 $(PYTHON) setup.py build_ext --inplace
	CBALLS=$(CURDIR)/$(EXEC) PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=$(CURDIR):$(PYTHONPATH) \
		$(PYTHON) tests/make_tests/test_octree_balls4_no_smoothing.py

.PHONY: test-octree-balls4-edge test-octree-balls4-mpi
test-octree-balls4-edge: $(EXEC) cyballs-static-lib
	CBALLS_STATIC_LIBRARY_READY=1 $(PYTHON) setup.py build_ext --inplace
	PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=$(CURDIR):$(CURDIR)/python:$(PYTHONPATH) \
		$(PYTHON) tests/make_tests/test_octree_balls4_edge.py --cballs $(CURDIR)/$(EXEC) --cython

test-octree-balls4-mpi: $(EXEC)
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tests/make_tests/test_octree_balls4_edge.py \
		--cballs $(CURDIR)/$(EXEC) --mpi-command "$(if $(MPIEXEC),$(MPIEXEC),mpiexec) -n 2"

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

.PHONY: test-shear
test-shear: cyballs
	PYTHONDONTWRITEBYTECODE=1 $(PYTHON) tests/make_tests/test_shear_octree_omp.py

test-sanitizer-smoke: $(EXEC) test-cell-production
	cd tests && bash ./make_tests/run_test_sanitizer_smoke
#
#E to test cBalls under different profiles
#


ifeq ($(CLASSLIBON),1)
cyballs: cyballs-static-lib $(PYTHON_FILES) $(MAKEFILE_LIST)
	export CC=$(CC); export CBALLS_LIB=$(EXEC); export CBALLS_CPPFLAGS='$(OPT2) $(INCLUDES)'; export CBALLS_STATIC_LIBRARY_READY=1; \
	output=$$($(PYTHON) -m pip install . --no-build-isolation 2>&1); \
	status=$$?; \
	echo "$$output"; \
	if [ $$status -ne 0 ]; then \
	    if echo "$$output" | grep -q "ERROR: Cannot uninstall"; then \
	        site_packages=$$($(PYTHON) -c "import sysconfig; print(sysconfig.get_paths()['purelib'])" || $(PYTHON) -c "import site; print(site.getsitepackages()[0])"); \
	        echo "Cleaning up previous installation in: $$site_packages"; \
	        rm -rf "$$site_packages"/cyballs*.egg-info; \
	        rm -rf "$$site_packages"/cyballs*.dist-info; \
	        $(PYTHON) -m pip install .; \
	    else \
	        exit $$status; \
	    fi; \
	fi
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
	rm -f $(MDIR)/cyballs*.so $(MDIR)/cyballs*.pyd $(MDIR)/cyballs*.dylib
	rm -rf $(MDIR)/python/build
	rm -rf cyballs.egg-info
	rm -rf build;
	rm -rf dist;
