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
	runtime_context.o memory_catalog.o common_histogram.o smooth_pivots.o mpi_runtime.o scalar_window_io.o

PYTHON_FILES = python/cyballs.pyx python/resource_api.pxi python/qualification_api.pxi setup.py python/ccyballs.pxd.in
CBALLS_OBJECTS = $(sort $(OBJS) $(TOOLS) $(SOURCE) $(EXTERNAL) $(EXTERNALCXX))
CBALLS_COMPILE_CONFIG = $(WRKDIR)/compile-settings.json
CBALLS_LINK_CONFIG = $(WRKDIR)/link-settings.json
export CBALLS_OBJECTS VENDORED_SOURCE_PATTERNS

all: $(EXEC) lib$(EXEC).a cyballs

# Reader headers may be dropped by later legacy HEADERFILES assignments.
main.o: $(MDIR)/include/cli_validation.h $(MDIR)/include/ia_tutor_cli.h

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
	@printf '__CBALLS_GSL_CFLAGS__=%s\n' '$(GSL_CFLAGS)'
	@printf '__CBALLS_GSL_LDFLAGS__=%s\n' '$(GSL_LDFLAGS)'
	@printf '__CBALLS_GSL_LIBS__=%s\n' '$(GSL_LIBS)'
	@printf '__CBALLS_CFITSIO_CFLAGS__=%s\n' '$(CFITSIO_INCL)'
	@printf '__CBALLS_CFITSIO_LDFLAGS__=%s\n' '$(CFITSIO_LDFLAGS)'
	@printf '__CBALLS_CFITSIO_LIBS__=%s\n' '$(CFITSIO_LIBS)'
	@printf '__CBALLS_CFITSIO_RPATH__=%s\n' '$(CFITSIO_RPATH)'
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


# Test, regression, and benchmark targets.
include $(MACHINES_DIR)Makefile_tests


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

# Resolve settings before stamping both native and Cython compilation units.
CBALLS_BUILD_HEADER = $(WRKDIR)/cballs_build_fingerprint.h
CBALLS_FINGERPRINT_VARIABLES = $(CBALLS_MAKE_INFO_VARIABLES) OBJS TOOLS SOURCE EXTERNALCXX EXEC
export CBALLS_FINGERPRINT_VARIABLES OBJS TOOLS SOURCE EXTERNALCXX EXEC
.PHONY: cballs-fingerprint-force print-build-fingerprint
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

# Do not treat .d files as profile Makefiles in the fingerprint prerequisites.
-include $(wildcard $(addprefix $(WRKDIR)/,$(CBALLS_OBJECTS:.o=.d)))

# Generation is explicit; every build rejects a stale checked-in catalogue.
.PHONY: check-capabilities generate-capabilities
check-capabilities:
	@$(PYTHON) scripts/generate_capabilities.py --check
generate-capabilities:
	@$(PYTHON) scripts/generate_capabilities.py
$(CBALLS_OBJECTS): | check-capabilities
