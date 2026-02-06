# ----- MAKE FILE FOR cballs CODE -----
# Mario A. Rodriguez-Meza, Ciudad de Mexico, 27.04.2023
#
#
# Nothing to do in this file. Make your settings in Makefile_settings file only
#
MACHINES_DIR = ./
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

OBJS = main.o cballsio.o cballs.o startrun.o testdata.o treeload.o \
	cballsutils.o search.o

PYTHON_FILES = python/cballys.pyx python/setup.py python/ccballys.pxd

all: $(EXEC) lib$(EXEC).a cballys

lib$(EXEC).a: $(OBJS) $(EXTERNAL)
	$(AR)  $@ $(addprefix build/, $(OBJS) $(TOOLS) $(SOURCE) $(EXTERNAL) $(EXTERNALCXX))

$(EXEC): $(OBJS) $(EXTERNAL) $(MAIN)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $(EXEC) $(addprefix build/,$(notdir $^)) $(MLIBS) $(FITSIOLIBS)


cballys: libcballs.a python/cballys.pyx python/ccballys.pxd
	export CC=$(CC); output=$$($(PYTHON) -m pip install . 2>&1); \
    echo "$$output"; \
    if echo "$$output" | grep -q "ERROR: Cannot uninstall"; then \
        site_packages=$$($(PYTHON) -c "import distutils.sysconfig; print(distutils.sysconfig.get_python_lib())" || $(PYTHON) -c "import site; print(site.getsitepackages()[0])") && \
        echo "Cleaning up previous installation in: $$site_packages" && \
        rm -rf $$site_packages/cballys* && \
        $(PYTHON) -m pip install .; \
    fi

.PHONY : clean
clean: .base
	rm -rf $(WRKDIR);
	rm -f $(EXEC)
	rm -f lib$(EXEC).a
	rm -f $(MDIR)/python/cballys.c
	rm -rf $(MDIR)/python/build

