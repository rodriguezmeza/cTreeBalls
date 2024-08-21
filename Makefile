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

all: $(EXEC) lib$(EXEC).a
# With cpython:
#all: $(EXEC) lib$(EXEC).a $(EXEC)y

lib$(EXEC).a: $(OBJS) $(EXTERNAL)
	$(AR)  $@ $(addprefix build/, $(OBJS) $(TOOLS) $(SOURCE) $(EXTERNAL) $(EXTERNALCXX))

#B To work with cpp files...
ifeq ($(NR3ON),1)
include $(MACHINES_DIR)addons/nr3/Makefile_nr3
endif
#E

$(EXEC): $(OBJS) $(EXTERNAL) $(MAIN)
	$(CC) $(OPTFLAG) $(OMPFLAG) $(LDFLAG) -o $(EXEC) $(addprefix build/,$(notdir $^)) $(MLIBS)


$(EXEC)y: lib$(EXEC).a python/$(EXEC)y.pyx python/c$(EXEC)y.pxd
ifdef OMPFLAG
	cp python/setup.py python/autosetup.py
else
	grep -v "lgomp lgsl lgslcblas" python/setup.py > python/autosetup.py
endif
	cd python; export CC=$(CC); $(PYTHON) autosetup.py install || $(PYTHON) autosetup.py install --user
	rm python/autosetup.py


.PHONY : clean
clean: .base
	rm -rf $(WRKDIR);
	rm -f $(EXEC)
	rm -f lib$(EXEC).a

#	rm -f $(MDIR)/python/$(EXEC)y.c
#	rm -rf $(MDIR)/python/build
