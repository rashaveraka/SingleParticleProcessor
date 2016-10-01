# Module.mk for pyroot module
# Copyright (c) 2004 Rene Brun and Fons Rademakers
#
# Authors: Pere Mato, Wim Lavrijsen, 22/4/2004

MODNAME      := pyroot
MODDIR       := $(ROOT_SRCDIR)/bindings/$(MODNAME)
MODDIRS      := $(MODDIR)/src
MODDIRI      := $(MODDIR)/inc

PYROOTDIR    := $(MODDIR)
PYROOTDIRS   := $(PYROOTDIR)/src
PYROOTDIRI   := $(PYROOTDIR)/inc

##### python64 #####
ifeq ($(ARCH),macosx64)
ifeq ($(MACOSX_MINOR),5)
PYTHON64S    := $(MODDIRS)/python64.c
PYTHON64O    := $(call stripsrc,$(PYTHON64S:.c=.o))
PYTHON64     := bin/python64
PYTHON64DEP  := $(PYTHON64O:.o=.d)
endif
endif

##### libPyROOT #####
PYROOTL      := $(MODDIRI)/LinkDef.h
PYROOTDS     := $(call stripsrc,$(MODDIRS)/G__PyROOT.cxx)
PYROOTDO     := $(PYROOTDS:.cxx=.o)
PYROOTDH     := $(PYROOTDS:.cxx=.h)

PYROOTH      := $(filter-out $(MODDIRI)/LinkDef%,$(wildcard $(MODDIRI)/*.h))
PYROOTS      := $(filter-out $(MODDIRS)/G__%,$(wildcard $(MODDIRS)/*.cxx))
PYROOTO      := $(call stripsrc,$(PYROOTS:.cxx=.o))

PYROOTDEP    := $(PYROOTO:.o=.d) $(PYROOTDO:.o=.d)

PYROOTLIB    := $(LPATH)/libPyROOT.$(SOEXT)
ifeq ($(ARCH),win32)
PYROOTPYD    := bin/$(notdir $(PYROOTLIB:.$(SOEXT)=.pyd))
endif
PYROOTMAP    := $(PYROOTLIB:.$(SOEXT)=.rootmap)

ROOTPYS      := $(wildcard $(MODDIR)/*.py)
ROOTAASS     := $(wildcard $(MODDIR)/ROOTaaS/* $(MODDIR)/ROOTaaS/*/* $(MODDIR)/ROOTaaS/*/*/*)
# Above includes ROOTaaS/config which is a directory; filter those out.
# Problem: $(dir $(ROOTAASS)) gives ROOTaaS/config/ thus patsubst %/, %
ROOTAASS     := $(filter-out $(sort $(patsubst %/,%,$(dir $(ROOTAASS)))),$(ROOTAASS))

ifeq ($(ARCH),win32)
ROOTPY       := $(subst $(MODDIR),bin,$(ROOTPYS))
ROOTAAS      := $(subst $(MODDIR),bin,$(ROOTAASS))
bin/%.py: $(MODDIR)/%.py
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cp $< $@
bin/ROOTaaS/%: $(MODDIR)/ROOTaaS/%
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cp -R $< $@
else
ROOTPY       := $(subst $(MODDIR),$(LPATH),$(ROOTPYS))
ROOTAAS      := $(subst $(MODDIR),$(LPATH),$(ROOTAASS))
$(LPATH)/%.py: $(MODDIR)/%.py
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cp $< $@
$(LPATH)/ROOTaaS/%: $(MODDIR)/ROOTaaS/%
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cp -R $< $@
endif
ROOTPYC      := $(ROOTPY:.py=.pyc)
ROOTPYO      := $(ROOTPY:.py=.pyo)
ROOTAASC     := $(ROOTAAS:.py=.pyc)
ROOTAASO     := $(ROOTAAS:.py=.pyo)

# used in the main Makefile
ALLHDRS      += $(patsubst $(MODDIRI)/%.h,include/%.h,$(PYROOTH))
ALLLIBS      += $(PYROOTLIB)
ALLMAPS      += $(PYROOTMAP)

ALLEXECS     += $(PYTHON64)
INCLUDEFILES += $(PYTHON64DEP)

# include all dependency files
INCLUDEFILES += $(PYROOTDEP)

##### local rules #####
.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME)

include/%.h:    $(PYROOTDIRI)/%.h
		cp $< $@

%.pyc: %.py;    python -c 'import py_compile; py_compile.compile( "$<" )'
%.pyo: %.py;    python -O -c 'import py_compile; py_compile.compile( "$<" )'

$(PYROOTLIB):   $(PYROOTO) $(PYROOTDO) $(ROOTPY) $(ROOTPYC) $(ROOTPYO) \
                $(ROOTLIBSDEP) $(PYTHONLIBDEP) \
                $(ROOTAAS) $(ROOTAASC) $(ROOTAASO)

		@$(MAKELIB) $(PLATFORM) $(LD) "$(LDFLAGS)" \
		  "$(SOFLAGS)" libPyROOT.$(SOEXT) $@ \
		  "$(PYROOTO) $(PYROOTDO)" \
		  "$(ROOTULIBS) $(RPATH) $(ROOTLIBS) $(PYROOTLIBEXTRA) \
		   $(PYTHONLIBDIR) $(PYTHONLIB) $(PYTHONLIBFLAGS)"
ifeq ($(ARCH),win32)
	link -dll -nologo -IGNORE:4001 -machine:ix86 -export:initlibPyROOT \
	lib/libPyROOT.lib -nodefaultlib kernel32.lib msvcrt.lib \
	-out:$(PYROOTPYD)
	@(if [ -f $(PYROOTPYD).manifest ]; then \
		mt -nologo -manifest $(PYROOTPYD).manifest \
			-outputresource\:$(PYROOTPYD)\;2 ; \
		rm -f $(PYROOTPYD).manifest ; \
	fi)
	@rm -f bin/libPyROOT.lib
	@rm -f bin/libPyROOT.exp
endif

$(call pcmrule,PYROOT)
	$(noop)

$(PYROOTDS):    $(PYROOTH) $(PYROOTL) $(ROOTCLINGEXE) $(call pcmdep,PYROOT)
		$(MAKEDIR)
		@echo "Generating dictionary $@..."
		$(ROOTCLINGSTAGE2) -f $@ $(call dictModule,PYROOT) -c -writeEmptyRootPCM $(PYROOTH) $(PYROOTL)

$(PYROOTMAP):   $(PYROOTH) $(PYROOTL) $(ROOTCLINGEXE) $(call pcmdep,PYROOT)
		$(MAKEDIR)
		@echo "Generating rootmap $@..."
		$(ROOTCLINGSTAGE2) -r $(PYROOTDS) $(call dictModule,PYROOT) -c $(PYROOTH) $(PYROOTL)

$(PYTHON64):    $(PYTHON64O)
		$(CC) $(LDFLAGS) -o $@ $(PYTHON64O) \
		   $(PYTHONLIBDIR) $(PYTHONLIB)

all-$(MODNAME): $(PYROOTLIB) $(PYTHON64)

clean-$(MODNAME):
		@rm -f $(PYROOTO) $(PYROOTDO) $(PYTHON64O)

clean::         clean-$(MODNAME)

distclean-$(MODNAME): clean-$(MODNAME)
		@rm -f $(PYROOTDEP) $(PYROOTDS) $(PYROOTDH) $(PYROOTLIB) \
		   $(ROOTPY) $(ROOTPYC) $(ROOTPYO) $(PYROOTMAP) \
		   $(PYROOTPYD) $(PYTHON64DEP) $(PYTHON64)
		@rm -rf $(LPATH)/ROOTaaS bin/ROOTaaS

distclean::     distclean-$(MODNAME)

##### extra rules ######
$(PYROOTO): CXXFLAGS += $(PYTHONINCDIR:%=-I%)
$(PYTHON64O): CFLAGS += $(PYTHONINCDIR:%=-I%)
ifeq ($(GCC_MAJOR),4)
$(PYROOTO): CXXFLAGS += -fno-strict-aliasing
endif
ifneq ($(CLANG_MAJOR)$(GCC_MAJOR),)
# Building with clang or GCC
$(PYROOTO) $(PYTHON64O) $(PYROOTDO): CXXFLAGS += -Wno-error=format 
endif

ifneq ($(CLANG_MAJOR),)
# Building with clang 
$(PYROOTO) $(PYTHON64O) $(PYROOTDO): CXXFLAGS += -Wno-ignored-attributes
endif
