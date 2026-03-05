include ../../Makefile.common

ROOT = ../..
COMETSEARCH = $(ROOT)/CometSearch
STATIC ?= 0
WITH_MZMLB ?= 0
HDF5_DIR ?=
HDF5_LIB_C ?= $(HDF5_DIR)/lib/libhdf5.a
HDF5_LIB_CPP ?= $(HDF5_DIR)/lib/libhdf5_cpp.a
HDF5_EXTRA_LIBS ?= -ldl -lz -lm -lpthread
MSTOOLKIT_MAKE_ARGS :=

ifeq ($(STATIC),1)
   override CXX := x86_64-conda-linux-gnu-g++
   override CC := x86_64-conda-linux-gnu-gcc
endif

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
   override CXXFLAGS += -O3 -std=c++14 -fpermissive -Wall -Wextra -Wno-char-subscripts -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D__LINUX__ -D_NOSQLITE \
      -I$(ROOT)/$(MSTOOLKIT)/include -I$(ROOT)/$(MSTOOLKIT)/extern/expat-2.2.9/lib -I$(ROOT)/$(MSTOOLKIT)/extern/zlib-1.2.11 \
      -I$(ROOT)/CometSearch -I$(ROOT)/$(ASCOREPRO)/include
else
   override CXXFLAGS += -O3 -std=c++14 -fpermissive -Wall -Wextra -Wno-char-subscripts -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D__LINUX__ -D_NOSQLITE \
      -I$(ROOT)/$(MSTOOLKIT)/include -I$(ROOT)/$(MSTOOLKIT)/extern/expat-2.2.9/lib -I$(ROOT)/$(MSTOOLKIT)/extern/zlib-1.2.11 \
      -I$(ROOT)/CometSearch -I$(ROOT)/$(ASCOREPRO)/include
endif

EXECNAME = cometplus
WORKER_EXECNAME = cometplus_prefilter_worker
OBJDIR = obj
OBJS = $(OBJDIR)/CometPlus.o
OBJS += $(OBJDIR)/NovelModeUtils.o
OBJS += $(OBJDIR)/NovelModeUtilsMerge.o
OBJS += $(OBJDIR)/NovelModeUtilsIdx.o
OBJS += $(OBJDIR)/NovelModeUtilsPrefilter.o
OBJS += $(OBJDIR)/CometPlusInputFiles.o
OBJS += $(OBJDIR)/CometPlusProcessCmd.o
OBJS += $(OBJDIR)/CometPlusRuntimeUtils.o
OBJS += $(OBJDIR)/CometPlusNovelWorkflow.o
OBJS += $(OBJDIR)/CometPlusPrefilterWorkflow.o
OBJS += $(OBJDIR)/CometPlusPrintParams.o
OBJS += $(OBJDIR)/CometPlusParamHelp.o
OBJS += $(OBJDIR)/CometPlusParamParser.o
WORKER_OBJS = $(OBJDIR)/PrefilterWorker.o
WORKER_OBJS += $(OBJDIR)/NovelModeUtils.o
WORKER_OBJS += $(OBJDIR)/NovelModeUtilsMerge.o
WORKER_OBJS += $(OBJDIR)/NovelModeUtilsIdx.o
WORKER_OBJS += $(OBJDIR)/NovelModeUtilsPrefilter.o
DEPS_STAMP = $(OBJDIR)/.deps_built

LIBPATHS = -L$(ROOT)/$(MSTOOLKIT) -L$(ROOT)/CometSearch -L$(ROOT)/$(ASCOREPRO)
LIBS = -lcometsearch -lmstoolkit -lmstoolkitextern -lascorepro -lm -lpthread
LDFLAGS :=

ifeq ($(WITH_MZMLB),1)
ifeq ($(strip $(HDF5_DIR)),)
$(error WITH_MZMLB=1 requires HDF5_DIR=/absolute/path/to/hdf5-static)
endif
ifeq ($(wildcard $(HDF5_LIB_C)),)
$(error WITH_MZMLB=1 requires static HDF5 C library: $(HDF5_LIB_C))
endif
ifeq ($(wildcard $(HDF5_LIB_CPP)),)
$(error WITH_MZMLB=1 requires static HDF5 C++ library: $(HDF5_LIB_CPP))
endif
MSTOOLKIT_MAKE_ARGS += HDF5=1 HDF5_DIR=$(HDF5_DIR) HDF5_LIB_C=$(HDF5_LIB_C) HDF5_LIB_CPP=$(HDF5_LIB_CPP)
LIBS += $(HDF5_EXTRA_LIBS)
endif

ifdef MSYSTEM
   LIBS += -lws2_32
endif

ifeq ($(STATIC),1)
   OBJS += $(OBJDIR)/StaticCompat.o
   WORKER_OBJS += $(OBJDIR)/StaticCompat.o
   LDFLAGS += -static
endif

all: $(EXECNAME) $(WORKER_EXECNAME)

static:
	$(MAKE) STATIC=1 all

$(DEPS_STAMP): | $(OBJDIR)
ifeq ($(WITH_MZMLB),1)
	$(MAKE) -C $(ROOT)/$(MSTOOLKIT) CC=$(CC) CXX=$(CXX) $(MSTOOLKIT_MAKE_ARGS) clean
endif
	$(MAKE) -C $(ROOT)/$(MSTOOLKIT) CC=$(CC) CXX=$(CXX) $(MSTOOLKIT_MAKE_ARGS) all
	$(MAKE) -C $(ROOT)/$(ASCOREPRO) CC=$(CC) CXX=$(CXX) all
	$(MAKE) -C $(COMETSEARCH) CC=$(CC) CXX=$(CXX) all
	@touch $(DEPS_STAMP)

$(EXECNAME): $(OBJS) $(DEPS_STAMP)
ifeq ($(UNAME_S),Darwin)
	${CXX} $(OBJS) -headerpad_max_install_names -o ${EXECNAME} $(LDFLAGS) $(CXXFLAGS) $(LIBPATHS) $(LIBS)
else
	${CXX} $(OBJS) -o ${EXECNAME} $(LDFLAGS) $(CXXFLAGS) $(LIBPATHS) $(LIBS)
endif

$(WORKER_EXECNAME): $(WORKER_OBJS) $(DEPS_STAMP)
ifeq ($(UNAME_S),Darwin)
	${CXX} $(WORKER_OBJS) -headerpad_max_install_names -o ${WORKER_EXECNAME} $(LDFLAGS) $(CXXFLAGS) $(LIBPATHS) $(LIBS)
else
	${CXX} $(WORKER_OBJS) -o ${WORKER_EXECNAME} $(LDFLAGS) $(CXXFLAGS) $(LIBPATHS) $(LIBS)
endif

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/CometPlus.o: CometPlus/CometPlus.cpp | $(OBJDIR)
	${CXX} ${CXXFLAGS} CometPlus/CometPlus.cpp -c -o $@

$(OBJDIR)/NovelModeUtils.o: CometPlus/NovelModeUtils.cpp | $(OBJDIR)
	${CXX} ${CXXFLAGS} CometPlus/NovelModeUtils.cpp -c -o $@

$(OBJDIR)/NovelModeUtilsMerge.o: CometPlus/NovelModeUtilsMerge.cpp | $(OBJDIR)
	${CXX} ${CXXFLAGS} CometPlus/NovelModeUtilsMerge.cpp -c -o $@

$(OBJDIR)/NovelModeUtilsIdx.o: CometPlus/NovelModeUtilsIdx.cpp | $(OBJDIR)
	${CXX} ${CXXFLAGS} CometPlus/NovelModeUtilsIdx.cpp -c -o $@

$(OBJDIR)/NovelModeUtilsPrefilter.o: CometPlus/NovelModeUtilsPrefilter.cpp | $(OBJDIR)
	${CXX} ${CXXFLAGS} CometPlus/NovelModeUtilsPrefilter.cpp -c -o $@

$(OBJDIR)/CometPlusInputFiles.o: CometPlus/CometPlusInputFiles.cpp | $(OBJDIR)
	${CXX} ${CXXFLAGS} CometPlus/CometPlusInputFiles.cpp -c -o $@

$(OBJDIR)/CometPlusProcessCmd.o: CometPlus/CometPlusProcessCmd.cpp | $(OBJDIR)
	${CXX} ${CXXFLAGS} CometPlus/CometPlusProcessCmd.cpp -c -o $@

$(OBJDIR)/CometPlusRuntimeUtils.o: CometPlus/CometPlusRuntimeUtils.cpp | $(OBJDIR)
	${CXX} ${CXXFLAGS} CometPlus/CometPlusRuntimeUtils.cpp -c -o $@

$(OBJDIR)/CometPlusNovelWorkflow.o: CometPlus/CometPlusNovelWorkflow.cpp | $(OBJDIR)
	${CXX} ${CXXFLAGS} CometPlus/CometPlusNovelWorkflow.cpp -c -o $@

$(OBJDIR)/CometPlusPrefilterWorkflow.o: CometPlus/CometPlusPrefilterWorkflow.cpp | $(OBJDIR)
	${CXX} ${CXXFLAGS} CometPlus/CometPlusPrefilterWorkflow.cpp -c -o $@

$(OBJDIR)/CometPlusPrintParams.o: CometPlus/CometPlusPrintParams.cpp | $(OBJDIR)
	${CXX} ${CXXFLAGS} CometPlus/CometPlusPrintParams.cpp -c -o $@

$(OBJDIR)/CometPlusParamHelp.o: CometPlus/CometPlusParamHelp.cpp | $(OBJDIR)
	${CXX} ${CXXFLAGS} CometPlus/CometPlusParamHelp.cpp -c -o $@

$(OBJDIR)/CometPlusParamParser.o: CometPlus/CometPlusParamParser.cpp | $(OBJDIR)
	${CXX} ${CXXFLAGS} CometPlus/CometPlusParamParser.cpp -c -o $@

$(OBJDIR)/PrefilterWorker.o: CometPlus/PrefilterWorker.cpp | $(OBJDIR)
	${CXX} ${CXXFLAGS} CometPlus/PrefilterWorker.cpp -c -o $@

$(OBJDIR)/StaticCompat.o: CometPlus/StaticCompat.cpp | $(OBJDIR)
	${CC} -O2 -c CometPlus/StaticCompat.cpp -o $@

clean:
	rm -rf $(OBJDIR) $(EXECNAME) $(WORKER_EXECNAME)
