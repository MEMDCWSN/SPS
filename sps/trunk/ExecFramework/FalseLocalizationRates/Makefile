XERCES_ROOT := ../../../src/xerces-c-3.1.1
ROOT := $(shell pwd)
#SPS_VERSION := $(shell date +%y%m%d)
LIBCORE_DIR := ../../libCore
LIBCORE_SRC := $(LIBCORE_DIR)/libcore.a
#SPS_VERSION := $(shell svnversion .)

CC := g++
#CPPFLAGS := -DSPS_VERSION=$(SPS_VERSION) -I$(XERCES_ROOT)/src -I.. -I../.. -I../../FalseLocalizationRates
CPPFLAGS := -DSPS_VERSION=$(SPS_VERSION) -I$(LIBCORE_DIR) -I$(XERCES_ROOT)/src -I.. -I../.. -I../../FalseLocalizationRates

type=debug
ifeq ($(type),debug)
  CPPFLAGS += -O0 -g
else
  CPPFLAGS += -O3 -ffunction-sections -fdata-sections -Wno-deprecated -Wno-write-strings
endif

# 'make type=32'
ifeq ($(type),32)
  LDFLAGS += -march=i686 -m32
  CPPFLAGS += -march=i686 -m32
  PREFIX := ../../../install/linux-g++

  # gcc 4
  ifeq ($(shell gcc -dumpversion | cut -d. -f 1), 4)
    LDFLAGS += -Wl,--gc-sections
  endif
endif

# 'make'
ifndef type
  # gcc 4
  ifeq ($(shell gcc -dumpversion | cut -d. -f 1), 4)
    LDFLAGS += -Wl,--gc-sections
    CPPFLAGS += -mtune=generic

  # gcc 3
  else
    CPPFLAGS += -march=$(shell uname -m)
  endif
endif

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
  ifneq (,$(findstring x86_64,$(shell gcc -dumpmachine)))
    PREFIX ?= ../../../install/linux-g++-64
  else
    PREFIX ?= ../../../install/linux-g++
  endif

  LDFLAGS += -L. -L.. -L$(PREFIX)/lib 

# Windows
else
  ifneq (,$(findstring x86_64,$(shell gcc -dumpmachine)))
    PREFIX ?= ../../../install/win32-g++-64
  else
    PREFIX ?= ../../../install/win32-g++
  endif

  LDFLAGS += -L. -L../.. -L$(LIBCORE_DIR) -L../../FalseLocalizationRates -L$(PREFIX)/lib -L$(PREFIX)/bin/cygwin
endif

LIBRARIES = -lpthread -lspecnets -lflr -lcore

CPPLIB_SRC = \
  ../CommandLineParser.cpp \
  ../ExecBase.cpp \
  ExecModuleFactoryLP.cpp \
  ExecPairSpectra.cpp \
  ExecFlr.cpp \
  ExecLpSolver.cpp \
  FileLoadingUtils.cpp \
  ../FileUtils.cpp \
  ../Filters.cpp \
  ../ParameterList.cpp \
  ../ParallelExecution.cpp \
  ../ParallelPbsExecution.cpp \
  ../ParallelSgeExecution.cpp \
  ../ParallelThreadedExecution.cpp \
  ../ExecMergeConvert.cpp \
  ../SgeGridMonitor.cpp \
  ../../Specific.cpp \
  ../../StatusFile.cpp \
  
CPPLIB_HDR = \
  ../CommandLineParser.h \
  ../ExecBase.h \
  ExecModuleFactoryLP.h \
  ExecPairSpectra.h \
  ExecFlr.h \
  ExecLpSolver.h \
  FileLoadingUtils.h \
  ../FileUtils.h \
  ../Filters.h \
  ../ParameterList.h \
  ../ParallelExecution.h \
  ../ParallelPbsExecution.h \
  ../ParallelSgeExecution.h \
  ../ParallelThreadedExecution.h \
    ../ExecMergeConvert.h \
  ../SgeGridMonitor.h \
  ../../Specific.h \
  ../../StatusFile.h \

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

EXECUTABLES = \
	main_execmodule \
	main_flr

.PHONY: subdirs $(SUBDIRS) clean all

all: $(EXECUTABLES) main_execmodule

clean:
	rm -f *.o $(EXECUTABLES) 

.o.o: $(CPPLIB_SRC) $(CPPLIB_HDR)
	$(CC) $(CPPFLAGS) $<

main_execmodule: main_execmoduleLP.o $(CPPLIB_OBJ) $(CPPLIB_HDR)
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LIBRARIES)

main_flr: main_flr.o $(CPPLIB_OBJ) $(CPPLIB_HDR)
	$(CC) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LIBRARIES)

# dependencies
depends : $(CPPLIB_SRC) $(CPPLIB_HDR)
	@$(CC) $(CPPFLAGS) -MM $^ > depends

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
-include depends
endif

