ROOT := $(shell pwd)
XERCES_ROOT := ../../../src/xerces-c-3.1.1
GNUPLOT_ROOT := ../../../src/gnuplot-4.2.5
SPS_VERSION := $(shell svnversion .)
INSTALL_DIR := ../../../install
REPORT_COMMON_FILES := ../common
COMMON_FILES := ../..
LIBCORE_DIR := $(COMMON_FILES)/libCore
LIBCORE_SRC := $(LIBCORE_DIR)/libcore.a

LIBSPSREPORTS_DIR := ../common
LIBSPSREPORTS_SRC := $(LIBSPSREPORTS_DIR)/libSpsReports.a


CC := g++
# CC := i686-w64-mingw32-g++

CPPFLAGS := -fPIC -DSPS_VERSION=$(SPS_VERSION) -I$(XERCES_ROOT)/src -I. -I$(REPORT_COMMON_FILES) -I$(COMMON_FILES)/include -I$(COMMON_FILES) -I$(COMMON_FILES)/ExecFramework -I$(LIBCORE_DIR) -pthread

LDFLAGS :=


# 'make type=32'
ifeq ($(type),32)
  LDFLAGS += -march=i686 -m32
  CPPFLAGS += -march=i686 -m32
  PREFIX := $(INSTALL_DIR)/linux-g++

  # gcc 4
  ifeq ($(shell gcc -dumpversion | cut -d. -f 1), 4)
    LDFLAGS += -Wl,--gc-sections
  endif
endif

# 'make type=static'
ifeq ($(type),static)
  LDFLAGS += -static
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
    PREFIX ?= $(INSTALL_DIR)/linux-g++-64
  else
    PREFIX ?= $(INSTALL_DIR)/linux-g++
  endif

  LDFLAGS += -L$(PREFIX)/lib

# Windows
else
  ifneq (,$(findstring x86_64,$(shell gcc -dumpmachine)))
    PREFIX ?= $(INSTALL_DIR)/win32-g++
  else
    PREFIX ?= $(INSTALL_DIR)/win32-g++
  endif

  LDFLAGS += -L$(PREFIX)/lib -L$(PREFIX)/bin/cygwin
endif




# type=debug
ifeq ($(type),debug)
  CPPFLAGS += -O0 -g
else
  CPPFLAGS += -O3 -ffunction-sections -fdata-sections -Wno-deprecated -Wno-write-strings
endif



CONT_SRC = \
	ContigSpectrum.cpp \
	PlotContig.cpp \
	AlignplotInterface.cpp \
	$(REPORT_COMMON_FILES)/Tokenizer.cpp

CONT_HDR = \
	ContigSpectrum.h \
	PlotContig.h \
	alignplot.h \
	AlignplotInterface.h \
	$(REPORT_COMMON_FILES)/Tokenizer.h


CPPLIB_SRC = \
  $(COMMON_FILES)/ExecFramework/CommandLineParser.cpp \
  $(COMMON_FILES)/ExecFramework/ParameterList.cpp  \
  $(COMMON_FILES)/alignment_scoring.cpp \
  $(COMMON_FILES)/base64.cpp \
  $(COMMON_FILES)/db_fasta.cpp \
  $(COMMON_FILES)/abruijn.cpp \
  $(COMMON_FILES)/ms1.cpp \
	$(COMMON_FILES)/SpectrumPair.cpp \
	$(COMMON_FILES)/SpectrumPairSet.cpp \
	$(COMMON_FILES)/mzxml.cpp

CPPLIB_HDR = \
  $(COMMON_FILES)/ExecFramework/CommandLineParser.h \
  $(COMMON_FILES)/ExecFramework/ParameterList.h  \
  $(COMMON_FILES)/alignment_scoring.h \
  $(COMMON_FILES)/base64.h \
  $(COMMON_FILES)/db_fasta.h \
  $(COMMON_FILES)/abruijn.h \
  $(COMMON_FILES)/ms1.h \
	$(COMMON_FILES)/SpectrumPair.h \
	$(COMMON_FILES)/SpectrumPairSet.h \
	$(COMMON_FILES)/mzxml.h



CONT_OBJ=$(CONT_SRC:.cpp=.o)

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)



EXECUTABLES = \
	alignplot



all: $(EXECUTABLES)
	# strip  $(EXECUTABLES)


clean:
	rm -f *.o $(EXECUTABLES) depends *.gch *~ > /dev/null 2>&1 | :


.cpp.o: $(CPPLIB_SRC)
	$(CC) $(CPPFLAGS) $< -c -o $@

.o.o: $(CPPLIB_SRC) $(CPPLIB_HDR)
	$(CC) $(CPPFLAGS) $<


$(GNUPLOT_ROOT)/Makefile:
	cd $(GNUPLOT_ROOT); ./configure --without-readline --without-lisp-files --without-tutorial --without-row-help --without-x --disable-wxwidgets; cd $(ROOT)



# Linux

ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))

$(XERCES_ROOT)/Makefile:
	chmod +x $(XERCES_ROOT)/configure
	chmod +x $(XERCES_ROOT)/config/pretty-make
	cd $(XERCES_ROOT); ./configure --disable-pretty-make --disable-network --disable-threads --without-curl --without-icu --disable-shared ; cd $(ROOT)

$(PREFIX)/lib/libxerces-c.a: $(XERCES_ROOT)/Makefile
	$(MAKE) -C $(XERCES_ROOT)/src LDFLAGS="$(LDFLAGS)" CPPFLAGS="$(CPPFLAGS)"
	mv $(XERCES_ROOT)/src/.libs/libxerces-c.a $(PREFIX)/lib

alignplot: $(PREFIX)/lib/libxerces-c.a alignplot.o $(CONT_OBJ) $(CPPLIB_OBJ) $(LIBCORE_SRC) $(LIBSPSREPORTS_SRC)
	$(CC) $(CPPFLAGS) $(LDFLAGS) alignplot.o $(CONT_OBJ) $(CPPLIB_OBJ) $(LIBCORE_SRC) $(LIBSPSREPORTS_SRC) -o $@ -Wall -lxerces-c


# windows

else
ifneq (,$(findstring cygwin,$(shell gcc -dumpmachine)))

$(XERCES_ROOT)/Makefile:
	chmod +x $(XERCES_ROOT)/configure
	chmod +x $(XERCES_ROOT)/config/pretty-make
	cd $(XERCES_ROOT); ./configure LDFLAGS=-no-undefined --disable-pretty-make --disable-network --disable-threads --without-curl --without-icu --disable-shared ; cd $(ROOT)
	# cd $(XERCES_ROOT); ./configure --disable-pretty-make --disable-network --disable-threads --without-curl --without-icu --host=i686-w64-mingw32 --enable-static --disable-shared LDFLAGS=-no-undefined CC=$(CC) CXX=$(CC) ; cd $(ROOT)

$(PREFIX)/bin/cygwin/libxerces-c.a: $(XERCES_ROOT)/Makefile
	$(MAKE) -C $(XERCES_ROOT)/src LDFLAGS="$(LDFLAGS)" CPPFLAGS="$(CPPFLAGS)"
	mv $(XERCES_ROOT)/src/.libs/libxerces-c.a $(PREFIX)/bin/cygwin

$(PREFIX)/bin/gnuplot.dll: $(GNUPLOT_ROOT)/Makefile
	$(MAKE) -C $(GNUPLOT_ROOT)/src gnuplot.exe LDFLAGS="$(LDFLAGS) -shared -Wl,-o,gnuplot.dll" CFLAGS="$(CPPFLAGS) -fPIC -Dmain=gnuplot_main"
	mv $(GNUPLOT_ROOT)/src/gnuplot.dll $(PREFIX)/bin

alignplot: $(PREFIX)/bin/cygwin/libxerces-c.a alignplot.o $(CONT_OBJ) $(CPPLIB_OBJ) $(LIBCORE_SRC) $(LIBSPSREPORTS_SRC) # $(PREFIX)/bin/gnuplot.dll
	$(CC) $(CPPFLAGS) $(LDFLAGS) alignplot.o $(CONT_OBJ) $(CPPLIB_OBJ) $(LIBCORE_SRC) $(LIBSPSREPORTS_SRC) -o $@ -Wall -lxerces-c
	strip contplot


# end
endif
endif
