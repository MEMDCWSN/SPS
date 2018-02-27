################################################################################
# Definitions
################################################################################
BASE_DIR := ../../..
#-------------------------------------------------------------------------------
include $(BASE_DIR)/Makefiles/Makefile.defs

################################################################################
# Conditions
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.cond

################################################################################
# Compile rules - general
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.general

################################################################################
# Source files - SPECPLOT
################################################################################
SPEC_SRC = \
	PlotSpectrum.cpp \
	SpecplotInterface.cpp \
	SpecplotInterface2.cpp
#	SpectrumLabel.cpp \

SPEC_HDR = \
	PlotSpectrum.h \
	SpecplotInterface.h \
	SpecplotInterface2.h
#	SpectrumLabel.h \

################################################################################
# Source files - CPPLIB
################################################################################
CPPLIB_SRC = \
  $(COMMON_FILES)/CommandLineParser.cpp \
  $(COMMON_FILES)/ParameterList.cpp  \
  $(COMMON_FILES)/abruijn.cpp \
  $(COMMON_FILES)/ms1.cpp \
	$(COMMON_FILES)/projectionutils.cpp \
	$(COMMON_FILES)/SpecnetsGraph.cpp \
	$(COMMON_FILES)/SpectrumPair.cpp \
	$(COMMON_FILES)/SpectrumPairSet.cpp \
  $(COMMON_FILES)/mzxml.cpp \
  $(COMMON_FILES)/PWizInterface.cpp

CPPLIB_HDR = \
  $(COMMON_FILES)/CommandLineParser.h \
  $(COMMON_FILES)/ParameterList.h  \
  $(COMMON_FILES)/abruijn.h \
  $(COMMON_FILES)/ms1.h \
  $(COMMON_FILES)/mzxml.h \
  $(COMMON_FILES)/PWizInterface.h \
	$(COMMON_FILES)/projectionutils.h \
	$(COMMON_FILES)/SpecnetsGraph.h \
	$(COMMON_FILES)/SpectrumPair.h \
	$(COMMON_FILES)/SpectrumPairSet.h

################################################################################
# Targets
################################################################################
EXECUTABLES = \
	specplot

all: $(EXECUTABLES)

################################################################################
# Compile rules - Specific
################################################################################

SPEC_OBJ=$(SPEC_SRC:.cpp=.o)

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

specplot: $(PREFIX)/$(SUFFIX)/libxerces-c.a specplot.o $(SPEC_OBJ) $(CPPLIB_OBJ) $(LIBCORE_SRC) $(LIBSPS_SRC) $(LIBSPSREPORTS_SRC)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ -Wall -lxerces-c $(LIBS) $(PWIZ_LIBS)
#	$(CXX) $(LDFLAGS) specplot.o $(SPEC_OBJ) $(CPPLIB_OBJ) $(LIBCORE_SRC) $(LIBSPSREPORTS_SRC) -o $@ -Wl,-s -lxerces-c -Wl,-Bdynamic -lncurses -lgd  $(LIBS)

################################################################################
# Compile rules - cleanup
################################################################################
clean:
	rm -f *.o $(EXECUTABLES) depends *.gch *~ *.bak > /dev/null 2>&1 | :

################################################################################

