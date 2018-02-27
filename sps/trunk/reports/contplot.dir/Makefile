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
# Source files - CONTPLOT
################################################################################
CONT_SRC = \
	ContigSequence.cpp \
	ContigSpectrum.cpp \
	PlotContig.cpp \
	ContplotInterface.cpp \
	$(REPORT_COMMON_FILES)/Tokenizer.cpp

CONT_HDR = \
	ContigSequence.h \
	ContigSpectrum.h \
	PlotContig.h \
	contplot.h \
	ContplotInterface.h \
	$(REPORT_COMMON_FILES)/Tokenizer.h

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
	$(COMMON_FILES)/SpectrumPairSet.cpp

CPPLIB_HDR = \
  $(COMMON_FILES)/CommandLineParser.h \
  $(COMMON_FILES)/ParameterList.h  \
  $(COMMON_FILES)/abruijn.h \
  $(COMMON_FILES)/ms1.h \
	$(COMMON_FILES)/projectionutils.h \
	$(COMMON_FILES)/SpecnetsGraph.h \
	$(COMMON_FILES)/SpectrumPair.h \
	$(COMMON_FILES)/SpectrumPairSet.h

################################################################################
# Targets
################################################################################
EXECUTABLES = \
	contplot

all: $(EXECUTABLES)

CONT_OBJ=$(CONT_SRC:.cpp=.o)

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

contplot: contplot.o $(CONT_OBJ) $(CPPLIB_OBJ) $(LIBCORE_SRC) $(LIBSPS_SRC) $(LIBSPSREPORTS_SRC) # $(PREFIX)/bin/gnuplot.dll
	$(CXX) $(CPPFLAGS) $(LDFLAGS) contplot.o $(CONT_OBJ) $(CPPLIB_OBJ) $(LIBCORE_SRC) $(LIBSPS_SRC) $(LIBSPSREPORTS_SRC) -o $@ -Wall  $(LIBS)

################################################################################
# Compile rules - cleanup
################################################################################
clean:
	rm -f *.o $(EXECUTABLES) depends *.gch *~ *.bak > /dev/null 2>&1 | :

################################################################################

