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
# Source files - CLUSTER
################################################################################
PAIRS_SRC = \
  PairsInfo.cpp \
  PairsInfoInterface.cpp \
  $(COMMON_FILES)/projectionutils.cpp \
  $(COMMON_FILES)/SpectrumPair.cpp \
  $(COMMON_FILES)/SpecnetsGraph.cpp \
  $(COMMON_FILES)/SpectrumPairSet.cpp \
  $(COMMON_FILES)/CommandLineParser.cpp \
  $(COMMON_FILES)/ParameterList.cpp \
  $(COMMON_FILES)/Specific.cpp

PAIRS_HDR = \
  PairsInfo.h \
  PairsInfoInterface.h \
  $(COMMON_FILES)/projectionutils.h \
  $(COMMON_FILES)/SpectrumPair.h \
  $(COMMON_FILES)/SpecnetsGraph.h \
  $(COMMON_FILES)/SpectrumPairSet.h \
  $(COMMON_FILES)/CommandLineParser.h \
  $(COMMON_FILES)/ParameterList.h \
  $(COMMON_FILES)/Specific.h

################################################################################
# Targets
################################################################################
EXECUTABLES = \
	pairsinfo

all: $(EXECUTABLES)

################################################################################
# Compile rules - specific
################################################################################
PAIRS_OBJ=$(PAIRS_SRC:.cpp=.o)

pairsinfo: $(PAIRS_OBJ) $(LIBCORE_SRC) $(LIBSPSREPORTS_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -Wall $(LIBS)

################################################################################
# Compile rules - cleanup
################################################################################
clean:
	rm -f *.o $(EXECUTABLES) depends *.gch *~ *.bak > /dev/null 2>&1 | :

################################################################################
