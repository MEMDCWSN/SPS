################################################################################
# Definitions
################################################################################
BASE_DIR := ../..
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
# Source files - SPECCHECK
################################################################################
SPECCHECK_SRC = \
  specCheck.cpp \
  SpecCheckInterface.cpp \
  $(LIBCORE_DIR)/utils.cpp \
  $(COMMON_FILES)/base64.cpp \
  $(COMMON_FILES)/mzxml.cpp \
  $(COMMON_FILES)/ExecFramework/CommandLineParser.cpp \
  $(COMMON_FILES)/ExecFramework/ParameterList.cpp

SPECCHECK_HDR = \
  specCheck.h \
  SpecCheckInterface.h \
  $(LIBCORE_DIR)/utils.h \
  $(COMMON_FILES)/base64.h \
  $(COMMON_FILES)/mzxml.h \
  $(COMMON_FILES)/ExecFramework/CommandLineParser.h \
  $(COMMON_FILES)/ExecFramework/ParameterList.h

################################################################################
# Compile rules - specific
################################################################################
SPECCHECK_OBJ=$(SPECCHECK_SRC:.cpp=.o)

specCheck: $(SPECCHECK_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -Wall -lxerces-c $(LIBS)

################################################################################
# Compile rules - target
################################################################################
EXECUTABLES = \
	specCheck

all: $(EXECUTABLES)
	$(STRIP) $(EXECUTABLES)

################################################################################
# Compile rules - cleanup
################################################################################
clean:
	rm -f *.o $(EXECUTABLES) *.gch *~

################################################################################
