################################################################################
# Makefile for the reports library
#
# The reports library contain classes and methods designed specifically to
# generate work with the SPS framweork report generation module.
# The classes included are used in several report generation tools and
# programs, like specplot, contplot, spsReportServer and spsReports.

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
# Source files - CPPREPORTLIB
################################################################################
CPPREPORTLIB_SRC = \
	ClusterData.cpp \
	ReportTableGenerator.cpp \
	ReportTableBase.cpp \
	ReportTableHeader.cpp \
	ReportTableProtein.cpp \
	ReportTableProteinCoverage.cpp \
	ReportTableContig.cpp \
	ReportTableClusterConsensus.cpp \
	ReportTableInputSpectra.cpp \
	ReportConfig.cpp \
	spsFiles.cpp \
  ReportBase64.cpp \
	RendererBase.cpp \
	RendererGnu.cpp \
	PlotBase.cpp \
	Tokenizer.cpp \
	$(COMMON_FILES)/StatusFile.cpp \
	$(COMMON_FILES)/Specific.cpp

################################################################################
# Targets
################################################################################
all: $(LIBSPSREPORTS_SRC)

CPPREPORTLIB_OBJ=$(CPPREPORTLIB_SRC:.cpp=.o)

$(LIBSPSREPORTS_SRC): $(CPPREPORTLIB_OBJ)
	$(AR) rc $@ $^
	$(RANLIB) $@

$(LIBSPSREPORTS_SRC_DYNAMIC): $(CPPREPORTLIB_OBJ)
	$(CXX) -shared  -Wl,-rpath,$(ROOT) -Wl,-soname,$@.1 -o $@.1.0 $^
	ln -s $@.1.0 $@.1
	ln -s $@.1 $@

################################################################################
# Compile rules - cleanup
################################################################################
clean:
	rm -f *.o *.a depends *.gch *~ *.bak > /dev/null 2>&1 | :

################################################################################
