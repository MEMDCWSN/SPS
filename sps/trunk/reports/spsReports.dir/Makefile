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
# Source files - SPS
################################################################################
SPS_SRC = \
	ReportBase.cpp \
	ReportHeader.cpp \
	ReportProtein.cpp \
	ReportContig.cpp \
  ReportCluster.cpp \
	ReportInputSpectra.cpp \
	ReportRendererBase.cpp \
	ReportRendererHtml.cpp \
	ReportRendererHtmlClient.cpp \
	ReportRendererHtmlSingle.cpp \
	ReportRendererHtmlDynamic.cpp \
	ReportRendererPdf.cpp \
	ReportModuleBase.cpp \
	ReportModuleFactory.cpp \
	ReportModuleSpecplot.cpp \
	ReportModuleContplot.cpp \
	$(COMMON_FILES)/SpectrumAnnotStatistics.cpp \
	$(COMMON_FILES)/SpectrumAnnotParameterList.cpp \
	ReportInterface.cpp \
	spsReport.cpp

SPS_HRD= \
	ReportBase.h \
	ReportHeader.h \
	ReportProtein.h \
	ReportContig.h \
  ReportCluster.h \
	ReportInputSpectra.h \
	ReportRendererBase.h \
	ReportRendererHtml.h \
	ReportRendererHtmlClient.h \
	ReportRendererHtmlSingle.h \
	ReportRendererHtmlDynamic.h \
	ReportRendererPdf.h \
	ReportModuleBase.h \
	ReportModuleFactory.h \
	ReportModuleSpecplot.h \
	ReportModuleContplot.h \
	$(COMMON_FILES)/SpectrumAnnotStatistics.h \
	$(COMMON_FILES)/SpectrumAnnotParameterList.h \
	ReportInterface.h \
	spsReport.h


################################################################################
# Source files - CONTPLOT
################################################################################
CONT_SRC = \
	$(CONTPLOT_DIR)/ContigSequence.cpp \
	$(CONTPLOT_DIR)/ContigSpectrum.cpp \
	$(CONTPLOT_DIR)/PlotContig.cpp \
	$(CONTPLOT_DIR)/ContplotInterface.cpp

CONT_HDR = \
	$(CONTPLOT_DIR)/ContigSequence.h \
	$(CONTPLOT_DIR)/ContigSpectrum.h \
	$(CONTPLOT_DIR)/PlotContig.h \
	$(CONTPLOT_DIR)/contplot.h \
	$(CONTPLOT_DIR)/ContplotInterface.h

################################################################################
# Source files - SPECPLOT
################################################################################
SPEC_SRC = \
	$(SPECPLOT_DIR)/PlotSpectrum.cpp \
	$(SPECPLOT_DIR)/SpecplotInterface.cpp

SPEC_HDR = \
	$(SPECPLOT_DIR)/PlotSpectrum.h \
	$(SPECPLOT_DIR)/SpecplotInterface.h

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
	$(COMMON_FILES)/Specific.cpp \
  $(COMMON_FILES)/StatusFile.cpp

CPPLIB_HDR = \
  $(COMMON_FILES)/CommandLineParser.h \
  $(COMMON_FILES)/ParameterList.h  \
  $(COMMON_FILES)/abruijn.h \
  $(COMMON_FILES)/ms1.h \
	$(COMMON_FILES)/projectionutils.h \
	$(COMMON_FILES)/SpecnetsGraph.h \
	$(COMMON_FILES)/SpectrumPair.h \
	$(COMMON_FILES)/SpectrumPairSet.h \
	$(COMMON_FILES)/Specific.h \
  $(COMMON_FILES)/StatusFile.h

################################################################################
# Compile rules - target
################################################################################
EXECUTABLES = spsReport

all: $(EXECUTABLES)

################################################################################
# Compile rules - specific
################################################################################
SPEC_OBJ=$(SPEC_SRC:.cpp=.o)

SPS_OBJ=$(SPS_SRC:.cpp=.o)

CONT_OBJ=$(CONT_SRC:.cpp=.o)

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

$(EXECUTABLES): $(SPS_OBJ)

spsReport: $(PREFIX)/$(SUFFIX)/libz.a $(PREFIX)/$(SUFFIX)/libpng.a $(PREFIX)/$(SUFFIX)/libhpdf.a $(SPS_OBJ) $(SPEC_OBJ) $(CPPLIB_OBJ) $(CONT_OBJ) $(LIBCORE_SRC) $(LIBSPS_SRC) $(LIBSPSREPORTS_SRC)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $(SPS_OBJ) $(SPEC_OBJ) $(CPPLIB_OBJ) $(CONT_OBJ) $(LIBSPSREPORTS_SRC) $(LIBCORE_SRC) $(LIBSPS_SRC) -o $@ -lhpdf -lpng -lz -lm $(LIBS)

################################################################################
# Compile rules - cleanup
################################################################################
clean:
	rm -f *.o $(EXECUTABLES) depends *.gch *~ *.bak > /dev/null 2>&1 | :

################################################################################



