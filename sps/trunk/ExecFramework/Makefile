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
# External Libraries
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.libs

################################################################################
# Compile rules - general
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.general

################################################################################
# Libraries
################################################################################
LIBRARIES = -lspecnets -lhpdf -lpng -lz -lm -fopenmp -lxerces-c
#LIBRARIES = -lpthread -lspecnets -lxerces-c

LIBRARIES_SRC = \
	$(COMMON_FILES)/libspecnets.a \
	$(LIBSPSREPORTS_SRC) \
 	$(LIBCORE_SRC) \
 	$(LIBSPS_SRC)

################################################################################
# Source files - CPPLIB
################################################################################
CPPLIB_SRC = \
  AlignmentUtils.cpp \
  ExecAlignment.cpp \
  ExecAssembly.cpp \
  ExecBase.cpp \
  ExecCreateSpectralLibrary.cpp \
  ExecDeconvoluteMS2.cpp \
  ExecFilterAligns.cpp \
  ExecFilterSpectra.cpp \
  ExecFdrPeptide.cpp \
  ExecFilterContigPairs.cpp \
  ExecFilterPairs.cpp \
  ExecFilterStarPairs.cpp \
  ExecHomologyAssembly.cpp \
  ExecIntersectingSpecProb.cpp \
  ExecKLFilter.cpp \
  ExecGenerateTags.cpp \
	ExecClusterBins.cpp \
  ExecMainSpecnets.cpp \
  ExecMergeConvert.cpp \
  ExecMetaAssembly.cpp \
  ExecCCMSMetaSPSSpecnetsParams.cpp \
  ExecModuleFactory.cpp \
  ExecMsCluster.cpp \
  ExecNetworksPropagation.cpp \
  ExecTagFilterPairs.cpp \
  ExecReportSpsplot.cpp \
  ExecReportSPSStats.cpp \
  ExecParallelAssembly.cpp \
  ExecPrmScoring.cpp \
  ExecPrmClustering.cpp \
  ExecProjection.cpp \
  ExecProjectionStatistics.cpp \
  ExecProteosafeToSpecnetsParams.cpp \
  ExecMolecularParallelPairs.cpp \
  ExecQCSpectralPairs.cpp \
  ExecQCSpectrum.cpp \
  ExecQualityControl.cpp \
  ExecSpectralLibrarySearch.cpp \
  ExecSpectralLibrarySearchMolecular.cpp \
  ExecSpectralLibraryDecoyCreation.cpp \
  ExecSpectralLibrarySLGFTraining.cpp \
  FirstCorrectTags.cpp \
  ExecSpectralLibrarySLGFCreation.cpp \
  ExecSpectralLibrarySearchSLGF.cpp \
  ExecSpectraSyntheticSpectraScanExtraction.cpp \
  ExecSpectralLibraryPeakStatistics.cpp \
  ExecCCMSMetabolomicsSpecnetsParams.cpp \
  ExecSpectralLibraryBatchRemove.cpp \
  ExecSpectralProbability.cpp \
  ExecLibraryView.cpp \
  ExecSpecNetsPropagation.cpp \
  ExecSpectraExtraction.cpp \
  ExecSpectraExtractionSingle.cpp \
  ExecSpectraExtractionTable.cpp \
  ExecSpectraUpdateSingle.cpp \
  ExecStarGF.cpp \
  ExecCreateSpectralLibraryMolecular.cpp \
  ExecStatistics.cpp \
  ExecSvmStatistics.cpp \
  ExecTagSearch.cpp \
  ExecTrainMassError.cpp \
  ExecPenaltyGen.cpp \
  ExecSpecProtAlign.cpp \
  ExecSpecProtAlignTgtDecoy.cpp \
  ExecSpecTagGen.cpp \
  ExecProtProtAlign.cpp \
  FileUtils.cpp \
  Filters.cpp \
  ParallelExecution.cpp \
  ParallelPbsExecution.cpp \
  ParallelSgeExecution.cpp \
  ParallelThreadedExecution.cpp \
	$(REPORT_DIR)/ReportBase.cpp \
	$(REPORT_DIR)/ReportHeader.cpp \
	$(REPORT_DIR)/ReportProtein.cpp \
	$(REPORT_DIR)/ReportContig.cpp \
  $(REPORT_DIR)/ReportCluster.cpp \
	$(REPORT_DIR)/ReportInputSpectra.cpp \
	$(REPORT_DIR)/ReportRendererBase.cpp \
	$(REPORT_DIR)/ReportRendererHtml.cpp \
	$(REPORT_DIR)/ReportRendererHtmlClient.cpp \
	$(REPORT_DIR)/ReportRendererHtmlSingle.cpp \
	$(REPORT_DIR)/ReportRendererHtmlDynamic.cpp \
	$(REPORT_DIR)/ReportRendererPdf.cpp \
	$(REPORT_DIR)/ReportModuleBase.cpp \
	$(REPORT_DIR)/ReportModuleFactory.cpp \
	$(REPORT_DIR)/ReportModuleSpecplot.cpp \
	$(REPORT_DIR)/ReportModuleContplot.cpp \
	$(REPORT_DIR)/ReportInterface.cpp \
	$(SPECPLOT_DIR)/PlotSpectrum.cpp \
	$(SPECPLOT_DIR)/SpecplotInterface.cpp \
	$(CONTPLOT_DIR)/ContigSequence.cpp \
	$(CONTPLOT_DIR)/ContigSpectrum.cpp \
	$(CONTPLOT_DIR)/PlotContig.cpp \
	$(CONTPLOT_DIR)/ContplotInterface.cpp \
	$(LIBSPS_DIR)/AssembledPeak.cpp \
   $(LIBSPS_DIR)/AbruijnEdge.cpp \
   $(LIBSPS_DIR)/AbruijnGraph.cpp \
   $(LIBSPS_DIR)/AbruijnNode.cpp \
   $(LIBSPS_DIR)/SpectrumAlignment.cpp \
   $(LIBSPS_DIR)/SpectrumAlignmentSet.cpp \
   $(LIBSPS_DIR)/SpectrumAligner.cpp \
  SgeGridMonitor.cpp \
  $(COMMON_FILES)/Specific.cpp \
  $(COMMON_FILES)/StatusFile.cpp \
      $(COMMON_FILES)/mzxml.cpp \
  ExecGenoMS.cpp
#  SetMerger.cpp \
#  SpectralPairs.cpp \
#  SpectrumPair.cpp \
#  SpectrumPairSet.cpp \
# ExecReportProteinCoverage.cpp \
#   ExecGFNetwork.cpp \
#  ExecSpecNetworkEval.cpp \


CPPLIB_HDR = \
  AlignmentUtils.h \
  ExecAlignment.h \
  ExecAssembly.h \
  ExecBase.h \
  ExecCreateSpectralLibrary.h \
  ExecDeconvoluteMS2.h \
  ExecFilterAligns.h \
  ExecFdrPeptide.h \
  ExecFilterContigPairs.h \
  ExecFilterPairs.h \
  ExecFilterStarPairs.h \
  ExecFilterSpectra.h \
  ExecHomologyAssembly.h \
  ExecIntersectingSpecProb.h \
  ExecKLFilter.h \
  ExecGenerateTags.h \
	ExecClusterBins.h \
  ExecMainSpecnets.h \
  ExecMergeConvert.h \
  ExecMetaAssembly.h \
  ExecModuleFactory.h \
  ExecMsCluster.h \
  ExecNetworksPropagation.h \
  ExecTagFilterPairs.h \
  ExecReportSpsplot.h \
  ExecReportSPSStats.h \
  ExecParallelAssembly.h \
  ExecProjection.h \
  ExecProjectionStatistics.h \
  ExecProteosafeToSpecnetsParams.h \
  ExecQCSpectralPairs.h \
  ExecQCSpectrum.h \
  ExecQualityControl.h \
  ExecSpectralProbability.h \
  ExecSpectralLibrarySearch.h \
  ExecSpectralLibrarySearchMolecular.h \
  ExecSpectralLibraryDecoyCreation.h \
  ExecSpectralLibrarySLGFTraining.h \
  FirstCorrectTags.h \
  ExecSpectralLibrarySLGFCreation.h \
  ExecSpectralLibrarySearchSLGF.h \
  ExecSpectralLibraryBatchRemove.h \
  ExecCCMSMetabolomicsSpecnetsParams.h \
  ExecSpectraSyntheticSpectraScanExtraction.h \
  ExecSpectralLibraryPeakStatistics.h \
  ExecLibraryView.h \
  ExecSpectraExtraction.h \
  ExecSpectraUpdateSingle.h \
  ExecSpectraExtractionSingle.h \
  ExecSpectraExtractionTable.h \
  ExecStarGF.h \
  ExecCreateSpectralLibraryMolecular.h \
  ExecPrmScoring.h \
  ExecPrmClustering.h \
  ExecTagSearch.h \
  ExecTrainMassError.h \
  ExecPenaltyGen.h \
  ExecSpecProtAlign.h \
  ExecSpecProtAlignTgtDecoy.h \
  ExecSpecTagGen.cpp \
  ExecProtProtAlign.h \
  ExecSpecNetsPropagation.h \
  ExecStatistics.h \
  ExecSvmStatistics.h \
  FileUtils.h \
  Filters.h \
  ParallelExecution.h \
  ParallelPbsExecution.h \
  ParallelSgeExecution.h \
  ParallelThreadedExecution.h \
	$(REPORT_DIR)/ReportBase.h \
	$(REPORT_DIR)/ReportHeader.h \
	$(REPORT_DIR)/ReportProtein.h \
	$(REPORT_DIR)/ReportContig.h \
  $(REPORT_DIR)/ReportCluster.h \
	$(REPORT_DIR)/ReportInputSpectra.h \
	$(REPORT_DIR)/ReportRendererBase.h \
	$(REPORT_DIR)/ReportRendererHtml.h \
	$(REPORT_DIR)/ReportRendererHtmlClient.h \
	$(REPORT_DIR)/ReportRendererHtmlSingle.h \
	$(REPORT_DIR)/ReportRendererHtmlDynamic.h \
	$(REPORT_DIR)/ReportRendererPdf.h \
	$(REPORT_DIR)/ReportModuleBase.h \
	$(REPORT_DIR)/ReportModuleFactory.h \
	$(REPORT_DIR)/ReportModuleSpecplot.h \
	$(REPORT_DIR)/ReportModuleContplot.h \
	$(REPORT_DIR)/ReportInterface.h \
	$(SPECPLOT_DIR)/PlotSpectrum.h \
	$(SPECPLOT_DIR)/SpecplotInterface.h \
	$(CONTPLOT_DIR)/ContigSequence.h \
	$(CONTPLOT_DIR)/ContigSpectrum.h \
	$(CONTPLOT_DIR)/PlotContig.h \
	$(CONTPLOT_DIR)/contplot.h \
	$(CONTPLOT_DIR)/ContplotInterface.h \
  SgeGridMonitor.h \
  $(COMMON_FILES)/Specific.h \
  $(COMMON_FILES)/StatusFile.h \
  ExecGenoMS.h
#  SpectralPairs.h \
#  SetMerger.h \
#  SpectrumPair.h \
#  SpectrumPairSet.h \
#  ExecReportProteinCoverage.h \
#  ExecGFNetwork.h \
#  ExecSpecNetworkEval.h \

MAIN_SPECNETS_SRC = \
  main_specnets_helpers.cpp \
  main_specnets_performs.cpp \
  main_specnets_perform_assembly.cpp \
  main_specnets_perform_filteralign.cpp \
  main_specnets_perform_genoms.cpp  \
  main_specnets_perform_specprotalign.cpp

MAIN_SPECNETS_HDR = \
  main_specnets_helpers.h \
  main_specnets_performs.h \
  main_specnets_perform_assembly.h \
  main_specnets_perform_filteralign.h \
  main_specnets_perform_genoms.h \
  main_specnets_perform_specprotalign.h

################################################################################
# Compile rules - cleanup
################################################################################
clean:
	rm -f *.o $(EXECUTABLES) main_specnets main_execmodule main_grid_execution depends *.gch *~ *.bak

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

MAIN_SPECNETS_OBJ=$(MAIN_SPECNETS_SRC:.cpp=.o)

################################################################################
# Compile rules - specific
################################################################################
.PHONY: subdirs $(SUBDIRS) clean all

EXECUTABLES = \
	main_specnets \
	main_execmodule \
	main_grid_execution \
	main_compute_fdr \
	main_diffspecs

ifeq ($(type),debug)
all: $(EXECUTABLES)
else
all: $(EXECUTABLES)
	$(STRIP) $(EXECUTABLES)
endif

main_specnets: main_specnets.o $(MAIN_SPECNETS_OBJ) $(CPPLIB_OBJ) $(LIBRARIES_SRC)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LIBRARIES) $(LIBS)

main_execmodule: main_execmodule.o $(CPPLIB_OBJ) $(LIBRARIES_SRC)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LIBRARIES) $(LIBS)

main_grid_execution: main_grid_execution.o $(CPPLIB_OBJ) $(LIBRARIES_SRC)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LIBRARIES) $(LIBS)

main_compute_fdr: main_compute_fdr.o $(CPPLIB_OBJ) $(LIBRARIES_SRC)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LIBRARIES) $(LIBS)

main_diffspecs: main_diffspecs.o $(LIBRARIES_SRC)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LIBRARIES) $(LIBS)



################################################################################
# Compile rules - convert
################################################################################
# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
#convert: convert.o mzxml.o $(CPPLIB_OBJ)
#	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c-3.0 $(LIBS)

# Windows
else
#convert: convert.o mzxml.o $(CPPLIB_OBJ)
#	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c-3-0 $(LIBS)
endif

################################################################################
# dependencies
################################################################################
depends : $(CPPLIB_SRC)  $(CPPLIB_HDR)
	@$(CXX) $(CPPFLAGS) -MM $^ > depends

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
-include depends
endif

################################################################################
