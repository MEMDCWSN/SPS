################################################################################
# Definitions
################################################################################
BASE_DIR := ..
MASTER := YES
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
# Gnuplot
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.gnuplot

################################################################################
# Compile rules - general
################################################################################
include $(BASE_DIR)/Makefiles/Makefile.general

################################################################################
# Source files - CPPLIB
################################################################################
CPPLIB_SRC = \
	abruijn.cpp \
	alignment_modmut.cpp \
	AlignmentPenaltyBased.cpp \
	CommandLineParser.cpp \
	clusters.cpp \
	cosine_similarity.cpp \
	denovo.cpp \
	FdrPeptide.cpp \
	filters.cpp \
	graph.cpp \
	inspect_parse.cpp \
	MappedContig.cpp \
	MappedContigStatTable.cpp \
   MappedContigSetTable.cpp \
	MappedSpecnets.cpp \
	MappedSpectrum.cpp \
	MappedSPSStatTable.cpp \
	ms1.cpp \
	msn.cpp \
	projectionutils.cpp \
	ParameterList.cpp \
	PenaltyMatrix.cpp \
	PepnovoTags.cpp \
	NetworkedPSM.cpp \
	PeptideSpectrumMatchNetwork.cpp \
	SetMerger.cpp \
	SpectralPairs.cpp \
	spectrum_window_filter.cpp \
	SpectrumPair.cpp \
	SpectrumAnnotParameterList.cpp \
	SpectrumAnnotStatistics.cpp \
	SpectrumLibrary.cpp \
	SpectrumPairSet.cpp \
	SpecnetsGraph.cpp \
	SvmScaleParameterList.cpp \
	SpectralLibrary.cpp \
	SvmModel.cpp \
	tags.cpp \
	ion.cpp \
#	dekel_align.cpp \
#	hash.cpp \
#	batch.cpp \

################################################################################
# Source files - META_SPS
################################################################################
META_SPS_SRC = MetaSPS/prm_alignment.cpp \
		MetaSPS/PRMSpecEdge.cpp \
		MetaSPS/PairedSpecSet.cpp \
		MetaSPS/CombineContigs.cpp \
		MetaSPS/Contig.cpp \
		MetaSPS/ContigNetwork.cpp

################################################################################
# Libraries - specnets
################################################################################
STATICLIB = libspecnets.a

DYNAMICLIB = libspecnets.so

################################################################################
# Executables & scripts
################################################################################
EXECUTABLES = \
	convert
#	MetaSPS/makeClustersFromAbinfo \
#	mergeConvert \
#	filterpairs \
#	starden \
#	filterstarpairs \
#	tagsearch \
#	protid \
#	matchma \
#	mergefilter \
#	homglue \
#	pathproj \
#	starden2 \

SCRIPTS = \
	run_csps.sh \
	run_specnets.sh \
	grid.sh \
	grid_csps_start.sh \
	grid_sps_resume.sh \
	grid_sps_start.sh

################################################################################
# Compile rules - install
################################################################################
install: all
	mkdir -p $(PREFIX)/bin
	mkdir -p $(PREFIX)/share/gnuplot/4.2/PostScript
	cp -f $(EXECUTABLES) $(SCRIPTS) $(PREFIX)/bin
	cp -f $(GNUPLOT_ROOT)/term/PostScript/prologue.ps $(PREFIX)/share/gnuplot/4.2/PostScript
ifneq (,$(shell which mcc 2>/dev/null))
	cp -f csps specnets $(PREFIX)/bin
endif

################################################################################
# Compile rules - cleanup
################################################################################
clean:
	$(info *****************************************************)
	$(info ******************** cleaning up ********************)
	$(info *****************************************************)
	@yes | rm -rf *.o *.a MetaSPS/*.o $(BIN_DIR)/* PepNovo_bin convert csps specnets MsCluster_bin pspectrum pcontig depends *.gch *~ *.bak; cd libSPS; make clean $(MAKE_OUTPUT); cd .. ; cd libCore; make clean $(MAKE_OUTPUT) ; cd ../reports; make clean $(MAKE_OUTPUT); cd ../ExecFramework;make clean $(MAKE_OUTPUT) ;cd .. ; cd utils ; make clean $(MAKE_OUTPUT) ; cd ..

distclean: clean
	$(info ***************************************************************)
	$(info ******************** cleaning up libraries ********************)
	$(info ***************************************************************)
	cd $(SRC_DIR); ./clean.sh ; cd $(ROOT)
	cd $(PREFIX)/$(SUFFIX) ; rm *.a $(MAKE_OUTPUT) ; cd $(ROOT)

#	$(MAKE) -C $(XERCES_ROOT) distclean
#	$(MAKE) -C $(GNUPLOT_ROOT) distclean
#	$(MAKE) -C $(PEPNOVO_DIR) clean
#	$(MAKE) -C $(MSCLUSTER_DIR) clean
#	$(MAKE) -C $(CLUSTAL_W_DIR) clean
#	$(MAKE) -C $(PDFHARU_ROOT) distclean
#	$(MAKE) -C $(LIBZ_ROOT) distclean
#	$(MAKE) -C $(LIBPNG_ROOT) distclean
#	cd $(SRC_DIR); ./clean.sh ; cd $(ROOT)
#	cd $(PREFIX)/$(SUFFIX) ; rm *.a ; cd $(ROOT)

################################################################################
# Compile rules - specific
################################################################################
CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

META_SPS_OBJ=$(META_SPS_SRC:.cpp=.o)

.PHONY: subdirs $(SUBDIRS) clean all


ifeq ($(type),debug)
all: $(SUBDIRS) ExternalLibs clustal_wTarget pepnovoTarget msclusterTarget gnuplotExe libSPSTarget libCoreTarget $(STATICLIB) ReportsTarget ExecTarget utilsTarget qcTarget $(EXECUTABLES)
	cp $(CLUSTAL_W_DIR)/clustalw $(BIN_DIR)
	cp $(MSCLUSTER_DIR)/MSCluster_bin $(BIN_DIR)/MsCluster_bin
	cp $(PEPNOVO_DIR)/PepNovo_bin $(BIN_variaDIR)
	cp $(EXEC_DIR)/main_specnets $(BIN_DIR)
	cp $(EXEC_DIR)/main_execmodule $(BIN_DIR)
	cp $(EXEC_DIR)/main_diffspecs $(BIN_DIR)
	cp $(EXEC_DIR)/main_grid_execution $(BIN_DIR)
	cp $(EXEC_DIR)/main_compute_fdr $(BIN_DIR)
	cp $(UTILS_DIR)/main_analyze_align $(BIN_DIR)
	cp $(UTILS_DIR)/main_compute_variants $(BIN_DIR)
	cp $(UTILS_DIR)/main_comparepsms $(BIN_DIR)
	cp $(UTILS_DIR)/main_convertpsm $(BIN_DIR)
	cp $(UTILS_DIR)/main_findtaglen $(BIN_DIR)
	cp $(UTILS_DIR)/main_makeptmtable $(BIN_DIR)
	cp $(UTILS_DIR)/main_scorecontig $(BIN_DIR)
	cp $(UTILS_DIR)/main_specdump $(BIN_DIR)
	cp $(UTILS_DIR)/util_parsimony $(BIN_DIR)
	cp $(REPORTS_DIR)/spsReports.dir/spsReport $(BIN_DIR)
	cp $(REPORTS_DIR)/spsReportServer.dir/spsReportServer $(BIN_DIR)
	cp $(REPORTS_DIR)/specplot.dir/specplot $(BIN_DIR)
	cp $(REPORTS_DIR)/contplot.dir/contplot $(BIN_DIR)
	cp $(REPORTS_DIR)/uu64Encode.dir/uu64Encode $(BIN_DIR)
	cp $(REPORTS_DIR)/clusterInfo.dir/clusterinfo $(BIN_DIR)
	cp $(REPORTS_DIR)/pairsinfo.dir/pairsinfo $(BIN_DIR)
	cp $(REPORTS_DIR)/VeraBd.ttf $(BIN_DIR)
	cp $(REPORTS_DIR)/prologue.ps $(BIN_DIR)
	cp ./convert $(BIN_DIR)
	cp ./*.txt $(BIN_DIR)
	cp ./resources/*.txt $(BIN_DIR)
	cp -r -f ./resources $(BIN_DIR)
	cp -r -f lib $(BIN_DIR)
	cp -r -f libs $(BIN_DIR)
else
all: $(SUBDIRS) ExternalLibs clustal_wTarget pepnovoTarget msclusterTarget gnuplotExe libSPSTarget libCoreTarget $(STATICLIB) ReportsTarget ExecTarget utilsTarget qcTarget $(EXECUTABLES)
	cp $(CLUSTAL_W_DIR)/clustalw $(BIN_DIR)
	cp $(MSCLUSTER_DIR)/MSCluster_bin $(BIN_DIR)/MsCluster_bin
	cp $(PEPNOVO_DIR)/PepNovo_bin $(BIN_DIR)
	cp $(EXEC_DIR)/main_specnets $(BIN_DIR)
	cp $(EXEC_DIR)/main_execmodule $(BIN_DIR)
	cp $(EXEC_DIR)/main_diffspecs $(BIN_DIR)
	cp $(EXEC_DIR)/main_grid_execution $(BIN_DIR)
	cp $(EXEC_DIR)/main_compute_fdr $(BIN_DIR)
	cp $(UTILS_DIR)/main_analyze_align $(BIN_DIR)
	cp $(UTILS_DIR)/main_compute_variants $(BIN_DIR)
	cp $(UTILS_DIR)/main_comparepsms $(BIN_DIR)
	cp $(UTILS_DIR)/main_convertpsm $(BIN_DIR)
	cp $(UTILS_DIR)/main_findtaglen $(BIN_DIR)
	cp $(UTILS_DIR)/main_makeptmtable $(BIN_DIR)
	cp $(UTILS_DIR)/main_scorecontig $(BIN_DIR)
	cp $(UTILS_DIR)/main_specdump $(BIN_DIR)
	cp $(UTILS_DIR)/util_parsimony $(BIN_DIR)
	cp $(REPORTS_DIR)/spsReports.dir/spsReport $(BIN_DIR)
	cp $(REPORTS_DIR)/spsReportServer.dir/spsReportServer $(BIN_DIR)
	cp $(REPORTS_DIR)/specplot.dir/specplot $(BIN_DIR)
	cp $(REPORTS_DIR)/contplot.dir/contplot $(BIN_DIR)
	cp $(REPORTS_DIR)/uu64Encode.dir/uu64Encode $(BIN_DIR)
	cp $(REPORTS_DIR)/clusterInfo.dir/clusterinfo $(BIN_DIR)
	cp $(REPORTS_DIR)/pairsinfo.dir/pairsinfo $(BIN_DIR)
	cp $(REPORTS_DIR)/VeraBd.ttf $(BIN_DIR)
	cp $(REPORTS_DIR)/prologue.ps $(BIN_DIR)
	cp ./convert $(BIN_DIR)
	cp ./*.txt $(BIN_DIR)
	cp ./resources/*.txt $(BIN_DIR)
	cp -r -f ./resources $(BIN_DIR)
	cp -r -f lib $(BIN_DIR)
	cp -r -f libs $(BIN_DIR)
	#rsync -aq --exclude=*/.svn* ./resources $(BIN_DIR)
	$(STRIP) $(BIN_DIR)/clustalw $(BIN_DIR)/MsCluster_bin $(BIN_DIR)/main_specnets $(BIN_DIR)/main_execmodule $(BIN_DIR)/main_diffspecs $(BIN_DIR)/main_grid_execution $(BIN_DIR)/main_analyze_align $(BIN_DIR)/main_comparepsms $(BIN_DIR)/main_convertpsm $(BIN_DIR)/main_findtaglen $(BIN_DIR)/main_scorecontig $(BIN_DIR)/main_specdump $(BIN_DIR)/util_parsimony $(BIN_DIR)/spsReport $(BIN_DIR)/spsReportServer $(BIN_DIR)/specplot $(BIN_DIR)/contplot $(BIN_DIR)/uu64Encode $(BIN_DIR)/clusterinfo $(BIN_DIR)/pairsinfo $(BIN_DIR)/convert
endif


#subdirs: $(LIBCORE_SRC) $(STATICLIB) $(SUBDIRS)

SUBDIRS =

$(SUBDIRS):
	$(MAKE) -C $@

################################################################################
# Compile rules - libraries
################################################################################
ExecTarget: $(STATICLIB)
	$(info )
	$(info ****************************************************************)
	$(info ******************** Building ExecFramework ********************)
	$(info ****************************************************************)
	$(info )
	$(MAKE) -j4 -C $(EXEC_DIR) $(MAKE_OUTPUT)

ReportsTarget: $(STATICLIB)
	$(info )
	$(info **********************************************************)
	$(info ******************** Building Reports ********************)
	$(info **********************************************************)
	$(info )
	$(MAKE) -C $(REPORTS_DIR) $(MAKE_OUTPUT)

libSPSTarget:
	$(info )
	$(info **************************************************************)
	$(info ******************** Building SPS library ********************)
	$(info **************************************************************)
	$(info )
	$(MAKE) -j4 -C $(LIBSPS_DIR) $(MAKE_OUTPUT)

libCoreTarget:
	$(info )
	$(info ***************************************************************)
	$(info ******************** Building Core library ********************)
	$(info ***************************************************************)
	$(info )
	$(MAKE) -j4 -C $(LIBCORE_DIR) $(MAKE_OUTPUT)
	$(info ***************************************************************)


$(STATICLIB): $(CPPLIB_OBJ) $(META_SPS_OBJ) $(LIBCORE_SRC)
	$(AR) rc $@ $^
	$(RANLIB) $@

$(LIBSPS_SRC):
	$(MAKE) -C $(LIBSPS_DIR) $(MAKE_OUTPUT)

$(LIBCORE_SRC):
	$(MAKE) -C $(LIBCORE_DIR) $(MAKE_OUTPUT)


$(DYNAMICLIB): $(CPPLIB_OBJ)
	$(CXX) -shared  -Wl,-rpath,$(ROOT) -Wl,-soname,$@.1 -o $@.1.0 $^
	ln -s $@.1.0 $@.1
	ln -s $@.1 $@

################################################################################
# Compile rules - external tools
################################################################################
pepnovoTarget:
	$(info )
	$(info **********************************************************)
	$(info ******************** Building PepNovo ********************)
	$(info **********************************************************)
	$(info )
	$(MAKE) -C $(PEPNOVO_DIR) $(MAKE_OUTPUT)

msclusterTarget:
	$(info )
	$(info ************************************************************)
	$(info ******************** Building MsCluster ********************)
	$(info ************************************************************)
	$(info )
	$(MAKE) -C $(MSCLUSTER_DIR) $(MAKE_OUTPUT)

clustal_wTarget:
	$(info )
	$(info ***********************************************************)
	$(info ******************** Building ClustalW ********************)
	$(info ***********************************************************)
	$(info )
	$(MAKE) -C $(CLUSTAL_W_DIR) $(MAKE_OUTPUT)

################################################################################
# Compile rules - utils directory
################################################################################
utilsTarget:
	$(info )
	$(info ********************************************************)
	$(info ******************** Building utils ********************)
	$(info ********************************************************)
	$(info )
	$(MAKE) -j4 -C $(UTILS_DIR) $(MAKE_OUTPUT)

################################################################################
# Compile rules - qc directory
################################################################################
qcTarget:
	$(MAKE) -j4 -C $(QC_TOOLS) $(MAKE_OUTPUT)


################################################################################
# Compile rules - test tools
################################################################################
test_cs: test_cs.cpp $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(CPPFLAGS) $^ -o $@ $(LIBS) $(PWIZ_LIBS)

test: test.cpp $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS) $(PWIZ_LIBS) -pthread

test_nb: test_nb.cpp $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS) $(PWIZ_LIBS) -pthread

convert: convert2.o mzxml.o PWizInterface.o Specific2.o $(CPPLIB_OBJ) $(LIBCORE_SRC) $(LIBSPS_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -static -lxerces-c -lcore $(LIBS) $(PWIZ_LIBS) -pthread  -lgomp $(MAKE_OUTPUT)

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))

#convert: convert.o mzxml.o $(CPPLIB_OBJ) $(LIBCORE_SRC) $(LIBSPS_SRC)
#	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS)

#convert: convert2.o mzxml.o PWizInterface.o $(CPPLIB_OBJ) $(LIBCORE_SRC) $(LIBSPS_SRC)
#	$(CXX) $(LDFLAGS) $^ -o $@ -static -lxerces-c $(LIBS) $(PWIZ_LIBS)

MetaSPS/makeClustersFromAbinfo: MetaSPS/makeClustersFromAbinfo.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS)

mergeConvert: mergeConvert.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS)

# Windows
else
#convert: convert.o libCore/spectrum.o libCore/IsoEnvelope.o mzxml.o Base64.o alignment_scoring.o libCore/label.o libCore/Logger.o ms1.o libCore/utils.o libCore/aminoacid.o libCore/PeptideSpectrumMatch.o libCore/PeptideSpectrumMatchSet.o libCore/DelimitedTextReader.o libCore/SpecSet.o libCore/spectrum_scoring.o libCore/inputParams.o libCore/mzrange.o denovo.o clusters.o PepnovoTags.o NetworkedPSM.o abruijn.o graph.o SpectrumPair.o SetMerger.o SpectrumPairSet.o $(PREFIX)/bin/mingw/libxerces-c.a
#	$(CC) $(LDFLAGS) $^ -o $@ -l ws2_32 -lxerces-c
pklbin_to_mgf: pklbin_to_mgf.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS)

#convert: convert.o $(CPPLIB_OBJ) $(LIBCORE_SRC)  $(LIBSPS_SRC)
#	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS)

mergeConvert: mergeConvert.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS)
endif

filterpairs: filterpairs.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@  $(LIBS)

starden: starden.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@  $(LIBS)

filterstarpairs: filterstarpairs.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@  $(LIBS)

tagsearch: tagsearch.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

pathproj: pathproj.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

findconnected: findconnected.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

msmod: msmod.o FloatMods.o FloatMods.h $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

#-------------------------------------------------------------------------------
# SPS modules
masab: masab.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

#-------------------------------------------------------------------------------
# CSPS modules
protid: protid.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)
matchma: matchma.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

matchme: matchme.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

mergefilter: mergefilter.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@  $(LIBS)

homglue: homglue.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

#-------------------------------------------------------------------------------
# NRP-sequencing modules
cycseq: cycseq.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

#-------------------------------------------------------------------------------
#misc
cosine_similarity: cosine_similarity_test.o cosine_similarity.o $(CPPLIB_OBJ) $(LIBCORE_SRC)
	$(CXX) $(LDFLAGS) $^ -o $@ -lxerces-c $(LIBS)


################################################################################
# dependencies
################################################################################
depends : $(CPPLIB_SRC) spsstats.cpp dbg.cpp
	$(CXX) $(CPPFLAGS) -MM $^ > depends

# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
-include depends
endif

################################################################################
