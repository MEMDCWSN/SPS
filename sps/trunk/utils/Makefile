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
# Libraries
################################################################################
LIBRARIES = -lpthread -lspecnets -lcore -lSpsReports -lxerces-c -lhpdf -lpng -lz -lm

################################################################################
# Targets
################################################################################
EXECUTABLES = \
	main_analyze_align \
	main_compute_variants \
	main_comparepsms \
	main_comparevariants \
	main_convertpsm \
	main_convert_psm_to_proteosafe \
	main_createdecoydb \
	main_create_psm_specific_penalties \
	main_dump_pairs_histogram \
  main_filtercomponents \
  main_filterpsms \
  main_filter_psms_by_network \
	main_find_supported_mods \
	main_findtaglen \
	main_havetag \
	main_makeptmtable \
	main_makedb_from_psms \
	util_parsimony \
  main_removepsmspecs \
  main_split_db \
	main_scorecontig \
	main_specdump \
	main_subtract_psms
	#main_havetag

.PHONY: subdirs $(SUBDIRS) clean all

ifeq ($(type),debug)
all: $(EXECUTABLES)
else
all: $(EXECUTABLES)
	$(STRIP) $(EXECUTABLES)
endif

################################################################################
# Compile rules: Cleanup
################################################################################
clean:
	rm -f *.o $(EXECUTABLES)

################################################################################
# Compile rules: Executables
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

main_analyze_align: main_analyze_align.cpp  $(LIBRARIES_SRC)
	$(CXX) main_analyze_align.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_analyze_align $(LIBS) $(PWIZ_LIBS)

main_compute_variants: main_compute_variants.cpp  $(LIBRARIES_SRC)
	$(CXX) main_compute_variants.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_compute_variants $(LIBS) $(PWIZ_LIBS)

main_comparepsms: main_comparepsms.cpp  $(LIBRARIES_SRC)
	$(CXX) main_comparepsms.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_comparepsms $(LIBS) $(PWIZ_LIBS)

main_comparevariants: main_comparevariants.cpp  $(LIBRARIES_SRC)
	$(CXX) main_comparevariants.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_comparevariants $(LIBS) $(PWIZ_LIBS)
	
main_convertpsm: main_convertpsm.cpp  $(LIBRARIES_SRC)
	$(CXX) main_convertpsm.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_convertpsm $(LIBS) $(PWIZ_LIBS)

main_convert_psm_to_proteosafe: main_convert_psm_to_proteosafe.cpp  $(LIBRARIES_SRC)
	$(CXX) main_convert_psm_to_proteosafe.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_convert_psm_to_proteosafe $(LIBS) $(PWIZ_LIBS)

main_createdecoydb: main_createdecoydb.cpp  $(LIBRARIES_SRC)
	$(CXX) main_createdecoydb.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_createdecoydb $(LIBS) $(PWIZ_LIBS)

main_create_psm_specific_penalties: main_create_psm_specific_penalties.cpp  $(LIBRARIES_SRC)
	$(CXX) main_create_psm_specific_penalties.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_create_psm_specific_penalties $(LIBS) $(PWIZ_LIBS)

main_dump_pairs_histogram: main_dump_pairs_histogram.cpp  $(LIBRARIES_SRC)
	$(CXX) main_dump_pairs_histogram.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_dump_pairs_histogram $(LIBS)

main_filtercomponents: main_filtercomponents.cpp  $(LIBRARIES_SRC)
	$(CXX) main_filtercomponents.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_filtercomponents $(LIBS) $(PWIZ_LIBS)

main_filterpsms: main_filterpsms.cpp  $(LIBRARIES_SRC)
	$(CXX) main_filterpsms.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_filterpsms $(LIBS) $(PWIZ_LIBS)
	
main_filter_psms_by_network: main_filterpsms.cpp  $(LIBRARIES_SRC)
	$(CXX) main_filter_psms_by_network.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_filter_psms_by_network $(LIBS) $(PWIZ_LIBS)
	
main_find_supported_mods: main_find_supported_mods.cpp  $(LIBRARIES_SRC)
	$(CXX) main_find_supported_mods.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_find_supported_mods $(LIBS) $(PWIZ_LIBS)	

main_findtaglen: main_findtaglen.cpp  $(LIBRARIES_SRC)
	$(CXX) main_findtaglen.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_findtaglen $(LIBS) $(PWIZ_LIBS)

main_havetag: main_havetag.cpp  $(LIBRARIES_SRC)
	$(CXX) main_havetag.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_havetag $(LIBS) $(PWIZ_LIBS)

main_makedb_from_psms: main_makedb_from_psms.cpp  $(LIBRARIES_SRC)
	$(CXX) main_makedb_from_psms.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_makedb_from_psms $(LIBS) $(PWIZ_LIBS)

main_makeptmtable: main_makeptmtable.cpp  $(LIBRARIES_SRC)
	$(CXX) main_makeptmtable.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_makeptmtable $(LIBS) $(PWIZ_LIBS)

main_removepsmspecs: main_removepsmspecs.cpp  $(LIBRARIES_SRC)
	$(CXX) main_removepsmspecs.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_removepsmspecs $(LIBS) $(PWIZ_LIBS)

main_scorecontig: main_scorecontig.cpp  $(LIBRARIES_SRC)
	$(CXX) main_scorecontig.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_scorecontig $(LIBS) $(PWIZ_LIBS)

main_specdump: main_specdump.cpp  $(LIBRARIES_SRC)
	$(CXX) main_specdump.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_specdump $(LIBS) $(PWIZ_LIBS)

main_subtract_psms: main_subtract_psms.cpp  $(LIBRARIES_SRC)
	$(CXX) main_subtract_psms.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_subtract_psms $(LIBS) $(PWIZ_LIBS)

main_split_db: main_split_db.cpp  $(LIBRARIES_SRC)
	$(CXX) main_split_db.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o main_split_db $(LIBS) $(PWIZ_LIBS)

util_parsimony: util_parsimony.cpp  $(LIBRARIES_SRC)
	$(CXX) util_parsimony.cpp $(CPPFLAGS) $(LDFLAGS) $(LIBRARIES) -o util_parsimony $(LIBS) $(PWIZ_LIBS)

################################################################################
# Compile rules: Dependencies
################################################################################
# Linux
ifneq (,$(findstring linux,$(shell gcc -dumpmachine)))
-include depends
endif

################################################################################
