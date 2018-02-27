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
# Source files - CPPLIB
################################################################################
CPPLIB_SRC = \
	aminoacid.cpp \
	alignment_scoring.cpp \
	BaseGraph.cpp \
	Cluster.cpp \
	ClusterSet.cpp \
	db_fasta.cpp \
	DeconvSpectrum.cpp \
	DelimitedTextReader.cpp \
	Edge.cpp \
	FilterableGraph.cpp \
	GFTable.cpp \
	GFTableIntersection.cpp \
	inputParams.cpp \
	IsoEnvelope.cpp \
	IonMass.cpp \
	label.cpp \
	Logger.cpp \
	mzrange.cpp \
	MassErrorModel.cpp \
	MathUtils.cpp \
	Node.cpp \
	OutputTable.cpp \
	PeptideSpectrumMatch.cpp \
	PeptideSpectrumMatchSet.cpp \
	SpecSet.cpp \
	spectrum_scoring.cpp \
	spectrum.cpp \
	base64.cpp \
	utils.cpp

################################################################################
# Compile rules: specific
################################################################################
CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

################################################################################
# Target
################################################################################
all: $(LIBCORE_SRC)

$(LIBCORE_SRC): $(CPPLIB_OBJ)
	$(AR) rc $@ $^
	$(RANLIB) $@

$(LIBCORE_SRC_DYNAMIC): $(CPPLIB_OBJ)
	$(CXX) -shared  -Wl,-rpath,$(ROOT) -Wl,-soname,$@.1 -o $@.1.0 $^
	ln -s $@.1.0 $@.1
	ln -s $@.1 $@

################################################################################
# Compile rules: cleanup
################################################################################
clean:
	rm -f *.o *.a depends *~ *.bak

################################################################################
