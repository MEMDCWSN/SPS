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


CPPLIB_SRC = \
	FalseLocalizationRates.cpp \
	LinearEquation.cpp \
	Timer.cpp \
	../cosine_similarity.cpp 

CPPLIB_OBJ=$(CPPLIB_SRC:.cpp=.o)

STATICLIB = libflr.a

DYNAMICLIB = libflr.so


clean: 
	rm -f *.o $(STATICLIB)

$(STATICLIB): $(CPPLIB_OBJ) $(LIBCORE_SRC)
	ar rcs $@ $^


all: $(STATICLIB)

#$(DYNAMICLIB): $(CPPLIB_OBJ) 
#	$(CC) -shared  -Wl,-rpath,$(ROOT) -Wl,-soname,$@.1 -o $@.1.0 $^
#	ln -s $@.1.0 $@.1
#	ln -s $@.1 $@