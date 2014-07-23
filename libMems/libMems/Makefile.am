
if DEBUG
D_CXXFLAGS = -Wall -g -D__GNDEBUG__
endif
OPTIMIZATION = -O2 -funroll-loops -fomit-frame-pointer -ftree-vectorize 
AM_CFLAGS = $(OPTIMIZATION) @DEPS_CFLAGS@ -DUSE_POSIX_AIO @OPENMP_CFLAGS@
AM_CXXFLAGS = $(OPTIMIZATION) @DEPS_CFLAGS@ @BOOST_CPPFLAGS@ $(D_CXXFLAGS) @EXTRA_CXX_FLAGS@ @OPENMP_CXXFLAGS@ 
AM_LDFLAGS = $(OPTIMIZATION)

AM_CXXFLAGS += -DBOOST_FILESYSTEM_VERSION=2

LIBMEMS_H = \
RepeatHash.h      MatchHashEntry.h \
DNAFileSML.h     MemorySML.h         MatchProjectionAdapter.h \
DNAMemorySML.h  MatchFinder.h           SortedMerList.h IntervalList.h \
FileSML.h      gnAlignedSequences.h  Interval.h        \
MemHash.h      AbstractMatch.h    SlotAllocator.h \
Aligner.h   Match.h     MatchList.h Matrix.h NumericMatrix.h \
Islands.h   MaskedMemHash.h   SeedMasks.h GappedAlignment.h \
MuscleInterface.h GappedAligner.h PhyloTree.h SparseAbstractMatch.h \
DenseAbstractMatch.h RepeatMatch.h UngappedLocalAlignment.h \
AbstractGappedAlignment.h CompactGappedAlignment.h HybridAbstractMatch.h \
twister.h SubstitutionMatrix.h RepeatMatchList.h \
Backbone.h ProgressiveAligner.h PairwiseMatchAdapter.h PairwiseMatchFinder.h \
SeedOccurrenceList.h TreeUtilities.h SuperInterval.h GreedyBreakpointElimination.h \
LCB.h DistanceMatrix.h Scoring.h configuration.h Memory.h Files.h gnRAWSequence.h

HOMOLOGYHMM_H = HomologyHMM/homology.h HomologyHMM/dptables.h HomologyHMM/algebras.h HomologyHMM/parameters.h

DMSML_H = \
dmSML/asyncio.h dmSML/alinuxaio.h dmSML/aPOSIXaio.h \
dmSML/alibc.h dmSML/awin32aio.h dmSML/buffer.h \
dmSML/util.h dmSML/sorting.h dmSML/dmsort.h \
dmSML/timing.h dmSML/sml.h

LIBMEMS_SRC = \
RepeatHash.cpp       \
DNAFileSML.cpp       MatchFinder.cpp       \
DNAMemorySML.cpp     MemorySML.cpp        SortedMerList.cpp \
FileSML.cpp          MemHash.cpp          MatchHashEntry.cpp \
Interval.cpp	     IntervalList.cpp     twister.c \
gnAlignedSequences.cpp                     \
MatchList.cpp        Aligner.cpp \
Islands.cpp          MaskedMemHash.cpp    GappedAlignment.cpp \
MuscleInterface.cpp  PhyloTree.cpp         \
RepeatMatchList.cpp  RepeatMatch.cpp \
Backbone.cpp	PairwiseMatchFinder.cpp	ProgressiveAligner.cpp \
SuperInterval.cpp	GreedyBreakpointElimination.cpp

HOMOLOGYHMM_SRC = \
HomologyHMM/algebras.cc HomologyHMM/homology.cc HomologyHMM/homologymain.cc

DMSML_SRC = \
dmSML/asyncio.c dmSML/alinuxaio.c dmSML/aPOSIXaio.c \
dmSML/alibc.c dmSML/awin32aio.c dmSML/buffer.c \
dmSML/util.c dmSML/sorting.c dmSML/dmsort.c \
dmSML/timing.c dmSML/sml.c


libmems_includedir=$(includedir)/$(GENERIC_LIBRARY_NAME)-$(GENERIC_API_VERSION)/$(GENERIC_LIBRARY_NAME)
libmems_include_HEADERS = $(LIBMEMS_H) 

# build libraries with gcc (no suffix)
lib_LTLIBRARIES = libMems-1.6.la
libMems_1_6_la_SOURCES = $(LIBMEMS_SRC) $(HOMOLOGYHMM_SRC) $(DMSML_SRC)

libMems_1_6_la_LDFLAGS= -version-info $(GENERIC_LIBRARY_VERSION)  

homologyhmm_includedir=$(includedir)/$(GENERIC_LIBRARY_NAME)-$(GENERIC_API_VERSION)/$(GENERIC_LIBRARY_NAME)/HomologyHMM
homologyhmm_include_HEADERS = $(HOMOLOGYHMM_H)

dmsml_includedir=$(includedir)/$(GENERIC_LIBRARY_NAME)-$(GENERIC_API_VERSION)/$(GENERIC_LIBRARY_NAME)/dmSML
dmsml_include_HEADERS = $(DMSML_H)

EXTRA_DIST = \
HomologyHMM/homology.xml 

#EXTRA_PROGRAMS = TestSML TestSMLstatic


#TestSMLstatic_SOURCES = TestSML.cpp
#TestSMLstatic_INCLUDES = -I$(top_srcdir)/include/ `wx-config --cxxflags`
#TestSMLstatic_LDFLAGS = -static $(top_builddir)/libMems/libMems.a `wx-config --static --libs`




