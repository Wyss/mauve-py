
if DEBUG
D_CXXFLAGS = -Wall -g -DCOMMAND_LINE -D__GNDEBUG__
endif

AM_CXXFLAGS = $(D_CXXFLAGS) @OPENMP_CFLAGS@

LIBGENOME_H = \
gnABISource.h      gnDataBaseSource.h  gnLocation.h         gnSetup.h \
gnBaseFeature.h    gnDebug.h           gnMultiSpec.h        gnSourceFactory.h \
gnBaseFilter.h     gnDefs.h            gnSourceHeader.h     gnPosSpecificTranslator.h \
gnBaseHeader.h     gnFASSource.h       gnSourceQualifier.h \
gnBaseQualifier.h  gnSourceSpec.h      gnStringQualifier.h\
gnBaseSource.h     gnFeature.h         gnProteinSequence.h  gnStringHeader.h \
gnBaseSpec.h       gnFileContig.h      gnRAWSource.h \
gnClone.h          gnFileSource.h      gnRNASequence.h      gnStringSpec.h \
gnCompare.h        gnFilter.h          gnSEQSource.h        gnStringTools.h \
gnContigSpec.h     gnFragmentSpec.h    gnTranslator.h \
gnDNASequence.h    gnGBKSource.h       gnVersion.h   gnFastTranslator.h \
gnGenomeSpec.h      gnSequence.h  gnException.h gnExceptionCode.h \
IntervalSequenceTree.h OmpGuard.h


LIBGENOME_SRC = \
gnFilter.cpp gnCompare.cpp gnTranslator.cpp gnSourceFactory.cpp \
gnFeature.cpp gnSourceQualifier.cpp gnFileContig.cpp gnSourceSpec.cpp \
gnFileSource.cpp gnFASSource.cpp \
gnFragmentSpec.cpp gnStringSpec.cpp gnGBKSource.cpp gnStringTools.cpp \
gnGenomeSpec.cpp  gnLocation.cpp  \
gnRAWSource.cpp gnBaseFeature.cpp gnSEQSource.cpp gnSequence.cpp \
gnContigSpec.cpp gnSourceHeader.cpp gnException.cpp \
gnFastTranslator.cpp gnPosSpecificTranslator.cpp gnDefs.cpp

library_includedir=$(includedir)/$(GENERIC_LIBRARY_NAME)-$(GENERIC_API_VERSION)/$(GENERIC_LIBRARY_NAME)

library_include_HEADERS = $(LIBGENOME_H)

lib_LTLIBRARIES = libGenome-1.3.la
libGenome_1_3_la_SOURCES = $(LIBGENOME_SRC)
libGenome_1_3_la_LDFLAGS= -version-info $(GENERIC_LIBRARY_VERSION) 



INCLUDES = -I$(top_srcdir)

check_PROGRAMS = TestRevComp
TestRevComp_SOURCES = TestRevComp.cpp
TestRevComp_LDADD = $(top_builddir)/src/libGenome.a

EXTRA_PROGRAMS = test_o_matic testSource testSourceFactory testSourceSeq testTests testgnSequence coordMapper
test_o_matic_SOURCES = test-o-matic.cpp
test_o_matic_LDADD = $(top_builddir)/src/libGenome.a
test_o_matic_LDFLAGS = @STATIC_FLAG@
testSource_SOURCES = testSource.cpp
testSource_LDADD = $(top_builddir)/src/libGenome.a
testSourceFactory_SOURCES = testSourceFactory.cpp
testSourceFactory_LDADD = $(top_builddir)/src/libGenome.a
testSourceSeq_SOURCES = testSourceSeq.cpp
testSourceSeq_LDADD = $(top_builddir)/src/libGenome.a
testTests_SOURCES = testTests.cpp
testTests_LDADD = $(top_builddir)/src/libGenome.a
testgnSequence_SOURCES = testgnSequence.cpp
testgnSequence_LDADD = $(top_builddir)/src/libGenome.a
coordMapper_SOURCES = coordMapper.cpp
coordMapper_LDADD = $(top_builddir)/src/libGenome.a
