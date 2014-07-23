
if DEBUG
D_CXXFLAGS = -Wall -g -DCOMMAND_LINE -D__GNDEBUG__ 
else
P_CXXFLAGS = -DCOMMAND_LINE
endif
OPTIMIZATION = -O2 -funroll-loops -fomit-frame-pointer -ftree-vectorize 
AM_CXXFLAGS = $(OPTIMIZATION) $(D_CXXFLAGS) $(P_CXXFLAGS) 
AM_LDFLAGS = @STATIC_FLAG@
LIBRARY_CL = $(DEPS_LIBS) 
#AM_LDADD = $(DEPS_LIBS) 
INCLUDES = @DEPS_CFLAGS@ 

bin_PROGRAMS = mauveAligner mauveStatic scoreAlignment \
uniqueMerCount toRawSequence \
mfa2xmfa addUnalignedIntervals \
toMultiFastA getAlignmentWindows uniquifyTrees \
toGrimmFormat mauveToXMFA \
stripGapColumns progressiveMauve progressiveMauveStatic \
extractBCITrees createBackboneMFA \
repeatoire alignmentProjector stripSubsetLCBs \
projectAndStrip makeBadgerMatrix randomGeneSample getOrthologList \
bbFilter bbAnalyze backbone_global_to_local xmfa2maf coordinateTranslate

EXTRA_PROGRAMS = bbBreakOnGenes mauveMpatrol mauveEfence toGBKsequence \
multiToRawSequence unalign makeMc4Matrix multiEVD evd scoreALU \
calculateBackboneCoverage2 sortContigs countInPlaceInversions gappiness \
joinAlignmentFiles extractBackbone2 pairCompare \
calculateCoverage calculateBackboneCoverage extractBackbone transposeCoordinates

mauveAligner_SOURCES = mauveAligner.cpp mauveAligner.h
mauveAligner_LDFLAGS = $(OPTIMIZATION)  
mauveAligner_LDADD = $(DEPS_LIBS) 
EXTRA_mauveAligner_SOURCES = getopt.c getopt.h getopt1.c
mauveAligner_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

mauveMpatrol_SOURCES = mauveAligner.cpp
mauveMpatrol_LDFLAGS = -lmpatrol -lbfd -liberty $(LIBRARY_CL) `wx-config --libs`
EXTRA_mauveMpatrol_SOURCES = getopt.c getopt.h getopt1.c
mauveMpatrol_DEPENDENCIES = @GETOPT_LONG_SYSTEM@



mauveEfence_SOURCES = mauveAligner.cpp
mauveEfence_LDADD = -lefence $(DEPS_LIBS)
EXTRA_mauveEfence_SOURCES = getopt.c getopt.h getopt1.c
mauveEfence_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

mauveStatic_SOURCES = mauveAligner.cpp
mauveStatic_LDADD = $(LIBRARY_CL)
EXTRA_mauveStatic_SOURCES = getopt.c getopt.h getopt1.c
mauveStatic_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

#mauveAligner4_SOURCES = mauveAligner.cpp
#mauveAligner4_LDADD = $(LIBRARY_CL)
#mauveAligner4_CXXFLAGS = MAX_SEQ_COUNT=4
#EXTRA_mauveAligner4_SOURCES = getopt.c getopt.h getopt1.c
#mauveAligner4_DEPENDENCIES = @GETOPT_LONG_SYSTEM@


calculateCoverage_SOURCES = calculateCoverage.cpp
calculateCoverage_LDADD = $(LIBRARY_CL)
EXTRA_calculateCoverage_SOURCES = getopt.c getopt.h getopt1.c
calculateCoverage_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

scoreAlignment_SOURCES = scoreAlignment.cpp
scoreAlignment_LDADD = $(LIBRARY_CL)
EXTRA_scoreAlignment_SOURCES = getopt.c getopt.h getopt1.c
scoreAlignment_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

extractBackbone_SOURCES = extractBackbone.cpp
extractBackbone_LDADD = $(LIBRARY_CL)
EXTRA_extractBackbone_SOURCES = getopt.c getopt.h getopt1.c
extractBackbone_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

toRawSequence_SOURCES = toRawSequence.cpp
toRawSequence_LDADD = $(LIBRARY_CL)
EXTRA_toRawSequence_SOURCES = getopt.c getopt.h getopt1.c
toRawSequence_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

transposeCoordinates_SOURCES = transposeCoordinates.cpp
transposeCoordinates_LDADD = $(LIBRARY_CL)
EXTRA_transposeCoordinates_SOURCES = getopt.c getopt.h getopt1.c
transposeCoordinates_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

uniqueMerCount_SOURCES = uniqueMerCount.cpp
uniqueMerCount_LDADD = $(LIBRARY_CL)
EXTRA_uniqueMerCount_SOURCES = getopt.c getopt.h getopt1.c
uniqueMerCount_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

calculateBackboneCoverage_SOURCES = calculateBackboneCoverage.cpp
calculateBackboneCoverage_LDADD = $(LIBRARY_CL)
EXTRA_calculateBackboneCoverage_SOURCES = getopt.c getopt.h getopt1.c
calculateBackboneCoverage_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

mfa2xmfa_SOURCES = mfa2xmfa.cpp
mfa2xmfa_LDADD = $(LIBRARY_CL)
EXTRA_mfa2xmfa_SOURCES = getopt.c getopt.h getopt1.c
mfa2xmfa_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

xmfa2maf_SOURCES = xmfa2maf.cpp
xmfa2maf_LDADD = $(LIBRARY_CL)
EXTRA_xmfa2maf_SOURCES = getopt.c getopt.h getopt1.c
xmfa2maf_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

toGBKsequence_SOURCES = toGBKsequence.cpp
toGBKsequence_LDADD = $(LIBRARY_CL)
EXTRA_toGBKsequence_SOURCES = getopt.c getopt.h getopt1.c
toGBKsequence_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

multiToRawSequence_SOURCES = multiToRawSequence.cpp
multiToRawSequence_LDADD = $(LIBRARY_CL)
EXTRA_multiToRawSequence_SOURCES = getopt.c getopt.h getopt1.c
multiToRawSequence_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

unalign_SOURCES = unalign.cpp
unalign_LDADD = $(LIBRARY_CL)
EXTRA_unalign_SOURCES = getopt.c getopt.h getopt1.c
unalign_DEPENDENCIES = @GETOPT_LONG_SYSTEM@

addUnalignedIntervals_SOURCES = addUnalignedIntervals.cpp
addUnalignedIntervals_LDADD = $(LIBRARY_CL)

toMultiFastA_SOURCES = toMultiFastA.cpp
toMultiFastA_LDADD = $(LIBRARY_CL)

getAlignmentWindows_SOURCES = getAlignmentWindows.cpp
getAlignmentWindows_LDADD = $(LIBRARY_CL)

extractBackbone2_SOURCES = extractBackbone2.cpp
extractBackbone2_LDADD = $(LIBRARY_CL)

uniquifyTrees_SOURCES = uniquifyTrees.cpp
uniquifyTrees_LDADD = $(LIBRARY_CL)


countInPlaceInversions_SOURCES = countInPlaceInversions.cpp
countInPlaceInversions_LDADD = $(LIBRARY_CL)

toGrimmFormat_SOURCES = toGrimmFormat.cpp
toGrimmFormat_LDADD = $(LIBRARY_CL)

joinAlignmentFiles_SOURCES = joinAlignmentFiles.cpp
joinAlignmentFiles_LDADD = $(LIBRARY_CL)

mauveToXMFA_SOURCES = mauveToXMFA.cpp
mauveToXMFA_LDADD = $(LIBRARY_CL)

stripGapColumns_SOURCES = stripGapColumns.cpp
stripGapColumns_LDADD = $(LIBRARY_CL)

gappiness_SOURCES = gappiness.cpp
gappiness_LDADD = $(LIBRARY_CL)


progressiveMauve_SOURCES = progressiveMauve.cpp UniqueMatchFinder.h UniqueMatchFinder.cpp 
progressiveMauve_LDFLAGS = 
progressiveMauve_LDADD = $(LIBRARY_CL)

progressiveMauveStatic_SOURCES = progressiveMauve.cpp UniqueMatchFinder.h UniqueMatchFinder.cpp 
progressiveMauveStatic_LDADD = $(LIBRARY_CL)

sortContigs_SOURCES = sortContigs.cpp
sortContigs_LDADD = $(LIBRARY_CL)

extractBCITrees_SOURCES = extractBCITrees.cpp
extractBCITrees_LDADD = $(LIBRARY_CL)

calculateBackboneCoverage2_SOURCES = calculateBackboneCoverage2.cpp
calculateBackboneCoverage2_LDADD = $(LIBRARY_CL)

createBackboneMFA_SOURCES = createBackboneMFA.cpp
createBackboneMFA_LDADD = $(LIBRARY_CL)

pairCompare_SOURCES = pairCompare.cpp
pairCompare_LDADD = $(LIBRARY_CL)

repeatoire_SOURCES = repeatoire.cpp MatchRecord.h SeedMatchEnumerator.h
repeatoire_LDADD = $(LIBRARY_CL) 

scoreALU_SOURCES = scoreALU.cpp
scoreALU_LDADD = $(LIBRARY_CL) 

evd_SOURCES = evd.cpp
evd_LDADD = $(LIBRARY_CL)

alignmentProjector_SOURCES = alignmentProjector.cpp
alignmentProjector_LDADD = $(LIBRARY_CL)

stripSubsetLCBs_SOURCES = stripSubsetLCBs.cpp
stripSubsetLCBs_LDADD = $(LIBRARY_CL)

projectAndStrip_SOURCES = projectAndStrip.cpp
projectAndStrip_LDADD = $(LIBRARY_CL)

makeBadgerMatrix_SOURCES = makeBadgerMatrix.cpp
makeBadgerMatrix_LDADD = $(LIBRARY_CL)

multiEVD_SOURCES = multiEVD.cpp
multiEVD_LDADD = $(LIBRARY_CL)

randomGeneSample_SOURCES = randomGeneSample.cpp
randomGeneSample_LDADD = $(LIBRARY_CL)

getOrthologList_SOURCES = getOrthologList.cpp
getOrthologList_LDADD = $(LIBRARY_CL)

bbFilter_SOURCES = bbFilter.cpp
bbFilter_LDADD = $(LIBRARY_CL)

bbAnalyze_SOURCES = bbAnalyze.cpp
bbAnalyze_LDADD = $(LIBRARY_CL)

makeMc4Matrix_SOURCES = makeMc4Matrix.cpp
makeMc4Matrix_LDADD = $(LIBRARY_CL)

bbBreakOnGenes_SOURCES = bbBreakOnGenes.cpp
bbBreakOnGenes_LDADD = $(LIBRARY_CL)

backbone_global_to_local_SOURCES = backbone_global_to_local.cpp
backbone_global_to_local_LDADD = $(LIBRARY_CL)

coordinateTranslate_SOURCES = coordinateTranslate.cpp
coordinateTranslate_LDADD = $(LIBRARY_CL)

