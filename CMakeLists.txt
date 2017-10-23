cmake_minimum_required(VERSION 2.8)

include(GNUInstallDirs)

project(inbix)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g -fopenmp -fno-var-tracking-assignments -Wno-write-strings -std=gnu++0x -D_FILE_OFFSET_BITS=64 -Dfopen64=fopen -DSYS=UNIX -DUNIX -DWITH_LAPACK=1 -DLAPCK -DWITH_ZLIB=1 -DWITH_WEBCHECK=0 -DSKIP -DWITH_R_PLUGINS=1")

add_executable(inbix
options.cpp parse.cpp lapackf.cpp zfstream.cpp zed.cpp whap.cpp webcheck.cpp trio.cpp tinput.cpp tdt.cpp tag.cpp step.cpp sockets.cpp simul.cpp sharing.cpp setscreen.cpp segment.cpp r.cpp Rconnection.cpp qualscores.cpp qfam.cpp proxy.cpp profile.cpp prephap.cpp poo.cpp plink.cpp phase.cpp perm.cpp pdriver.cpp nonfounderphasing.cpp nlist.cpp multiple.cpp multi.cpp mishap.cpp mh.cpp metaem.cpp metaanal.cpp merge.cpp lookup.cpp lookup2.cpp locus.cpp linput.cpp legacy.cpp informative.cpp impute.cpp hotel.cpp homozyg.cpp haplowindow.cpp haploTDT.cpp haploQTL.cpp haplohelper.cpp haploCC.cpp hapglm.cpp gxe.cpp gvar.cpp greport.cpp glm.cpp genome.cpp genogroup.cpp genoerr.cpp genepi.cpp genedrop.cpp flip.cpp fisher.cpp filters.cpp epi.cpp em.cpp dosage.cpp dfam.cpp dcdflib.cpp crandom.cpp cnvqt.cpp cnv.cpp clumpld.cpp cfamily.cpp bmerge.cpp blox.cpp binput.cpp assoc.cpp annot.cpp mds.cpp cluster.cpp input.cpp linear.cpp model.cpp logistic.cpp output.cpp sets.cpp elf.cpp idhelp.cpp BirdseedData.cpp ChiSquared.cpp Statistics.cpp DgeData.cpp regain.cpp stats.cpp TreeProbability.cpp TreeClassification.cpp DataFloat.cpp DataChar.cpp SNReliefF.cpp RReliefF.cpp DistanceMetrics.cpp TreeSurvival.cpp Deseq.cpp Edger.cpp Insilico.cpp Data.cpp ReliefSeqController.cpp InteractionNetwork.cpp CentralityRanker.cpp ArmadilloFuncs.cpp EvaporativeCooling.cpp AttributeRanker.cpp ReliefFSeq.cpp ForestSurvival.cpp ForestProbability.cpp Forest.cpp ArgumentHandler.cpp utility.cpp TreeRegression.cpp ForestClassification.cpp ForestRegression.cpp EpistasisEQtl.cpp Tree.cpp helper.cpp ReliefF.cpp DataDouble.cpp PlinkInternalsDataset.cpp DatasetInstance.cpp Dataset.cpp PlinkInternalsDatasetInstance.cpp RandomForest.cpp EvaporativeCoolingPrivacy.cpp DcVar.cpp inbix.cpp 
)

# Define headers for this library. PUBLIC headers are used for
# compiling the library, and will be added to consumers' build
# paths.
target_include_directories(inbix PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# library dependencies, all available through packages in Debian/Ubuntu
target_link_libraries(inbix
    dl
    crypt
    pthread
    gomp
    z
    boost_program_options
    gsl
    gslcblas
    armadillo
    lapack
    blas)

# 'make install' to the correct locations (provided by GNUInstallDirs).
INSTALL(TARGETS inbix RUNTIME DESTINATION bin)
