# -----------------------------------------------------------------
#   Makefile for inbix - In Silico Bioinformatics
#   
#   Compilation options
#       R plugins                   WITH_R_PLUGINS
#       Web-based version check     WITH_WEBCHECK
#       Ensure 32-bit binary        FORCE_32BIT 
#       (Ignored)                   WITH_ZLIB
#       Link to LAPACK              WITH_LAPACK
#       Force dynamic linking       FORCE_DYNAMIC
#
# ---------------------------------------------------------------------

# Set this variable to either UNIX, MAC or WIN
SYS = UNIX

# Leave blank after "=" to disable; put "= 1" to enable
WITH_R_PLUGINS = 1
WITH_WEBCHECK = 1
FORCE_32BIT = 
WITH_ZLIB = 1
WITH_LAPACK = 
#FORCE_DYNAMIC = $(INBIX_FORCE_DYNAMIC)
FORCE_DYNAMIC = 1

# Put C++ compiler here; Windows has it's own specific version
CXX_UNIX = g++
CXX_WIN = g++.exe

# Any other compiler flags here ( -Wall, -g, etc)
# optimized mode
# uncomment below for debug mode
CXXFLAGS += -g -I. -D_FILE_OFFSET_BITS=64 -Dfopen64=fopen \
  -Wno-write-strings -std=c++11
# uncomment below for optimized mode
#CXXFLAGS += -O3 -I. -D_FILE_OFFSET_BITS=64 -Dfopen64=fopen \
#  -DNDEBUG -Wno-write-strings -std=c++11
LDFLAGS += -L/usr/local/lib -lgsl -lboost_program_options

# Misc
LIB_LAPACK = /usr/lib/liblapack.so
LIB_BLAS = /usr/lib/libblas.so
LIB_IGRAPH = /usr/lib/libigraph.so

OUTPUT = inbix

# Some system specific flags
ifndef FORCE_DYNAMIC
 CXXFLAGS += -static
endif

ifeq ($(SYS),UNIX)
 CXXFLAGS += -DUNIX
 CXX = $(CXX_UNIX)
 ifndef FORCE_DYNAMIC
  CXXFLAGS += -static
  LIB += /usr/local/lib/libarmadillo.a
  LIB += /usr/lib/liblapack.a
  LIB += /usr/lib/libblas/libblas.a
  LIB += /usr/lib/gcc/x86_64-linux-gnu/4.4/libgfortran.a
 else
  LIB += -lgomp
  LIB += -larmadillo
  LIB += -llapack
  LIB += -lopenblas
  LIB += -lgsl
  LIB += -lgslcblas
 endif
endif

ifeq ($(SYS),MAC)
 CXXFLAGS += -DUNIX
 CXX = $(CXX_UNIX)
endif

ifdef FORCE_32BIT
 CXXFLAGS += -m32
endif

# Flags for web-based version check
ifdef WITH_WEBCHECK
ifeq ($(SYS),WIN)
 LIB += -lwsock32
endif
else
 CXXFLAGS += -DSKIP
endif

CXXFLAGS += -fopenmp

SRC = inbix.cpp plink.cpp options.cpp input.cpp binput.cpp tinput.cpp \
genome.cpp helper.cpp stats.cpp filters.cpp locus.cpp multi.cpp crandom.cpp \
cluster.cpp mds.cpp output.cpp informative.cpp assoc.cpp epi.cpp  \
prephap.cpp phase.cpp trio.cpp tdt.cpp sharing.cpp genepi.cpp sets.cpp	\
perm.cpp mh.cpp genedrop.cpp gxe.cpp merge.cpp hotel.cpp multiple.cpp \
haploCC.cpp haploTDT.cpp poo.cpp webcheck.cpp qfam.cpp linear.cpp \
bmerge.cpp parse.cpp mishap.cpp legacy.cpp homozyg.cpp segment.cpp  \
model.cpp logistic.cpp glm.cpp dcdflib.cpp elf.cpp dfam.cpp fisher.cpp	\
linput.cpp sockets.cpp lookup.cpp proxy.cpp pdriver.cpp haploQTL.cpp  \
haplohelper.cpp haplowindow.cpp genogroup.cpp nonfounderphasing.cpp \
clumpld.cpp genoerr.cpp em.cpp impute.cpp metaem.cpp profile.cpp  \
nlist.cpp whap.cpp simul.cpp gvar.cpp cnv.cpp step.cpp greport.cpp  \
flip.cpp qualscores.cpp cnvqt.cpp cfamily.cpp setscreen.cpp idhelp.cpp	\
tag.cpp hapglm.cpp lookup2.cpp blox.cpp zed.cpp dosage.cpp annot.cpp  \
metaanal.cpp ArmadilloFuncs.cpp regain.cpp InteractionNetwork.cpp \
CentralityRanker.cpp EpistasisEQtl.cpp Dataset.cpp AttributeRanker.cpp	\
DistanceMetrics.cpp PlinkInternalsDataset.cpp RReliefF.cpp ReliefF.cpp	\
ReliefFSeq.cpp SNReliefF.cpp DatasetInstance.cpp Insilico.cpp DgeData.cpp  \
PlinkInternalsDatasetInstance.cpp EvaporativeCooling.cpp Deseq.cpp Edger.cpp \
BirdseedData.cpp ChiSquared.cpp Statistics.cpp ReliefSeqController.cpp \
ArgumentHandler.cpp DataChar.cpp Data.cpp DataDouble.cpp DataFloat.cpp \
ForestClassification.cpp Forest.cpp ForestProbability.cpp RandomForest.cpp \
ForestRegression.cpp ForestSurvival.cpp TreeClassification.cpp \
Tree.cpp TreeProbability.cpp TreeRegression.cpp TreeSurvival.cpp utility.cpp \
EvaporativeCoolingPrivacy.cpp

HDR = plink.h options.h helper.h stats.h crandom.h sets.h phase.h \
perm.h model.h linear.h logistic.h dcdflib.h ipmpar.h cdflib.h	\
fisher.h sockets.h haplowindow.h genogroup.h clumpld.h nlist.h whap.h \
gvar.h cnv.h cfamily.h idhelp.h zed.h StringUtils.h \
ArmadilloFuncs.h regain.h InteractionNetwork.h CentralityRanker.h \
EpistasisEQtl.h Dataset.h AttributeRanker.h DistanceMetrics.h \
PlinkInternalsDataset.h RReliefF.h ReliefF.h ReliefFSeq.h SNReliefF.h  \
DatasetInstance.h Insilico.h BestN.h DgeData.h BirdseedData.h \
ChiSquared.h Statistics.h PlinkInternalsDatasetInstance.h \
ReliefSeqController.h EvaporativeCooling.h Deseq.h Edger.h \
ArgumentHandler.h DataDouble.h Data.h Forest.h ForestRegression.h \
RandomForest.h globals.h Tree.h TreeRegression.h utility.h DataChar.h \
DataFloat.h ForestClassification.h ForestProbability.h ForestSurvival.h \
TreeClassification.h TreeProbability.h TreeSurvival.h version.h \
EvaporativeCoolingPrivacy.h

ifdef WITH_R_PLUGINS
CXXFLAGS += -DWITH_R_PLUGINS
HDR += sisocks.h Rsrv.h Rconnection.h config.h
SRC += r.cpp Rconnection.cpp
ifeq ($(SYS),MAC)
LIB += -ldl
endif
ifeq ($(SYS),UNIX)
LIB += -ldl -lcrypt -pthread
endif
endif

ifdef WITH_ZLIB
CXXFLAGS += -DWITH_ZLIB
HDR += zfstream.h
SRC += zfstream.cpp
LIB += -lz
endif

ifdef WITH_LAPACK
CXXFLAGS += -DWITH_LAPACK
HDR += lapackf.h
SRC += lapackf.cpp
LIB += $(LIB_LAPACK) 
endif

OBJ = $(SRC:.cpp=.o)

all : $(OUTPUT) 

$(OUTPUT) :
	$(CXX) $(LDFLAGS) -o $(OUTPUT) $(OBJ) $(LIB) 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(CXXFLAGS) -c $*.cpp

.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:

clean:
	rm -f *.o *~ inbix

install:
	cp inbix /usr/local/bin

