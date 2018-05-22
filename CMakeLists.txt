cmake_minimum_required(VERSION 2.7)

project(inbix)

enable_language(CXX Fortran)

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  message(STATUS "Mac OS X system detected, setting compiler variables")
  set(CMAKE_C_COMPILER "clang" CACHE STRING "C compiler" FORCE)
  set(CMAKE_CXX_COMPILER "clang++" CACHE STRING "C++ compiler" FORCE)
  set(CMAKE_CXX_STANDARD 11)
  set(INBIX_CXX_DEFINES "-D_FILE_OFFSET_BITS=64 -Dfopen64=fopen -DSYS=MAC -DMAC -DWITH_LAPACK=1 -DLAPACK -DWITH_ZLIB=1 -DWITH_WEBCHECK=0 -DSKIP -DWITH_R_PLUGINS=1") 
  # turn off the nit-picky clang compiler warnings in old PLINK code
  set(INBIX_WARNINGS_OFF "-Wno-return-type -Wno-header-guard -Wno-writable-strings -Wno-user-defined-warnings")
  set(OMP_FLAG_AppleClang "-Xclang -fopenmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INBIX_CXX_FLAGS} ${INBIX_CXX_DEFINES} ${OMP_FLAG_AppleClang} ${INBIX_WARNINGS_OFF}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread -L /usr/local/lib -liomp5")
else(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  message(STATUS "Linux system detected, setting compiler variables")
  set(CMAKE_C_COMPILER gcc CACHE STRING "C compiler" FORCE)
  set(CMAKE_CXX_COMPILER g++ CACHE STRING "C++ compiler" FORCE)
# enable either debugging mode with embedded symbols
  # or optimized with no debugging symbols
  # set(BUILD_DEBUG OFF CACHE BOOL "" FORCE)
  set(BUILD_DEBUG ON CACHE BOOL "" FORCE)
  if(BUILD_DEBUG)
    set(INBIX_CXX_FLAGS "-g -fno-var-tracking-assignments -Wno-write-strings -std=c++11")
  else()
    set(INBIX_CXX_FLAGS "-O2 -fno-var-tracking-assignments -Wno-write-strings -std=c++11")
  endif()
  set(INBIX_CXX_DEFINES "-D_FILE_OFFSET_BITS=64 -Dfopen64=fopen -Dunix -DSYS=UNIX -DUNIX -DWITH_LAPACK=1 -DLAPCK -DWITH_ZLIB=1 -DWITH_WEBCHECK=0 -DSKIP -DWITH_R_PLUGINS=1") 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INBIX_CXX_FLAGS} ${INBIX_CXX_DEFINES}")
endif()

# build static
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
#set(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++ -pthread")
#set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
# build shared
set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)

add_executable(inbix
options.cpp parse.cpp lapackf.cpp zfstream.cpp zed.cpp whap.cpp webcheck.cpp
trio.cpp tinput.cpp tdt.cpp tag.cpp step.cpp sockets.cpp simul.cpp sharing.cpp
setscreen.cpp segment.cpp r.cpp Rconnection.cpp qualscores.cpp qfam.cpp
proxy.cpp profile.cpp prephap.cpp poo.cpp plink.cpp phase.cpp perm.cpp
pdriver.cpp nonfounderphasing.cpp nlist.cpp multiple.cpp multi.cpp mishap.cpp
mh.cpp metaem.cpp metaanal.cpp merge.cpp lookup.cpp lookup2.cpp locus.cpp
linput.cpp legacy.cpp informative.cpp impute.cpp hotel.cpp homozyg.cpp
haplowindow.cpp haploTDT.cpp haploQTL.cpp haplohelper.cpp haploCC.cpp
hapglm.cpp gxe.cpp gvar.cpp greport.cpp glm.cpp genome.cpp genogroup.cpp
genoerr.cpp genepi.cpp genedrop.cpp flip.cpp fisher.cpp filters.cpp epi.cpp
em.cpp dosage.cpp dfam.cpp dcdflib.cpp crandom.cpp cnvqt.cpp cnv.cpp
clumpld.cpp cfamily.cpp bmerge.cpp blox.cpp binput.cpp assoc.cpp annot.cpp
mds.cpp cluster.cpp input.cpp linear.cpp model.cpp logistic.cpp output.cpp
sets.cpp elf.cpp idhelp.cpp BirdseedData.cpp ChiSquared.cpp Statistics.cpp
DgeData.cpp regain.cpp stats.cpp TreeProbability.cpp TreeClassification.cpp
DataFloat.cpp DataChar.cpp SNReliefF.cpp RReliefF.cpp DistanceMetrics.cpp
TreeSurvival.cpp Deseq.cpp Edger.cpp Insilico.cpp Data.cpp
ReliefSeqController.cpp InteractionNetwork.cpp CentralityRanker.cpp
ArmadilloFuncs.cpp EvaporativeCooling.cpp AttributeRanker.cpp ReliefFSeq.cpp
ForestSurvival.cpp ForestProbability.cpp Forest.cpp ArgumentHandler.cpp
utility.cpp TreeRegression.cpp ForestClassification.cpp ForestRegression.cpp
EpistasisEQtl.cpp Tree.cpp helper.cpp ReliefF.cpp DataDouble.cpp
PlinkInternalsDataset.cpp DatasetInstance.cpp Dataset.cpp
PlinkInternalsDatasetInstance.cpp RandomForest.cpp
EvaporativeCoolingPrivacy.cpp DcVar.cpp inbix.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

find_package(Threads REQUIRED)
if(Threads_FOUND)
  message(STATUS "Threads package loaded")
  target_link_libraries(inbix ${CMAKE_THREAD_LIBS_INIT})
else()
  message(FATAL_ERROR "Could not find the 'Threads' package")
endif()

OPTION (USE_OpenMP "Use OpenMP" ON)
IF(USE_OpenMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    message(STATUS "OpenMP found")
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    message(STATUS "New make flags for OpenMP: ${CMAKE_CXX_FLAGS}")
  else()
    message(STATUS "OpenMP not found, flags might not be necessary")
  ENDIF()
ENDIF()

find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
  target_link_libraries(inbix ${LAPACK_LIBRARIES})
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
  message(STATUS "LAPACK package loaded")
else()
  message(FATAL_ERROR "Could not find the 'LAPACK' package")
endif()

find_package(BLAS REQUIRED)
if(BLAS_FOUND)
  target_link_libraries(inbix ${BLAS_LIBRARIES})
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
  message(STATUS "BLAS package loaded")
else()
  message(FATAL_ERROR "Could not find the 'BLAS' package")
endif()

find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIRS})
  target_link_libraries(inbix ${ZLIB_LIBRARIES})
  message(STATUS "ZLIB package loaded")
else()
  message(FATAL_ERROR "Could not find the 'ZLIB' package")
endif()

find_package(GSL REQUIRED)
if(GSL_FOUND)
  include_directories(${GSL_INCLUDE_DIRS})
  target_link_libraries(inbix ${GSL_LIBRARIES})
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${GSL_CLBAS_LIBRARY} ${GSL_LIBRARY}")
  message(STATUS "GSL package loaded")
else()
  message(FATAL_ERROR "Could not find the 'GSL' package")
endif()

find_package(Armadillo REQUIRED)
if(ARMADILLO_FOUND)
  include_directories(${ARMADILLO_INCLUDE_DIRS})
  target_link_libraries(inbix ${ARMADILLO_LIBRARIES})
  message(STATUS "Armadillo package loaded")
else()
  message(FATAL_ERROR "Could not find the 'Armadillo' package")
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  message(STATUS "Mac OS X system detected")
  # find_package(GCRYPT REQUIRED)
  # is this next 'if' necessary, or do anything?
  find_library(GCRYPT_LIBRARY NAMES gcrypt)
  if(GCRYPT-NOTFOUND)
    message(FATAL_ERROR "Could not find the 'gcrypt' library")
  else()
  message(STATUS "Mac OS X, using gcrypt")
    target_link_libraries(inbix ${GCRYPT_LIBRARY})
  endif()
elseif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  message(STATUS "Linux system detected")
  message(STATUS "Linux, using crypt")
  find_library(CRYPT_LIBRARY NAMES crypt)
  target_link_libraries(inbix ${CRYPT_LIBRARY})
else(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  messge(STATUS "Windows is not supported")
endif()

INSTALL(TARGETS inbix RUNTIME DESTINATION bin)
