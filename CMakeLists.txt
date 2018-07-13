cmake_minimum_required(VERSION 3.7)

project(inbix)

enable_language(CXX Fortran)

set(CMAKE_CXX_STANDARD 11)
set(INBIX_WARNINGS_OFF "-Wno-return-type -Wno-header-guard -Wno-writable-strings -Wno-user-defined-warnings -Wno-write-strings ")
set(INBIX_CXX_DEFINES_COMMON "-D_FILE_OFFSET_BITS=64 -Dfopen64=fopen -DWITH_LAPACK=1 -DLAPACK -DWITH_ZLIB=1 -DWITH_WEBCHECK=0 -DSKIP -DWITH_R_PLUGINS=1") 
set(INBIX_CXX_DEFINES_LINUX "${INBIX_CXX_DEFINES_COMMON} -Dunix -DSYS=UNIX -DUNIX ")
set(INBIX_CXX_DEFINES_MAC "${INBIX_CXX_DEFINES_COMMON} -DSYS=MAC -DMAC ")

if(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  message(STATUS "Mac OS X system detected, setting compiler variables")
  set(CMAKE_C_COMPILER "clang" CACHE STRING "C compiler" FORCE)
  set(CMAKE_CXX_COMPILER "clang++" CACHE STRING "C++ compiler" FORCE)
  # turn off the nit-picky clang compiler warnings in old PLINK code
  set(OMP_FLAG_AppleClang "-Xclang -fopenmp")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${INBIX_CXX_FLAGS} ${INBIX_CXX_DEFINES_MAC} ${OMP_FLAG_AppleClang} ${INBIX_WARNINGS_OFF}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L /usr/local/lib -liomp5")
else(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  message(STATUS "Linux system detected, setting compiler variables")
  set(CMAKE_C_COMPILER gcc CACHE STRING "C compiler" FORCE)
  set(CMAKE_CXX_COMPILER g++ CACHE STRING "C++ compiler" FORCE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-var-tracking-assignments ${INBIX_CXX_FLAGS} ${INBIX_CXX_DEFINES_LINUX} ${INBIX_WARNINGS_OFF}")
endif()

# enable either debug or optimized compiler options
set(BUILD_DEBUG ON CACHE BOOL "" FORCE)
if(BUILD_DEBUG)
  # debugging mode with embedded symbols
  set(INBIX_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -g")
else()
  # or optimized with no debugging symbols
  set(INBIX_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
endif()

# build static
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static")
#set(CMAKE_EXE_LINKER_FLAGS "-static -static-libgcc -static-libstdc++ -pthread")
#set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
# build shared
set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)

# source code location
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(src)
file(GLOB_RECURSE SRC_LIST src/*.cpp src/*.h)

# Executable
add_executable(${PROJECT_NAME} ${SRC_LIST})

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

# -----------------------------------------------------------------------------
# From: https://vicrucann.github.io/tutorials/quick-cmake-doxygen/
# first we can indicate the documentation build as an option and 
# set it to OFF by default
option(BUILD_DOC "Build documentation" OFF)

# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_BINARY_DIR}/../docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # note the option ALL which allows to build the docs together with 
    # the application
    add_custom_target( doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen needs to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
# -----------------------------------------------------------------------------

INSTALL(TARGETS inbix RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
INSTALL(DIRECTORY man DESTINATION ${CMAKE_INSTALL_PREFIX}/man)
