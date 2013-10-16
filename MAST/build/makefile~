#PLATFORM is either mac, machome or linux

SOURCE_DIR=./

# Binary type
ifeq (${BUILD64},0)
	BINTYPE=32bit
else
	BINTYPE=64bit
endif

ifeq (${PLATFORM},Darwin)
    GCC_PRE = g++
    GCC_COMP = clang++
    PYTHON = python
    ARPACK_LIB_DIR=/opt/local/lib
    ARPACK_LIB=-larpack
    LAPACK_BLAS_LIBS=-llapack -lblas 
    LAPACK_BLAS_LIB_DIR=/usr/lib
    GFORT_SYSTEM_LIB_DIR=/opt/local/lib/gcc48
    GFORT_SYSTEM_LIBS=-lgfortran -lgomp
    PLPLOT_INCLUDE_DIR=/opt/local/include
    HAVE_PLPLOT_LIBS=1
    PLPLOT_LIB_DIR=/opt/local/lib
    PLPLOT_LIBS=-lplplotcxxd -lplplotd
    METIS_INCLUDE_DIR=/opt/local/include/metis
    METIS_LIB_DIR=/opt/local/lib
    METIS_LIBS=-lmetis
    TBB_LIB_DIR=/opt/local/tbb41_20121003oss/lib
    TBB_INCLUDE_DIR=/opt/local/tbb41_20121003oss/include
    TBB_LIBS=-ltbb -ltbbmalloc
    BOOST_DIR=/opt/local/include
    FESYSTEM_DIR=../../../FESystem/FESystem/
    LIBMESH_DIR=../../../libmesh
    TRILINOS_DIR=/usr/local/trilinos
    VTK_INCLUDE_DIR=/opt/local/include/vtk-5.10
    MPI_INCLUDE_DIR=/opt/local/include/openmpi
    MPI_LIB_DIR=/opt/local/lib/openmpi
    GCC_FLAGS_DEBUG = -g3 -O0 -Wall -DDEBUG -felide-constructors 
    GCC_FLAGS_OPT = -DNDEBUG -O2 -felide-constructors -Qunused-arguments 
    GCC_FLAGS_PRO = -DNDEBUG -O2 -g -felide-constructors -Qunused-arguments 
endif


ifeq (${PLATFORM},Linux)
    GCC_PRE = g++
    GCC_COMP = g++
    PYTHON = python
    ARPACK_LIB_DIR=/usr/include 
    ARPACK_LIB=
    LAPACK_BLAS_LIBS=
    LAPACK_BLAS_LIB_DIR=/usr/lib
    GFORT_SYSTEM_LIB_DIR=/usr/include
    GFORT_SYSTEM_LIBS=
    PLPLOT_INCLUDE_DIR=/usr/include
    HAVE_PLPLOT_LIBS=0
    PLPLOT_LIB_DIR=
    PLPLOT_LIBS=
    METIS_INCLUDE_DIR=/home/bhatiam/usr/include
    METIS_LIB_DIR=/home/bhatiam/usr/lib
    METIS_LIBS=
    HDF_LIB_DIR=/apps/global/oss/hdf5/1.8.0/lib
    TBB_LIB_DIR=/usr/local/intel/2013-u1/tbb/lib/intel64
    TBB_INCLUDE_DIR=/usr/local/intel/2013-u1/tbb/include
    TBB_LIBS=-ltbb -ltbbmalloc
    BOOST_DIR=/usr/include
    FESYSTEM_DIR=../../../FESystem/FESystem/
    LIBMESH_DIR=../../../libmesh
    TRILINOS_DIR=
    VTK_LIB_DIR=/home/bhatiam/usr/vtk/lib/vtk-5.10
    VTK_LIBS=-lvtkIO -lvtkCommon -lvtkFiltering
    VTK_INCLUDE_DIR=/home/bhatiam/usr/vtk/include/vtk-5.10
    MPI_INCLUDE_DIR=/usr/local/openmpi/openmpi-1.6.3/x86_64/ib/gcc/include
    MPI_LIB_DIR=/usr/local/openmpi/openmpi-1.6.3/x86_64/ib/gcc/lib
    GCC_FLAGS_DEBUG = -g3 -O0 -Wall -std=c++0x
    GCC_FLAGS_OPT = -DNDEBUG -O2 -felide-constructors -std=c++0x
    GCC_FLAGS_PRO = -DNDEBUG -O2 -g -felide-constructors -std=c++0x
endif




D      = cd
ECHO    = echo
LN      = ln
LNFLAGS = -s
MAKE    = make
RM      = rm
RMFLAGS = -f
SHELL   = /bin/sh
ARFLAGS = rv
RANLIB   = ranlib
AR = ar


SRC_SUBDIR=../src
BUILD_DIR=${PWD}

EXTERNAL_INCLUDE_DIRS+=-I${FESYSTEM_DIR}/src
EXTERNAL_INCLUDE_DIRS+=-I${BOOST_DIR}
EXTERNAL_INCLUDE_DIRS+=-I${TBB_INCLUDE_DIR}
EXTERNAL_INCLUDE_DIRS+=-I${METIS_INCLUDE_DIR}
EXTERNAL_INCLUDE_DIRS+=-I${SRC_SUBDIR}
EXTERNAL_INCLUDE_DIRS+=-I${LIBMESH_DIR}/${SCALAR}/include
EXTERNAL_INCLUDE_DIRS+=-I${TRILINOS_DIR}/include
EXTERNAL_INCLUDE_DIRS+=-I${MPI_INCLUDE_DIR}
EXTERNAL_INCLUDE_DIRS+=-I${VTK_INCLUDE_DIR}
ifeq (${HAVE_PLPLOT_LIBS},1)
    EXTERNAL_INCLUDE_DIRS+=-I${PLPLOT_INCLUDE_DIR}
endif

EXTERNAL_LIB_DIRS=-L${FESYSTEM_DIR}/build
EXTERNAL_LIB_DIRS+=-L${SOURCE_DIR}
EXTERNAL_LIB_DIRS+=-L${LAPACK_BLAS_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${ARPACK_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${METIS_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${TBB_LIB_DIR} -Wl,-rpath,${TBB_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${GFORT_SYSTEM_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${MPI_LIB_DIR}
EXTERNAL_LIB_DIRS+=-Wl,-rpath,${VTK_LIB_DIR}
ifeq (${HAVE_PLPLOT_LIBS},1)
    EXTERNAL_LIB_DIRS+=-L${PLPLOT_LIB_DIR}
endif
FESYSTEM_LIB_DIR=-L${FESYSTEM_DIR}/build
LIBMESH_LIB_DIR=-L${LIBMESH_DIR}/${SCALAR}/lib

EXTERNAL_LIBS+=${ARPACK_LIB}
EXTERNAL_LIBS+=${LAPACK_BLAS_LIBS} 
EXTERNAL_LIBS+=${METIS_LIBS}
EXTERNAL_LIBS+=${TBB_LIBS}
#EXTERNAL_LIBS+=${GFORT_SYSTEM_LIBS}
ifeq (${HAVE_PLPLOT_LIBS},1)
    EXTERNAL_LIBS+=${PLPLOT_LIBS}
endif
FESYSTEM_LIBS=-lfesystem-${METHOD}
LIBMESH_LIBS=-lmesh_${METHOD}
MPI_LIBS=-lmpi_cxx -lmpi


ifeq (${HAVE_PLPLOT_LIBS},1)
    COMPILATION_DEFINES=-DHAVE_PLPLOT
else
    COMPILATION_DEFINES=
endif


FESYSTEM_LIB_SOURCES = $(wildcard ${SRC_SUBDIR}/*.cpp ${SRC_SUBDIR}/*/*.cpp ${SRC_SUBDIR}/*/*/*.cpp)
DEP_FILES= $(patsubst %.cpp,%-cpp.dep,${FESYSTEM_LIB_SOURCES})

all: mast-${METHOD}-${SCALAR}

mast-dbg-real: $(patsubst %.cpp,%-${METHOD}-${SCALAR}-cpp.o,${FESYSTEM_LIB_SOURCES}) 
	${GCC_COMP} ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} ${LIBMESH_LIB_DIR} $^ ${FESYSTEM_LIBS}  ${LIBMESH_LIBS}  ${EXTERNAL_LIB_DIRS} ${EXTERNAL_LIBS} ${MPI_LIBS} -o $@

mast-opt-real: $(patsubst %.cpp,%-${METHOD}-${SCALAR}-cpp.o,${FESYSTEM_LIB_SOURCES}) 
	${GCC_COMP} ${GCC_FLAGS_OPT} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} ${LIBMESH_LIB_DIR} $^ ${FESYSTEM_LIBS}  ${LIBMESH_LIBS}  ${EXTERNAL_LIB_DIRS} ${EXTERNAL_LIBS} ${MPI_LIBS} -o $@

mast-oprof-real: $(patsubst %.cpp,%-${METHOD}-${SCALAR}-cpp.o,${FESYSTEM_LIB_SOURCES}) 
	${GCC_COMP} ${GCC_FLAGS_PRO} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} ${LIBMESH_LIB_DIR} $^ ${FESYSTEM_LIBS} ${LIBMESH_LIBS}  ${EXTERNAL_LIB_DIRS} ${EXTERNAL_LIBS} ${MPI_LIBS} -o $@

mast-dbg-complex: $(patsubst %.cpp,%-${METHOD}-${SCALAR}-cpp.o,${FESYSTEM_LIB_SOURCES}) 
	${GCC_COMP} ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} ${LIBMESH_LIB_DIR} $^ ${FESYSTEM_LIBS} ${LIBMESH_LIBS}  ${EXTERNAL_LIB_DIRS} ${EXTERNAL_LIBS} ${MPI_LIBS} -o $@

mast-opt-complex: $(patsubst %.cpp,%-${METHOD}-${SCALAR}-cpp.o,${FESYSTEM_LIB_SOURCES}) 
	${GCC_COMP} ${GCC_FLAGS_OPT} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} ${LIBMESH_LIB_DIR} $^ ${FESYSTEM_LIBS} ${LIBMESH_LIBS}  ${EXTERNAL_LIB_DIRS} ${EXTERNAL_LIBS} ${MPI_LIBS} -o $@

mast-oprof-complex: $(patsubst %.cpp,%-${METHOD}-${SCALAR}-cpp.o,${FESYSTEM_LIB_SOURCES}) 
	${GCC_COMP} ${GCC_FLAGS_PRO} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} ${LIBMESH_LIB_DIR} $^ ${FESYSTEM_LIBS}  ${LIBMESH_LIBS}  ${EXTERNAL_LIB_DIRS} ${EXTERNAL_LIBS} ${MPI_LIBS} -o $@

%-cpp.dep: %.cpp 
	${GCC_PRE} -M  ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp > $@
	${PYTHON} ${BUILD_DIR}/modify_dep_file.py $@ ${*D}

%-dbg-cpp.o: %.cpp
	${GCC_COMP} -c ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-opt-cpp.o: %.cpp
	${GCC_COMP} -c ${GCC_FLAGS_OPT} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-dbg-real-cpp.o: %.cpp
	${GCC_COMP} -c ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-opt-real-cpp.o: %.cpp
	${GCC_COMP} -c ${GCC_FLAGS_OPT} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-oprof-real-cpp.o: %.cpp
	${GCC_COMP} -c ${GCC_FLAGS_PRO} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-dbg-complex-cpp.o: %.cpp
	${GCC_COMP} -c ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-opt-complex-cpp.o: %.cpp
	${GCC_COMP} -c ${GCC_FLAGS_OPT} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-oprof-complex-cpp.o: %.cpp
	${GCC_COMP} -c ${GCC_FLAGS_PRO} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@


clean:
	rm -rf $(patsubst %.cpp,%-${METHOD}-${SCALAR}-cpp.o,${FESYSTEM_LIB_SOURCES})  mast-${METHOD}-${SCALAR}
	rm -f ${DEP_FILES}

-include ${DEP_FILES}

