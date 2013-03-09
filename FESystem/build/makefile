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
    LAPACK_BLAS_LIBS=-llapack -lblas -latlas
    LAPACK_BLAS_LIB_DIR=/usr/lib
    GFORT_SYSTEM_LIB_DIR=/opt/local/lib/gcc47
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
    LIBMESH_DIR=../../../libmesh
    TRILINOS_DIR=/usr/local/trilinos
    VTK_INCLUDE_DIR=/opt/local/include/vtk-5.10
    MPI_INCLUDE_DIR=/opt/local/include/openmpi
    MPI_LIB_DIR=/opt/local/lib/openmpi
endif

ifeq (${PLATFORM},Darwinhome)
    GCC_PRE = g++
    GCC_COMP = clang++
    PYTHON = python
    ARPACK_LIB_DIR=/opt/local/lib/
    ARPACK_LIB=-larpack
    LAPACK_BLAS_LIBS=-llapack -lblas -latlas
    LAPACK_BLAS_LIB_DIR=/usr/lib
    GFORT_SYSTEM_LIB_DIR=/opt/local/lib/gcc47
    GFORT_SYSTEM_LIBS=-lgfortran -lgomp
    HAVE_PLPLOT_LIBS=1
    PLPLOT_INCLUDE_DIR=/opt/local/include
    PLPLOT_LIB_DIR=/opt/local/lib
    PLPLOT_LIBS=-lplplotcxxd -lplplotd
    METIS_INCLUDE_DIR=/opt/local/include/metis
    METIS_LIB_DIR=/opt/local/lib
    METIS_LIBS=-lmetis
    TBB_LIB_DIR=/opt/local/tbb41_20121003oss/lib
    TBB_INCLUDE_DIR=/opt/local/tbb41_20121003oss/include
    TBB_LIBS=-ltbb -ltbbmalloc
    PLATFORM=Darwin
    BOOST_DIR=/opt/local/include
    LIBMESH_DIR=../../../libmesh
    TRILINOS_DIR=/usr/local/trilinos
    VTK_INCLUDE_DIR=/opt/local/include/vtk-5.10
    MPI_INCLUDE_DIR=/opt/local/include/openmpi
    MPI_LIB_DIR=/opt/local/lib/openmpi
endif

ifeq (${PLATFORM},Linux)
    GCC_PRE = g++
    GCC_COMP = g++
    PYTHON = python
    ARPACK_LIB_DIR=/home/bhatiam/usr/lib/
    ARPACK_LIB=-larpack
    LAPACK_BLAS_LIBS=-llapack -lblas
    LAPACK_BLAS_LIB_DIR=/home/bhatiam/usr/lib/
    GFORT_SYSTEM_LIB_DIR=/apps/global/oss/TOOLCHAIN/Linux/apps/gcc-4.7.0/lib/
    GFORT_SYSTEM_LIBS=-lgfortran -lgomp
    HAVE_PLPLOT_LIBS=0
    METIS_INCLUDE_DIR=/home/bhatiam/usr/metis/include
    METIS_LIB_DIR=/home/bhatiam/usr/metis/lib
    METIS_LIBS=-lmetis
    TBB_LIB_DIR=/home/bhatiam/usr/tbb/tbb41_20121003oss/lib/intel64/cc4.1.0_libc2.4_kernel2.6.16.21
    TBB_INCLUDE_DIR=/home/bhatiam/usr/tbb/tbb41_20121003oss/include
    TBB_LIBS=-ltbb -ltbbmalloc
    BOOST_DIR=/home/bhatia/usr/include
endif

ifeq (${PLATFORM},script)
    ARPACK_LIB_DIR=${BASEDIR}/install/${KERNEL}-${BINTYPE}/lib
    LAPACK_BLAS_LIB_DIR=${BASEDIR}/install/${KERNEL}-${BINTYPE}/lib
    GFORT_SYSTEM_LIB_DIR=${GCCHOME}/lib
    GFORT_SYSTEM_LIBS=-lgfortran -lgomp
    HAVE_PLPLOT_LIBS=0
    PLATFORM=linux
endif


GCC_FLAGS_DEBUG = -g3 -O0 -Wall 
GCC_FLAGS_OPT = -DNDEBUG -O2 -felide-constructors -Qunused-arguments   
#GCC_FLAGS_OPT = -g0 -Os 

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
APP_SUBDIR=../applications
BUILD_DIR=${PWD}

EXTERNAL_INCLUDE_DIRS+=-I${BOOST_DIR}
EXTERNAL_INCLUDE_DIRS+=-I${TBB_INCLUDE_DIR}
EXTERNAL_INCLUDE_DIRS+=-I${METIS_INCLUDE_DIR}
EXTERNAL_INCLUDE_DIRS+=-I${SRC_SUBDIR}
EXTERNAL_INCLUDE_DIRS+=-I${SRC_SUBDIR}/../examples
EXTERNAL_INCLUDE_DIRS+=-I${SRC_SUBDIR}/../applications
EXTERNAL_INCLUDE_DIRS+=-I${LIBMESH_DIR}/${SCALAR}/include
EXTERNAL_INCLUDE_DIRS+=-I${TRILINOS_DIR}/include
EXTERNAL_INCLUDE_DIRS+=-I${MPI_INCLUDE_DIR}
EXTERNAL_INCLUDE_DIRS+=-I${VTK_INCLUDE_DIR}
ifeq (${HAVE_PLPLOT_LIBS},1)
    EXTERNAL_INCLUDE_DIRS+=-I${PLPLOT_INCLUDE_DIR}
endif

EXTERNAL_LIB_DIRS=-L${SOURCE_DIR}
EXTERNAL_LIB_DIRS+=-L${ARPACK_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${LAPACK_BLAS_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${METIS_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${TBB_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${GFORT_SYSTEM_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${MPI_LIB_DIR}
ifeq (${HAVE_PLPLOT_LIBS},1)
    EXTERNAL_LIB_DIRS+=-L${PLPLOT_LIB_DIR}
endif
LIBMESH_LIB_DIR=-L${LIBMESH_DIR}/${SCALAR}/lib

EXTERNAL_LIBS=-lfesystem-${METHOD}
EXTERNAL_LIBS+=${ARPACK_LIB}
EXTERNAL_LIBS+=${LAPACK_BLAS_LIBS} 
EXTERNAL_LIBS+=${METIS_LIBS}
EXTERNAL_LIBS+=${TBB_LIBS}
EXTERNAL_LIBS+=${GFORT_SYSTEM_LIBS}
ifeq (${HAVE_PLPLOT_LIBS},1)
    EXTERNAL_LIBS+=${PLPLOT_LIBS}
endif
LIBMESH_LIBS=-lmesh_${METHOD}
MPI_LIBS=-lmpi_cxx -lmpi


ifeq (${HAVE_PLPLOT_LIBS},1)
    COMPILATION_DEFINES=-DHAVE_PLPLOT
else
    COMPILATION_DEFINES=
endif


FESYSTEM_LIB_SOURCES = $(wildcard ${SRC_SUBDIR}/*.cpp ${SRC_SUBDIR}/*/*.cpp ${SRC_SUBDIR}/*/*/*.cpp)
FESYSTEM_APP_SOURCES = $(wildcard ${APP_SUBDIR}/*.cpp ${APP_SUBDIR}/*/*.cpp ${APP_SUBDIR}/*/*/*.cpp)
FESYSTEM_EXAMPLE_SOURCES = $(wildcard ${SRC_SUBDIR}/../examples/*.cpp)
DEP_FILES= $(patsubst %.cpp,%-cpp.dep,${FESYSTEM_LIB_SOURCES} ${FESYSTEM_APP_SOURCES} ${FESYSTEM_EXAMPLE_SOURCES})

all: fesystem-${METHOD} euler-${METHOD}-${SCALAR}

euler-dbg-real: $(patsubst %.cpp,%-${METHOD}-${SCALAR}-cpp.o,${FESYSTEM_APP_SOURCES}) 
	${GCC_COMP} ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} ${LIBMESH_LIB_DIR}  ${LIBMESH_LIBS}  ${EXTERNAL_LIB_DIRS} $^ ${MPI_LIBS} -o $@

euler-opt-real: $(patsubst %.cpp,%-${METHOD}-${SCALAR}-cpp.o,${FESYSTEM_APP_SOURCES}) 
	${GCC_COMP} ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} ${LIBMESH_LIB_DIR}  ${LIBMESH_LIBS}  ${EXTERNAL_LIB_DIRS} $^ ${MPI_LIBS} -o $@

euler-dbg-complex: $(patsubst %.cpp,%-${METHOD}-${SCALAR}-cpp.o,${FESYSTEM_APP_SOURCES}) 
	${GCC_COMP} ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} ${LIBMESH_LIB_DIR}  ${LIBMESH_LIBS}  ${EXTERNAL_LIB_DIRS} $^ ${MPI_LIBS} -o $@

euler-opt-complex: $(patsubst %.cpp,%-${METHOD}-${SCALAR}-cpp.o,${FESYSTEM_APP_SOURCES}) 
	${GCC_COMP} ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} ${LIBMESH_LIB_DIR}  ${LIBMESH_LIBS}  ${EXTERNAL_LIB_DIRS} $^ ${MPI_LIBS} -o $@

fesystem-dbg: $(patsubst %.cpp,%-${METHOD}-cpp.o,${FESYSTEM_EXAMPLE_SOURCES})  libfesystem-dbg.a 
	${GCC_COMP} ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} ${LIBMESH_LIB_DIR}  ${LIBMESH_LIBS}  ${EXTERNAL_LIB_DIRS} $^ ${MPI_LIBS} -o $@

fesystem-opt: $(patsubst %.cpp,%-${METHOD}-${SCALAR}-cpp.o,${FESYSTEM_EXAMPLE_SOURCES})  libfesystem-opt.a 
	${GCC_COMP} ${GCC_FLAGS_OPT} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} ${EXTERNAL_LIB_DIRS} $^ ${EXTERNAL_LIBS} ${MPI_LIBS} -o $@

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

%-dbg-complex-cpp.o: %.cpp
	${GCC_COMP} -c ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-opt-complex-cpp.o: %.cpp
	${GCC_COMP} -c ${GCC_FLAGS_OPT} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@


libfesystem-${METHOD}.a: $(patsubst %.cpp,%-${METHOD}-cpp.o,${FESYSTEM_LIB_SOURCES}) 
	ar -rv $@ $^ 

clean:
	rm -rf $(patsubst %.cpp,%-${METHOD}-cpp.o,${FESYSTEM_EXAMPLE_SOURCES}) $(patsubst %.cpp,%-${METHOD}-cpp.o,${FESYSTEM_LIB_SOURCES})  libfesystem-${METHOD}.a  fesystem-${METHOD}
	rm -rf $(patsubst %.cpp,%-${METHOD}-${SCALAR}-cpp.o,${FESYSTEM_APP_SOURCES})  euler-${METHOD}-${SCALAR}

clean-deps:
	rm -f ${DEP_FILES}

-include ${DEP_FILES}

