#PLATFORM is either mac, machome or linux

SOURCE_DIR=./

# Binary type
ifeq (${BUILD64},0)
	BINTYPE=32bit
else
	BINTYPE=64bit
endif

ifeq (${PLATFORM},Darwin)
    ARPACK_LIB_DIR=/opt/local/lib
    LAPACK_BLAS_LIB_DIR=/usr/lib
    GFORT_SYSTEM_LIB_DIR=/opt/local/lib/gcc47
    GFORT_SYSTEM_LIBS=-lgfortran -lgomp
    PLPLOT_INCLUDE_DIR=/opt/local/include
    HAVE_PLPLOT_LIBS=1
    PLPLOT_LIB_DIR=/opt/local/lib
    PLPLOT_LIBS=-lplplotcxxd -lplplotd
    METIS_LIB_DIR=/opt/local/lib
    METIS_LIBS=-lmetis
    BOOST_DIR=/opt/local/include
endif

ifeq (${PLATFORM},Darwinhome)
    ARPACK_LIB_DIR=/opt/local/lib/
    LAPACK_BLAS_LIB_DIR=/usr/lib
    GFORT_SYSTEM_LIB_DIR=/opt/local/lib/gcc47
    GFORT_SYSTEM_LIBS=-lgfortran -lgomp
    HAVE_PLPLOT_LIBS=1
    PLPLOT_INCLUDE_DIR=/opt/local/include
    PLPLOT_LIB_DIR=/opt/local/lib
    PLPLOT_LIBS=-lplplotcxxd -lplplotd
    METIS_LIB_DIR=/opt/local/lib
    METIS_LIBS=-lmetis
    PLATFORM=Darwin
    BOOST_DIR=/opt/local/include
endif

ifeq (${PLATFORM},Linux)
    ARPACK_LIB_DIR=/home/bhatiam/usr/lib/
    LAPACK_BLAS_LIB_DIR=/home/bhatiam/usr/lib/
    GFORT_SYSTEM_LIB_DIR=/apps/global/oss/TOOLCHAIN/Linux/apps/gcc-4.7.0/lib/
    GFORT_SYSTEM_LIBS=-lgfortran -lgomp
    HAVE_PLPLOT_LIBS=0
    METIS_LIB_DIR=/opt/local/lib
    METIS_LIBS=-lmetis
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

GCC = g++
PYTHON = python

GCC_FLAGS_DEBUG = -g3 -O0 -Wall 
GCC_FLAGS_OPT = -g3 -O3 -ftree-vectorize 

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

EXTERNAL_INCLUDE_DIRS+=-I${BOOST_DIR}
EXTERNAL_INCLUDE_DIRS+=-I${SRC_SUBDIR}
EXTERNAL_INCLUDE_DIRS+=-I${SRC_SUBDIR}/../examples
ifeq (${HAVE_PLPLOT_LIBS},1)
    EXTERNAL_INCLUDE_DIRS+=-I${PLPLOT_INCLUDE_DIR}
endif

EXTERNAL_LIB_DIRS=-L${SOURCE_DIR}
EXTERNAL_LIB_DIRS+=-L${ARPACK_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${LAPACK_BLAS_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${METIS_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${GFORT_SYSTEM_LIB_DIR}
ifeq (${HAVE_PLPLOT_LIBS},1)
    EXTERNAL_LIB_DIRS+=-L${PLPLOT_LIB_DIR}
endif

EXTERNAL_LIBS=-lfesystem-${METHOD}
EXTERNAL_LIBS+=-larpack -llapack -lblas -latlas
EXTERNAL_LIBS+=${METIS_LIBS}
EXTERNAL_LIBS+=${GFORT_SYSTEM_LIBS}
ifeq (${HAVE_PLPLOT_LIBS},1)
    EXTERNAL_LIBS+=${PLPLOT_LIBS}
endif


ifeq (${HAVE_PLPLOT_LIBS},1)
    COMPILATION_DEFINES=-DHAVE_PLPLOT
else
    COMPILATION_DEFINES=
endif


SOURCES = $(wildcard ${SRC_SUBDIR}/*.cpp ${SRC_SUBDIR}/*/*.cpp ${SRC_SUBDIR}/*/*/*.cpp)
DEP_FILES= $(patsubst %.cpp,%-cpp.dep,${SOURCES})
EXAMPLE_SOURCES = $(wildcard ${SRC_SUBDIR}/../examples/*.cpp)
EXAMPLE_DEP_FILES= $(patsubst %.cpp,%-cpp.dep,${EXAMPLE_SOURCES})

all: fesystem-meshtest-${METHOD}

fesystem-meshtest-dbg: $(patsubst %.cpp,%-${METHOD}-cpp.o,${EXAMPLE_SOURCES})  libfesystem-dbg.a 
	${GCC} ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} ${EXTERNAL_LIB_DIRS} $^ ${EXTERNAL_LIBS} -o $@

fesystem-meshtest-opt: $(patsubst %.cpp,%-${METHOD}-cpp.o,${EXAMPLE_SOURCES})  libfesystem-opt.a 
	${GCC} ${GCC_FLAGS_OPT} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} ${EXTERNAL_LIB_DIRS} $^ ${EXTERNAL_LIBS} -o $@

%-cpp.dep: %.cpp 
	${GCC} -M  ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp > $@
	${PYTHON} ${BUILD_DIR}/modify_dep_file.py $@ ${*D}

%-dbg-cpp.o: %.cpp
	${GCC} -c ${GCC_FLAGS_DEBUG} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@

%-opt-cpp.o: %.cpp
	${GCC} -c ${GCC_FLAGS_OPT} ${COMPILATION_DEFINES} ${EXTERNAL_INCLUDE_DIRS} $*.cpp -o $@


libfesystem-${METHOD}.a: $(patsubst %.cpp,%-${METHOD}-cpp.o,${SOURCES}) 
	ar -rv $@ $^ 

clean:
	rm -rf ${DEP_FILES} $(patsubst %.cpp,%-${METHOD}-cpp.o,${EXAMPLE_SOURCES}) $(patsubst %.cpp,%-${METHOD}-cpp.o,${SOURCES})  libfesystem-${METHOD}.a  fesystem-meshtest-${METHOD}*

clean-deps:
	rm -f ${DEP_FILES}

-include ${DEP_FILES}

