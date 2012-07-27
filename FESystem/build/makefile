#either mac, machome or linux
PLATFORM=mac


SOURCE_DIR=./

ifeq (${PLATFORM},mac)
    ARPACK_LIB_DIR=/Users/bhatiam/Documents/usr/lib/
    LAPACK_BLAS_LIB_DIR=/usr/lib
    GFORT_SYSTEM_LIB_DIR=/opt/local/lib/gcc47/i386
    GFORT_SYSTEM_LIBS=-lgfortran -lgomp
    PLPLOT_INCLUDE_DIR=/opt/local/include
    PLPLOT_LIB_DIR=/opt/local/lib
    PLPLOT_LIBS=-lplplotcxxd -lplplotd
endif

ifeq (${PLATFORM},machome)
    ARPACK_LIB_DIR=/Users/manav/Documents/codes/numerical_libs/arpack
    LAPACK_BLAS_LIB_DIR=/usr/lib
    GFORT_SYSTEM_LIB_DIR=/opt/local/lib/gcc47/i386
    GFORT_SYSTEM_LIBS=-lgfortran -lgomp
    PLPLOT_INCLUDE_DIR=/opt/local/include
    PLPLOT_LIB_DIR=/opt/local/lib
    PLPLOT_LIBS=-lplplotcxxd -lplplotd
    PLATFORM=mac
endif

ifeq (${PLATFORM},linux)
    ARPACK_LIB_DIR=/home/bhatiam/usr/lib/
    LAPACK_BLAS_LIB_DIR=/home/bhatiam/usr/lib/
    GFORT_SYSTEM_LIB_DIR=/apps/global/oss/TOOLCHAIN/Linux/apps/gcc-4.7.0/lib/
    GFORT_SYSTEM_LIBS=-lgfortran -lgomp
    PLPLOT_INCLUDE_DIR=/home/bhatiam/usr/plplot/include
    PLPLOT_LIB_DIR=/home/bhatiam/usr/plplot/lib
    PLPLOT_LIBS=-lplplotcxxd -lplplotd
endif


GCC = g++
PYTHON = python

GCC_FLAGS_DEBUG = -g3 -O0 -Wall -m32
GCC_FLAGS_OPT = -g3 -O3 -m32 -ftree-vectorize 

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

EXTERNAL_INCLUDE_DIRS+=-I${SRC_SUBDIR}
EXTERNAL_INCLUDE_DIRS+=-I${FESYSTEM_DIR}/src
EXTERNAL_INCLUDE_DIRS+=-I${PLPLOT_INCLUDE_DIR}

EXTERNAL_LIB_DIRS=-L${SOURCE_DIR}
EXTERNAL_LIB_DIRS+=-L${ARPACK_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${LAPACK_BLAS_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${PLPLOT_LIB_DIR}
EXTERNAL_LIB_DIRS+=-L${GFORT_SYSTEM_LIB_DIR}


EXTERNAL_LIBS=-lfesystem-${METHOD}
EXTERNAL_LIBS+=-larpack -llapack -lblas
EXTERNAL_LIBS+=${GFORT_SYSTEM_LIBS}
EXTERNAL_LIBS+=${PLPLOT_LIBS}


SOURCES = $(wildcard ${SRC_SUBDIR}/*.cpp ${SRC_SUBDIR}/*/*.cpp ${SRC_SUBDIR}/*/*/*.cpp)
DEP_FILES= $(patsubst %.cpp,%-cpp.dep,${SOURCES})

all: fesystem-meshtest-${METHOD}

fesystem-meshtest-dbg: ../examples/MeshTest.cpp  libfesystem-dbg.a 
	${GCC} ${GCC_FLAGS_DEBUG} ${EXTERNAL_INCLUDE_DIRS} ${EXTERNAL_LIB_DIRS} $^ ${EXTERNAL_LIBS} -o $@

fesystem-meshtest-opt: ../examples/MeshTest.cpp  libfesystem-opt.a 
	${GCC} ${GCC_FLAGS_OPT} ${EXTERNAL_INCLUDE_DIRS} ${EXTERNAL_LIB_DIRS} $^ ${EXTERNAL_LIBS} -o $@

%-cpp.dep: %.cpp 
	${GCC} -M -I${EXTERNAL_INCLUDE_DIRS} -I${SRC_SUBDIR}  $*.cpp > $@
	${PYTHON} ${BUILD_DIR}/modify_dep_file.py $@ ${*D}

%-dbg-cpp.o: %.cpp
	${GCC} -c ${GCC_FLAGS_DEBUG} ${EXTERNAL_INCLUDE_DIRS} -I${SRC_SUBDIR}  $*.cpp -o $@

%-opt-cpp.o: %.cpp
	${GCC} -c ${GCC_FLAGS_OPT} ${EXTERNAL_INCLUDE_DIRS} -I${SRC_SUBDIR}  $*.cpp -o $@


libfesystem-${METHOD}.a: $(patsubst %.cpp,%-${METHOD}-cpp.o,${SOURCES}) 
	ar -rv $@ $^ 

clean:
	rm -rf ${DEP_FILES} $(patsubst %.cpp,%-${METHOD}-cpp.o,${SOURCES})  libfesystem-${METHOD}.a  fesystem-meshtest-${METHOD}*

clean-deps:
	rm -f ${DEP_FILES}

-include ${DEP_FILES}

