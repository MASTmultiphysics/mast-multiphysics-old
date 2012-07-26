PYTHON = python
#GCC = g++-mp-4.7
GCC = g++

#GCC_FLAGS_DEBUG = -g3 -O0 -Wall -fpermissive
#GCC_FLAGS_OPT = -g0 -O3 -ftree-vectorize
GCC_FLAGS_DEBUG = -g3 -O0 -Wall -arch i386 -fpermissive
GCC_FLAGS_OPT = -g0 -O3 -arch i386 -ftree-vectorize

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



BUILD_DIR = ${PWD}
SRC_SUBDIR =../src
SOURCES = $(wildcard ${SRC_SUBDIR}/*.cpp ${SRC_SUBDIR}/*/*.cpp ${SRC_SUBDIR}/*/*/*.cpp)
DEP_FILES= $(patsubst %.cpp,%-cpp.dep,${SOURCES})

LAPACK_INCLUDE_DIR=
LAPACK_LIBDIR=/usr/lib

EXTERNAL_INCLUDE_DIRS=-I/opt/local/include
#EXTERNAL_INCLUDE_DIRS+=-I${LAPACK_INCLUDE_DIR}
EXTERNAL_INCLUDE_DIRS+=-I../src
EXTERNAL_INCLUDE_DIRS+=-I/home/bhatiam/usr/plplot32/include
EXTERNAL_INCLUDE_DIRS+=-I/opt/local/include/plplot
EXTERNAL_LIB_DIRS=-L${LAPACK_LIBDIR}
EXTERNAL_LIB_DIRS+=-L${BUILD_DIR}
EXTERNAL_LIB_DIRS+=-L/home/bhatiam/usr/plplot32/lib
EXTERNAL_LIB_DIRS+=-L/usr/lib
EXTERNAL_LIB_DIRS+=-L/Users/bhatiam/Documents/usr/lib
EXTERNAL_LIB_DIRS+=-L/opt/local/lib
EXTERNAL_LIB_DIRS+=-L/opt/local/lib/gcc47/i386/
EXTERNAL_LIBS=-larpack
EXTERNAL_LIBS+=-llapack
EXTERNAL_LIBS+=-lblas
EXTERNAL_LIBS+=-lplplotcxxd
EXTERNAL_LIBS+=-lplplotd
EXTERNAL_LIBS+=-lgfortran


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
	rm -f ${DEP_FILES} $(patsubst %.cpp,%-${METHOD}-cpp.o,${SOURCES})  libfesystem-${METHOD}.a  fesystem-meshtest-${METHOD}

clean-deps:
	rm -f ${DEP_FILES}

-include ${DEP_FILES}

