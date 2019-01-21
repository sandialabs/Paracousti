include Makefile_local
ifndef SRCHOME
ifdef ROOTDIR
SRCHOME=${ROOTDIR}/ParAcousti
else
SRCHOME=~/src/ParAcousti
endif
endif
ifndef FC
FC=gfortran
endif
ifndef CC
CC=mpicc
endif
ifndef CXX
CXX=mpic++
endif
ifndef MPICXX
MPICXX=mpic++
endif
ifndef OBJDIR
OBJDIR=.
endif
ifndef BINDIR
BINDIR=~/src/ParAcousti/acoustic
endif

COMMON_FLAGS=-w -O3 -mtune=native ${COMMON_FLAGS_LOCAL}

C_FLAGS=  ${COMMON_FLAGS} --param large-function-growth=1000 --param inline-unit-growth=500 --param max-inline-insns-single=800 \
 -DUSE_WARNING=0 -DUSE_MASS_PASS=1 \
 -DMAX_NUM_SPACE_COEFFS=2 \
 -DCOMPILE_DATE="\"`date`\"" -DUSE_MPI_SEND=1 -DUSE_VAMPIR=0 \
 -DUSE_IMMEDIATE_MPI_SEND=1 -DPARALLEL_PROGRAM=1 \
 -fexcess-precision=fast -ftree-vectorize -ffast-math -ftree-vectorizer-verbose=1 ${C_FLAGS_LOCAL}

F_FLAGS=${COMMON_FLAGS} 

LIBS= -L${SRCHOME}/lib -lnetcdf -lm -lc -lstdc++ ${LIBS_LOCAL}

UTIL=${SRCHOME}/utils
SGFD=${SRCHOME}/sgfd
UTILOBJ=${UTIL}/${OBJDIR}
SGFDOBJ=${SGFD}/${OBJDIR}

CFILES=${SGFD}/sgfd_util.c ${UTIL}/io_procs.c ${UTIL}/nstdutil_c.c \
${UTIL}/mpi_procs.c

CCFILES=parallel_acousti.cc acousticBoundary.cc ${SGFD}/sgfd.cc ${UTIL}/selector.cc \
${SGFD}/extra_output.cc ${SGFD}/movingSgfdSources.cc \
${SGFD}/sgfdCommandParser.cc ${SGFD}/sgfdLoops.cc ${SGFD}/sgfdReceiverGroups.cc \
${SGFD}/sgfdReceivers.cc ${SGFD}/sgfdReceiverNetwork.cc ${SGFD}/sgfdSources.cc \
${UTIL}/commandParser.cc \
${SGFD}/sgfdBoundary.cc ${UTIL}/readfile.cc ${UTIL}/nstdutil.cc ${UTIL}/xtrautil.cc \
updateAcousticCPMLFull.cc acousticCPMLNameSpace.cc fixed_acoustic.cc \
updateAcousticSpongeFull.cc acousticSpongeNameSpace.cc \
updateAcousticMPMLFull.cc acousticMPMLNameSpace.cc \
updateAcousticAttenMPMLFull.cc acousticAttenMPMLNameSpace.cc \
acoustic_control.cc paracousti.cc

FORTFILES=

C_OBJ_FILES=$(CFILES:.c=.o)
CC_OBJ_FILES=$(CCFILES:.cc=.o)
FORT_OBJ_FILES=$(FORTFILES:.f=.o)

C_D_FILES=$(CFILES:.c=.d)
CC_D_FILES=$(CCFILES:.cc=.d)

OBJECTS=${C_OBJ_FILES} ${CC_OBJ_FILES} ${FORT_OBJ_FILES} 

INCLUDES= -I${INCHOME} -I${SRCHOME}/include -I${UTIL} -I${SGFD} -I/usr/local/include -I/usr/include ${INCLUDES_LOCAL}

main:
	make Acoustic

Acoustic: Makefile
Acoustic: ${OBJECTS}
	$(MPICXX) ${C_FLAGS} \
		-o ${BINDIR}/ParAcousti \
		${OBJECTS} ${LIBS}

$(C_D_FILES): %.d: %.c
	$(SHELL) -ec '${CC} -MM ${CFLAGS} $(INCLUDES) -MT $(patsubst %.c,%.o,$<) $< \
		| sed '\''s|\($*\)\.o[ :]*|\1.o $@ : |g'\'' > $@; \
		[ -s $@ ] || rm -f $@'

$(CC_D_FILES): %.d: %.cc
	$(SHELL) -ec '${CXX} -MM ${CFLAGS} $(INCLUDES) -MT $(patsubst %.cc,%.o,$<) $< \
		| sed '\''s|\($*\)\.o[ :]*|\1.o $@ : |g'\'' > $@; \
		[ -s $@ ] || rm -f $@'

include $(C_D_FILES) $(CC_D_FILES)

setup:
	mkdir -p ${OBJDIR}
	mkdir -p ${BINDIR}
	mkdir -p ../sgfd/${OBJDIR}
	mkdir -p ../utils/${OBJDIR}

cleanall:
	/bin/rm -f  ${BINDIR}/ParAcousti \
		${OBJECTS} ${C_D_FILES} ${CC_D_FILES}

clean:
	/bin/rm -f ${OBJECTS} ${C_D_FILES} ${CC_D_FILES}

$(C_OBJ_FILES): %.o: %.c
	$(CC) $(C_FLAGS) -c $< $(INCLUDES) -o $@

$(CC_OBJ_FILES): %.o: %.cc
	$(CXX) $(C_FLAGS) -c $< $(INCLUDES) -o $@

$(FORT_OBJ_FILES): %.o: %.f
	$(FC) $(F_FLAGS) -c $< -o $@

