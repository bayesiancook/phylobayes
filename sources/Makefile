#To compile with gsl, uncomment the following line
# USE_GSL=-DGSL

#To build the executables for mac osX, comment the following line
#STATIC=-static


ifdef USE_GSL
GSLLIB=-lgsl -lgslcblas
else
GSLLIB=
endif

LIBS=$(GSLLIB)

CC=g++
SYSLIB=
INCLUDES=
CPPFLAGS= -Wall -O3 $(USE_GSL) $(INCLUDES)
LDFLAGS=$(STATIC) 
INSTALL_DIR=
INSTALL_LIB=
SRCS=Clock.cpp Random.cpp mystring.cpp mystream.cpp PolyNode.cpp Chrono.cpp TaxaParameters.cpp Bipartition.cpp BipartitionList.cpp  TreeList.cpp Consensus.cpp Tree.cpp SubMatrix.cpp MCParameters.cpp PhyloBayes.cpp PhyloProbs.cpp Move.cpp Misc.cpp  Pruning.cpp Chain.cpp Sample.cpp PhyloEM.cpp correl.cpp correlation.cpp linalg.cpp
OBJS=$(patsubst %.cpp,%.o,$(patsubst %.c,%.o,$(SRCS)))
ALL_SRCS:=$(wildcard *.cpp *.c)
ALL_OBJS=$(patsubst %.cpp,%.o,$(patsubst %.c,%.o,$(ALL_SRCS)))

PROGSDIR=../data
ALL=pb readpb readdiv tracecomp bpcomp readcv cvrep sumcv ppred tree2ps bf stopafter subsample ancestral 
#ALL=makeconcat makerandomconcat taxasubsample pb stoppb readpb readdiv tracecomp bpcomp readcv cvrep sumcv ppred tree2ps bf stopafter subsample ancestral
PROGS=$(addprefix $(PROGSDIR)/, $(ALL))

# If we are on a windows platform, executables are .exe files
# [S.Bigaret] we assume here that we are compiling w/ MinGW, not Cygwin
EXEEXT=
uname := $(shell uname -s 2>/dev/null)
ifeq (MINGW,$(findstring MINGW,$(uname)))
	EXEEXT=.exe
	PROGS=$(addsuffix $(EXEEXT), $(ALL))
endif

.PHONY: all clean
all: $(PROGS)

# Rules for generate the dependencies automatically

%.d: %.cpp
	@echo "Generating dependencies for $<..."; \
	 set -e; rm -f $@; $(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	 sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; rm -f $@.$$$$

# Rules to create .o files from .cpp files
%.o: %.cpp %.d
	$(CC) -c $(CPPFLAGS) $*.cpp

# Include the dependencies unless the request was to clean eveything up
ifneq ($(MAKECMDGOALS),clean)
-include $(ALL_OBJS:.o=.d)
endif

# Targets

subdata$(EXEEXT): $(PROGSDIR)/subdata$(EXEEXT)
$(PROGSDIR)/subdata$(EXEEXT): SubData.o $(OBJS)
	$(CC) SubData.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

makeconcat$(EXEEXT): $(PROGSDIR)/makeconcat$(EXEEXT)
$(PROGSDIR)/makeconcat$(EXEEXT): MakeConcat.o $(OBJS)
	$(CC) MakeConcat.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

makerandomconcat$(EXEEXT): $(PROGSDIR)/makerandomconcat$(EXEEXT)
$(PROGSDIR)/makerandomconcat$(EXEEXT): MakeRandomConcat.o $(OBJS)
	$(CC) MakeRandomConcat.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

taxasubsample$(EXEEXT): $(PROGSDIR)/taxasubsample$(EXEEXT)
$(PROGSDIR)/taxasubsample$(EXEEXT): TaxaSubSample.o $(OBJS)
	$(CC) TaxaSubSample.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

subsample$(EXEEXT): $(PROGSDIR)/subsample$(EXEEXT)
$(PROGSDIR)/subsample$(EXEEXT): SubSample.o $(OBJS)
	$(CC) SubSample.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

stopafter$(EXEEXT): $(PROGSDIR)/stopafter$(EXEEXT)
$(PROGSDIR)/stopafter$(EXEEXT): StopAfter.o $(OBJS)
	$(CC) StopAfter.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

sumcv$(EXEEXT): $(PROGSDIR)/sumcv$(EXEEXT)
$(PROGSDIR)/sumcv$(EXEEXT): SumCV.o $(OBJS)
	$(CC) SumCV.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readcv$(EXEEXT): $(PROGSDIR)/readcv$(EXEEXT)
$(PROGSDIR)/readcv$(EXEEXT): ReadCV.o $(OBJS)
	$(CC) ReadCV.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

cvrep$(EXEEXT): $(PROGSDIR)/cvrep$(EXEEXT)
$(PROGSDIR)/cvrep$(EXEEXT): CVrep.o $(OBJS)
	$(CC) CVrep.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

ancestral$(EXEEXT): $(PROGSDIR)/ancestral$(EXEEXT)
$(PROGSDIR)/ancestral$(EXEEXT): SampleSub.o $(OBJS)
	$(CC) SampleSub.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

stoppb$(EXEEXT): $(PROGSDIR)/stoppb$(EXEEXT)
$(PROGSDIR)/stoppb$(EXEEXT): StopChain.o $(OBJS)
	$(CC) StopChain.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

ppred$(EXEEXT): $(PROGSDIR)/ppred$(EXEEXT)
$(PROGSDIR)/ppred$(EXEEXT):	PostPred.o $(OBJS)
	$(CC) PostPred.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

pb$(EXEEXT): $(PROGSDIR)/pb$(EXEEXT)
$(PROGSDIR)/pb$(EXEEXT): PB.o $(OBJS)
	$(CC) PB.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

bf$(EXEEXT): $(PROGSDIR)/bf$(EXEEXT)
$(PROGSDIR)/bf$(EXEEXT): BF.o $(OBJS)
	$(CC) BF.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

tree2ps$(EXEEXT): $(PROGSDIR)/tree2ps$(EXEEXT)
$(PROGSDIR)/tree2ps$(EXEEXT): TreeToPS.o $(OBJS)
	$(CC) TreeToPS.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

tracecomp$(EXEEXT): $(PROGSDIR)/tracecomp$(EXEEXT)
$(PROGSDIR)/tracecomp$(EXEEXT): TraceComp.o $(OBJS)
	$(CC) TraceComp.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

bpcomp$(EXEEXT): $(PROGSDIR)/bpcomp$(EXEEXT)
$(PROGSDIR)/bpcomp$(EXEEXT): BPCompare.o $(OBJS)
	$(CC) BPCompare.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readpb$(EXEEXT): $(PROGSDIR)/readpb$(EXEEXT)
$(PROGSDIR)/readpb$(EXEEXT):	ReadSample.o $(OBJS)
	$(CC) ReadSample.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readdiv$(EXEEXT): $(PROGSDIR)/readdiv$(EXEEXT)
$(PROGSDIR)/readdiv$(EXEEXT): ReadDiv.o $(OBJS)
	$(CC) ReadDiv.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

readcov$(EXEEXT): $(PROGSDIR)/readcov$(EXEEXT)
$(PROGSDIR)/readcov$(EXEEXT): ReadCov.o $(OBJS)
	$(CC) ReadCov.o $(OBJS) $(LDFLAGS) $(LIBS) -o $@

clean:
	-rm -f *.o *.d
	-rm -f $(PROGS)
#../data/stop ../data/pb ../data/read ../data/bpcomp ../data/tracecomp ../data/tree2ps

