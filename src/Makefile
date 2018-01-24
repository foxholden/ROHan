
CXX      = g++ -std=c++11 -fpermissive #-g # -pg
LIBGAB   = ../lib/libgab/

CXXFLAGS = -Wall -lm -O3 -lz -I../lib/libgab/ -I../lib/libgab//gzstream/ -I../lib/bamtools/include/ -I../lib/bamtools/src/  -c
LDFLAGS  = -lgsl -lgslcblas -lz  -lpthread 
LDLIBS   = ../lib/bamtools/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o ../lib/bamtools/build/src/api/libbamtools.a


all: ../lib/bamtools/build/src/api/libbamtools.a ../lib/libgab/utils.o  rohan refPlusSegToGeno sitesPlusGeno ref2DNAprof #testComp testCompCont

%.o: %.cpp #../lib/libgab/utils.o  
	${CXX} ${CXXFLAGS} $^ -o $@

#testComp:	testComp.o 
#	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)
#testCompCont:	testCompCont.o ../lib/libgab/utils.o  
#	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


refPlusSegToGeno:	refPlusSegToGeno.o ../lib/libgab/utils.o   ../lib/libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

sitesPlusGeno:	sitesPlusGeno.o ../lib/libgab/utils.o   ../lib/libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

rohan:       rohan.o Hmm.o HmmState.o DataChunk.o DataToWrite.o GenoResults.o PositionResult.o miscfunc.o ../lib/bamtools/build/src/api/libbamtools.a GenomicWindows.o GenomicRange.o  ../lib/libgab/FastQObj.o ../lib/libgab/FastQParser.o ../lib/libgab/utils.o  ../lib/libgab/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

../lib/libgab/utils.h:
	rm -rfv ../lib/libgab/
	mkdir ../lib/libgab/
	git clone --recursive https://github.com/grenaud/libgab.git ../lib/libgab/

ref2DNAprof.o:	ref2DNAprof.cpp
	${CXX} ${CXXFLAGS} ref2DNAprof.cpp


ref2DNAprof:	ref2DNAprof.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)



../lib/libgab/utils.o: ../lib/bamtools/build/src/api/libbamtools.a  ../lib/libgab/utils.h
	make -C ../lib/libgab/

../lib/bamtools/src/api/BamAlignment.h:
	rm -rf ../lib/bamtools/
	mkdir ../lib/bamtools/
	git clone --recursive https://github.com/pezmaster31/bamtools.git ../lib/bamtools/

../lib/bamtools/build/build/src/api/libbamtools.a: ../lib/bamtools/src/api/BamAlignment.h
	cd ../lib/bamtools/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..

clean :
	rm -f *.o testComp  testCompCont  refPlusSegToGeno sitesPlusGeno  testPrior  rohan  ref2DNAprof


