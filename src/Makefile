
CXX      = g++ #-g -pg
LIBGAB   = ../lib/libgab/

CXXFLAGS = -Wall -lm -O3 -lz -I../lib/libgab/ -I../lib/libgab//gzstream/ -I../lib/bamtools/include/ -I../lib/bamtools/src/  -c
LDFLAGS  = -lz  -lpthread 
LDLIBS   = ../lib/bamtools/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o ../lib/bamtools/lib/libbamtools.a


all: testComp testCompCont rohan refPlusSegToGeno sitesPlusGeno ref2DNAprof

%.o: %.cpp ../lib/libgab/utils.o  
	${CXX} ${CXXFLAGS} $^ -o $@

testComp:	testComp.o 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


testCompCont:	testCompCont.o ../lib/libgab/utils.o  
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


refPlusSegToGeno:	refPlusSegToGeno.o ../lib/libgab/utils.o   ../lib/libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

sitesPlusGeno:	sitesPlusGeno.o ../lib/libgab/utils.o   ../lib/libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

rohan:       rohan.o DataChunk.o DataToWrite.o GenoResults.o PositionResult.o miscfunc.o ../lib/bamtools/lib/libbamtools.a GenomicWindows.o GenomicRange.o ../lib/bamtools/lib/libbamtools.a ../lib/libgab/FastQObj.o ../lib/libgab/FastQParser.o ../lib/libgab/utils.o  ../lib/libgab/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

../lib/libgab/utils.h:
	rm -rfv ../lib/libgab/
	mkdir ../lib/libgab/
	git clone --recursive https://github.com/grenaud/libgab.git ../lib/libgab/

ref2DNAprof.o:	ref2DNAprof.cpp
	${CXX} ${CXXFLAGS} ref2DNAprof.cpp


ref2DNAprof:	ref2DNAprof.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)



../lib/libgab/utils.o: ../lib/bamtools/lib/libbamtools.a  ../lib/libgab/utils.h
	make -C ../lib/libgab/

../lib/bamtools/src/api/BamAlignment.h:
	rm -rf ../lib/bamtools/
	mkdir ../lib/bamtools/
	git clone --recursive https://github.com/pezmaster31/bamtools.git ../lib/bamtools/

../lib/bamtools/lib/libbamtools.a: ../lib/bamtools/src/api/BamAlignment.h
	cd ../lib/bamtools/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..

clean :
	rm -f testComp.o testComp testCompCont.o testCompCont refPlusSegToGeno.o refPlusSegToGeno sitesPlusGeno.o sitesPlusGeno testPrior.o testPrior rohan.o rohan ref2DNAprof.o ref2DNAprof


