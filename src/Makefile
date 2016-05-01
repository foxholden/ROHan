
CXX      = g++ -g
LIBGAB   = ../lib/libgab/

CXXFLAGS = -Wall -lm -O3 -lz -I../lib/libgab/ -I../lib/libgab//gzstream/ -I../lib/bamtools/include/ -I../lib/bamtools/src/  -c
LDFLAGS  = -lz  -lpthread 
LDLIBS   = ../lib/bamtools/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o ../lib/bamtools/lib/libbamtools.a


all: testComp testCompCont rohan refPlusSegToGeno sitesPlusGeno

%.o: %.cpp
	${CXX} ${CXXFLAGS} $^ -o $@


testComp.o:	testComp.cpp
	${CXX} ${CXXFLAGS} testComp.cpp

testComp:	testComp.o ../lib/libgab/utils.o  
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

testCompCont.o:	testCompCont.cpp
	${CXX} ${CXXFLAGS} testCompCont.cpp

testCompCont:	testCompCont.o ../lib/libgab/utils.o  
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

refPlusSegToGeno.o:	refPlusSegToGeno.cpp
	${CXX} ${CXXFLAGS} refPlusSegToGeno.cpp

refPlusSegToGeno:	refPlusSegToGeno.o ../lib/libgab/utils.o   ../lib/libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

sitesPlusGeno.o:	sitesPlusGeno.cpp
	${CXX} ${CXXFLAGS} sitesPlusGeno.cpp

sitesPlusGeno:	sitesPlusGeno.o ../lib/libgab/utils.o   ../lib/libgab/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

rohan.o:     rohan.cpp
	${CXX} ${CXXFLAGS} rohan.cpp

rohan:       rohan.o DataChunk.o DataToWrite.o GenoResults.o PositionResult.o miscfunc.o ../lib/bamtools/lib/libbamtools.a GenomicWindows.o GenomicRange.o ../lib/bamtools/lib/libbamtools.a ../lib/libgab/FastQObj.o ../lib/libgab/FastQParser.o ../lib/libgab/utils.o  ../lib/libgab/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

../lib/libgab/utils.h:
	rm -rfv ../lib/libgab/
	mkdir ../lib/libgab/
	git clone --recursive https://github.com/grenaud/libgab.git ../lib/libgab/


../lib/libgab/utils.o: ../lib/bamtools/lib/libbamtools.a  ../lib/libgab/utils.h
	make -C ../lib/libgab/

../lib/bamtools/src/api/BamAlignment.h:
	rm -rf ../lib/bamtools/
	mkdir ../lib/bamtools/
	git clone --recursive https://github.com/pezmaster31/bamtools.git ../lib/bamtools/

../lib/bamtools/lib/libbamtools.a: ../lib/bamtools/src/api/BamAlignment.h
	cd ../lib/bamtools/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..

clean :
	rm -f testComp.o testComp testCompCont.o testCompCont refPlusSegToGeno.o refPlusSegToGeno sitesPlusGeno.o sitesPlusGeno testPrior.o testPrior rohan.o rohan

