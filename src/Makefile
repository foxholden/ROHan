
CXX      = g++
LIBGAB   = /home/gabriel/lib/

CXXFLAGS = -Wall -lm -O3 -lz -I${LIBGAB} -I${LIBGAB}/gzstream/ -I../bamtools/include/ -I../bamtools/src/  -c
LDFLAGS  = -lz  -lpthread 
LDLIBS   = ../bamtools/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o ../bamtools/lib/libbamtools.a


all: testComp testCompCont testPrior rohan

%.o: %.cpp
	${CXX} ${CXXFLAGS} $^ -o $@

testPrior.o:	testPrior.cpp
	${CXX} ${CXXFLAGS} testPrior.cpp

testPrior:	testPrior.o ${LIBGAB}utils.o  
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

testComp.o:	testComp.cpp
	${CXX} ${CXXFLAGS} testComp.cpp

testComp:	testComp.o ${LIBGAB}utils.o  
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

testCompCont.o:	testCompCont.cpp
	${CXX} ${CXXFLAGS} testCompCont.cpp

testCompCont:	testCompCont.o ${LIBGAB}utils.o  
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

rohan.o:     rohan.cpp
	${CXX} ${CXXFLAGS} rohan.cpp

rohan:       rohan.o miscfunc.o ../bamtools/lib/libbamtools.a GenomicWindows.o GenomicRange.o ../bamtools/lib/libbamtools.a ${LIBGAB}FastQObj.o ${LIBGAB}FastQParser.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

../bamtools/src/api/BamAlignment.h:
	rm -rf ../bamtools/
	git clone --recursive https://github.com/pezmaster31/bamtools.git ../

../bamtools/lib/libbamtools.a: ../bamtools/src/api/BamAlignment.h
	cd ../bamtools/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..

clean :
	rm -f testComp.o testComp testCompCont.o testCompCont testPrior.o testPrior rohan.o rohan

