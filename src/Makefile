
CXX      = g++ -std=c++11 -g -pg #-fpermissive #-g # -pg
LIBGAB   = ../lib/libgab/

CXXFLAGS = -Wall -lm -O3 -lz -I../lib/libgab/ -I../lib/libgab//gzstream/ -I../lib/samtools/ -I../lib/htslib/  -I../lib/libharu/include/   -c
LDFLAGS  = -lgsl -lgslcblas -llzma -lbz2 -lz  -lpthread 
LDLIBS   = ../lib/samtools/bedidx.o ../lib/htslib/libhts.a ../lib/samtools/libbam.a ../lib/samtools/libst.a 


all:  ../lib/libgab/utils.o ../lib/libharu/src/.libs/libhpdf.a ../lib/htslib/libhts.a ../lib/samtools/bedidx.o rohan refPlusSegToGeno sitesPlusGeno ref2DNAprof #testComp testCompCont

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

rohan:       rohan.o PdfWriter.o Hmm.o HmmState.o DataChunk.o DataToWrite.o GenoResults.o PositionResult.o miscfunc.o  GenomicWindows.o GenomicRange.o  ../lib/libgab/FastQObj.o ../lib/libgab/FastQParser.o ../lib/libharu/src/.libs/libhpdf.a ../lib/libgab/utils.o  ../lib/libgab/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


ref2DNAprof.o:	ref2DNAprof.cpp
	${CXX} ${CXXFLAGS} ref2DNAprof.cpp


ref2DNAprof:	ref2DNAprof.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


../lib/htslib/libhts.a:  ../lib/htslib/hts_internal.h
	cd ../lib/htslib/ && make && cd ../..

../lib/htslib/hts_internal.h:
	rm -rf ../lib/htslib/
        git clone --recursive https://github.com/samtools/htslib.git ../lib/htslib/


../lib/samtools/bedidx.o:  ../lib/samtools/sam.h
	cd ../lib/samtools/ && make && cd ../..

../lib/samtools/sam.h:
	rm -rf ../lib/samtools/
        git clone --recursive https://github.com/samtools/samtools.git ../lib/samtools/

../lib/libgab/utils.o:  ../lib/libgab/utils.h
	cd ../lib/libgab/ &&  make utils.o && cd ../..

../lib/libgab/utils.h:
	rm -rfv ../lib/libgab/
	mkdir ../lib/libgab/
	git clone --depth 1 https://github.com/grenaud/libgab.git ../lib/libgab/

../lib/libharu/include/hpdf.h:
	rm -rf ../lib/libharu/
	mkdir ../lib/libharu/
	git clone --depth 1 https://github.com/libharu/libharu.git ../lib/libharu/

../lib/libharu/src/.libs/libhpdf.a: ../lib/libharu/include/hpdf.h
	cd ../lib/libharu/ && ./buildconf.sh && ./configure && make

clean :
	rm -f *.o testComp  testCompCont  refPlusSegToGeno sitesPlusGeno  testPrior  rohan  ref2DNAprof


