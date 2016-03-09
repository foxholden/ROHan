
CXX      = g++   
LIBGAB   = /home/gabriel/lib/

CXXFLAGS = -Wall -lm -O3 -lz -I${LIBGAB} -I${LIBGAB}/gzstream/ -c
LDFLAGS  = -lz


all: testComp testCompCont 

testComp.o:	testComp.cpp
	${CXX} ${CXXFLAGS} testComp.cpp

testComp:	testComp.o ${LIBGAB}utils.o  
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

testCompCont.o:	testCompCont.cpp
	${CXX} ${CXXFLAGS} testCompCont.cpp

testCompCont:	testCompCont.o ${LIBGAB}utils.o  
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

clean :
	rm -f testComp.o testComp testCompCont.o testCompCont

