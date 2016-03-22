/*
 * DataChunk
 * Date: Mar-21-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef DataChunk_h
#define DataChunk_h

#include "GenomicWindows.h"

using namespace std;

class DataChunk{
private:
    
    public:
    //vector<BamAlignment>  dataToProcess;    
    GenomicRange rangeGen;
    int rank;
    
    DataChunk();
    DataChunk(const DataChunk & other);
    ~DataChunk();
    DataChunk & operator= (const DataChunk & other);
};

class CompareDataChunk {
public:
    bool operator() ( DataChunk * cd1, DataChunk * cd2)  {
        //comparison code here
	return ( cd1->rank > cd2->rank );
    }
};


#endif
