/*
 * DataToWrite
 * Date: Mar-22-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef DataToWrite_h
#define DataToWrite_h

using namespace std;

#include <vector>


#include "GenomicRange.h"
#include "PositionResult.h"

class DataToWrite{
 private:
    
 public:
    vector<PositionResult *>  * vecPositionResults;
    GenomicRange rangeGen;
    int rank;
    hResults hetEstResults;

    DataToWrite();
    DataToWrite(const DataToWrite & other);
    ~DataToWrite();
    DataToWrite & operator= (const DataToWrite & other);
};


class CompareDataToWrite {
public:
    bool operator() ( DataToWrite * cd1, DataToWrite * cd2)  {
        //comparison code here
	return ( cd1->rank > cd2->rank );
    }
};


#endif
