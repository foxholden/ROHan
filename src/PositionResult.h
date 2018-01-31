/*
 * PositionResult
 * Date: Mar-21-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef PositionResult_h
#define PositionResult_h



#include <string>

#include <api/BamConstants.h>
#include <api/BamMultiReader.h>

#include "utils.h"

using namespace std;
using namespace BamTools;

class PositionResult{
 private:
    
 public:

    /* int          refID; */
    unsigned int pos ;
    char         refB;
    /* char         altB; */
    int avgMQ;
    int dp;
    int         baseC[4];
    /* int          refC; */
    /* int          altC; */
    // 0  1  2  3  4  5  6  7  8  9
    //AA,AC,AG,AT,CC,CG,CT,GG,GT,TT
    long double  ll[10];
    
    int         mostLikelyGenoIdx;
    int         mostLikelyGenoHetIdx;

    long double gq;
    
    /* long double  lqual; */
    /* long double  llCov; */
    /* int          geno; */
    /* char         genoS [2]; */

    
    PositionResult();
    PositionResult(const PositionResult & other);
    ~PositionResult();
    string toString(const RefVector * references,const int & refID) const;

    PositionResult & operator= (const PositionResult & other);
    //    friend ostream & operator<<(ostream & os, const PositionResult & ct);
    pair<char,char> hetIndex2Bases() const;
    char homoIndex2Base() const;
    int base2HomoIndex(const char b) const;
    int bases2hetIndex(const char c1,const char c2) const;
};



#endif
