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
    /* int         baseC[4]; */
    /* int          refC; */
    /* int          altC; */

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
    pair<char,char> hetIndex2Bases();
    char homoIndex2Base();
};



#endif
