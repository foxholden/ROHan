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

    int          refID;
    unsigned int pos ;
    char         refB;
    char         altB;
    int          refC;
    int          altC;

    long double  rrll;
    long double  rall;
    long double  aall;

    long double  lqual;
    long double  llCov;
    int          geno;

    
    PositionResult();
    PositionResult(const PositionResult & other);
    ~PositionResult();
    string toString(const RefVector  references) const;

    PositionResult & operator= (const PositionResult & other);
    //    friend ostream & operator<<(ostream & os, const PositionResult & ct);
};



#endif
