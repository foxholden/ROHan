/*
 * miscfunc
 * Date: Jun-08-2014 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef miscfunc_h
#define miscfunc_h

#include <stdlib.h>
#include <vector>
#include <set>
#include <string>
#include <gzstream.h>


#include "utils.h"
#include "GenomicRange.h"

/* #define DEBUGCOV */
using namespace std;

static vector<long double> lnFactVec;

#define HMMCODEMIN 0
#define HMMCODEMID 1
#define HMMCODEMAX 2


#define MIN2(a,b) (((a)<(b))?(a):(b))
#define MAX2(a,b) (((a)>(b))?(a):(b))

#define MIN3(a,b,c) MIN2(MIN2(a,b),c)
#define MAX3(a,b,c) MAX2(MAX2(a,b),c)



typedef struct { 
    long double p;
    int idx;
} emission;


typedef struct { 
    vector<int> seq;
    long double llik;
 } hmmpath;



typedef struct { 
    vector< vector<long double > > m;
    long double llik;
 } fbreturnVal;

 typedef struct{

    long double hAvg;
    long double pAvg;
    long double sAvg;
    
    long double hMin;    
    long double pMin;    
    long double sMin;
	
    long double hMax;
    long double pMax;
    long double sMax;
    
    uint64_t rohSegments   ;
    uint64_t nonrohSegments;
    uint64_t unsureSegments;

    double avgLengthROHSegments;

    fbreturnVal postprob;
 } hmmRes;


typedef struct { 
    int maxReadLength;
    bool isPe;
 } rgInfo;

typedef struct { 
    long double  h;
    long double  hLow;
    long double  hHigh;
    long double  errb;   
    unsigned int sites;
    bool         hasConverged;
 } hResults;


typedef struct { 
    long double s[12];
 } substitutionRates;


//  model->obs
//  0  A->A 
//  1  A->C 
//  2  A->G 
//  3  A->T 
//  4  C->A 
//  5  C->C 
//  6  C->G 
//  7  C->T 
//  8  G->A 
//  9  G->C 
//  10 G->G 
//  11 G->T 
//  12 T->A 
//  13 T->C 
//  14 T->G 
//  15 T->T 

typedef struct { 
    long double s[16];
} probSubstition;

typedef struct { 
    long double gl[16];
} babdlikelihood;

typedef struct { 
    long double p[4][4];
} diNucleotideProb;


//frequency of A,C,G,T
typedef struct { 
    long double f[4];
 } alleleFrequency;

typedef struct { 
    char b;//observed base
    int bq;//base quality phred scale
    int mpq;//mapping qual
    long double probDDc;
    long double probADc;
    long double probDAc;
    long double probAAc;

    long double probDDe;
    long double probADe;
    long double probDAe;
    long double probAAe;

 } singleBase;


typedef struct { 
    string id;
    int leftCoord;//base quality phred scale
    int rightCoord;
 } region;


typedef struct { 
    string chr;
    unsigned int coord;

    char ancAllele;
    char derAllele;
    
    vector<singleBase> sitesBAMd;
    vector<singleBase> sitesBAMa;
    vector<double>     freqDerived;
 } singleSite;


//To store consensus information
typedef struct { 
    long double perror[4];
    long double phred[4];
    long double perrorC[4];
    long double phredC[4];

    char ref;
    char consensus; //for the endogenous
} PHREDgeno;


typedef struct { 

    char ref;
    char base; //for the endogenous
    int  pos;
    double q;

    double aprob;
    double cprob;
    double gprob;
    double tprob;

} logRecord;
//#define DEBUGSINGLEREAD

//To store a single read
typedef struct { 
    uint8_t base;
    uint8_t qual;
    uint8_t mapq;       
    uint8_t pos5p;       
    uint8_t lengthF;             
    /* int dist5p; */
    /* int dist3p; */
#ifdef DEBUGSINGLEREAD
    //string name; //TO REMOVE
#endif
    bool isrv;
} singleRead;

//To store contamination likelihood
typedef struct { 
    string filename;    
    long double contaminationRate;
    long double logLike;
    //    long double logLike;
} contaminationEstimate;


typedef struct { 
    vector<singleRead> readsVec;
    /* long double mapqAvg; */
    //int cov;
    char refBase;
    //int refID;
    int posAlign;
    int avgMQ;
    int          baseC[4];
    bool skipPosition;
} positionInformation;


typedef struct { 
    long double likeBaseNoindel[4];
    long double likeBaseNoindelCont[4][4];
    long double likeBaseNoindelNoBoundary[4];           //when we do not consider bases at the ends of reads
    long double likeBaseNoindelContNoBoundary[4][4];    //when we do not consider bases at the ends of reads

    int  covPerBase[4];
    int  covPerBaseNoBoundary[4];

    long double mapqAvg;
    
    int numDel;
    long double llikDeletion;
    long double llikNoDeletion;

    long double llikDeletionBoth;
    long double llikDeletionCont;
    long double llikDeletionEndo;
    long double llikDeletionNone;

    set<string> allInserts;

    vector<string> insertionRight;
    map< pair<string,string> , long double> insertion2loglikeEndoCont; //key is (endo ins,cont ins) to log likelihood

    map<string,int> insertion2count;
    map<string,long double> insertion2loglike;
    //map<string,long double> insertion2loglikeCont;

    int cov;
} singlePosInfo;



void readNucSubstitionRatesFreq(const string filename,     vector<substitutionRates> & subVec);
void readNucSubstitionFreq(     const string filename,     vector<probSubstition>    & subVec);
void readDNABaseFreq(const string filename, alleleFrequency & dnaDefaultFreq);

void readIlluminaError(const string errFile,probSubstition & illuminaErrorsProb);


inline long double lnfact (unsigned int n){
    long double toreturn=0.0;
    for(unsigned int i=1;i<=n;i++){
	toreturn +=  logl( (long double)i );
    }
    return toreturn;
}

inline void initMiscFuncVariables(){
    for(double i=0;i<50;i++){
    	lnFactVec.push_back(lnfact(i) );
     }    
}

//taken from gsl/randist/multinomial.c
inline long double gsl_ran_multinomial_lnpdf (const size_t K, const long double p[], const unsigned int n[]){
  size_t k;
  unsigned int N      =   0;
  long double log_pdf = 0.0;
  long double norm    = 0.0;

  for (k = 0; k < K; k++){
      N += n[k];
  }
  
  for (k = 0; k < K; k++){
      norm += p[k];
  }
  
  log_pdf = lnfact (N);

  for (k = 0; k < K; k++){
      log_pdf -= lnfact (n[k]);
  }

  for (k = 0; k < K; k++){
      log_pdf += logl (p[k] / norm) * n[k];
  }

  return log_pdf;
}

inline long double gsl_ran_multinomial_pdf (const size_t K,const long double p[], const unsigned int n[]){
    return expl (gsl_ran_multinomial_lnpdf (K, p, n));
}


void populatedCoverateVector(      const string  programName       ,  vector<long double> * cov2ProbSite, long double rateForPoissonCov, int maxcov);
void populatedCoverateVectorSingle(const string  directoryProgram  ,  vector<long double> * cov2ProbSite, long double rateForPoissonCov, int maxcov);
/* void readMTConsensus(const string consensusFile,map<int, PHREDgeno> & pos2phredgeno,int & sizeGenome,vector<int> & posOfIndels); */
/* void readMTAlleleFreq(const string freqFile,	map<int, alleleFrequency> & pos2allelefreq); */

//! Simple method to read a bed file 
/*!
 *
 * This method does not check for the records being ordered
  \param filetoread : String with the full path to the file to read
  \return           : Return(stack) the  value of a vector of GenomicRange objects
  \sa  readBEDSortedfile()
*/

inline vector<GenomicRange> readBEDfile(string filetoread){
    vector<GenomicRange> toReturn;
    ifstream myFile;
    myFile.open(filetoread.c_str(), ios::in);
    string line;

    if (myFile.is_open()){
	while ( getline (myFile,line)){	    
	    string chrName;
	    uint64_t startCoord;
	    uint64_t endCoord;
	    vector<string> temp=allTokens(line,'\t');
	    if(temp.size() < 3){
		cerr << "Error in readBEDfile(): the following line does not have at least 3 fields"<<line<<endl;
		exit(1);		
	    }

	    chrName     = destringify<string>(temp[0]);
	    startCoord  = destringify<uint64_t>(temp[1])+1; //the left coordinate is zero based
	    endCoord    = destringify<uint64_t>(temp[2]);
	    GenomicRange toadd (chrName,startCoord,endCoord);
	    toReturn.push_back(toadd);
	}
	myFile.close();
    }else{
	cerr << "Error in readBEDfile(): Unable to open file "<<filetoread<<endl;
	exit(1);
    }

    return toReturn;
}

/* long double computeLL(const int                   al1Current    , */
/* 		      const int                   al2Current    ,		       */
/* 		      const vector<int>         & obsBase       , */
/* 		      const vector<long double> & probDeam1to2  , //rate of deamination from al1 to al2 */
/* 		      const vector<long double> & probDeam2to1  , //rate of deamination from al2 to al1 */
/* 		      const vector<int>         & obsQual       , */
/* 		      const long double           contRate      , */
/* 		      const int                   alContCurrent , */
/* 		      const vector<long double> & mismappingProb); */

// Returns log( exp(x)+exp(y) ), but does so without causing
// overflow or loss of precision.
/* inline double oplusl( long double x, long double y ){ */
/*     return x > y  */
/*         ? x + log1pl( expl(y-x ) )  */
/*         : y + log1pl( expl(x-y ) ) ; */
/* } */

// Returns log10( pow(10,x)+pow(10,y) ), but does so without causing
// overflow or loss of precision.
/* template <typename T> */
/* inline T oplusInit(T x,T y ){ */
/*     if( x == 0 ){ //no initialized, a log = 0 should not exist */
/* 	return y; */
/*     } */

/*     return x > y  */
/*         ? x + log1p( pow( 10, y-x ) ) / log(10) */
/*         : y + log1p( pow( 10, x-y ) ) / log(10) ; */
/* } */


#endif
