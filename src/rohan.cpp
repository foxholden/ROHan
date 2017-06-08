#include <iostream>
#include <fstream>
#include <queue>
#include <algorithm>   
#include <cfloat>   
//#include <random>

//TODO
// test robustness of deamination
// test estimate of h across various coverage ok at 10X10,000
// add mappability track



#include "api/internal/io/BgzfStream_p.h"
#include <api/BamConstants.h>
#include <api/BamMultiReader.h>
#include <utils/bamtools_fasta.h>
#include <utils/bamtools_options.h>
#include <utils/bamtools_pileup_engine.h>
#include <utils/bamtools_utilities.h>

// #include "api/BamMultiReader.h"
// #include "api/BamReader.h"
// #include "api/BamWriter.h"
// #include "api/BamAux.h"

#include "GenomicWindows.h"

#include "miscfunc.h"
#include "utils.h"

using namespace std;
using namespace BamTools;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define DEBUGCOV
//#define DEBUGILLUMINAFREQ
//#define DEBUGINITSCORES
//#define DEBUGDEFAULTFREQ
//#define DEBUGDEAM
//#define DEBUGHCOMPUTE
//#define DEBUGCOMPUTELL
//#define HETVERBOSE
// #define COVERAGETVERBOSE
#define DUMPTRIALLELIC //hack to remove tri-allelic, we need to account for them

#define MAXLENGTHFRAGMENT 500     // maximal length for fragment
#define MAXMAPPINGQUAL    257     // maximal mapping quality, should be sufficient as mapping qualities are encoded using 8 bits
#define MAXCOV             50     // maximal coverage

char offsetQual=33;

long double likeMatch        [MAXMAPPINGQUAL];
long double likeMismatch     [MAXMAPPINGQUAL];

long double likeMatchProb    [MAXMAPPINGQUAL];
long double likeMismatchProb [MAXMAPPINGQUAL];

long double likeMatchMap        [MAXMAPPINGQUAL];
long double likeMismatchMap     [MAXMAPPINGQUAL];

long double likeMatchProbMap    [MAXMAPPINGQUAL];
long double likeMismatchProbMap [MAXMAPPINGQUAL];


vector< vector<long double> > binomVec (MAXCOV+1,vector<long double>(MAXCOV+1,0)) ;
unsigned int totalBasesSum;
unsigned int totalSitesSum;

//Substitution rates due to deamination
vector<probSubstition> sub5p;
vector<probSubstition> sub3p;
probSubstition defaultSubMatch;
probSubstition illuminaErrorsProb;

//default basepairs frequencies
alleleFrequency dnaDefaultBases;
vector<alleleFrequency> defaultDNA5p;
vector<alleleFrequency> defaultDNA3p;


vector<long double> cov2probPoisson;

long double contrate=0.0;
long double rateForPoissonCov;
long double pdfRateForPoissonCov;

// // Returns logl( expl(x)+expl(y) )
// inline long double oplusl(long double x, long double y ){
//     return x > y 
//         ? x + log1pl( expl( y-x ) )
//         : y + log1pl( expl( x-y ) )  ;
// }


long double randomPMismatch4Bases =  ( (long double)(3) ) /  ( (long double)(4) );  // 3/4
long double randomPMatch4Bases    =  1.0-randomPMismatch4Bases;                     // 1/4

long double randomPMmToABase      =  ( (long double)(1) ) /  ( (long double)(3) );  // 1/3, prob of a mismatch to given base

string genoIdx2Code [10] = {"AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"};

// 0	0	0
// 4	1	1
// 7	2	2
// 9	3	3
vector<int> genoPriority (10,0);


#include "DataChunk.h"
#include "DataToWrite.h"
#include "GenoResults.h"






















int    timeThreadSleep =    10;
int    timeSleepWrite  =    1;

bool      readDataDone = false;
unsigned int sizeChunk =  1000000;

string                                                             bamFileToOpen;
queue< DataChunk * >                                               queueDataToprocess;
queue< DataChunk * >                                               queueDataForCoverage;

priority_queue<DataToWrite *, vector<DataToWrite *>, CompareDataToWrite> queueDataTowrite;

pthread_mutex_t  mutexQueue   = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  mutexCounter = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  mutexRank    = PTHREAD_MUTEX_INITIALIZER;

//GLOBALLY accessed
map<unsigned int, int>       threadID2Rank;


//! Chunk of code to check if a certain thread call failed
/*!
  This block is calls by the pthread

*/				
#define checkResults(string, val) {             \
 if (val) {                                     \
     cerr<<"Failed with "<<val<<" at "<<string<<endl;	\
   exit(1);                                     \
 }                                              \
}
 

//#define debugInitscores


// typedef struct{
//     long double ll;
//     long double expal1; //expectation of # of allele 1
//     long double expal2; //expectation of # of allele 2

// } computeLLRes;


// //! A method to compute the pdf of poisson dist.
// /*!
//   This method computes the probability density function for a poisson dist.
// */
// long double pdfPoisson(const long double l,const long double k ) {
//     return expl(k*logl(l)-lgammal(k+1.0)-l);
// }


//! A method to initialize various probability scores to avoid recomputation
/*!
  This method is called by the main after capturing the arguments
*/
void initScores(){
    totalBasesSum=0;
    totalSitesSum=0;


    //homoz
    genoPriority[0] = 0;
    genoPriority[1] = 4;
    genoPriority[2] = 7;
    genoPriority[3] = 9;

    //hetero
    genoPriority[4] = 1;
    genoPriority[5] = 2;
    genoPriority[6] = 3;
    genoPriority[7] = 5;
    genoPriority[8] = 6;
    genoPriority[9] = 8;


    for(int i=0;i<2;i++){
        likeMatch[i]          = log1pl(    -randomPMatch4Bases );          
        likeMismatch[i]       = logl  (     randomPMismatch4Bases );

        likeMatchProb[i]              =    randomPMatch4Bases;    // 1/4
        likeMismatchProb[i]           =    randomPMismatch4Bases; //3/4
    }

    for(int i=2;i<MAXMAPPINGQUAL;i++){
        likeMatch[i]          = log1pl(    -pow(10.0,i/-10.0) );          
        likeMismatch[i]       = logl  (     pow(10.0,i/-10.0) );

        likeMatchProb[i]              = 1.0-pow(10.0,i/-10.0);
        likeMismatchProb[i]           =     pow(10.0,i/-10.0);
    }


    for(int i=0;i<MAXMAPPINGQUAL;i++){
        likeMatchMap[i]        = log1pl(    -pow(10.0,i/-10.0) );          
        likeMismatchMap[i]     = logl  (     pow(10.0,i/-10.0) );

        likeMatchProbMap[i]           = 1.0-pow(10.0,i/-10.0);
        likeMismatchProbMap[i]        =     pow(10.0,i/-10.0);
    }



    for(int i=1;i<=MAXCOV;i++){
	//cout<<i<<endl;

	for(int j=0;j<=i;j++){	    
	    binomVec[i][j] = ( logl(nChoosek(i,j))+logl(powl(0.5,i)) );	     
	}

    }

#ifdef DEBUGINITSCORES
    cout<<"q= "<<"QUAL"<<"\t"<<"likeMatch"<<"\t"<<"likeMismatch"<<"\t"<<"likeMatchProb"<<"\t"<<"likeMismatchProb"<<endl;
    for(int i=0;i<MAXMAPPINGQUAL;i++){
	cout<<"q= "<<i<<"\t"<<likeMatch[i]<<"\t"<<likeMismatch[i]<<"\t"<<likeMatchProb[i]<<"\t"<<likeMismatchProb[i]<<endl;
    }
    exit(1);
#endif

}//end initScores





//! A method to initialize the deamination probabilities
/*!
  This method is called by the main 
*/
void initDeamProbabilities(const string & deam5pfreqE,const string & deam3pfreqE){


    vector<substitutionRates> sub5pT;
    vector<substitutionRates> sub3pT;

    readNucSubstitionRatesFreq(deam5pfreqE,sub5pT);
    readNucSubstitionRatesFreq(deam3pfreqE,sub3pT);;

    //5'
    for(unsigned int i=0;i<sub5pT.size();i++){
	probSubstition toadd;
	for(int nuc1=0;nuc1<4;nuc1++){
	    double probIdentical=1.0;

	    for(int nuc2=0;nuc2<4;nuc2++){
		int nuc = nuc1*4+nuc2;
		if(nuc1==nuc2) continue;
		int ind2 = dimer2indexInt( nuc1,nuc2  );
		probIdentical = probIdentical-sub5pT[i].s[ ind2 ];
		toadd.s[nuc]  =               sub5pT[i].s[ ind2 ];
	    }
	    
	    if(probIdentical<0){
		cerr<<"Error with deamination profile, identity probability is less than 0"<<endl;
		exit(1);
	    }

	    int nuc = nuc1*4+nuc1;
	    toadd.s[nuc] = probIdentical;	    	    	    
	}
	sub5p.push_back(toadd);
    }

    //copying to the rest of the fragment the last position
    for(unsigned int i=(sub5p.size()-1);i<MAXLENGTHFRAGMENT;i++){     
	sub5p.push_back( sub5p[sub5p.size()-1] );
    }
    


    //3'
    for(unsigned int i=0;i<sub3pT.size();i++){
	probSubstition toadd;
	for(int nuc1=0;nuc1<4;nuc1++){
	    double probIdentical=1.0;

	    for(int nuc2=0;nuc2<4;nuc2++){
		int nuc = nuc1*4+nuc2;
		if(nuc1==nuc2) continue;
		int ind2 = dimer2indexInt( nuc1,nuc2  );
		probIdentical = probIdentical-sub3pT[i].s[ ind2 ];
		toadd.s[nuc]  =               sub3pT[i].s[ ind2 ];
	    }
	    
	    if(probIdentical<0){
		cerr<<"Error with deamination profile, identity probability is less than 0"<<endl;
		exit(1);
	    }

	    int nuc = nuc1*4+nuc1;
	    toadd.s[nuc] = probIdentical;	    	    	    
	}
	sub3p.push_back(toadd);
    }

    //copying to the rest of the fragment the last position
    for(unsigned int i=(sub3p.size()-1);i<MAXLENGTHFRAGMENT;i++){
	sub3p.push_back( sub3p[sub3p.size()-1] );
    }



#ifdef DEBUGDEAM
    cerr<<"-- 5' deamination rates --"<<endl;
    for(unsigned int i=0;i<sub5p.size();i++){
    	cerr<<"i="<<i<<" - ";
    	for(int nuc1=0;nuc1<4;nuc1++){
    	    for(int nuc2=0;nuc2<4;nuc2++){
    		int nuc = nuc1*4+nuc2;
    		//cout<<nuc1<<"\t"<<nuc2<<"\t"
    		cerr<<sub5p[i].s[nuc]<<" ";
    	    }
    	    cerr<<" - ";
    	}
    	cerr<<endl;
    }

    cerr<<"-- 3' deamination rates --"<<endl;
    for(unsigned int i=0;i<sub3p.size();i++){
    	cerr<<"i="<<i<<" - ";
    	for(int nuc1=0;nuc1<4;nuc1++){
    	    for(int nuc2=0;nuc2<4;nuc2++){
    		int nuc = nuc1*4+nuc2;
    		//cout<<nuc1<<"\t"<<nuc2<<"\t"
    		cerr<<sub3p[i].s[nuc]<<" ";
    	    }
    	    cerr<<" - ";
    	}
    	cerr<<endl;
    }
#endif




    //if no deamination, cannot have a mismatch
    for(int nuc1=0;nuc1<4;nuc1++){
	for(int nuc2=0;nuc2<4;nuc2++){
	    int nuc = nuc1*4+nuc2;
	    if(nuc1==nuc2)
		defaultSubMatch.s[ nuc ] = 1.0;	
	    else
		defaultSubMatch.s[ nuc ] = 0.0;	
	}
    }



}//end initDeamProbabilities







//! A method to initialize the probability of seeing a base if a fragment is mismapped
/*!
  This method is called by the main after reading the deamination profiles
  as deamination will bias this.
*/
void initDefaultBaseFreq(const string & dnafreqFile){

    readDNABaseFreq( dnafreqFile  ,  dnaDefaultBases );

    for(unsigned int i=0;i<MAXLENGTHFRAGMENT;i++){     
	//defaultDNA5p
	alleleFrequency dnaFreq5p_;
	alleleFrequency dnaFreq3p_;

	for(int b=0;b<4;b++){//original base
	    dnaFreq5p_.f[b]=0;
	    dnaFreq3p_.f[b]=0;
	}

	for(int b1=0;b1<4;b1++){//original base

	    for(int b2=0;b2<4;b2++){
		int b = b2*4+b1;
		//cout<<"ACGT"[b1]<<"\t"<<"ACGT"[b2]<<"\t"<<b<<"\t"<<dnaDefaultBases.f[b1]<<"\t"<<sub5p[i].s[b]<<"\t"<<sub3p[i].s[b]<<endl;
		dnaFreq5p_.f[b1] += dnaDefaultBases.f[b1]*sub5p[i].s[b];
		dnaFreq3p_.f[b1] += dnaDefaultBases.f[b1]*sub3p[i].s[b];
	    }
	}
	
	defaultDNA5p.push_back(dnaFreq5p_);
	defaultDNA3p.push_back(dnaFreq3p_);	
    }

#ifdef DEBUGDEFAULTFREQ

    cerr<<"-- Defautl base frequencies --"<<endl;
    cerr<<dnaDefaultBases.f[0]<<"\t"<<dnaDefaultBases.f[1]<<"\t"<<dnaDefaultBases.f[2]<<"\t"<<dnaDefaultBases.f[3]<<endl;

    cerr<<"-- 5' base frequencies --"<<endl;

    for(unsigned int i=0;i<MAXLENGTHFRAGMENT;i++){     
    	cerr<<"i="<<i<<" - ";
	cerr<<defaultDNA5p[i].f[0]<<"\t"<<defaultDNA5p[i].f[1]<<"\t"<<defaultDNA5p[i].f[2]<<"\t"<<defaultDNA5p[i].f[3]<<endl;
    }

    cerr<<"-- 3' base frequencies --"<<endl;

    for(unsigned int i=0;i<MAXLENGTHFRAGMENT;i++){     
    	cerr<<"i="<<i<<" - ";
	cerr<<defaultDNA3p[i].f[0]<<"\t"<<defaultDNA3p[i].f[1]<<"\t"<<defaultDNA3p[i].f[2]<<"\t"<<defaultDNA3p[i].f[3]<<endl;
    }


#endif


}//end initDeamProbabilities





//! A method to initialize the likelihood of observing data to avoid recomputation
/*!
  This method is called by the main after calling initDefaultBaseFreq
*/
void initLikelihoodScores(){

    //we need 4 theoretical bases X 4 observed bases X MAXMAPPINGQUAL mapping quality X MAXLENGTHFRAGMENT positions on the fragment X base quality

    vector<diNucleotideProb> baseQualEffect;
    for(int i=0;i<MAXMAPPINGQUAL;i++){
	//we know the probability of P(b1->b2) after deamination is 
    }
    
}




//! A method to compute the prob. of seeing an observed base given an original allele
/*!
  This method is called by the computeLL to compute the prob. of seeing observed base ob 
  with error rate q, a of deamination probDeam and original allele al. This method is called 
  for each base.
*/
inline long double computeBaseAl2Obs(const int al,
				     const int ob,
				     const int q,
				     const probSubstition * probDeam      , //rate of deamination
				     const bool isRev,
				     const long double mismappingProb){
    int dinucal2al=-1;

    if( isRev ){		    
	dinucal2al =     complementInt(al)*4+complementInt(al);     //genotype is 1, observed is obsBase
    }else{
	dinucal2al =                   al*4+               al;       //genotype is 1, observed is obsBase
    }

    if(probDeam->s[dinucal2al] == 1.0){ //simple case, no deamination
	int dinucIndexal2ob;


	if( isRev ){
	    dinucIndexal2ob =     complementInt(al)*4+complementInt(ob);     //genotype is al, observed is ob
	}else{
	    dinucIndexal2ob =                   al *4+              ob;      //genotype is al, observed is ob
	}
                                                            
	
	return ( (1.0-mismappingProb)*(
				       likeMatchProb[    q ]*defaultSubMatch.s[   dinucIndexal2ob]
				       + 
				       likeMismatchProb[ q ]*illuminaErrorsProb.s[dinucIndexal2ob] 
				       )
		 +
		 mismappingProb*randomPMatch4Bases);

    }else{
	
	long double sumProbToReturn = 0.0;
	for(int alpostdeam=0;alpostdeam<4;alpostdeam++){
	    int dinucal2ald = -1;
	    if( isRev ){		    
		dinucal2ald =     complementInt(al)*4+complementInt(alpostdeam);     //genotype is 1, observed is obsBase
	    }else{
		dinucal2ald =                   al*4+               alpostdeam;       //genotype is 1, observed is obsBase
	    }
	    

	    if(probDeam->s[dinucal2ald] > 0.0){ //has deamination       	
		int dinucald2ob = -1;
		if( isRev ){		    
		    dinucald2ob =     complementInt(alpostdeam)*4+complementInt(ob);     //genotype is 1, observed is obsBase
		}else{
		    dinucald2ob =                   alpostdeam*4+               ob;       //genotype is 1, observed is obsBase
		}

		sumProbToReturn += probDeam->s[dinucal2ald]*( (1.0-mismappingProb)*(
										    likeMatchProb[    q ]*defaultSubMatch.s[   dinucald2ob]
										    + 
										    likeMismatchProb[ q ]*illuminaErrorsProb.s[dinucald2ob] 
										    )
							      +
							      mismappingProb*randomPMatch4Bases);
	    }
	}

	return sumProbToReturn;
    }

    return -1;
}//end computeBaseAl2Obs



//! A method to compute the prob. of seeing observed ob for a given genotype al1Current and al1Current
/*!
  This method is called by class heteroComputerVisitor : public PileupVisitor, computes the probability of
  observing the set of bases obsBase with qualities obsQual. 
*/
inline long double computeLL(const int                           al1Current    ,
			     const int                           al2Current    ,		      
			     const vector<int>                 & obsBase       ,
			     const vector<int>                 & obsQual       ,
			     // const vector<long double> & probDeam1to2  , //rate of deamination from al1 to al2
			     // const vector<long double> & probDeam2to1  , //rate of deamination from al2 to al1
			     const vector<probSubstition *>    & probDeam      , //rate of deamination from al2 to al1
			     const vector<bool>                & isRev         ,//false = plus strand, true = reverse		      
			     const long double                   contRate      ,
			     const int                           alCCurrent ,
			     const vector<long double>         & mismappingProb
			     ){
    //computeLLRes toreturn;
    long double llik=0;
    long double llik1=0;
    long double llik2=0;
    long double llikC=0;
    
#ifdef DEBUGCOMPUTELL
    cout<<endl<<"al1="<<al1Current<<"\tal2="<<al2Current<<endl;
#endif

    for(int i=0;i<int(obsBase.size());i++){//iterating over each observed base

#ifdef DEBUGCOMPUTELL
	cout<<i<<"\tob="<<obsBase[i]<<"\ta1="<<al1Current<<"\ta2="<<al2Current<<"\talc"<<alCCurrent<<endl;
#endif


	llikC = computeBaseAl2Obs(alCCurrent  ,
				  obsBase[i]  ,
				  obsQual[i]  ,
				  &defaultSubMatch,//no deamination for the contaminant
				  isRev[i]    ,
				  mismappingProb[i]);
	



	
	long double llikAl1t=0;
	long double llikAl2t=0;
	
	


	
	llikAl1t = computeBaseAl2Obs(al1Current   ,
				     obsBase[i]   ,
				     obsQual[i]   ,
				     probDeam[i]  , //deamination for the endogenous
				     isRev[i]     ,
				     mismappingProb[i]);



	llikAl2t = computeBaseAl2Obs(al2Current   ,
				     obsBase[i]   ,
				     obsQual[i]   ,
				     probDeam[i] , //deamination for the endogenous
				     isRev[i]     ,
				     mismappingProb[i]);



#ifdef DEBUGCOMPUTELL
	// cout<<i<<"\tdeam1to2="<<(probSubDeam1toObs)<<endl;
	// cout<<i<<"\tdeam2to1="<<(probSubDeam2toObs)<<endl;

	cout<<i<<"\tlm="<<likeMatchProb[obsQual[i]]<<"\tlmm="<<likeMismatchProb[obsQual[i]]<<"\tlmism="<<mismappingProb[i]<<endl;
	cout<<i<<"\tla1="<<llikAl1t <<"\tla2="<<llikAl2t<<endl;
#endif
	// exit(1);

	long double llikTE   = (llikAl1t + llikAl2t)/2.0 ;   //endogenous likelihood  0.5*al1 + 0.5*al2



	long double llikT   = (1.0-contRate)*(  llikTE ) + (contRate)*llikC  ;	
	llik               += logl(llikT);

 	llik1=oplusInitnatl( llik1, logl(llikAl1t) );
	llik2=oplusInitnatl( llik2, logl(llikAl2t) );

#ifdef DEBUGCOMPUTELL
	cout<<i<<"\tLLa1="<<llikAl1t <<"\tLLa2="<<llikAl2t<<"\tLLC="<<llikC<<"\tLLE="<<llikTE<<"\tLLT="<<llikT<<"\tlogLLT="<<logl(llikT)<<endl;
#endif

    }
	//cout<<"CC\t"<<llikCC<<"\t"<<llikCC1<<"\t"<<llikCC2<<endl;    
#ifdef DEBUGCOMPUTELL
    cout<<"ll1="<<llik1<<"\tll2="<<llik2<<"\tp1="<<expl(llik1)<<"\tp2="<<expl(llik2)<<endl;
#endif

    //BEGIN BINOMIAL COEFFICIENT
    long double expal1  = ((long double)(obsBase.size())) * ( expl(llik1) / expl(oplusnatl(llik1,llik2))) ;
    int expal_f = floorl(expal1);
    int expal_c =  ceill(expal1);

    long double expal_fb= binomVec[int(obsBase.size())][expal_f];
    long double expal_cb= binomVec[int(obsBase.size())][expal_c];

    long double expal_fr = 1-(expal1-expal_f);
    long double expal_cr = 1-(expal_c-expal1);
    //END BINOMIAL COEFFICIENT
    // if(expal_f == expal_c){
    // 	expal_fr = 1.0;
    // 	expal_cr = 0.0;
    // }

#ifdef DEBUGCOMPUTELL
    cout<<"expal1="<<expal1<<endl;
    cout<<"expal_f="<<expal_f<<"\texpal_fb="<<expal_fb<<"\texpal_fr="<<expal_fr<<endl;
    cout<<"expal_c="<<expal_c<<"\texpal_cb="<<expal_cb<<"\texpal_cr="<<expal_cr<<endl;
#endif

    //long double expal1=roundl(  sizeAr * ( expl(llik1) / llik1+llik2)) );
    //long double expal2  = sizeAr-expal1;
    long double binomE2 = expal_fr*expal_fb+expal_cr*expal_cb;
    
#ifdef DEBUGCOMPUTELL
    cout<<"a1="<<al1Current<<" a2="<<al2Current<<"\tll="<<llik<<"\tE[1]="<<expal1<<"\tE[2]="<<(double(obsBase.size())-expal1)<<endl;
    cout<<"b="<<binomE2<<"\tb="<<expl(binomE2)<<"\tb+l="<<(binomE2+llik)<<endl;
#endif
    
    // toreturn.ll     = (binomE2+llik);
    // toreturn.expal1 = expal1;
    // toreturn.expal2 = expal2;
    
    return (binomE2+llik);
} //end computeLL


class coverageComputeVisitor : public PileupVisitor {
  
public:
    coverageComputeVisitor(const RefVector& references,unsigned int leftCoord, unsigned int rightCoord)
	: PileupVisitor()
	, m_references(references)
	, m_leftCoord(leftCoord)
	, m_rightCoord(rightCoord)
    { 
	totalBases=0;
	totalSites=0;
	
    }
    ~coverageComputeVisitor(void) {}
  
    // PileupVisitor interface implementation

    
    void Visit(const PileupPosition& pileupData) {   
	//bool foundOneFragment=false;
	if(pileupData.Position < int(m_leftCoord)   || 
	   pileupData.Position > int(m_rightCoord) ){
	    return ;
	}
	//cout<<m_leftCoord<<"\t"<<m_rightCoord<<"\t"<<pileupData.Position<<endl;

	for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){
	    if( pileupData.PileupAlignments[i].IsCurrentDeletion &&
	    	pileupData.PileupAlignments[i].IsNextInsertion ){
	    	continue;
	    }
	    //foundOneFragment=true;		  
	    totalBases++;
	}

	//if(foundOneFragment)
	totalSites++;
    }
    
    unsigned int getTotalBases() const{
    	return totalBases;
    }

    unsigned int getTotalSites() const{
    	return totalSites;
    }

private:
    RefVector m_references;
    //Fasta * m_fastaReference;
    unsigned int totalBases;
    unsigned int totalSites;
    unsigned int m_leftCoord;
    unsigned int m_rightCoord;
    
};//end coverageComputeVisitor




class heteroComputerVisitor : public PileupVisitor {
  
public:
    heteroComputerVisitor(const RefVector& references, 
			  const int refID,
			  const unsigned int leftCoord,
			  const unsigned int rightCoord,
			  const long double contRate,
			  vector<PositionResult *> * dataToWriteOut)
	: PileupVisitor()
	, m_references(references)
	, m_refID(refID)
	, m_leftCoord(leftCoord)
	, m_rightCoord(rightCoord)
	, m_contRate(contRate)
	, m_dataToWriteOut( dataToWriteOut)
    { 
    }
    ~heteroComputerVisitor(void) { }
  
    // PileupVisitor interface implementation

    
    void Visit(const PileupPosition& pileupData) {   


	if(pileupData.Position < int(m_leftCoord)   || 
	   pileupData.Position > int(m_rightCoord) ){
	    return ;
	}

	int                 totalBases=0 ;
	int                 counterB  [4];
	//long double         llBaseDeam[4];
	vector<int>              obsBase      ;
	vector<int>              obsQual      ;
	vector<long double>      mmProb       ; //mismapping probability
	vector<probSubstition *> substitutionRatesPerRead;
	vector<bool>             isRevVec;
	// vector<probSubstition *> substitutionRatesPerRead    (pileupData.PileupAlignments.size(),NULL );
	// vector<bool>             isRevVec                    (pileupData.PileupAlignments.size(),false);
	
	//reserve
	//cout<<"init1 "<<obsBase.size()<<"\t"<<substitutionRatesPerRead.size()<<"\t"<<pileupData.PileupAlignments.size()<<endl;
	//TODO Debub resize
	// obsBase.resize(pileupData.PileupAlignments.size());
	// obsQual.resize(pileupData.PileupAlignments.size());
	// mmProb.resize(pileupData.PileupAlignments.size());
	// substitutionRatesPerRead.resize(pileupData.PileupAlignments.size());
	// isRevVec.resize(pileupData.PileupAlignments.size());
	// cout<<"init2 "<<obsBase.size()<<"\t"<<substitutionRatesPerRead.size()<<"\t"<<pileupData.PileupAlignments.size()<<endl;


	unsigned int                posAlign = pileupData.Position+1;

	for(unsigned int i=0;i<4;i++){
	    counterB[i]   = 0;
	    //llBaseDeam[i] = 0.0;
	}
	//vector<bool> includeFragment (pileupData.PileupAlignments.size(),false);

	for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){


	    if( pileupData.PileupAlignments[i].IsCurrentDeletion   ||
	    	pileupData.PileupAlignments[i].IsNextInsertion     ||
	    	pileupData.PileupAlignments[i].IsNextDeletion      ||
		(pileupData.PileupAlignments[i].DeletionLength>0)  ||
		(pileupData.PileupAlignments[i].InsertionLength>0) ){		
		//includeFragment was initialized as false
	    	continue;
	    }


	    if(i>=MAXCOV){
		break;
	    }

	    char  b   =     pileupData.PileupAlignments[i].Alignment.QueryBases[ pileupData.PileupAlignments[i].PositionInAlignment ];
	    if(!isResolvedDNA(b)){ 
		continue; 
	    }//avoid Ns
	    int bIndex = baseResolved2int(b);
	    int   q    = int(pileupData.PileupAlignments[i].Alignment.Qualities[  pileupData.PileupAlignments[i].PositionInAlignment ]-offsetQual); 
	    int   m    = int(pileupData.PileupAlignments[i].Alignment.MapQuality);
	  

	    // BEGIN DEAMINATION COMPUTATION
            //zero base distance to the 5p/3p end
            int dist5p=-1;
            int dist3p=-1;

	    bool isRev = pileupData.PileupAlignments[i].Alignment.IsReverseStrand();
	    //isRevVec[i]=isRev;
	    
            if( isRev ){
                dist5p = pileupData.PileupAlignments[i].Alignment.QueryBases.size() - pileupData.PileupAlignments[i].PositionInAlignment-1;
                dist3p = pileupData.PileupAlignments[i].PositionInAlignment;
            }else{
                dist5p = pileupData.PileupAlignments[i].PositionInAlignment;
                dist3p = pileupData.PileupAlignments[i].Alignment.QueryBases.size() - pileupData.PileupAlignments[i].PositionInAlignment-1;
            }
	    
	    // cout<<"5p="<<dist5p<<" 3p="<<dist3p<<endl;

            // probSubstition * probSubMatchToUseEndo = &defaultSubMatch ;
            // probSubstition * probSubMatchToUseCont = &defaultSubMatch ;
            probSubstition * probSubMatchToUseEndo = &defaultSubMatch ;
            // substitutionRates * probSubMatchToUseCont = &defaultSubMatch ;
	    

            if(dist5p <= (int(sub5p.size()) -1)){ //out of range
                probSubMatchToUseEndo = &sub5p[  dist5p ];                      
            }else{
		//TODO what to do in between
		//probSubMatchToUseEndo = &sub5p[  int(sub5p.size()) -1 ];                      
	    }

            if(dist3p <= (int(sub3p.size()) -1)){ // out of range
                probSubMatchToUseEndo = &sub3p[  dist3p ];
            }else{
		//TODO what to do in between
		//probSubMatchToUseEndo = &sub3p[ int(sub3p.size()) -1 ];
	    }

            //we have substitution probabilities for both... take the closest
            if(dist5p <= (int(sub5p.size()) -1) &&
               dist3p <= (int(sub3p.size()) -1) ){
                    
                if(dist5p < dist3p){
                    probSubMatchToUseEndo = &sub5p[  dist5p ];
                }else{
                    probSubMatchToUseEndo = &sub3p[  dist3p ];
                }
                    
            }
	    //cout<<"add= "<<i<<"\t"<<bIndex<<"\t"<<probSubMatchToUseEndo<<"\t"<<obsBase.size()<<"\t"<<substitutionRatesPerRead.size()<<endl;
	    //substitutionRatesPerRead[i] = probSubMatchToUseEndo;

  
	    //cout<<"pos "<<posAlign<<" "<<bIndex<<" "<<b<<" "<<pileupData.PileupAlignments[i].Alignment.Name<<endl;
	    counterB[ bIndex ]++;
	    totalBases++;

	    obsBase.push_back(             bIndex   );
	    obsQual.push_back(                  q   );
	    mmProb.push_back(  likeMismatchProbMap[m]  );
	    substitutionRatesPerRead.push_back( probSubMatchToUseEndo );
	    isRevVec.push_back(isRev);
	    //	    includeFragment[i]=true;

	}//END FOR EACH READ

	if(totalBases == 0) 
	    return;

	PositionResult * prToAdd=new PositionResult();
	prToAdd->refID = m_refID;
	prToAdd->pos   = posAlign;

	int counterUnique=0;
	for(int i=0;i<4;i++){
	    prToAdd->baseC[i] = counterB[i];
	    if(counterB[i]!=0) 
		counterUnique++;

	}


	if(m_contRate == 0){
	    int alc=randomInt(0,3);//random cont base
	    int genoIdx=0;
	    for(int al1=0;al1<4;al1++){
		for(int al2=0;al2<4;al2++){
		    if(al2<al1)
			continue;

		    //cout<<"ACGT"[al1]<<"\t"<<"ACGT"[al2]<<"\t"<<counterB[al1]<<"\t"<<counterB[al2]<<endl;
		    
		    
		    //TO REMOVE
		    // if(obsBase.size() != substitutionRatesPerRead.size() ){
		    // 	cerr<<"problem2 "<<obsBase.size()<<" "<<substitutionRatesPerRead.size()<<endl;
		    // 	cerr<<prToAdd->refID<<"\t"<<prToAdd->pos<<endl;
		    // 	exit(1);
		    // }

		    // for(unsigned int i=0;i<substitutionRatesPerRead.size();i++){
		    // 	cout<<"ll "<<i<<"\t"<<obsBase[i]<<"\t"<<substitutionRatesPerRead[i]<<endl;
		    // }


		    long double ll=computeLL(al1         ,//al1
					     al2         ,//al2		      		  
					     obsBase     ,
					     obsQual     ,
					     substitutionRatesPerRead,
					     isRevVec    ,
					     m_contRate  ,
					     alc         ,//cont
					     mmProb      );


		    //cout<<"ACGT"[al1]<<"\t"<<"ACGT"[al2]<<"\t"<<ll<<endl;
		    prToAdd->ll[genoIdx]=ll;
		    genoIdx++;
		}
	    }

	}else{
	    for(int genoIdx=0;genoIdx<10;genoIdx++)
		prToAdd->ll[genoIdx] = 0.0;
	    long double factorScale = logl( ( (long double)1) / ((long double)counterUnique) );

	    for(int alc=0;alc<4;alc++){//for each contaminant
		if(counterB[alc]==0) //just concentrate on bases that are there
		    continue;

				
		int genoIdx=0;
		for(int al1=0;al1<4;al1++){
		    for(int al2=0;al2<4;al2++){
			if(al2<al1)
			    continue;
			
			
			
			long double ll=computeLL(al1         ,//al1
						 al2         ,//al2		      		  
						 obsBase     ,
						 obsQual     ,
						 substitutionRatesPerRead,
						 isRevVec    ,
						 m_contRate  ,
						 alc         ,//cont
						 mmProb      );
						
			prToAdd->ll[genoIdx]=oplusInitnatl(prToAdd->ll[genoIdx],ll+factorScale);

#ifdef DEBUGCOMPUTELL
			cout<<"ACGT"[al1]<<"\t"<<"ACGT"[al2]<<"\t""\t"<<"ACGT"[alc]<<"\t"<<ll<<"\t"<<factorScale<<"\t"<<prToAdd->ll[genoIdx]<<endl<<"--------------------"<<endl;
#endif
			genoIdx++;
		    }
		}

	    }
	}	
	

	





	//TODO

	// char refB="ACGT"[ref];
	// char altB="ACGT"[alt];
	// PositionResult * prToAdd=new PositionResult();
	// prToAdd->refID = m_refID;
	// prToAdd->pos   = posAlign;
	// prToAdd->refB  = refB;
	// prToAdd->altB  = altB;
	// prToAdd->refC  = counterB[ref];
	// prToAdd->altC  = counterB[alt];



       


	//normalize
	long double llEvidence=prToAdd->ll[0] ;
	for(int genoIdx=1;genoIdx<10;genoIdx++){
	    llEvidence = oplusnatl(llEvidence,prToAdd->ll[genoIdx] );
	    //cout<<genoIdx<<"\t"<<prToAdd->ll[genoIdx]<<"\t"<<llEvidence<<endl;
	}

	for(int genoIdx=0;genoIdx<10;genoIdx++){
	    prToAdd->ll[genoIdx] = prToAdd->ll[genoIdx] - llEvidence;
	    //cout<<genoIdx<<"\t"<<prToAdd->ll[genoIdx]<<"\t"<<llEvidence<<endl;
	}

	
	//	exit(1);
	//TODO
       // prToAdd->rrll  = oplusnatl( rrllCr+logl(0.5), rrllCa+logl(0.5));
       // prToAdd->rall  = oplusnatl( rallCr+logl(0.5), rallCa+logl(0.5));
       // prToAdd->aall  = oplusnatl( aallCr+logl(0.5), aallCa+logl(0.5));
	vector<long double> arrLL (10,0.0); 

	
	//for(int g=0;g<10;g++){
	for(int genoIdx=0;genoIdx<10;genoIdx++){
	    arrLL[genoIdx] = prToAdd->ll[genoIdx];
	}
	// cout<<"-----"<<endl;
	
	// arrLL[0]       = prToAdd->rrll;
	// arrLL[1]       = prToAdd->rall;
	// arrLL[2]       = prToAdd->aall;
	sort (arrLL.begin(), arrLL.end() , greater<long double>());
	
	// for(int genoIdx=0;genoIdx<10;genoIdx++){
	//     cout<<genoIdx<<"\t"<<arrLL[genoIdx]<<endl;
	// }

	prToAdd->lqual = (arrLL[1]-arrLL[0]);//log ratio of most likely to second most likely

	// //1st most likely = 0
	// //2nd most likely = 1


	

	for(int g=0;g<10;g++){
	    int genoIdx = genoPriority[g];
	    //cout<<genoIdx<<"\t"<<prToAdd->ll[genoIdx]<<"\t"<<arrLL[0]<<endl;
	    if( prToAdd->ll[genoIdx] == arrLL[0]){
		if(genoIdx == 0 ||
		   genoIdx == 4 ||
		   genoIdx == 7 ||
		   genoIdx == 9 ){
		    prToAdd->geno = 0;     //rr
		}else{
		    prToAdd->geno = 1;     //ra
		}

		prToAdd->genoS[0] = genoIdx2Code[genoIdx][0];
		prToAdd->genoS[1] = genoIdx2Code[genoIdx][1];
		//cout<<genoIdx<<"\t"<<prToAdd->geno<<"\t"<<prToAdd->genoS[0]<<prToAdd->genoS[1]<<endl;
		break;
	    }
	}
	// cout<<prToAdd->geno<<endl;


	
	//prToAdd->lqual = -1;
       // if(arrLL[2]     == prToAdd->rrll){
       // 	   prToAdd->geno = 0;     //rr
       // }else{
       // 	   if(arrLL[2] == prToAdd->rall){
       // 	       prToAdd->geno = 1; //ra
       // 	   }else{
       // 	       prToAdd->geno = 2; //aa
       // 	   }
       // }

	//prToAdd->llCov = logl( pdfPoisson( (long double)totalBases, rateForPoissonCov)/pdfRateForPoissonCov );
	//prToAdd->llCov = logl( pdfPoisson( (long double)totalBases, rateForPoissonCov)  );

	m_dataToWriteOut->push_back(prToAdd);
       
    }//end Visit()
    

private:
    RefVector m_references;
    //Fasta * m_fastaReference;
    // unsigned int totalBases;
    // unsigned int totalSites;
    int          m_refID;
    unsigned int m_leftCoord;
    unsigned int m_rightCoord;
    long double  m_contRate;
    vector<PositionResult *> * m_dataToWriteOut;
};//heteroComputerVisitor










//! Method called for each thread
/*!
  

*/				
void *mainHeteroComputationThread(void * argc){

    int   rc;
#ifdef HETVERBOSE    
    int rankThread=0;
#endif
    rc = pthread_mutex_lock(&mutexRank);
    checkResults("pthread_mutex_lock()\n", rc);

    threadID2Rank[*(int *)pthread_self()]  = threadID2Rank.size()+1;

#ifdef HETVERBOSE    
    rankThread = threadID2Rank[*(int *)pthread_self()];
#endif
    
    rc = pthread_mutex_unlock(&mutexRank);
    checkResults("pthread_mutex_unlock()\n", rc);

 checkqueue:    
    // stackIndex=-1;
    //check stack

    
    rc = pthread_mutex_lock(&mutexQueue);
    checkResults("pthread_mutex_lock()\n", rc);


    bool foundData=false;
    
#ifdef HETVERBOSE
    cerr<<"Thread #"<<rankThread <<" started and is requesting data"<<endl;
#endif

    DataChunk    * currentChunk;


    if(!queueDataToprocess.empty()){    
 	foundData=true;
 	currentChunk = queueDataToprocess.front();
 	queueDataToprocess.pop();
#ifdef HETVERBOSE
 	cerr<<"Thread #"<<rankThread<<" is reading "<<currentChunk->rank<<endl;
#endif
	//cout<<"rank "<< &(currentChunk->dataToProcess) <<endl;
    }

    
  

    if(!foundData){
	rc = pthread_mutex_unlock(&mutexQueue);
	checkResults("pthread_mutex_unlock()\n", rc);


	if(readDataDone){
#ifdef HETVERBOSE
	    cerr<<"Thread #"<<rankThread<<" is done"<<endl;
#endif
	    return NULL;	
	}else{
#ifdef HETVERBOSE
	    cerr<<"Thread #"<<rankThread<<" sleeping for "<<timeThreadSleep<<endl;
#endif
	    sleep(timeThreadSleep);
	    goto checkqueue;
	}

    }else{
	//release stack
	rc = pthread_mutex_unlock(&mutexQueue);
	checkResults("pthread_mutex_unlock()\n", rc);
    }

    //////////////////////////////////////////////////////////////
    //                BEGIN COMPUTATION                         //
    //////////////////////////////////////////////////////////////

    //cout<<currentChunk->rangeGen<<endl;

#ifdef HETVERBOSE
    cerr<<"Thread #"<<rankThread<<" is reading "<<currentChunk->rangeGen<<endl;
#endif

    //sleep(10);


    


    BamReader reader;
    if ( !reader.Open(bamFileToOpen) ) {
	cerr << "Could not open input BAM file:" << bamFileToOpen <<endl;
    	exit(1);
    }

    reader.LocateIndex();

    if(!reader.HasIndex()){
    	cerr << "The BAM file: " << bamFileToOpen <<" does not have an index"<<endl;
    	exit(1);
    }

    // retrieve reference data
    const RefVector  references = reader.GetReferenceData();
    const int        refID      = reader.GetReferenceID( currentChunk->rangeGen.getChrName() );
    
    // cerr<<"Thread #"<<rankThread<<" refID "<<refID<<endl;    
    // cerr<<"Thread #"<<rankThread<<" "<<references[0].RefName<<endl;

    BamRegion bregion (refID, 
		       currentChunk->rangeGen.getStartCoord(), 
		       refID, 
		       currentChunk->rangeGen.getEndCoord()   );

    bool setRegionRes=reader.SetRegion( bregion   );

    // cerr<<"Thread #"<<rankThread<<" "<<references[0].RefName<<"\t"<<setRegionRes<<endl;

    if( refID==-1 ||
       !setRegionRes){
    	cerr << "Heterozygous computation: could not set region "<<currentChunk->rangeGen<<" for BAM file:" << bamFileToOpen <<" "<<currentChunk->rangeGen.getStartCoord()<<" "<< currentChunk->rangeGen.getEndCoord()<< endl;
    	exit(1);
    }

    DataToWrite  * dataToWrite = new DataToWrite();

    dataToWrite->rangeGen      =  currentChunk->rangeGen;
    dataToWrite->rank          =  currentChunk->rank;

    //dataToWrite->dataToWriteOut=new vector<PositionResult *>();
    heteroComputerVisitor* cv = new heteroComputerVisitor(references,
							  refID,
							  currentChunk->rangeGen.getStartCoord(), 
							  currentChunk->rangeGen.getEndCoord()  ,
							  contrate,
							  dataToWrite->vecPositionResults);

    

    PileupEngine pileup;
    pileup.AddVisitor(cv);

    BamAlignment al;
    //unsigned int numReads=0;
    while ( reader.GetNextAlignment(al) ) {
        pileup.AddAlignment(al);
    }

    //clean up
    pileup.Flush();
    reader.Close();
    //fastaReference.Close();
    
    //cerr<<"Thread #"<<rankThread <<" "<<cv->getTotalBases()<<"\t"<<cv->getTotalSites()<<"\t"<<double(cv->getTotalBases())/double(cv->getTotalSites())<<endl;

    delete cv;

	
#ifdef HETVERBOSE
    cerr<<"Thread #"<<rankThread <<" is done with computations"<<endl;
#endif

    //////////////////////////////////////////////////////////////
    //                END   COMPUTATION                         //
    //////////////////////////////////////////////////////////////
	

    //COUNTERS
    rc = pthread_mutex_lock(&mutexCounter);
    checkResults("pthread_mutex_lock()\n", rc);
    

    //TODO: EACH THREAD WRITE TO A TEMP FILE 
    //THEN COMBINED
    // std::ifstream if_a("a.txt", std::ios_base::binary);
    // std::ifstream if_b("b.txt", std::ios_base::binary);
    // std::ofstream of_c("c.txt", std::ios_base::binary);    
    // of_c << if_a.rdbuf() << if_b.rdbuf();

    queueDataTowrite.push(dataToWrite);

    rc = pthread_mutex_unlock(&mutexCounter);
    checkResults("pthread_mutex_unlock()\n", rc);

#ifdef HETVERBOSE
    cerr<<"Thread #"<<rankThread <<" is re-starting"<<endl;
#endif
    goto checkqueue;	   


    

#ifdef HETVERBOSE    
    cerr<<"Thread "<<rankThread<<" ended "<<endl;
#endif

    return NULL;

}// end mainHeteroComputationThread






queue< DataChunk * >  randomSubQueue(const queue< DataChunk * > queueDataToSubsample,unsigned int sizeToReturn){

    if( sizeToReturn > queueDataToSubsample.size()){
	cerr<<"Cannot subsample the queue to the size required"<<endl;
	exit(1);
    }

    queue< DataChunk *  >  toReturn=queueDataToSubsample;
    vector< DataChunk * >  myvectortemp;

    while(!toReturn.empty()){    
	DataChunk *dc = toReturn.front();
 	toReturn.pop();
	myvectortemp.push_back(dc);
    }

    random_shuffle ( myvectortemp.begin(), myvectortemp.end() );
    
    for(unsigned int i=0;i<sizeToReturn;i++){
	toReturn.push(myvectortemp[i]);
    }

    return toReturn;
} // end randomSubQueue


//TODO: GC bias for coverage?
void *mainCoverageComputationThread(void * argc){

    int   rc;
#ifdef COVERAGETVERBOSE    
    int rankThread=0;
#endif
    rc = pthread_mutex_lock(&mutexRank);
    checkResults("pthread_mutex_lock()\n", rc);

    threadID2Rank[*(int *)pthread_self()]  = threadID2Rank.size()+1;
#ifdef COVERAGETVERBOSE    
    rankThread = threadID2Rank[*(int *)pthread_self()];
#endif
    
    rc = pthread_mutex_unlock(&mutexRank);
    checkResults("pthread_mutex_unlock()\n", rc);

 checkqueue:    
    // stackIndex=-1;
    //check stack

    
    rc = pthread_mutex_lock(&mutexQueue);
    checkResults("pthread_mutex_lock()\n", rc);


    bool foundData=false;
#ifdef COVERAGETVERBOSE
    cerr<<"Thread coverage #"<<rankThread <<" started and is requesting data"<<endl;
#endif

    DataChunk * currentChunk;


    if(!queueDataForCoverage.empty()){    
 	foundData=true;
 	currentChunk = queueDataForCoverage.front();
 	queueDataForCoverage.pop();
#ifdef COVERAGETVERBOSE
 	cerr<<"Thread #"<<rankThread<<" is reading "<<currentChunk->rank<<endl;
#endif
	//cout<<"rank "<< &(currentChunk->dataToProcess) <<endl;
    }

    
  

    if(!foundData){
	rc = pthread_mutex_unlock(&mutexQueue);
	checkResults("pthread_mutex_unlock()\n", rc);


	if(readDataDone){
#ifdef COVERAGETVERBOSE
	    cerr<<"Thread #"<<rankThread<<" is done"<<endl;
#endif
	    return NULL;	
	}else{
#ifdef COVERAGETVERBOSE
	    cerr<<"Thread #"<<rankThread<<" sleeping for "<<timeThreadSleep<<endl;
#endif
	    sleep(timeThreadSleep);
	    goto checkqueue;
	}

    }else{
	//release stack
	rc = pthread_mutex_unlock(&mutexQueue);
	checkResults("pthread_mutex_unlock()\n", rc);
    }

    //////////////////////////////////////////////////////////////
    //                BEGIN COMPUTATION                         //
    //////////////////////////////////////////////////////////////

    //cout<<currentChunk->rangeGen<<endl;
#ifdef COVERAGETVERBOSE
    cerr<<"Thread #"<<rankThread<<" is reading "<<currentChunk->rangeGen<<endl;
#endif
    //sleep(10);


    


    BamReader reader;
    if ( !reader.Open(bamFileToOpen) ) {
	cerr << "Could not open input BAM file:" << bamFileToOpen <<endl;
    	exit(1);
    }

    reader.LocateIndex();

    if(!reader.HasIndex()){
    	cerr << "The BAM file: " << bamFileToOpen <<" does not have an index"<<endl;
    	exit(1);
    }

    // retrieve reference data
    const RefVector  references = reader.GetReferenceData();
    const int        refID      = reader.GetReferenceID( currentChunk->rangeGen.getChrName() );


#ifdef COVERAGETVERBOSE    
    cerr<<"Thread #"<<rankThread<<" refID "<<refID<<" "<<currentChunk->rangeGen.getStartCoord()<<" "<<currentChunk->rangeGen.getEndCoord()<<endl;    
#endif
    // cerr<<"Thread #"<<rankThread<<" "<<references[0].RefName<<endl;


    BamRegion bregion (refID, 
    		       int(currentChunk->rangeGen.getStartCoord()), 
    		       refID, 
    		       int(currentChunk->rangeGen.getEndCoord() )  );

    bool setRegionRes = reader.SetRegion( bregion   );


    if( refID==-1 ||
       !setRegionRes){	
    	cerr << "Coverage computation: could not set region "<<currentChunk->rangeGen<<" for BAM file:" << bamFileToOpen <<" "<<currentChunk->rangeGen.getStartCoord()<<" "<< currentChunk->rangeGen.getEndCoord()<< "\trefID:"<<refID<<"\tset region fail?:"<<booleanAsString(setRegionRes)<<endl;
    	exit(1);
    }

   
    coverageComputeVisitor* cv = new coverageComputeVisitor(references,
							    currentChunk->rangeGen.getStartCoord(),
							    currentChunk->rangeGen.getEndCoord());
    PileupEngine pileup;
    pileup.AddVisitor(cv);

    BamAlignment al;
    //unsigned int numReads=0;
    while ( reader.GetNextAlignment(al) ) {
        pileup.AddAlignment(al);
    }

    
    //clean up
    pileup.Flush();
    reader.Close();
    //fastaReference.Close();
#ifdef COVERAGETVERBOSE    
    cerr<<"Thread #"<<rankThread <<" "<<cv->getTotalBases()<<"\t"<<cv->getTotalSites()<<"\t"<<double(cv->getTotalBases())/double(cv->getTotalSites())<<endl;
#endif

    unsigned int totalBasesL=cv->getTotalBases();
    unsigned int totalSitesL=cv->getTotalSites();


    delete cv;
	
#ifdef COVERAGETVERBOSE    
    cerr<<"Thread #"<<rankThread <<" is done with computations"<<endl;
#endif
    //////////////////////////////////////////////////////////////
    //                END   COMPUTATION                         //
    //////////////////////////////////////////////////////////////
	

    //COUNTERS
    rc = pthread_mutex_lock(&mutexCounter);
    checkResults("pthread_mutex_lock()\n", rc);
    
//queueDataTowrite.push(currentChunk);
    totalBasesSum+=totalBasesL;
    totalSitesSum+=totalSitesL;

    rc = pthread_mutex_unlock(&mutexCounter);
    checkResults("pthread_mutex_unlock()\n", rc);

#ifdef COVERAGETVERBOSE    
    cerr<<"Thread #"<<rankThread <<" is re-starting"<<endl;
#endif

    goto checkqueue;	   


    

#ifdef COVERAGETVERBOSE        
    cerr<<"Thread "<<rankThread<<" ended "<<endl;
#endif

    return NULL;

} // mainCoverageComputationThread





//! Main method
/*!
  The main:
    calls initScores(), 
    captures the arguments
    reads the deamination profiles
*/

int main (int argc, char *argv[]) {
    setlocale(LC_ALL, "POSIX");



    ////////////////////////////////////
    // BEGIN Parsing arguments        //
    ////////////////////////////////////

    string cwdProg=getCWD(argv[0]);    
    string deam5pfreqE  = getFullPath(cwdProg+"../deaminationProfile/none.prof");
    string deam3pfreqE  = getFullPath(cwdProg+"../deaminationProfile/none.prof");
    string illuminafreq = getFullPath(cwdProg+"../illuminaProf/null.prof");
    string dnafreqFile  = getFullPath(cwdProg+"../DNAprof/default_30_20_20_30.freq");

    // cout<<deam5pfreqE<<endl;
    // cout<<deam3pfreqE<<endl;

    // return 1;

    //no contaminant deamination for now
    // string deam5pfreqC = getCWD(argv[0])+"deaminationProfile/none.prof";
    // string deam3pfreqC = getCWD(argv[0])+"deaminationProfile/none.prof";

    // vector<substitutionRates>    deam5PsubE;
    // vector<substitutionRates>    deam3PsubE;
    // vector<substitutionRates>    deam5PsubC;
    // vector<substitutionRates>    deam3PsubC;


    int    numberOfThreads   = 1;
    string outFileSiteLL;
    bool   outFileSiteLLFlag=false;

    string sampleName        = "sample";
    bool   useVCFoutput      = false;

    string genoFileAsInput    ="";
    bool   genoFileAsInputFlag=false;

    
    const string usage=string("\nThis program co-estimates heterozygosity rates and runs of homozygosity\n\n\t"+
                              string(argv[0])+                        
                              " [options] [fasta file] [bam file]  "+"\n\n"+
			      "\twhere:\n"+
			      "\t\t[fasta file]\t\tThe fasta file used for alignement\n"
			      "\t\t[bam file]\t\tThe aligned and indexed BAM file\n"+
			      "\n\n"
			      
                              "\n\tI/O options:\n"+
			      "\t\t"+"-o"+"\t"+"--out"  + "\t\t"   +    "[outfile]" +"\t\t"+"Output per-site likelihoods in BGZIP (default: none)"+"\n"+
			      "\t\t"+""  +"\t"+"--name" + "\t\t"   +    "[name]"    +"\t\t\t"+"Sample name (default: "+sampleName+")"+"\n"+
			      "\t\t"+""  +""+"--vcf"    + "\t\t\t" +    ""          +"\t\t\t"+"Use VCF as output format (default: "+booleanAsString(useVCFoutput)+")"+"\n"+
			      //"\t\t"+""+"\t"+"--ingeno"  + "\t\t"   +    "[infile]" +"\t\t"+"Read likelihoods in BGZIP and start comp. from there (default: none)"+"\n"+
			      "\n\tComputation options:\n"+
                              "\t\t"+"-t"+"\t"+""       +"\t\t"    +    "[threads]" +"\t\t"+"Number of threads to use (default: "+stringify(numberOfThreads)+")"+"\n"+
                              "\t\t"+""  +""+"--phred64"+"\t\t\t"  +    ""          +"\t\t"+"Use PHRED 64 as the offset for QC scores (default : PHRED33)"+"\n"+
			      "\t\t"+""  +""+"--size"       +"\t\t\t"    + "[window size]" +"\t\t"+"Size of windows in bp  (default: "+stringify(sizeChunk)+")"+"\n"+	      



                              "\n\tSample options:\n"+
                              "\t\t"+""  +""+"--cont"  +"\t\t\t"    +  "[cont rate:0-1]" +"\t\t"+"Present-day human contamination rate (default: "+stringify(contrate)+")"+"\n"+
                              // "\t\t"+"--phred64" +"\t\t\t\t"+"Use PHRED 64 as the offset for QC scores (default : PHRED33)"+"\n"+
			      
			      
                              "\n\tDeamination and error options:\n"+                                   
                              "\t\t"+""  +""+"--deam5p\t\t"+"[.prof file]" +"\t\t"+"5p deamination frequency for the endogenous\n\t\t\t\t\t\t\t\t(default: "+deam5pfreqE+")"+"\n"+
                              "\t\t"+""  +""+"--deam3p\t\t"+"[.prof file]" +"\t\t"+"3p deamination frequency for the endogenous\n\t\t\t\t\t\t\t\t(default: "+deam3pfreqE+")"+"\n"+
                              // "\t\t"+"-deam5pc [.prof file]" +"\t\t"+"5p deamination frequency for the contaminant (default: "+deam5pfreqC+")"+"\n"+
                              // "\t\t"+"-deam3pc [.prof file]" +"\t\t"+"3p deamination frequency for the contaminant (default: "+deam3pfreqC+")"+"\n"+			      
			      "\t\t"+""  +""+"--err\t\t\t"    +"[.prof file]"+"\t\t"    +" Illumina error profile (default: "+illuminafreq+")"+"\n"+
			      "\t\t"+""  +""+"--base\t\t\t"   +"[.freq file]"+"\t\t"    +" Frequency of DNA bases in the genome (default: "+dnafreqFile+")"+"\n"+

                              "");


    if( (argc== 1) ||
        (argc== 2 && string(argv[1]) == "-h") ||
        (argc== 2 && string(argv[1]) == "-help") ||
        (argc== 2 && string(argv[1]) == "--help") ){
        cout<<usage<<endl;
        return 1;
    }

    int lastOpt=1;

    for(int i=1;i<(argc);i++){ 

	//cout<<i<<"\t"<<string(argv[i])<<endl;

        if(string(argv[i])[0] != '-' ){
            lastOpt=i;
            break;
        }


	if( string(argv[i]) == "--ingeno" ){
	    genoFileAsInput     = string(argv[i+1]);
	    genoFileAsInputFlag = true;
            i++;
            continue;
	}
	
        if( string(argv[i]) == "--size"  ){
	    sizeChunk=destringify<unsigned int>(argv[i+1]);
            i++;
            continue;
        }
	
        if( string(argv[i]) == "--vcf"  ){
	    useVCFoutput=true;
            continue;
        }

        if( string(argv[i]) == "--name"  ){
            sampleName=string(argv[i+1]);
            i++;
            continue;
        }


        if( string(argv[i]) == "-t"  ){
            numberOfThreads=destringify<int>(argv[i+1]);
            i++;
            continue;
        }

        if( string(argv[i]) == "--cont"  ){
            contrate=destringify<long double>(argv[i+1]);
            i++;
            continue;
        }


        if( string(argv[i]) == "-o"    ||
	    string(argv[i]) == "--out" ){
            outFileSiteLL=string(argv[i+1]);
	    outFileSiteLLFlag=true;
            i++;
            continue;
        }


        if(string(argv[i]) == "--phred64"  ){
            offsetQual=64;
            continue;
        }

        if(string(argv[i]) == "--deam5p"  ){
            deam5pfreqE=string(argv[i+1]);
            i++;
            continue;
        }

        if(string(argv[i]) == "--deam3p"  ){
            deam3pfreqE=string(argv[i+1]);
            i++;
            continue;
        }

	if(string(argv[i]) == "--err"  ){
	    illuminafreq=string(argv[i+1]);
	    i++;
	    continue;
	}
	
	if(string(argv[i]) == "--base"  ){
	    dnafreqFile=string(argv[i+1]);
	    i++;
	    continue;
	}


	cerr<<"Error: unknown option "<<string(argv[i])<<endl;
	return 1;
    }


    string fastaFile         = string(argv[lastOpt]);
    bamFileToOpen            = string(argv[lastOpt+1]);
    string fastaIndex        = fastaFile+".fai";

    if( !isFile(fastaFile) ){
	cerr<<"The fasta file "<<fastaFile<<" does not exists"<<endl;
	return 1;	
    }

    if( !isFile(fastaIndex) ){
	cerr<<"The fasta file "<<fastaFile<<"  does not have an index: "<<fastaIndex<<endl;
	return 1;	
    }


    if(outFileSiteLLFlag)
	if( !strEndsWith(outFileSiteLL,".gz")){
	    cerr<<"The output file "<<outFileSiteLL<<" must end with .gz"<<endl;
	    return 1;	
	}




    //Testing BAM file
    BamReader reader;
    if ( !reader.Open(bamFileToOpen) ) {
	cerr << "Could not open input BAM file:" << bamFileToOpen <<endl;
    	exit(1);
    }

    reader.LocateIndex();

    if(!reader.HasIndex()){
    	cerr << "The BAM file: " << bamFileToOpen <<" does not have an index"<<endl;
    	exit(1);
    }

    // retrieve reference data
    const RefVector  references = reader.GetReferenceData();


    reader.Close();

    ////////////////////////////////////
    //   END Parsing arguments        //
    ////////////////////////////////////

















    ////////////////////////////////////
    // BEGIN Initializing scores      //
    ////////////////////////////////////
    initScores();
    ////////////////////////////////////
    //    END Initializing scores     //
    ////////////////////////////////////




       


    ///////////////////////////////////////
    // BEGIN read Illumina error profile //
    ///////////////////////////////////////    


    readIlluminaError(illuminafreq,illuminaErrorsProb);

#ifdef DEBUGILLUMINAFREQ
    for(unsigned int i=0;i<16;i++){
    	cout<<i<<"\t"<<illuminaErrorsProb.s[i]<<endl;	
    }
    //return 1;
#endif
    ///////////////////////////////////////
    // END read Illumina error profile //
    ///////////////////////////////////////    



    ////////////////////////////////////////////////////////////////////////
    //
    // BEGIN DEAMINATION PROFILE
    //
    ////////////////////////////////////////////////////////////////////////
    

    initDeamProbabilities(deam5pfreqE,deam3pfreqE);
    
    // return 1;
    ////////////////////////////////////////////////////////////////////////
    //
    // END  DEAMINATION PROFILE
    //
    ////////////////////////////////////////////////////////////////////////



    // for(int q=0;q<=63;q++){
    // 	long double mmp=0.000000000001;



    // 	int al1Current=1;
    // 	long double sum=0;
    // 	cout<<"q="<<q<<"\t";
    // 	for(int obs=0;obs<4;obs++){

    // 	    long double t=computeBaseAl2Obs(al1Current  ,
    // 					    obs         ,
    // 					    q           ,
    // 					    //&defaultSubMatch,
    // 					    &(sub3p[0]) ,
    // 					    false       ,
    // 					    mmp         );
    // 	    sum+=t;
    // 	    cout<<obs<<"\t"<<t<<"\t";
		

    // 	    // int dinucIndex1toObs =     al1Current               *4+              obs;      
    // 	    // //long double probSubDeam1toObs              = defaultSubMatch.s[dinucIndex1toObs];
    // 	    // long double probSubDeam1toObs              = sub3p[0].s[dinucIndex1toObs];
	    


    // 	}
    // 	cout<<"\t"<<sum<<endl;

    // }




    ////////////////////////////////////////////////////////////////////////
    //
    // BEGIN DNA BASE FREQUENCY
    //
    ////////////////////////////////////////////////////////////////////////

    initDefaultBaseFreq(dnafreqFile);

    
    //return 1;
    ////////////////////////////////////////////////////////////////////////
    //
    // END DNA BASE FREQUENCY
    //
    ////////////////////////////////////////////////////////////////////////




    // alleleFrequency dnaDefaultBases;
    // vector<alleleFrequency> defaultDNA5p;
    // vector<alleleFrequency> defaultDNA3p;

    // for(unsigned int i=0;i<sub5p.size();i++){
    // 	//defaultDNA5p
    // 	alleleFrequency toadd;
    // 	for(int b=0;b<4;b++){
	    
    // 	}
    // 	for(unsigned int i=0;i<sub5p.size();i++){
    // 	}
    // }
    // vector<alleleFrequency> defaultDNA5p;
    // vector<alleleFrequency> defaultDNA3p;



    ////////////////////////////////////////////////////////////////////////
    //
    // BEGIN COMPUTE P[OBSERVED BASE|THEORITICAL BASE]
    //
    ////////////////////////////////////////////////////////////////////////

    initLikelihoodScores();

    ////////////////////////////////////////////////////////////////////////
    //
    // END COMPUTE P[OBSERVED BASE|THEORITICAL BASE]
    //
    ////////////////////////////////////////////////////////////////////////





    int    bpToExtract       = sizeChunk;
    
    pthread_t             thread[numberOfThreads];
    int                   rc=0;



    GenomicWindows     rw  (fastaIndex,false);
    //TODO add genomic ranges in queue
    vector<GenomicRange> v = rw.getGenomicWindows(bpToExtract,0);
    if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
    
    unsigned int      rank=0;
    int          lastRank=-1;
    unsigned int sizeGenome=0;

    for(unsigned int i=0;i<v.size();i++){
	//cout<<v[i]<<endl;
	DataChunk * currentChunk = new DataChunk();
	
	currentChunk->rangeGen  = v[i];
	currentChunk->rank      = rank;
	lastRank                = rank;
	//sizeGenome             += v[i].getLength();
 
	queueDataToprocess.push(currentChunk);
	rank++;
    }
    readDataDone=true;
    //    return 1;

    vector<GenoResults *> vectorGenoResults;

    if(!genoFileAsInputFlag){

	/////////////////////////////
	// BEGIN  Compute coverage //
	/////////////////////////////
	int bpToComputeCoverage = 1000000;
	int genomicRegionsToUse = bpToComputeCoverage/bpToExtract;
	if( genomicRegionsToUse > int(queueDataToprocess.size())){
	    genomicRegionsToUse = int(queueDataToprocess.size());
	}



	queueDataForCoverage = randomSubQueue( queueDataToprocess,genomicRegionsToUse);


	pthread_mutex_init(&mutexQueue,   NULL);
	pthread_mutex_init(&mutexCounter, NULL);
	pthread_mutex_init(&mutexRank ,   NULL);

	for(int i=0;i<numberOfThreads;i++){
	    rc = pthread_create(&thread[i], NULL, mainCoverageComputationThread, NULL);
	    checkResults("pthread_create()\n", rc);
	}

	cerr<<"Creating threads for coverage calculation"<<endl;


	while(queueDataForCoverage.size()!=0){
	    cerr<<getDateString()<<" "<<getTimeString()<<" # of slices left to process: "<<queueDataForCoverage.size()<<"/"<<queueDataToprocess.size()<<endl;
	    sleep(timeThreadSleep);
	}
    
	//waiting for threads to finish
	for (int i=0; i <numberOfThreads; ++i) {	
	    rc = pthread_join(thread[i], NULL);
	    checkResults("pthread_join()\n", rc);
	}
	cerr<<"Coverage computations are done"<<endl;
	pthread_mutex_destroy(&mutexRank);
	pthread_mutex_destroy(&mutexQueue);
	pthread_mutex_destroy(&mutexCounter);

	//    cout<<"Final" <<" "<<totalBasesSum<<"\t"<<totalSitesSum<<"\t"<<double(totalBasesSum)/double(totalSitesSum)<<endl;
	//    pthread_exit(NULL);
    
	rateForPoissonCov    = ((long double)totalBasesSum)/((long double)totalSitesSum);
	long double rateForPoissonCovFloor = floorl(rateForPoissonCov);
	long double rateForPoissonCovCeil  =  ceill(rateForPoissonCov);

#ifdef DEBUGCOV
	cout<<"rateForPoissonCovFloor "<<rateForPoissonCovFloor<<endl;
	cout<<"rateForPoissonCovCeil "<<rateForPoissonCovCeil<<endl;
#endif

	cov2probPoisson = vector<long double> (MAXCOV,0);

	for(int i=0;i<=MAXCOV;i++){
	    long double florPPMF=poisson_pmfl( (long double)(i), rateForPoissonCovFloor);
	    long double ceilPPMF=poisson_pmfl( (long double)(i), rateForPoissonCovCeil);
	
	
	    cov2probPoisson[i] = (1.0-(rateForPoissonCov-rateForPoissonCovFloor))*florPPMF  + (1.0-(rateForPoissonCovCeil-rateForPoissonCov))*ceilPPMF;

#ifdef DEBUGCOV
	    cerr<<"florPPMF "<<i<<" = "<<florPPMF<<" w= "<<(1.0-(rateForPoissonCov-rateForPoissonCovFloor))<<endl;
	    cerr<<"ceilPPMF "<<i<<" = "<<ceilPPMF<<" w= "<<(1.0-(rateForPoissonCovCeil-rateForPoissonCov))<<endl;
	    cerr<<"cov2probPoisson["<<i<<"] = "<<cov2probPoisson[i]<<endl;
#endif

	}



	cerr<<"Results\tbp="<<totalBasesSum<<"\tsites="<<totalSitesSum<<"\tlambda="<<double(totalBasesSum)/double(totalSitesSum)<<endl;
	// for(int i=0;i<20;i++){
	// 	cout<<i<<"\t"<<pdfPoisson( (long double)i, rateForPoissonCov)/pdfPoisson( rateForPoissonCov, rateForPoissonCov)<<endl;
    // }


    // for(int i=0;i<100;i++){
    // 	cout<<i<<"\t"<<pdfPoisson( (long double)i, 20)/pdfPoisson( 20, 20)<<endl;
    // }

    return 1;


    ////////////////////////////
    // END   Compute coverage //
    ////////////////////////////



















    // doneReading=true;    

    ///////////////////////
    //  Compute hetero   //
    ///////////////////////
    cerr<<"Creating threads for heterozygosity calculation"<<endl;
    pthread_mutex_init(&mutexQueue,   NULL);
    pthread_mutex_init(&mutexCounter, NULL);
    pthread_mutex_init(&mutexRank ,   NULL);

    for(int i=0;i<numberOfThreads;i++){
	rc = pthread_create(&thread[i], NULL, mainHeteroComputationThread, NULL);
	checkResults("pthread_create()\n", rc);
    }
    // 	//threads are running here

    // unsigned int originalSize = queueDataToprocess.size();
    // while(queueDataToprocess.size()!=0){
    // 	cout<<"# of slices left to process: "<<queueDataToprocess.size()<<"/"<<originalSize<<endl;
    // 	sleep(timeThreadSleep);
    // }


    ///////////////////
    //Writing data out/
    ///////////////////


    Internal::BgzfStream  bgzipWriter;

    if(outFileSiteLLFlag){
	bgzipWriter.Open(outFileSiteLL, IBamIODevice::WriteOnly);
	if(!bgzipWriter.IsOpen()){
	    cerr<<"Cannot open file "<<outFileSiteLL<<" in bgzip writer"<<endl;
	    return 1;
	}
    

	string headerOutFile;
	if(useVCFoutput){
	    headerOutFile="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+sampleName+"\n";	
	}else{
	    headerOutFile="#CHROM\tPOS\tA\tC\tG\tT\tGENO\tGENOS\tQualL\tCovL\tAA\tAC\tAG\tAT\tCC\tCG\tCT\tGG\tGT\tTT\n";	
	}

	bgzipWriter.Write(headerOutFile.c_str(),headerOutFile.size());
    }

    bool wroteEverything=false;
    int lastWrittenChunk=-1;   
           
    while(!wroteEverything){

	//threads are running here
	rc = pthread_mutex_lock(&mutexCounter);
	checkResults("pthread_mutex_lock()\n", rc);
	
	bool wroteData=false;
	if(!queueDataTowrite.empty()){
	
	    DataToWrite *  dataToWrite= queueDataTowrite.top();

	    if( lastWrittenChunk == (dataToWrite->rank-1) ){ 	    //correct order
		queueDataTowrite.pop();
		rc = pthread_mutex_unlock(&mutexCounter);
		checkResults("pthread_mutex_unlock()\n", rc);

		cerr<<getDateString()<<" "<<getTimeString()<<" writing chunk#"<<dataToWrite->rank<<" with "<<dataToWrite->vecPositionResults->size()<<" records"<<endl;
		sizeGenome+=dataToWrite->vecPositionResults->size();
		    
		string strToWrite="";
		for(unsigned int i=0;i<dataToWrite->vecPositionResults->size();i++){
		    //cout<<i<<"\t"<<dataToWrite->vecPositionResults->at(i)->toString(references)<<endl;
		    strToWrite += dataToWrite->vecPositionResults->at(i)->toString(references);
		    //cout<<i<<"\t"<<dataToWrite->vecPositionResults->at(i)->toString(references);
		    if( (i%500) == 499){
			if(outFileSiteLLFlag){
			    bgzipWriter.Write(strToWrite.c_str(), strToWrite.size());
			}
			strToWrite="";
		    }
		    
		    GenoResults * toadd =  new GenoResults( dataToWrite->vecPositionResults->at(i) );
		    vectorGenoResults.push_back(toadd);
		}

		if(!strToWrite.empty()){
		    if(outFileSiteLLFlag){
			if(outFileSiteLLFlag){ bgzipWriter.Write(strToWrite.c_str(), strToWrite.size()); }
		    }
		}
		    
		
		wroteData=true;		
		lastWrittenChunk=dataToWrite->rank;
		
		if(dataToWrite->rank == lastRank)
		    wroteEverything=true;	
		delete dataToWrite;
	    }else{
		//do nothing, we have to wait for the chunk with the right rank
		rc = pthread_mutex_unlock(&mutexCounter);
		checkResults("pthread_mutex_unlock()\n", rc);
	
	    }

	}else{//end if queue not empty
	    rc = pthread_mutex_unlock(&mutexCounter);
	    checkResults("pthread_mutex_unlock()\n", rc);
	}

	if(!wroteData)
	    sleep(timeSleepWrite);
    }

    ///////////////////////
    //end Writing data out/
    ///////////////////////

    //return 1;

    // 	rc = pthread_mutex_lock(&mutexCounter);
    // 	checkResults("pthread_mutex_lock()\n", rc);
    // 	bool wroteData=false;
    // 	if(!queueDataTowrite.empty()){
	
    // 	    DataChunk *  dataToWrite= queueDataTowrite.top();

    // 	    if( lastWrittenChunk == (dataToWrite->rank-1) ){
    // 		queueDataTowrite.pop();
    // 		cout<<"writing "<<dataToWrite->rank<<endl;
    // 		//writing dataToWrite
    // 		for(unsigned int i=0;i<dataToWrite->dataToProcess.size();i++){
    // 		    outfile<< dataToWrite->dataToProcess[i].ids 
    // 			   << endl  
    // 			   << dataToWrite->dataToProcess[i].seqs << endl
    // 			   << "+" <<endl 
    // 			   << dataToWrite->dataToProcess[i].qual << endl;
    // 		}
		
    // 		wroteData=true;		
    // 		lastWrittenChunk=dataToWrite->rank;
		
    // 		if(dataToWrite->rank == lastRank)
    // 		    wroteEverything=true;	
    // 		delete dataToWrite;
    // 	    }else{
		
    // 	    }

    // 	}

    // 	rc = pthread_mutex_unlock(&mutexCounter);
    // 	checkResults("pthread_mutex_unlock()\n", rc);

    // 	if(!wroteData)
    // 	    sleep(timeThreadSleep);
    // }

    //waiting for threads to finish    
    for (int i=0; i <numberOfThreads; ++i) {
	rc = pthread_join(thread[i], NULL);
	checkResults("pthread_join()\n", rc);
    }

    cerr<<"Heterozygosity computations are done"<<endl;
    pthread_mutex_destroy(&mutexRank);
    pthread_mutex_destroy(&mutexQueue);
    pthread_mutex_destroy(&mutexCounter);


    if(outFileSiteLLFlag){
	bgzipWriter.Close();
    }

    }else{//if we suply the geno file as input

	string    lineGENOL;
	igzstream myFileGENOL;
	cerr<<"Reading genotypes: "<<genoFileAsInput<<" ...";
	myFileGENOL.open(genoFileAsInput.c_str(), ios::in);

	if (myFileGENOL.good()){
	    getline (myFileGENOL,lineGENOL);//header
	    while ( getline (myFileGENOL,lineGENOL)){
		//cout<<lineGENOL<<endl;
		GenoResults * toadd =  new GenoResults( lineGENOL );
		vectorGenoResults.push_back(toadd);
		sizeGenome++;
	    }
	    myFileGENOL.close();
	}else{
	    cerr << "Error: unable to open file with genotype likelihoods: "<<genoFileAsInput<<endl;
	    return 1;
	}
	cerr<<"..done"<<endl;

    }

    //THORFINN: END OF GENOTYPE LIKELIHOOD COMPUTATIONS


    //    return 1;
    //////////////////////////////////
    //                              //
    // BEGIN COMPUTE HETERO RATE    //
    //                              //
    //////////////////////////////////

    long double randomLog = -1.0*log(1000);
    vector<GenoResults *> vectorGenoResultsT;

    for(unsigned int i=0;i<vectorGenoResults.size();i++){
	// if(i == 796868){
	//cout<<*(vectorGenoResults[i])<<endl;
	// }
	long double sumProbHomoz=0;
	long double sumProbHeter=0;
	long double mostLikeHomozy =  -1.0*numeric_limits<long double>::infinity();
	long double mostLikeHetero =  -1.0*numeric_limits<long double>::infinity();



	for(int g=0;g<4;g++){
	    sumProbHomoz = oplusInitnatl( sumProbHomoz , vectorGenoResults[i]->ll[ genoPriority[g] ] );//+logl( ( 1/( (long double)4))) );

	    if( mostLikeHomozy < vectorGenoResults[i]->ll[ genoPriority[g] ]){
		mostLikeHomozy = vectorGenoResults[i]->ll[ genoPriority[g] ];
	    }

	}
	sumProbHomoz = sumProbHomoz + logl( ((long double)1)/((long double)4) );

	for(int g=4;g<10;g++){
	    sumProbHeter = oplusInitnatl( sumProbHeter , vectorGenoResults[i]->ll[ genoPriority[g] ] );//+logl( ( 1/( (long double)6))) );

	    if( mostLikeHetero  < vectorGenoResults[i]->ll[ genoPriority[g] ]){
		mostLikeHetero  = vectorGenoResults[i]->ll[ genoPriority[g] ];
	    }

	}
	sumProbHeter = sumProbHeter + logl( ((long double)1)/((long double)6) );
	//  = oplusnatl( vectorGenoResults[i]->rrll , vectorGenoResults[i]->aall );
	// long double sumProbHeter = oplusnatl( vectorGenoResults[i]->rall , vectorGenoResults[i]->rall );//twice
	
	long double minSumProb   = MIN( sumProbHomoz, sumProbHeter);

	long double sumProbAll   = oplusnatl( sumProbHomoz , sumProbHeter );
	

	long double confidence = minSumProb-sumProbAll;

	//confidence = MIN(0,confidence+randomLog);


	//confidence=1-expl(confidence);
	long double covCorrect = (logl(2)) * vectorGenoResults[i]->cov; //multiply by 2 = log(2)
	long double sumProbHoHe = oplusnatl( sumProbHomoz,sumProbHeter);
	long double lqual=-1.0*numeric_limits<long double>::infinity();

	if(mostLikeHomozy>mostLikeHetero){//if is homozygous, check if scaling the het changes
	    // long double mostLikeHeteroC         = MIN( mostLikeHetero+covCorrect , mostLikeHomozy);//takes care of coverage =1,2	   
	    lqual                               = (mostLikeHetero-mostLikeHomozy)/1.0;
	}else{//most likely is hetero
	    
	    //lqual                               = (mostLikeHomozy-mostLikeHetero)/2.0;
	    lqual                               = (mostLikeHomozy-mostLikeHetero)/1.0;
	    //lqual                               = (mostLikeHomozy-mostLikeHetero)/1.8;
	}
	vectorGenoResults[i]->expectedH     = expl(  sumProbHeter - sumProbHoHe );	

	//lqual = lqual*expl(vectorGenoResults[i]->llCov);
	

	//TODO check if covCorrect mitigates low cov problem
	
	
	//vectorGenoResults[i]->probAccurate =  (  (1.0-expl(vectorGenoResults[i]->lqual) ) ) * ( expl(vectorGenoResults[i]->llCov)  );	
	//vectorGenoResults[i]->probAccurate =  1.0-expl( (  lqual ) * ( expl(vectorGenoResults[i]->llCov)  ) );
	vectorGenoResults[i]->probAccurate =  1.0-expl(  lqual ) ;
	
	//cout<<i<<"\t"<<mostLikeHomozy<<"\t"<<mostLikeHetero<<"\t"<<sumProbHomoz<<"\t"<<sumProbHeter<<"\t"<<sumProbHoHe<<"\t"<<vectorGenoResults[i]->cov<<"\t"<<covCorrect<<"\th=\t"<<vectorGenoResults[i]->expectedH<<"\t"<<lqual<<"\t"<<vectorGenoResults[i]->lqual<<"\t"<<expl(vectorGenoResults[i]->llCov)<<"\t"<<vectorGenoResults[i]->probAccurate<<endl;

	//cout<<vectorGenoResults[i]->lqual <<"\t"<< randomLog<<endl;

	//	if( lqual < randomLog){
	vectorGenoResultsT.push_back( vectorGenoResults[i] );
	// } else {
	//     delete vectorGenoResults[i];
	// }
	//scale for coverage 2xhomo
	//cout<<"passed"<<endl;
	
	
	//cout<<fixed<<i<<"\th="<<vectorGenoResults[i]->expectedH<<"\t"<<"\t"<<sumProbHomoz<<"\t"<<sumProbAll<<"\thom\t"<<sumProbHomoz<<"\tHet\t"<<sumProbHeter<<endl;
	// if(i == 796868){
	//cout<<fixed<<i<<"\th=\t"<<vectorGenoResults[i]->expectedH<<"\tp[q]=\t"<<(1.0-expl(vectorGenoResults[i]->lqual ) ) <<"\tcov=\t"<<expl(vectorGenoResults[i]->llCov)<<"\tconf\t"<<confidence<<"\tP[acc]\t"<<vectorGenoResults[i]->probAccurate<<"\t"<<"\t"<<sumProbAll<<"\thom="<<sumProbHomoz<<"\tHet="<<sumProbHeter<<"\t"<<randomLog<<endl;
	//     return 1;
	// }

    }

    //return 1;

    vectorGenoResults = vectorGenoResultsT;
    // for(unsigned int i=0;i<vectorGenoResults.size();i++){

    // 	// if(i == 796868){
    // 	//cout<<*(vectorGenoResults[i])<<endl;
    // 	// }
    // 	long double sumProbHomoz=0;
    // 	long double sumProbHeter=0;
    // 	long double mostLikeHomozy =  -1.0*numeric_limits<long double>::infinity();
    // 	long double mostLikeHetero =  -1.0*numeric_limits<long double>::infinity();



    // 	for(int g=0;g<4;g++){
    // 	    sumProbHomoz = oplusInitnatl( sumProbHomoz , vectorGenoResults[i]->ll[ genoPriority[g] ]+logl( ( 1/( (long double)4))) );

    // 	    if( mostLikeHomozy < vectorGenoResults[i]->ll[ genoPriority[g] ]){
    // 		mostLikeHomozy = vectorGenoResults[i]->ll[ genoPriority[g] ];
    // 	    }

    // 	}


    // 	for(int g=4;g<10;g++){
    // 	    sumProbHeter = oplusInitnatl( sumProbHeter , vectorGenoResults[i]->ll[ genoPriority[g] ]+logl( ( 1/( (long double)6))) );

    // 	    if( mostLikeHetero  < vectorGenoResults[i]->ll[ genoPriority[g] ]){
    // 		mostLikeHetero  = vectorGenoResults[i]->ll[ genoPriority[g] ];
    // 	    }

    // 	}
    // 	//  = oplusnatl( vectorGenoResults[i]->rrll , vectorGenoResults[i]->aall );
    // 	// long double sumProbHeter = oplusnatl( vectorGenoResults[i]->rall , vectorGenoResults[i]->rall );//twice
	
    // 	long double minSumProb   = MIN( sumProbHomoz, sumProbHeter);

    // 	long double sumProbAll   = oplusnatl( sumProbHomoz , sumProbHeter );
	

    // 	long double confidence = minSumProb-sumProbAll;

    // 	//confidence = MIN(0,confidence+randomLog);


    // 	//confidence=1-expl(confidence);
    // 	long double covCorrect = (logl(2)) * vectorGenoResults[i]->cov; //multiply by 2 = log(2)
    // 	long double sumProbHoHe;
    // 	long double lqual=-1.0*numeric_limits<long double>::infinity();
    // 	if(mostLikeHomozy>mostLikeHetero){//if is homozygous, check if scaling the het changes
    // 	    long double mostLikeHeteroC         = MIN( mostLikeHetero+covCorrect , mostLikeHomozy);//takes care of coverage =1,2
    // 	    sumProbHoHe                         = oplusnatl( mostLikeHomozy,mostLikeHeteroC);
    // 	    vectorGenoResults[i]->expectedH     = expl(  mostLikeHeteroC - sumProbHoHe );
    // 	    lqual                               = (mostLikeHetero-mostLikeHomozy)/1.0;
    // 	}else{//most likely is hetero
    // 	    sumProbHoHe                         = oplusnatl( mostLikeHomozy,mostLikeHetero);
    // 	    vectorGenoResults[i]->expectedH     = expl(  mostLikeHetero - sumProbHoHe );
    // 	    lqual                               = (mostLikeHomozy-mostLikeHetero)/2.0;
    // 	}

	
	
    // 	cout<<mostLikeHomozy<<"\t"<<mostLikeHetero<<"\t"<<sumProbHoHe<<"\t"<<vectorGenoResults[i]->cov<<"\t"<<covCorrect<<"\t"<<vectorGenoResults[i]->expectedH<<"\t"<<lqual<<"\t"<<vectorGenoResults[i]->lqual<<"\t"<<expl(vectorGenoResults[i]->llCov)<<"\t"<<vectorGenoResults[i]->probAccurate<<endl;
    // }


#ifdef DEBUGHCOMPUTE
    long double h         = 1/ (long double)(randomInt(1000,20000));
    long double he        = 1/ (long double)(randomInt(1000,20000));

    long double hW        = 1/ (long double)(randomInt(1000,20000));
    long double hWe       = 1/ (long double)(randomInt(1000,20000));

    long double h_old     = 10000;
    long double hAccu     = 0.00000005;

    long double lambdaH   = 0.0000000001;
    long double lambdaHW  = 0.0000000001;

    int iterationsMax     = 10000;
    int iterationsGrad    = 1;

    while( iterationsGrad<iterationsMax){

	// if(h>=1){
	//     h=1-espilon;
	// }
	long double probNull = 0.000001;
	long double ll    = 0.0;
	long double llP   = 0.0;
	long double llPP  = 0.0;

	long double llW   = 0.0;
	long double llPW  = 0.0;
	long double llPPW = 0.0;


	for(unsigned int i=0;i<vectorGenoResults.size();i++){


	    long double llT;
	    long double llTP;
	    long double llTPP;

	    long double llTW;
	    long double llTPW;
	    long double llTPPW;


	    
	    long double pcorrect=MAX((vectorGenoResults[i]->probAccurate),probNull);
	    
	    llT  = expl(vectorGenoResults[i]->llCov)*logl( pcorrect
							   *
							   ( (1-h )*(1-vectorGenoResults[i]->expectedH) + h*vectorGenoResults[i]->expectedH )
							   +
							   (1-pcorrect)*
							   ( (1-hW)*(1-vectorGenoResults[i]->expectedH) + hW*vectorGenoResults[i]->expectedH )
							   );

	    
    	    llTP = 
	    	expl(vectorGenoResults[i]->llCov)*(  pcorrect*(2.0*vectorGenoResults[i]->expectedH-1) )
	    	/
	    	( pcorrect
	    	  *
	    	  ( (1-h)*(1-vectorGenoResults[i]->expectedH) + h*vectorGenoResults[i]->expectedH )
		  +
		  (1-pcorrect)*
		  ( (1-hW)*(1-vectorGenoResults[i]->expectedH) + hW*vectorGenoResults[i]->expectedH )
	    	  );


    	    llTPP = -1.0*
	    	expl(vectorGenoResults[i]->llCov)*(  powl( pcorrect,2.0) * powl((2.0*vectorGenoResults[i]->expectedH-1),2.0) )
	    	/
	    	powl( 
		     pcorrect
		     *
		     ( (1-h)*(1-vectorGenoResults[i]->expectedH) + h*vectorGenoResults[i]->expectedH )
		     +
		     (1-pcorrect)*
		     ( (1-hW)*(1-vectorGenoResults[i]->expectedH) + hW*vectorGenoResults[i]->expectedH )
		      ,2.0);



	    llTW  = llT;
	    
    	    llTPW = 
	    	expl(vectorGenoResults[i]->llCov)*(  (1-pcorrect)*(2.0*vectorGenoResults[i]->expectedH-1) )
	    	/
	    	( pcorrect
	    	  *
	    	  ( (1-h)*(1-vectorGenoResults[i]->expectedH) + h*vectorGenoResults[i]->expectedH )
		  +
		  (1-pcorrect)*
		  ( (1-hW)*(1-vectorGenoResults[i]->expectedH) + hW*vectorGenoResults[i]->expectedH )
	    	  );


    	    llTPPW = -1.0*
	    	expl(vectorGenoResults[i]->llCov)*(  powl( 1-pcorrect,2.0) * powl( (2.0*vectorGenoResults[i]->expectedH-1),2.0 ) )
	    	/
	    	powl( 
		     pcorrect
		     *
		     ( (1-h)*(1-vectorGenoResults[i]->expectedH) + h*vectorGenoResults[i]->expectedH )
		     +
		     (1-pcorrect)*
		     ( (1-hW)*(1-vectorGenoResults[i]->expectedH) + hW*vectorGenoResults[i]->expectedH )
		      ,2.0);





	    // cout<<i<<"\tp[acc]="<<vectorGenoResults[i]->probAccurate<<"\t"<<(1-vectorGenoResults[i]->probAccurate)<<"\t"<<MAX((1-vectorGenoResults[i]->probAccurate),probNull)<<"\th= "<<h<<"\t"<<vectorGenoResults[i]->expectedH<<"\tllT= "<<llT<<"\t"<<llTP<<"\t"<<llTPP<<"\tll= "<<ll<<"\t"<<llP<<"\t"<<llPP<<"\tllW= "<<llW<<"\t"<<llPW<<"\t"<<llPPW<<endl;

	    
	    ll   += llT;
	    llP  += llTP;
	    llPP += llTPP;


	    llW   += llTW;
	    llPW  += llTPW;
	    llPPW += llTPPW;

	}//for each item vector
	
	he  = 1.96/sqrtl( -1.0*llPP  );
	hWe = 1.96/sqrtl( -1.0*llPPW );

	cout<<fixed<<
	    "h= "<<h<<"\t"  <<(1-h)<<"\t"<<ll <<"\t"<<llP <<"\t"<<llPP <<"\t"<<he <<"\t"<<(h-he)  <<"\t"<<(h+he)<<"\t"<<
	    "hw= "<<hW<<"\t"<<(1-hW)<<"\t"<<llW<<"\t"<<llPW<<"\t"<<llPPW<<"\t"<<hWe<<"\t"<<(hW-hWe)<<"\t"<<(hW+hWe)<<
	    endl;
	//cout<<h<<"\t"<<(1-h)<<"\t"<<ll<<"\t"<<llP<<"\t"<<llPP<<endl;	


	long double hnew  = h        + lambdaH*llP;
	
	while( (hnew>=1) || (hnew<=0) ){
	    //cout<<setprecision(14)<<"lambdaH1\t"<<lambdaH<<"\t"<<h<<"\t"<<hnew<<"\t"<<llP<<endl;
	    lambdaH       = lambdaH/2.0;
	    hnew          = h        + lambdaH*llP;
	    //cout<<"lambdaH2\t"<<lambdaH<<"\t"<<h<<"\t"<<hnew<<"\t"<<llP<<endl;
	}
	
	long double hwNew = hW       + lambdaHW*llPW;
	while( (hwNew>=1) || (hwNew<=0) ){
	    //cout<<"lambdaHW1\t"<<lambdaHW<<"\t"<<hW<<"\t"<<hwNew<<"\t"<<llPW<<endl;
	    lambdaHW      = lambdaHW/2.0;
	    hwNew         = hW       + lambdaHW*llPW;
	    //cout<<"lambdaHW2\t"<<lambdaHW<<"\t"<<hW<<"\t"<<hwNew<<"\t"<<llPW<<endl;
	}
	
        h_old    = h;
	h        = hnew;
	hW       = hwNew;

	long double hDiff = fabsl(h_old - h);
        if(hDiff<hAccu){
            break;
        }
	

	iterationsGrad++;
    }//end loop of iterations


    long double hmin = (h-he);
    long double hmax = (h+he);
    if(hmin<0){
	hmax=hmax-hmin;
	hmin=0;
	//h=he/2;
    }


    // if(hmin<0){
    // 	hmin=0;
    // 	hmax=he;
    // 	h=he/2;
    // }

    cout<<"hetrate\t"<<h<<"\t"<<(hmin)<<"\t"<<(hmax)<<endl;
#endif    
    ////////////////////////////////
    //                            //
    //  end COMPUTE HETERO RATE   //
    //                            //
    ////////////////////////////////
        

    


    //////////////////////////////////
    //                              //
    // BEGIN HMM                    //
    //                              //
    //////////////////////////////////


    //////////////////////////////////
    //                              //
    //  END HMM                     //
    //                              //
    //////////////////////////////////



    for(unsigned int i=0;i<vectorGenoResults.size();i++){
	delete( vectorGenoResults[i] );
    }

    pthread_exit(NULL);

    
    return 0;
}

