#include <iostream>
#include <fstream>
#include <queue>
#include <algorithm>   
#include <cfloat>   
//#include <random>

//TODO
// code h estimate
//   precompute gl works
//   test at different h prior 0.01	-1.65138e+07
//                           0.000824	-1.65085e+07
//                           1e-08	-1.65141e+07


// HMM
// add mappability track?



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

#define PRECOMPUTELOG
//#define DEBUGCOV
//#define DEBUGILLUMINAFREQ
//#define DEBUGINITSCORES

// #define DEBUGINITLIKELIHOODSCORES
//#define DEBUGINITLIKELIHOODSCORES2

//#define DEBUGDEFAULTFREQ//print the default base frequency 
//#define DEBUGDEAM //to print deamination scores
//#define DEBUGHCOMPUTE
//#define DEBUGCOMPUTELLGENO
//#define DEBUGCOMPUTELL
//#define DEBUGCOMPUTELLEACHBASE



#define HETVERBOSE
//#define COVERAGETVERBOSE
#define DUMPTRIALLELIC //hack to remove tri-allelic, we need to account for them

#define MINLENGTHFRAGMENT     35      // mininam length for fragment
#define MAXLENGTHFRAGMENT     250     // maximal length for fragment


//#define MAXMAPPINGQUAL        257     // maximal mapping quality, should be sufficient as mapping qualities are encoded using 8 bits
#define MAXMAPPINGQUAL          37     // maximal mapping quality, should be sufficient as mapping qualities are encoded using 8 bits
#define MAXBASEQUAL             64      // maximal base quality score, greater qual score do not make a lot of sense

//#define MAXMAPPINGBASEQUAL    64      // maximal base quality score, should be sufficient as mapping qualities are encoded using 8 bits

#define MAXCOV             50     // maximal coverage

char offsetQual=33;

long double likeMatch           [MAXBASEQUAL+1];
long double likeMismatch        [MAXBASEQUAL+1];

long double likeMatchProb       [MAXBASEQUAL+1];
long double likeMismatchProb    [MAXBASEQUAL+1];

// long double likeMatchMap        [MAXMAPPINGQUAL+1];
// long double likeMismatchMap     [MAXMAPPINGQUAL+1];

long double likeMatchProbMap    [MAXMAPPINGQUAL+1];
long double likeMismatchProbMap [MAXMAPPINGQUAL+1];

long double TStoTVratio;

vector< vector<long double> > binomVec (MAXCOV+1,vector<long double>(MAXCOV+1,0)) ;
unsigned int totalBasesSum;
unsigned int totalSitesSum;

//Substitution rates due to deamination
vector<probSubstition> sub5p;
vector<probSubstition> sub3p;

vector<diNucleotideProb> sub5pDiNuc;
vector<diNucleotideProb> sub3pDiNuc;

probSubstition   defaultSubMatch;
diNucleotideProb defaultSubMatchMatrix;

probSubstition   illuminaErrorsProb;
diNucleotideProb illuminaErrorsProbMatrix;

//default basepairs frequencies
alleleFrequency         dnaDefaultBases;
vector<alleleFrequency> defaultDNA5p;
vector<alleleFrequency> defaultDNA3p;


// 1D: mapping quality 
// 2D: length of fragment
// 3D: base qual
// 4D: pos fragment
// 5D: 4X4 matrix
//vector< vector< vector< vector<diNucleotideProb> > > > mpq2Length2BaseQual2Pos2SubMatrix;

// // 1D: mapping quality 
// // 2D: pos away from 5'/3' end
// // 3D: base qual
// vector< vector< vector<diNucleotideProb> > >  mpq2Pos2BaseQual2SubMatrix5p;
// vector< vector< vector<diNucleotideProb> > >  mpq2Pos2BaseQual2SubMatrix3p;



typedef vector< vector<diNucleotideProb> > mpq2bsq2submatrix;

 // 1D: pos away from 5'/3' end
 // 2D: mapping quality 
 // 3D: base qual
vector< mpq2bsq2submatrix >  pos2mpq2BaseQual2SubMatrix5p;
vector< mpq2bsq2submatrix >  pos2mpq2BaseQual2SubMatrix3p;

// 1D: length of fragment
// 2D: pos fragment from the 5' end
vector< vector< mpq2bsq2submatrix * > > length2pos2mpq2bsq2submatrix;


// pointer to data structure
// 3D: mapping quality 
// 4D: base qual
// 5D: 4X4 matrix

//vector< 



//vector< vector< vector<diNucleotideProb> > > mqp;//FOR EACH FRAGMENT LENGTH



    // //dummy values
    // for(int L=0;L<MINLENGTHFRAGMENT;L++){ //for each fragment length
    // 	vector<diNucleotideProb> baseQualEffect;//FOR EACH QUAL SCORE
    // 	fragmentLengthEffect.push_back(baseQualEffect);
    // }

    // //int L=35;//fragment length    
    // for(int L=MINLENGTHFRAGMENT;L<MAXLENGTHFRAGMENT;L++){ //for each fragment length
    // 	cerr<<L<<endl;

    // 	vector<diNucleotideProb> baseQualEffect;//FOR EACH QUAL SCORE




vector<long double> * cov2probPoisson;

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

string genoIdx2Code10 [10] = {"AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"};
string babdIdx2Code16 [16] = {"AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"};
//                              0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
char  geno16_2_geno10 [16] = {  0 ,  1,   2,   3,   1,   4,   5,   6,   2,   5,   7 ,  8 ,  3 ,  6,   8,   9 };

 // genoPriority[0] = 0;//AA
 //    genoPriority[1] = 4; //CC
 //    genoPriority[2] = 7; //GG
 //    genoPriority[3] = 9; //TT

 //    //hetero
 //    genoPriority[4] = 1; //AC
 //    genoPriority[5] = 2; //AG
 //    genoPriority[6] = 3; //AT
 //    genoPriority[7] = 5; //CA
 //    genoPriority[8] = 6; //CG
 //    genoPriority[9] = 8; //CT

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

pthread_mutex_t  mutexQueue           = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  mutexCounter         = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t  mutexRank            = PTHREAD_MUTEX_INITIALIZER;
// pthread_mutex_t  mutexCoverageCounter = PTHREAD_MUTEX_INITIALIZER;


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


    //BASE QUALs, cannot go down to qual zero
    for(int i=0;i<2;i++){
        likeMatch[i]          = log1pl(    -randomPMatch4Bases );          
        likeMismatch[i]       = logl  (     randomPMismatch4Bases );

        likeMatchProb[i]              =    randomPMatch4Bases;    // 1/4
        likeMismatchProb[i]           =    randomPMismatch4Bases; //3/4
    }

    for(int i=2;i<=MAXBASEQUAL;i++){
        likeMatch[i]          = log1pl(    -powl(10.0,i/-10.0) );          
        likeMismatch[i]       = logl  (     powl(10.0,i/-10.0) );

        likeMatchProb[i]              = 1.0-powl(10.0,i/-10.0);
        likeMismatchProb[i]           =     powl(10.0,i/-10.0);
    }

    //mapping QUALs
    for(int i=0;i<=MAXMAPPINGQUAL;i++){
        // likeMatchMap[i]        = log1pl(    -powl(10.0,i/-10.0) );          
        // likeMismatchMap[i]     = logl  (     powl(10.0,i/-10.0) );

        likeMatchProbMap[i]           = 1.0-powl(10.0,i/-10.0);
        likeMismatchProbMap[i]        =     powl(10.0,i/-10.0);
    }


    for(int i=1;i<=MAXCOV;i++){
	//cout<<i<<endl;

	for(int j=0;j<=i;j++){	    
	    binomVec[i][j] = ( logl(nChoosek(i,j))+logl(powl(0.5,i)) );	     
	}

    }

#ifdef DEBUGINITSCORES
    cerr<<"qQUAL"<<"\t"<<"likeMatch"<<"\t"<<"likeMismatch"<<"\t"<<"likeMatchProb"<<"\t"<<"likeMismatchProb"<<endl;
    for(int i=0;i<=MAXBASEQUAL;i++){
	cerr<<"q= "<<i<<"\t"<<likeMatch[i]<<"\t"<<likeMismatch[i]<<"\t"<<likeMatchProb[i]<<"\t"<<likeMismatchProb[i]<<endl;
    }

    cerr<<"mQUAL"<<"\t"<<"likeMatchProb"<<"\t"<<"likeMismatchProb"<<endl;
    for(int i=0;i<=MAXMAPPINGQUAL;i++){
	//cerr<<"m= "<<i<<"\t"<<likeMatchMap[i]<<"\t"<<likeMismatchMap[i]<<"\t"<<likeMatchProbMap[i]<<"\t"<<likeMismatchProbMap[i]<<endl;
	cerr<<"m= "<<i<<"\t"<<likeMatchProbMap[i]<<"\t"<<likeMismatchProbMap[i]<<endl;

    }
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

    //copying to di-nucleotides
    for(unsigned int i=0;i<MAXLENGTHFRAGMENT;i++){     
	diNucleotideProb diNuctoadd;
	for(int nuc1=0;nuc1<4;nuc1++){
	    for(int nuc2=0;nuc2<4;nuc2++){
		int nuc = nuc1*4+nuc2;
		diNuctoadd.p[nuc1][nuc2]=sub5p[i].s[nuc];		
	    }	
	}
	sub5pDiNuc.push_back(diNuctoadd);
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

    //copying to di-nucleotides
    for(unsigned int i=0;i<MAXLENGTHFRAGMENT;i++){     
	diNucleotideProb diNuctoadd;
	for(int nuc1=0;nuc1<4;nuc1++){
	    for(int nuc2=0;nuc2<4;nuc2++){
		int nuc = nuc1*4+nuc2;
		diNuctoadd.p[nuc1][nuc2]=sub3p[i].s[nuc];		
	    }	
	}
	sub3pDiNuc.push_back(diNuctoadd);
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

    cerr<<"------------------"<<endl;

    for(unsigned int i=0;i<sub5p.size();i++){
    	cerr<<"i="<<i<<endl<<"\t";
    	for(int nuc1=0;nuc1<4;nuc1++)
	    cerr<<"ACGT"[nuc1]<<"\t";
	cerr<<endl;
    	for(int nuc1=0;nuc1<4;nuc1++){
	    cerr<<"ACGT"[nuc1]<<"\t";
    	    for(int nuc2=0;nuc2<4;nuc2++){
    		//int nuc = nuc1*4+nuc2;
    		//cout<<nuc1<<"\t"<<nuc2<<"\t"
    		cerr<<sub5pDiNuc[i].p[nuc1][nuc2]<<"\t";
    	    }
    	    cerr<<endl;
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


    cerr<<"------------------"<<endl;

    for(unsigned int i=0;i<sub3p.size();i++){
    	cerr<<"i="<<i<<endl<<"\t";
    	for(int nuc1=0;nuc1<4;nuc1++)
	    cerr<<"ACGT"[nuc1]<<"\t";
	cerr<<endl;
    	for(int nuc1=0;nuc1<4;nuc1++){
	    cerr<<"ACGT"[nuc1]<<"\t";
    	    for(int nuc2=0;nuc2<4;nuc2++){
    		//int nuc = nuc1*4+nuc2;
    		//cout<<nuc1<<"\t"<<nuc2<<"\t"
    		cerr<<sub3pDiNuc[i].p[nuc1][nuc2]<<"\t";
    	    }
    	    cerr<<endl;
    	}
    	cerr<<endl;
    }

#endif




    //if no deamination, cannot have a mismatch
    for(int b1=0;b1<4;b1++){
	for(int b2=0;b2<4;b2++){
	    int b = b1*4+b2;
	    if(b1==b2){
		defaultSubMatch.s[ b ]           = 1.0;	
		defaultSubMatchMatrix.p[b1][b2]  = 1.0;	
	    }else{
		defaultSubMatch.s[ b ]           = 0.0;	
		defaultSubMatchMatrix.p[b1][b2]  = 0.0;	
	    }
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

	for(int b1=0;b1<4;b1++){//    observed base

	    for(int b2=0;b2<4;b2++){ //original base base
		int b = b2*4+b1;
		//cerr<<"pos="<<i<<"\t"<<"ACGT"[b1]<<"\t"<<"ACGT"[b2]<<"\t"<<b<<"\t"<<dnaDefaultBases.f[b1]<<"\t"<<sub5p[i].s[b]<<"\t"<<sub3p[i].s[b]<<endl;
		dnaFreq5p_.f[b1] += dnaDefaultBases.f[b2]*sub5p[i].s[b];
		dnaFreq3p_.f[b1] += dnaDefaultBases.f[b2]*sub3p[i].s[b];
		//cerr<<"\t"<<dnaFreq5p_.f[b1]<<"\t"<<dnaFreq3p_.f[b1]<<endl;
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
	cerr<<defaultDNA5p[i].f[0]<<"\t"<<defaultDNA5p[i].f[1]<<"\t"<<defaultDNA5p[i].f[2]<<"\t"<<defaultDNA5p[i].f[3]<<"\t"<<(defaultDNA5p[i].f[0]+defaultDNA5p[i].f[1]+defaultDNA5p[i].f[2]+defaultDNA5p[i].f[3])<<endl;
    }

    cerr<<"-- 3' base frequencies --"<<endl;

    for(unsigned int i=0;i<MAXLENGTHFRAGMENT;i++){     
    	cerr<<"i="<<i<<" - ";
	cerr<<defaultDNA3p[i].f[0]<<"\t"<<defaultDNA3p[i].f[1]<<"\t"<<defaultDNA3p[i].f[2]<<"\t"<<defaultDNA3p[i].f[3]<<"\t"<<(defaultDNA3p[i].f[0]+defaultDNA3p[i].f[1]+defaultDNA3p[i].f[2]+defaultDNA3p[i].f[3])<<endl;
    }

    
#endif


}//end initDeamProbabilities





//! A method to initialize the likelihood of observing data to avoid recomputation
/*!
  This method is called by the main after calling initDefaultBaseFreq
*/
void initLikelihoodScores(){







//  // 1D: pos away from 5'/3' end
//  // 2D: mapping quality 
//  // 3D: base qual
// vector< vector< vector<diNucleotideProb> > >  pos2mpq2BaseQual2SubMatrix5p;
// vector< vector< vector<diNucleotideProb> > >  pos2mpq2BaseQual2SubMatrix3p;
    
    for(int l=0;l<( (MAXLENGTHFRAGMENT/2) +1);l++){//for each position

	mpq2bsq2submatrix  mpq2BaseQualSubMatrix5p;
	mpq2bsq2submatrix  mpq2BaseQualSubMatrix3p;
	
	for(int mq=0;mq<=MAXMAPPINGQUAL;mq++){ //for each mapping quality 




	    //we need 4 theoretical bases X 4 observed bases X MAXMAPPINGQUAL mapping quality X MAXLENGTHFRAGMENT positions on the fragment X base quality
	    vector<diNucleotideProb>  baseQual2SubMatrix5p;//FOR EACH QUAL SCORE
	    vector<diNucleotideProb>  baseQual2SubMatrix3p;//FOR EACH QUAL SCORE



	    for(int q=0;q<=MAXBASEQUAL;q++){ //for each base quality
		
		//vector<diNucleotideProb> pos2SubMatrix;//FOR EACH QUAL SCORE
		//vector< vector< vector<diNucleotideProb> > > mpq2Length2BaseQual2Pos2SubMatrix;//FOR EACH FRAGMENT LENGTH
		
		
		//we know the probability of P(b1->b2) after deamination is vector<diNucleotideProb> sub5pDiNuc; vector<diNucleotideProb> sub3pDiNuc;
		
		diNucleotideProb toAddForBaseQual5p;
		diNucleotideProb toAddForBaseQual3p;

		diNucleotideProb toAddForBaseQual5p_;
		diNucleotideProb toAddForBaseQual3p_;
		
		for(int bTheo=0;bTheo<4;bTheo++){
		    for(int bObs=0;bObs<4;bObs++){
			toAddForBaseQual5p_.p[bTheo][bObs]=0.0;
			toAddForBaseQual3p_.p[bTheo][bObs]=0.0;
		    }
		}
		    

		for(int bTheo=0;bTheo<4;bTheo++){               // each possible theoritical base
		    for(int bpstDeam=0;bpstDeam<4;bpstDeam++){	// each possible deaminated base		
			for(int bObs=0;bObs<4;bObs++){	        // each possible observed base

			    //cerr<<"tripleloop1\tl="<<l<<"\tmq="<<mq<<"\tq="<<q<<"\t"<<sub5pDiNuc[    l].p[bTheo][bpstDeam] <<"\t"<<likeMatchProb[q] <<"\t"<<   defaultSubMatchMatrix.p[bpstDeam][bObs]<<"\t"<<likeMismatchProb[q] <<"\t"<<illuminaErrorsProbMatrix.p[bpstDeam][bObs]<<"\t"<<toAddForBaseQual5p_.p[bTheo][bObs]<<"\t"<<toAddForBaseQual3p_.p[bTheo][bObs]<<"\t"<<sub5pDiNuc[    l].p[bTheo][bpstDeam]<<"*"<<(		    likeMatchProb[q]    *   defaultSubMatchMatrix.p[bpstDeam][bObs]				    +				    likeMismatchProb[q] * illuminaErrorsProbMatrix.p[bpstDeam][bObs]				)<<"\t"<<				sub3pDiNuc[    l].p[bTheo][bpstDeam]<<" * "<<(				    likeMatchProb[q]    *   defaultSubMatchMatrix.p[bpstDeam][bObs]				    +				    likeMismatchProb[q] * illuminaErrorsProbMatrix.p[bpstDeam][bObs]				)<<endl;

			    toAddForBaseQual5p_.p[bTheo][bObs] +=        
				sub5pDiNuc[    l].p[bTheo][bpstDeam] * (
				    likeMatchProb[q]    *   defaultSubMatchMatrix.p[bpstDeam][bObs]
				    +
				    likeMismatchProb[q] * illuminaErrorsProbMatrix.p[bpstDeam][bObs]
				);
			    
			    toAddForBaseQual3p_.p[bTheo][bObs] += 
				sub3pDiNuc[    l].p[bTheo][bpstDeam] * (
				    likeMatchProb[q]    *   defaultSubMatchMatrix.p[bpstDeam][bObs]
				    +
				    likeMismatchProb[q] * illuminaErrorsProbMatrix.p[bpstDeam][bObs]
				);

			    // cerr<<"tripleloop2\tl="<<l<<"\tmq="<<mq<<"\tq="<<q<<"\t"<<sub5pDiNuc[    l].p[bTheo][bpstDeam] <<"\t"<<likeMatchProb[q] <<"\t"<<   defaultSubMatchMatrix.p[bpstDeam][bObs]<<"\t"<<likeMismatchProb[q] <<"\t"<<illuminaErrorsProbMatrix.p[bpstDeam][bObs]<<"\t"<<toAddForBaseQual5p_.p[bTheo][bObs]<<"\t"<<toAddForBaseQual3p_.p[bTheo][bObs]<<endl;

			}//end each bObs
		    }//end each bpstDeam
		}//for each bTheo


		
		
		for(int bTheo=0;bTheo<4;bTheo++){               // each possible theoretical base
		    for(int bObs=0;bObs<4;bObs++){	        // each possible observed base
			
			//0.5 since this is the prior prob of having sampled a given chromosome
#ifdef PRECOMPUTELOG
			toAddForBaseQual5p.p[bTheo][bObs] = logl(0.5 * (likeMatchProbMap[mq]*toAddForBaseQual5p_.p[bTheo][bObs] + likeMismatchProbMap[mq] * defaultDNA5p[l].f[bObs]) );			
			toAddForBaseQual3p.p[bTheo][bObs] = logl(0.5 * (likeMatchProbMap[mq]*toAddForBaseQual3p_.p[bTheo][bObs] + likeMismatchProbMap[mq] * defaultDNA3p[l].f[bObs]) );

// 			if(mq==37 && q==38 && bObs==3 && bTheo == 0 && l==5 ){
			    
// 			    cerr<<setprecision(20)<<"doubleloop\tl="<<l<<"\tmq="<<mq<<"\tq="<<q<<"\t"<<toAddForBaseQual5p.p[bTheo][bObs]<<"\t=\t"<<likeMatchProbMap[mq]<<"*"<<toAddForBaseQual5p_.p[bTheo][bObs] <<"+"<< likeMismatchProbMap[mq]<<"*"<< defaultDNA5p[l].f[bObs]

// 				<<"="<<(likeMatchProbMap[mq]*toAddForBaseQual5p_.p[bTheo][bObs] + likeMismatchProbMap[mq] * defaultDNA5p[l].f[bObs] )
// 				<<"\t"<<logl(likeMatchProbMap[mq]*toAddForBaseQual5p_.p[bTheo][bObs] + likeMismatchProbMap[mq] * defaultDNA5p[l].f[bObs] )			
// 				<<"\t"<<logl(likeMatchProbMap[mq]*toAddForBaseQual5p_.p[bTheo][bObs] + likeMismatchProbMap[mq] * defaultDNA5p[l].f[bObs] )			
// 				<<"\t"<<toAddForBaseQual5p.p[bTheo][bObs]<<"\t"<<expl(toAddForBaseQual5p.p[bTheo][bObs])
// <<endl;
// 			}

#else
			// without logl
			toAddForBaseQual5p.p[bTheo][bObs] =  likeMatchProbMap[mq]*toAddForBaseQual5p_.p[bTheo][bObs] + likeMismatchProbMap[mq] * defaultDNA5p[l].f[bObs];			
			toAddForBaseQual3p.p[bTheo][bObs] =  likeMatchProbMap[mq]*toAddForBaseQual3p_.p[bTheo][bObs] + likeMismatchProbMap[mq] * defaultDNA3p[l].f[bObs];
			
			
			// if(mq==37 && q==38 && bObs==3 && bTheo == 0 && l==5 ){
			    
			//     cerr<<setprecision(20)<<"doubleloop\tl="<<l<<"\tmq="<<mq<<"\tq="<<q<<"\t"<<toAddForBaseQual5p.p[bTheo][bObs]<<"\t=\t"<<likeMatchProbMap[mq]<<"*"<<toAddForBaseQual5p_.p[bTheo][bObs] <<"+"<< likeMismatchProbMap[mq]<<"*"<< defaultDNA5p[l].f[bObs]<<"\t"<<toAddForBaseQual5p.p[bTheo][bObs]<<endl;
			// }
			
#endif
			// cerr<<"doubleloop\tl="<<l<<"\tmq="<<mq<<"\tq="<<q<<"\t"<<toAddForBaseQual5p.p[bTheo][bObs]<<"\t=\t"<<likeMatchProbMap[mq]<<"*"<<toAddForBaseQual5p_.p[bTheo][bObs] <<"+"<< likeMismatchProbMap[mq]<<"*"<< defaultDNA5p[l].f[bObs]<<endl;
			// cerr<<"doubleloop\tl="<<l<<"\tmq="<<mq<<"\tq="<<q<<"\t"<<toAddForBaseQual3p.p[bTheo][bObs]<<"\t=\t"<<likeMatchProbMap[mq]<<"*"<<toAddForBaseQual3p_.p[bTheo][bObs] <<"+"<< likeMismatchProbMap[mq]<<"*"<< defaultDNA3p[l].f[bObs]<<endl;
			
		    }//for each bObs
		}//for each bTheo, bpstDeam and bObs
		
		baseQual2SubMatrix5p.push_back(toAddForBaseQual5p);
		baseQual2SubMatrix3p.push_back(toAddForBaseQual3p);
		
		
#ifdef DEBUGINITLIKELIHOODSCORES
		cerr<<"ppos = "<<l<<" mq = "<<mq<<" ("<<likeMatchProbMap[mq]<<" "<<likeMismatchProbMap[mq]<<" )  q = "<<q<<endl;
		cerr<<"5'----------"<<endl;
		for(int nuc1=0;nuc1<4;nuc1++){
		    cerr<<"ACGT"[nuc1]<<"\t";		    
		    long double s=0;
		    for(int nuc2=0;nuc2<4;nuc2++){
			cerr<<toAddForBaseQual5p_.p[nuc1][nuc2]<<"\t";			
			s+=toAddForBaseQual5p_.p[nuc1][nuc2];
		    }
		    cerr<<"\t"<<s<<"\t";
		    s=0;
		    
		    for(int nuc2=0;nuc2<4;nuc2++){
			cerr<<defaultDNA5p[l].f[nuc2]<<"\t";			
			s+=defaultDNA5p[l].f[nuc2];
		    }

		    cerr<<"\t"<<s<<"\t";
		    s=0;
		    
		    for(int nuc2=0;nuc2<4;nuc2++){
			cerr<<toAddForBaseQual5p.p[nuc1][nuc2]<<"\t";
			s+=toAddForBaseQual5p.p[nuc1][nuc2];
		    }

		    cerr<<"\t"<<s<<endl;
		}
		cerr<<endl;
		cerr<<"3'----------"<<endl;
		for(int nuc1=0;nuc1<4;nuc1++){
		    cerr<<"ACGT"[nuc1]<<"\t";
		    long double s=0;

		    for(int nuc2=0;nuc2<4;nuc2++){
			cerr<<toAddForBaseQual3p_.p[nuc1][nuc2]<<"\t";
			s+=toAddForBaseQual3p_.p[nuc1][nuc2];
		    }
		    cerr<<"\t"<<s<<"\t";
		    s=0;
		    
		    for(int nuc2=0;nuc2<4;nuc2++){
			cerr<<defaultDNA3p[l].f[nuc2]<<"\t";
			s+=defaultDNA3p[l].f[nuc2];			
		    }
		    cerr<<"\t"<<s<<"\t";
		    s=0;
		    
		    for(int nuc2=0;nuc2<4;nuc2++){
			cerr<<toAddForBaseQual3p.p[nuc1][nuc2]<<"\t";
			s+=toAddForBaseQual3p.p[nuc1][nuc2];
		    }
		    cerr<<"\t"<<s<<endl;
		}
		cerr<<endl;
#endif
		    
	    }//for each base qual
	    


	    mpq2BaseQualSubMatrix5p.push_back(baseQual2SubMatrix5p);
	    mpq2BaseQualSubMatrix3p.push_back(baseQual2SubMatrix3p);
	    
	    
	}// for each mapping quality

	pos2mpq2BaseQual2SubMatrix5p.push_back(mpq2BaseQualSubMatrix5p);
	pos2mpq2BaseQual2SubMatrix3p.push_back(mpq2BaseQualSubMatrix3p);

	
    }//for each position
    


    //if(mq==37 && q==38 && bObs==0 && bTheo == 0 && l==5 ){
    //    cerr<<setprecision(20)<<"TEST1 "<<pos2mpq2BaseQual2SubMatrix5p[  12  ][37][38].p[0][3]<<"\t"<<pos2mpq2BaseQual2SubMatrix5p[  12  ][37][18].p[0][3]<<endl;
	//cerr<<"doubleloop\tl="<<l<<"\tmq="<<mq<<"\tq="<<q<<"\t"<<toAddForBaseQual5p.p[bTheo][bObs]<<"\t=\t"<<likeMatchProbMap[mq]<<"*"<<toAddForBaseQual5p_.p[bTheo][bObs] <<"+"<< likeMismatchProbMap[mq]<<"*"<< defaultDNA5p[l].f[bObs]<<"\t"<<toAddForBaseQual5p.p[bTheo][bObs]<<endl;
    //}
    

    //vector< vector< mpq2bsq2submatrix * > > length2pos2mpq2bsq2submatrix;
    //dummy values
    for(int L=0;L<MINLENGTHFRAGMENT;L++){//for each fragment length
	vector< mpq2bsq2submatrix * > vectorToAdd;	
	length2pos2mpq2bsq2submatrix.push_back(vectorToAdd);
    }

    for(int L=MINLENGTHFRAGMENT;L<=MAXLENGTHFRAGMENT;L++){//for each fragment length
	vector< mpq2bsq2submatrix * > vectorToAdd;
	for(int l=0;l<L;l++){//for each pos
	    //cerr<<l<<" "<<(L-l-1)<<" "<<L<<endl;	   
	    
	    if( l<(L/2) ){//use 5' substitutions
		vectorToAdd.push_back( &pos2mpq2BaseQual2SubMatrix5p[  l  ] );
		//cerr<< "5' size "<<pos2mpq2BaseQual2SubMatrix5p[  l  ].size() <<" sizemq0 "<<pos2mpq2BaseQual2SubMatrix5p[  l  ][0].size()<<endl;
	    }else{        //use 3' substitutions
		//cerr<< "3' size "<<pos2mpq2BaseQual2SubMatrix3p[  l  ].size() <<" sizemq0 "<<pos2mpq2BaseQual2SubMatrix3p[  l  ][0].size()<<endl;
		vectorToAdd.push_back( &pos2mpq2BaseQual2SubMatrix3p[L-l-1] );
	    }

	}
	length2pos2mpq2bsq2submatrix.push_back(vectorToAdd);
    }//for each fragment length
    
    // cerr<<setprecision(20)<<"TEST2 AT 38 "<<pos2mpq2BaseQual2SubMatrix5p[  12  ][37][38].p[0][3]<<"\tL2P="<<length2pos2mpq2bsq2submatrix[113][12]->at(37)[38].p[0][3]<<endl ;
    // cerr<<setprecision(20)<<"TEST2 AT 18 "<<pos2mpq2BaseQual2SubMatrix5p[  12  ][37][18].p[0][3]<<"\tL2P="<<length2pos2mpq2bsq2submatrix[113][12]->at(37)[18].p[0][3]<<endl ;

#ifdef DEBUGINITLIKELIHOODSCORES2

    for(int L=MINLENGTHFRAGMENT;L<=MAXLENGTHFRAGMENT;L++){//for each fragment length

	
	for(int l=0;l<L;l++){//for each pos
	    // cerr<<"pos "<<l<<"/"<<L<<endl;//<<"\t"<<length2pos2mpq2bsq2submatrix[L][l]->size()<<endl;
	    
	    for(int mq=0;mq<=MAXMAPPINGQUAL;mq++){ //for each mapping quality 
		//cerr<<"mq="<<mq<<"\t"<<endl;//length2pos2mpq2bsq2submatrix[L][l]->at(mq).size()<<endl;

		for(int q=0;q<=MAXBASEQUAL;q++){ //for each base quality
		    cerr<<"pos "<<l<<"/"<<L<<" mq="<<mq<<" bq="<<q<<endl;

		    //cerr<<" mq = "<<mq<<" ("<<likeMatchProbMap[mq]<<" "<<likeMismatchProbMap[mq]<<" )  q = "<<q<<endl;
		    for(int nuc1=0;nuc1<4;nuc1++){
			cerr<<"ACGT"[nuc1]<<"\t";		    
			long double s=0;//sum of probs

			for(int nuc2=0;nuc2<4;nuc2++){
			    cerr<<length2pos2mpq2bsq2submatrix[L][l]->at(mq)[q].p[nuc1][nuc2]<<"\t";
			    //length2pos2mpq2bsq2submatrix[L][l]->at(mq)
			    //toAddForBaseQual5p.p[nuc1][nuc2]<<"\t";
			    s+=length2pos2mpq2bsq2submatrix[L][l]->at(mq)[q].p[nuc1][nuc2];
			}

			cerr<<"\t"<<s<<endl;
		    }
		    cerr<<endl;
		    
		}//for each base qual	   		    

	    }// for each mapping quality


	    

	}
    }//for each fragment length

#endif    

    //cerr<<"done "<<endl;

    //exit(1);
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





inline void computeLL(vector<positionInformation> * piForGenomicWindow){

    cerr<<"computeLL "<<piForGenomicWindow->size()<<endl;

    long double h=0.00000001;
    //long double h=0.000824;
    //long double h=0.000735;
    //long double h=0.010000;
    diNucleotideProb priorGenotype;
    //compute prior genotype matrix
    
    for(int ba=0;ba<4;ba++){//ancestral base

	for(int bd=0;bd<4;bd++){//derived base
	    if(ba == bd){
		//priorGenotype.p[ba][bd]     = dnaDefaultBases.f[ba]  *    (1.0-h);
		priorGenotype.p[ba][bd]     = logl(dnaDefaultBases.f[ba])  +    logl(1.0-h);
		//cerr<<"ACGT"[ba]<<"\t"<<"ACGT"[bd]<<"\t"<<dnaDefaultBases.f[ba]<<" X "<<(1.0-h)<<endl;
	    }else{//mutation
		if( (ba%2)==(bd%2) ){//transition
		    //priorGenotype.p[ba][bd] = dnaDefaultBases.f[ba]  *  ( (h) * (TStoTVratio/(TStoTVratio+1.0)) );
		    priorGenotype.p[ba][bd] = logl(dnaDefaultBases.f[ba])  +   logl( (h) * (TStoTVratio/(TStoTVratio+1.0))     );
		    //cerr<<"ACGT"[ba]<<"\t"<<"ACGT"[bd]<<"\t"<<dnaDefaultBases.f[ba]<<" X "<< ( (h) * (TStoTVratio/(TStoTVratio+1.0)) )<<endl; 
		}else{
		    //priorGenotype.p[ba][bd] = dnaDefaultBases.f[ba]  * (( (h) * (        1.0/(TStoTVratio+1.0)) )/2.0);
		    priorGenotype.p[ba][bd] = logl(dnaDefaultBases.f[ba])  +   logl( (h) * (        1.0/(TStoTVratio+1.0)) /2.0);
		    //cerr<<"ACGT"[ba]<<"\t"<<"ACGT"[bd]<<"\t"<<dnaDefaultBases.f[ba]<<" X "<< ( (h) * (        1.0/(TStoTVratio+1.0)) )/2.0<<endl; 
		}
	    }
	}
    }
	
    // cerr<<setprecision(20)<<"TEST3 "<<"\tL2P="<<length2pos2mpq2bsq2submatrix[113][12]->at(37)[38].p[0][3]<<"\t"<<length2pos2mpq2bsq2submatrix[113][12]->at(37)[18].p[0][3]<<"\t"<<expl(length2pos2mpq2bsq2submatrix[113][12]->at(37)[38].p[0][3])<<"\t"<<expl(length2pos2mpq2bsq2submatrix[113][12]->at(37)[18].p[0][3])<<"\t"<<logl(0.5*length2pos2mpq2bsq2submatrix[113][12]->at(37)[38].p[0][3])<<"\t"<<logl(0.5*length2pos2mpq2bsq2submatrix[113][12]->at(37)[18].p[0][3])<<endl ;

#ifdef DEBUGCOMPUTELL

    cerr<<"h="<<h<<endl;
    cerr<<"\t";
    for(int ba=0;ba<4;ba++)
	cerr<<"ACGT"[ba]<<"\t";
    cerr<<endl;
    long double sumProb_=0.0;
    for(int ba=0;ba<4;ba++){
	cerr<<"ACGT"[ba]<<"\t";
	for(int bd=0;bd<4;bd++){
	    cerr<<expl(priorGenotype.p[ba][bd])<<"\t";
	    sumProb_+=expl(priorGenotype.p[ba][bd]);
	}
	cerr<<endl;
    }
    cerr<<endl;
    cerr<<"sum prior for each geno = "<<sumProb_<<endl;	   
#endif


    /////////////////////////////////////////
    //BEGIN pre-computing babdlikelihood
    /////////////////////////////////////////

    vector<babdlikelihood> vectorBaBdLikelihood;

    for(unsigned int p=0;p<piForGenomicWindow->size();p++){
	
	babdlikelihood bblikeForGivenPos;

	int         babdIdx          =0;

       	
	for(uint8_t ba=0;ba<4;ba++){//ancestral base
	    uint8_t ba_c = 3-ba;

	    for(uint8_t bd=0;bd<4;bd++){//derived base
		uint8_t bd_c = 3-bd;

				
		long double loglikelihoodForGivenBaBd          =0.0;


		
		for(unsigned int i=0;i<piForGenomicWindow->at(p).readsVec.size();i++){ //for each fragment at pos p
		    //Likelihood it comes from A
		    //char bObs=piForGenomicWindow->at(p).readsVec[i].base;
		    


#ifdef PRECOMPUTELOG
		    long double llA; //Likelihood it comes from A
		    long double llD; //Likelihood it comes from D

		    if(piForGenomicWindow->at(p).readsVec[i].isrv){
			llA = length2pos2mpq2bsq2submatrix
			    [piForGenomicWindow->at(p).readsVec[i].lengthF]
			    [piForGenomicWindow->at(p).readsVec[i].pos5p]->at( piForGenomicWindow->at(p).readsVec[i].mapq)
			    [piForGenomicWindow->at(p).readsVec[i].qual].p[ba_c][piForGenomicWindow->at(p).readsVec[i].base];

			llD = length2pos2mpq2bsq2submatrix
			    [piForGenomicWindow->at(p).readsVec[i].lengthF]
			    [piForGenomicWindow->at(p).readsVec[i].pos5p]->at( piForGenomicWindow->at(p).readsVec[i].mapq)
			    [piForGenomicWindow->at(p).readsVec[i].qual].p[bd_c][piForGenomicWindow->at(p).readsVec[i].base];


		    }else{

			llA = length2pos2mpq2bsq2submatrix
			    [piForGenomicWindow->at(p).readsVec[i].lengthF]
			    [piForGenomicWindow->at(p).readsVec[i].pos5p]->at( piForGenomicWindow->at(p).readsVec[i].mapq)
			    [piForGenomicWindow->at(p).readsVec[i].qual].p[ba][piForGenomicWindow->at(p).readsVec[i].base];

			llD = length2pos2mpq2bsq2submatrix
			    [piForGenomicWindow->at(p).readsVec[i].lengthF]
			    [piForGenomicWindow->at(p).readsVec[i].pos5p]->at( piForGenomicWindow->at(p).readsVec[i].mapq)
			    [piForGenomicWindow->at(p).readsVec[i].qual].p[bd][piForGenomicWindow->at(p).readsVec[i].base];


		    }
			



#else
		    
		    // //without precompute log
		    long double llA = logl(0.5* length2pos2mpq2bsq2submatrix
		    			   [piForGenomicWindow->at(p).readsVec[i].lengthF]
		    			   [piForGenomicWindow->at(p).readsVec[i].pos5p]->at( piForGenomicWindow->at(p).readsVec[i].mapq)
		    			   [piForGenomicWindow->at(p).readsVec[i].qual].p[ba][piForGenomicWindow->at(p).readsVec[i].base]);

		    //Likelihood it comes from D
		    long double llD = logl(0.5* length2pos2mpq2bsq2submatrix
		    			   [piForGenomicWindow->at(p).readsVec[i].lengthF]
		    			   [piForGenomicWindow->at(p).readsVec[i].pos5p]->at( piForGenomicWindow->at(p).readsVec[i].mapq)
		    			   [piForGenomicWindow->at(p).readsVec[i].qual].p[bd][piForGenomicWindow->at(p).readsVec[i].base]);
#endif

		    
			



		    
		    loglikelihoodForGivenBaBd += oplusnatl( llA, llD); //, adding probs of llA and llD, multiplying probabilities for each site, assuming independence, (\prod_{fragment} P(D|G))

		}//END  for each fragment at pos p

		//TODO, pre compute the log the prior
		//  product of (\prod_{fragment} P(D|G)) times the prior P(G) for the genotype
		
		bblikeForGivenPos.gl[babdIdx] = loglikelihoodForGivenBaBd;
		

		

		babdIdx++;
	    }//END for each derived base
	}//END for each ancestral base


	vectorBaBdLikelihood.push_back(  bblikeForGivenPos );
    }//END for each genomic position

    /////////////////////////////////////////
    //END pre-computing babdlikelihood
    /////////////////////////////////////////

    cerr<<"computeLL done computing babdlikelihood "<<piForGenomicWindow->size()<<endl;

    long double loglikelihoodForEveryPositionForEveryBaBd          =0.0;


    for(unsigned int p=0;p<piForGenomicWindow->size();p++){//every genomic position
	

#ifdef DEBUGCOMPUTELL

	//if(p>10000 && p<11000){
	    
	cerr<<p<<"\tpos="<<piForGenomicWindow->at(p).posAlign<<"\t"<<piForGenomicWindow->at(p).readsVec.size()<<endl;

	cerr<<"B\tQ\tMQ\t5p\tL"<<endl;

	for(unsigned int i=0;i<piForGenomicWindow->at(p).readsVec.size();i++){
	    cerr<<"ACGT"[piForGenomicWindow->at(p).readsVec[i].base]<<"\t"
		<<int(piForGenomicWindow->at(p).readsVec[i].qual)<<"\t"
		<<int(piForGenomicWindow->at(p).readsVec[i].mapq)<<"\t"
		<<int(piForGenomicWindow->at(p).readsVec[i].pos5p)<<"\t"
		<<int(piForGenomicWindow->at(p).readsVec[i].lengthF)<<"\t"
		<<piForGenomicWindow->at(p).readsVec[i].isrv<<"\t"
		//		<<piForGenomicWindow->at(p).readsVec[i].name<<"\t"
		<<endl;	    
	}
	//}    
#endif




	// typedef vector< vector<diNucleotideProb> > mpq2bsq2submatrix;
	
	
	//  // 2D: mapping quality 
	//  // 3D: base qual
	// vector< mpq2bsq2submatrix >  pos2mpq2BaseQual2SubMatrix5p;
	// vector< mpq2bsq2submatrix >  pos2mpq2BaseQual2SubMatrix3p;
	
	// // 1D: length of fragment
	// // 2D: pos fragment from the 5' end
	// vector< vector< mpq2bsq2submatrix * > > length2pos2mpq2bsq2submatrix;
	long double loglikelihoodForEveryBaBd          =0.0;
	vector<long double> vectorOfloglikelihoodForGivenBaBd     (16,0.0) ;
	vector<long double> vectorOfloglikelihoodForGivenGeno     (10,0.0) ;

	long double mostLikelyBaBd   =-1.0*numeric_limits<long double>::infinity();
	int         mostLikelyBaBdIdx=-1;
	int         babdIdx          =0;


	// 	for(unsigned int i=0;i<piForGenomicWindow->at(p).readsVec.size();i++){ //for each fragment at pos p
	// 		    //Likelihood it comes from A

	// // #ifdef DEBUGCOMPUTELLEACHBASE
	// // 	    //if(p>10000 && p<11000){
	// // 	    cerr<<"ACGT"[piForGenomicWindow->at(p).readsVec[i].base]<<"\t"
	// // 		<<"Q="<<int(piForGenomicWindow->at(p).readsVec[i].qual)<<"\t"
	// // 		<<"M="<<int(piForGenomicWindow->at(p).readsVec[i].mapq)<<"\t"
	// // 		<<"5="<<int(piForGenomicWindow->at(p).readsVec[i].pos5p)<<"\t"
	// // 		<<"L="<<int(piForGenomicWindow->at(p).readsVec[i].lengthF)<<"\t"
	// // 		<<"R="<<piForGenomicWindow->at(p).readsVec[i].isrv<<"\t"
	// // 		//		<<piForGenomicWindow->at(p).readsVec[i].name<<"\t"
	// // 		<<endl;
	// // 	    //cerr<<length2pos2mpq2bsq2submatrix[piForGenomicWindow->at(p).readsVec[i].lengthF][piForGenomicWindow->at(p).readsVec[i].pos5p]->size()<<endl; 
	// // 	    //}
	// // #endif
	// 	}
	





	for(uint8_t ba=0;ba<4;ba++){//ancestral base
	    uint8_t ba_c = 3-ba;

	    for(uint8_t bd=0;bd<4;bd++){//derived base
		uint8_t bd_c = 3-bd;

#ifdef DEBUGCOMPUTELLEACHBASE
		//if(p>10000 && p<11000)
		cerr<<endl<<"GENO:A="<<"ACGT"[ba]<<" ("<<"ACGT"[ba_c]<<") \tD="<<"ACGT"[bd]<<" ("<<"ACGT"[bd_c]<<")\tprior "<<expl(priorGenotype.p[ba][bd])<<endl;
		//cerr<<int(ba)<<"\t"<<int(ba_c)<<endl;
#endif
		
		long double loglikelihoodForGivenBaBdTimesPrior=0.0;
		long double loglikelihoodForGivenBaBd          =0.0;

			       	      
		//  product of (\prod_{fragment} P(D|G)) times the prior P(G) for the genotype
		//loglikelihoodForGivenBaBdTimesPrior = loglikelihoodForGivenBaBd + priorGenotype.p[ba][bd]; // (\prod_{fragment} P(D|G))*P(G)
		loglikelihoodForGivenBaBdTimesPrior = vectorBaBdLikelihood[p].gl[babdIdx] + priorGenotype.p[ba][bd]; // (\prod_{fragment} P(D|G))*P(G)
		
		
#ifdef DEBUGCOMPUTELLEACHBASE
		//if(p>10000 && p<11000)
		cerr<<"GENO:A="<<"ACGT"[ba]<<"\tD="<<"ACGT"[bd]<<"\tprior "<<expl(priorGenotype.p[ba][bd])<<"\tllForBaBD "<<loglikelihoodForGivenBaBd<<"\tllForBaBD*Prior "<<loglikelihoodForGivenBaBdTimesPrior<<"\tllForEveryBaBD "<<loglikelihoodForEveryBaBd<<"\tgenoLike "<<mostLikelyBaBd<<"\tmostLikeG "<<mostLikelyBaBdIdx<<endl;
#endif


		//adding probabilities for each genotype
		loglikelihoodForEveryBaBd  = oplusInitnatl( loglikelihoodForEveryBaBd , loglikelihoodForGivenBaBdTimesPrior); // \sum_{genotype} (\prod_{fragment} P(D|G))*P(G)

		vectorOfloglikelihoodForGivenBaBd[babdIdx] = loglikelihoodForGivenBaBdTimesPrior ;

		if(loglikelihoodForGivenBaBdTimesPrior>mostLikelyBaBd){
		    mostLikelyBaBd    = loglikelihoodForGivenBaBdTimesPrior;
		    mostLikelyBaBdIdx = babdIdx;
		}
		
#ifdef DEBUGCOMPUTELL
		//if(p>10000 && p<11000)
		cerr<<endl<<"---------------------------------------------------------------------------------"<<endl;
#endif
		babdIdx++;
	    }//END for each derived base
	}//END for each ancestral base



#ifdef DEBUGCOMPUTELL	
	//if(p>10000 && p<11000)
	cerr<<endl<<"---------------------------------------------------------------------------------"<<endl;
#endif


	//
	// BEGIN GENOTYPING
	//

	//Compute likelihood of all minus best BABD
	long double loglikelihoodForEveryBaBd_minusBest          =0.0;
	
	for(int g=0;g<16;g++){
	    if(g!= mostLikelyBaBdIdx)
		loglikelihoodForEveryBaBd_minusBest = oplusInitnatl( loglikelihoodForEveryBaBd_minusBest , vectorOfloglikelihoodForGivenBaBd[g] ); // \sum_{genotype} (\prod_{fragment} P(D|G))*P(G)
	    //cerr<<babdIdx2Code16[g]<<"\t"<<vectorOfloglikelihoodForGivenBaBd[g]<<endl;
	    //cerr<<vectorToString(vectorOfloglikelihoodForGivenBaBd,"\t")<<endl;
	}
	
	
	// if( (mostLikelyBaBdIdx !=  0) && 
	//     (mostLikelyBaBdIdx !=  5) && 
	//     (mostLikelyBaBdIdx != 10) && 
	//     (mostLikelyBaBdIdx != 15) ){
	//     //exit(1);
	//     cerr<<endl<<"---------------------------------------------------------------------------------"<<endl;
	//     cerr<<p<<"\tpos="<<piForGenomicWindow->at(p).posAlign<<"\t"<<piForGenomicWindow->at(p).readsVec.size()<<endl;

	//     for(unsigned int g=0;g<16;g++){
	// 	cerr<<babdIdx2Code16[g]<<"\t"<<vectorOfloglikelihoodForGivenBaBd[g]<<endl;
	// 	//cerr<<vectorToString(vectorOfloglikelihoodForGivenBaBd,"\t")<<endl;
	//     }

	//     cerr<<mostLikelyBaBd<<"\t"<<babdIdx2Code16[mostLikelyBaBdIdx]<<"\t"<<loglikelihoodForEveryBaBd_minusBest <<"\t"<< loglikelihoodForEveryBaBd<<"\t"<<( loglikelihoodForEveryBaBd_minusBest - loglikelihoodForEveryBaBd)<<endl;

	// } 

	//add to 10 genotypes
	// string babdIdx2Code10 [10] = {"AA","AC","AG","AT","CC","CG","CT","GG","GT","TT"};
	// string babdIdx2Code16 [16] = {"AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"};
	//                                 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
	//---------------------------------------------
	// BEGIN convert to 10 genotypes
	//---------------------------------------------
	//homo
	vectorOfloglikelihoodForGivenGeno[0] = vectorOfloglikelihoodForGivenBaBd[ 0];
	vectorOfloglikelihoodForGivenGeno[4] = vectorOfloglikelihoodForGivenBaBd[ 5];
	vectorOfloglikelihoodForGivenGeno[7] = vectorOfloglikelihoodForGivenBaBd[10];
	vectorOfloglikelihoodForGivenGeno[9] = vectorOfloglikelihoodForGivenBaBd[15];

	//hetero
	vectorOfloglikelihoodForGivenGeno[ 1] = oplusnatl(vectorOfloglikelihoodForGivenBaBd[ 1], vectorOfloglikelihoodForGivenBaBd[ 4]);
	vectorOfloglikelihoodForGivenGeno[ 2] = oplusnatl(vectorOfloglikelihoodForGivenBaBd[ 2], vectorOfloglikelihoodForGivenBaBd[ 8]);
	vectorOfloglikelihoodForGivenGeno[ 3] = oplusnatl(vectorOfloglikelihoodForGivenBaBd[ 3], vectorOfloglikelihoodForGivenBaBd[12]);
	vectorOfloglikelihoodForGivenGeno[ 5] = oplusnatl(vectorOfloglikelihoodForGivenBaBd[ 6], vectorOfloglikelihoodForGivenBaBd[ 9]);
	vectorOfloglikelihoodForGivenGeno[ 6] = oplusnatl(vectorOfloglikelihoodForGivenBaBd[ 7], vectorOfloglikelihoodForGivenBaBd[13]);
	vectorOfloglikelihoodForGivenGeno[ 8] = oplusnatl(vectorOfloglikelihoodForGivenBaBd[11], vectorOfloglikelihoodForGivenBaBd[14]);

	long double mostLikelyGeno   =-1.0*numeric_limits<long double>::infinity();
	int         mostLikelyGenoIdx=-1;
	long double loglikelihoodForEveryGeno          =0.0;
	long double loglikelihoodForEveryGeno_minusBest=0.0;

	for(unsigned int g=0;g<10;g++){	   
	    if(vectorOfloglikelihoodForGivenGeno[g] > mostLikelyGeno){
		mostLikelyGeno    = vectorOfloglikelihoodForGivenGeno[g];
		mostLikelyGenoIdx = g;
	    }
	}
	
	for(int g=0;g<10;g++){
	    if(g != mostLikelyGenoIdx)
		loglikelihoodForEveryGeno_minusBest = oplusInitnatl( loglikelihoodForEveryGeno_minusBest , vectorOfloglikelihoodForGivenGeno[g] ); // 
	    loglikelihoodForEveryGeno               = oplusInitnatl( loglikelihoodForEveryGeno           , vectorOfloglikelihoodForGivenGeno[g] ); // 
	}

	//---------------------------------------------
	// END convert to 10 genotypes
	//---------------------------------------------

	//	if(true){
	if( (mostLikelyGenoIdx !=  0) && 
	    (mostLikelyGenoIdx !=  4) && 
	    (mostLikelyGenoIdx !=  7) && 
	    (mostLikelyGenoIdx !=  9) ){

#ifdef DEBUGCOMPUTELLGENO
	    cerr<<endl<<"---------------------------------------------------------------------------------"<<endl;
	    cerr<<p<<"\tpos="<<piForGenomicWindow->at(p).posAlign<<"\t"<<piForGenomicWindow->at(p).readsVec.size()<<endl;
	    for(unsigned int i=0;i<piForGenomicWindow->at(p).readsVec.size();i++){ //for each fragment at pos p
	    
		cerr<<"ACGT"[piForGenomicWindow->at(p).readsVec[i].base]<<"\t"<<
		    "ACGT"[piForGenomicWindow->at(p).readsVec[i].isrv?(3-piForGenomicWindow->at(p).readsVec[i].base):piForGenomicWindow->at(p).readsVec[i].base]<<"\t"
		    <<"Q="<<int(piForGenomicWindow->at(p).readsVec[i].qual)<<"\t"
		    <<"M="<<int(piForGenomicWindow->at(p).readsVec[i].mapq)<<"\t"
		    <<"5="<<int(piForGenomicWindow->at(p).readsVec[i].pos5p)<<"\t"
		    <<"L="<<int(piForGenomicWindow->at(p).readsVec[i].lengthF)<<"\t"
		    <<"R="<<piForGenomicWindow->at(p).readsVec[i].isrv<<"\t"
		    //		<<piForGenomicWindow->at(p).readsVec[i].name<<"\t"
		    <<endl;		
	    }
	    for(unsigned int g=0;g<16;g++){
		cerr<<babdIdx2Code16[g]<<"\t"<<vectorOfloglikelihoodForGivenBaBd[g]<<endl;
	    }
	    cerr<<endl<<"---------------------------------------------------------------------------------"<<endl;
	    for(unsigned int g=0;g<10;g++){	       
		cerr<<genoIdx2Code10[g]<<"\t"<<vectorOfloglikelihoodForGivenGeno[g]<<endl;
	    }
	    cerr<<endl<<"---------------------------------------------------------------------------------"<<endl;
	    
	    cerr<<p<<"\tpos="<<piForGenomicWindow->at(p).posAlign<<"\tmostLikelyBaBd\t"<<mostLikelyBaBd<<"\t"<<babdIdx2Code16[mostLikelyBaBdIdx]<<"\t"<<loglikelihoodForEveryBaBd_minusBest <<"\t"<< loglikelihoodForEveryBaBd<<"\t"<<( loglikelihoodForEveryBaBd_minusBest - loglikelihoodForEveryBaBd)<<"\t"<<piForGenomicWindow->at(p).readsVec.size()<<endl;
	    cerr<<p<<"\tpos="<<piForGenomicWindow->at(p).posAlign<<"\tmostLikelyGeno\t"<<mostLikelyGeno<<"\t"<<genoIdx2Code10[mostLikelyGenoIdx]<<"\t"<<loglikelihoodForEveryGeno_minusBest <<"\t"<< loglikelihoodForEveryGeno<<"\t"<<( loglikelihoodForEveryGeno_minusBest - loglikelihoodForEveryGeno)<<"\t"<<piForGenomicWindow->at(p).readsVec.size()<<endl;
	    cerr<<endl<<"---------------------------------------------------------------------------------"<<endl;

#endif
	}//end for each geno


	//
	// END GENOTYPING
	//
	//exit(1);
	//     //exit(1);
	//     cerr<<endl<<"---------------------------------------------------------------------------------"<<endl;
	//     cerr<<p<<"\tpos="<<piForGenomicWindow->at(p).posAlign<<"\t"<<piForGenomicWindow->at(p).readsVec.size()<<endl;

	//     for(unsigned int g=0;g<16;g++){
	// 	cerr<<babdIdx2Code16[g]<<"\t"<<vectorOfloglikelihoodForGivenBaBd[g]<<endl;
	// 	//cerr<<vectorToString(vectorOfloglikelihoodForGivenBaBd,"\t")<<endl;
	//     }

	//     cerr<<mostLikelyBaBd<<"\t"<<babdIdx2Code16[mostLikelyBaBdIdx]<<"\t"<<loglikelihoodForEveryBaBd_minusBest <<"\t"<< loglikelihoodForEveryBaBd<<"\t"<<( loglikelihoodForEveryBaBd_minusBest - loglikelihoodForEveryBaBd)<<endl;

	// } 




	//product for each genomic position
	
	loglikelihoodForEveryPositionForEveryBaBd += loglikelihoodForEveryBaBd; // \prod_{site} \sum_{genotype} (\prod_{fragment} P(D|G))*P(G)
	
    }//END for each genomic position
    cout<<h<<"\t"<<loglikelihoodForEveryPositionForEveryBaBd<<endl;
    exit(1);
	
}

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
			  vector<PositionResult *> * dataToWriteOut,
			  const int threadID,
			  vector<positionInformation> * piForGenomicWindow)
	: PileupVisitor()
	, m_references(references)
	, m_refID(refID)
	, m_leftCoord(leftCoord)
	, m_rightCoord(rightCoord)
	, m_dataToWriteOut( dataToWriteOut)
	, m_threadID( threadID )
	, m_numberOfSites( 0 )
	, totalBases(0)
	, totalSites(0)   
	, m_piForGenomicWindow( piForGenomicWindow )
    { 
	//cerr<<"heteroComputerVisitor constructor"<<endl;
    }
    ~heteroComputerVisitor(void) { }
  
    // PileupVisitor interface implementation

    
    void Visit(const PileupPosition& pileupData) {   
	


	if(pileupData.Position < int(m_leftCoord)   || 
	   pileupData.Position > int(m_rightCoord) ){
	    return ;
	}

	if( (m_numberOfSites%50000)==0 &&
	    m_numberOfSites != 0 ){
	    cerr<<"Thread#"<<m_threadID<<" reading: "<<m_references[m_refID].RefName<<":"<<pileupData.Position<<" valid sites:\t"<<thousandSeparator(totalSites)<<endl;

	}

	m_numberOfSites++;


	// int                 totalBases=0 ;
	//int                 counterB  [4];
	//long double         llBaseDeam[4];


	// vector<int>              obsBase      ;
	// vector<int>              obsQual      ;
	// vector<int>              mmQual       ; //mismapping probability
	// vector<bool>             isRevVec;
	//vector<singleRead> singleReadToAdd;
	positionInformation piToAdd;
	piToAdd.posAlign     = pileupData.Position;
	piToAdd.skipPosition = false;


	unsigned int                posAlign = pileupData.Position+1;

	bool foundSites=false;
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
	    if(!isResolvedDNA(b)){ //avoid Ns
		continue; 
	    }

	    int bIndex = baseResolved2int(b);
	    int   q    = MIN( int(pileupData.PileupAlignments[i].Alignment.Qualities[  pileupData.PileupAlignments[i].PositionInAlignment ]-offsetQual), MAXBASEQUAL);
	    int   m    = MIN( int(pileupData.PileupAlignments[i].Alignment.MapQuality), MAXMAPPINGQUAL );
	    bool isRev = pileupData.PileupAlignments[i].Alignment.IsReverseStrand();
	    

	    totalBases++;
	    foundSites=true;

	    singleRead sr_;
	    sr_.base    = uint8_t(bIndex);
	    sr_.qual    = uint8_t(q);
	    sr_.mapq    = uint8_t(m);
	    sr_.lengthF = uint8_t(pileupData.PileupAlignments[i].Alignment.Length);

	    if(isRev){
		sr_.pos5p = uint8_t(  pileupData.PileupAlignments[i].Alignment.Length-pileupData.PileupAlignments[i].PositionInAlignment-1 ); 
		sr_.base = 3 - sr_.base;//complement
	    }else{
		sr_.pos5p = uint8_t(  pileupData.PileupAlignments[i].PositionInAlignment ); 
	    }
	    sr_.isrv=isRev;
	    //sr_.name=pileupData.PileupAlignments[i].Alignment.Name;//to remove

	    piToAdd.readsVec.push_back(sr_);
	    // obsBase.push_back( bIndex  );
	    // obsQual.push_back( q        );
	    // mmQual.push_back(  m        );
	    // isRevVec.push_back(isRev);

	    //mmProb.push_back(  likeMismatchProbMap[m]  );
	    // substitutionRatesPerRead.push_back( probSubMatchToUseEndo );
	    //	    includeFragment[i]=true;

	}//END FOR EACH READ

	if( foundSites ){
	    totalSites++;
	}
	
	m_piForGenomicWindow->push_back(piToAdd);

    }//end Visit()
    
    unsigned int getTotalSites() const{
	return totalSites;
    }

    unsigned int getTotalBases() const{
	return totalBases;
    }

private:
    RefVector m_references;
    //Fasta * m_fastaReference;
    unsigned int totalBases;
    unsigned int totalSites;

    int          m_threadID;
    int          m_refID;
    unsigned int m_leftCoord;
    unsigned int m_rightCoord;

    vector<PositionResult *> * m_dataToWriteOut;
    unsigned int m_numberOfSites;
    vector<positionInformation> * m_piForGenomicWindow;

};//heteroComputerVisitor









//! Method called for each thread
/*!
  

*/				
void *mainHeteroComputationThread(void * argc){
    cerr<<"mainHeteroComputationThread started"<<endl;

    int   rc;

#ifdef HETVERBOSE    
    int rankThread=0;
#endif

    cerr<<"mainHeteroComputationThread mutex1"<<endl;

    rc = pthread_mutex_lock(&mutexRank);
    checkResults("pthread_mutex_lock()\n", rc);

    cerr<<"mainHeteroComputationThread mutex21 ID="<<(*(int *)pthread_self())<<endl;    
    cerr<<"mainHeteroComputationThread mutex22 ID="<<(threadID2Rank.size()+1)<<endl;
    //cerr<<"mainHeteroComputationThread mutex2"<<endl;    

    threadID2Rank[*(int *)pthread_self()]  = threadID2Rank.size()+1;

    cerr<<"mainHeteroComputationThread mutex3"<<endl;
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
 	cerr<<"Thread #"<<rankThread<<" is reading chunk rank#"<<currentChunk->rank<<endl;
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
    cerr<<"Thread #"<<rankThread<<" is reading region "<<currentChunk->rangeGen<<endl;
#endif

    //sleep(10);


    


#ifdef HETVERBOSE
    cerr<<"Thread #"<<rankThread<<" is reading BAM file: "<<bamFileToOpen<<endl;
#endif

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

#ifdef HETVERBOSE
    cerr<<"Thread #"<<rankThread<<" refID   "<<refID<<endl;    
    cerr<<"Thread #"<<rankThread<<" refName "<<references[refID].RefName<<endl;
#endif

    BamRegion bregion (refID, 
		       currentChunk->rangeGen.getStartCoord(), 
		       refID, 
		       currentChunk->rangeGen.getEndCoord()   );

    bool setRegionRes=reader.SetRegion( bregion   );

#ifdef HETVERBOSE
    cerr<<"Thread #"<<rankThread<<" setting region: "<<references[refID].RefName<<":"<<currentChunk->rangeGen.getStartCoord()<<"-"<<currentChunk->rangeGen.getEndCoord()<<"\tsucces? "<<setRegionRes<<endl;
#endif

    if( refID==-1 ||
       !setRegionRes){
    	cerr << "Heterozygous computation: could not set region "<<currentChunk->rangeGen<<" for BAM file:" << bamFileToOpen <<" "<<currentChunk->rangeGen.getStartCoord()<<" "<< currentChunk->rangeGen.getEndCoord()<< endl;
    	exit(1);
    }

    DataToWrite  * dataToWrite = new DataToWrite();

    dataToWrite->rangeGen      =  currentChunk->rangeGen;
    dataToWrite->rank          =  currentChunk->rank;

    vector<positionInformation> * piForGenomicWindow = new vector<positionInformation> ();

    //dataToWrite->dataToWriteOut=new vector<PositionResult *>();
    heteroComputerVisitor* cv = new heteroComputerVisitor(references,
							  refID,
							  currentChunk->rangeGen.getStartCoord(), 
							  currentChunk->rangeGen.getEndCoord()  ,
							  dataToWrite->vecPositionResults,
							  rankThread,
							  piForGenomicWindow);

    

    PileupEngine pileup;
    pileup.AddVisitor(cv);

    BamAlignment al;
    //unsigned int numReads=0;
    while ( reader.GetNextAlignment(al) ) {
        //cout<<"mainHeteroComputationThread al.Name="<<al.Name<<endl;
	

	if(al.Length>=MINLENGTHFRAGMENT &&
	   al.Length<=MAXLENGTHFRAGMENT ){	   
	    pileup.AddAlignment(al);
	}
    }

    //clean up
    pileup.Flush();
    reader.Close();
    //fastaReference.Close();
    
    cerr<<"Thread #"<<rankThread <<" bases="<<thousandSeparator(cv->getTotalBases())<<"\tsites="<<thousandSeparator(cv->getTotalSites())<<"\t"<<double(cv->getTotalBases())/double(cv->getTotalSites())<<endl;






    delete cv;


    computeLL(piForGenomicWindow);
    delete piForGenomicWindow;
	

    //call computeLL here
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



queue< DataChunk * >  subFirstElemsQueue(const queue< DataChunk * > queueDataToSubsample,unsigned int sizeToReturn){

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

    //random_shuffle ( myvectortemp.begin(), myvectortemp.end() );
    
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
							    currentChunk->rangeGen.getEndCoord()   );
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

    double lambdaCov=0;
    bool   lambdaCovSpecified=false;

    TStoTVratio=2.1;
    // string genoFileAsInput    ="";
    // bool   genoFileAsInputFlag=false;

    
    const string usage=string("\nThis program co-estimates heterozygosity rates and runs of homozygosity\n\n\t"+
                              string(argv[0])+                        
                              " [options] [fasta file] [bam file]  "+"\n\n"+
			      "\twhere:\n"+
			      "\t\t\t\t\t[fasta file]\t\tThe fasta file used for alignement\n"
			      "\t\t\t\t\t[bam file]\t\tThe aligned and indexed BAM file\n"+
			      "\n\n"
			      
                              "\n\tI/O options:\n"+
			      "\t\t"+"-o"+","+"--out"  + "\t\t"   +    "[outfile]" +"\t\t"+"Output per-site likelihoods in BGZIP (default: none)"+"\n"+
			      "\t\t"+""  +""+"--name" + "\t\t\t"   +    "[name]"    +"\t\t\t"+"Sample name (default: "+sampleName+")"+"\n"+
			      "\t\t"+""  +""+"--vcf"    + "\t\t\t" +    ""          +"\t\t\t"+"Use VCF as output format (default: "+booleanAsString(useVCFoutput)+")"+"\n"+
			      //"\t\t"+""+"\t"+"--ingeno"  + "\t\t"   +    "[infile]" +"\t\t"+"Read likelihoods in BGZIP and start comp. from there (default: none)"+"\n"+
			      "\n\tComputation options:\n"+
                              "\t\t"+"-t"+"\t"+""       +"\t\t"    +    "[threads]" +"\t\t"+"Number of threads to use (default: "+stringify(numberOfThreads)+")"+"\n"+
                              "\t\t"+""  +""+"--phred64"+"\t\t\t"  +    ""          +"\t\t"+"Use PHRED 64 as the offset for QC scores (default : PHRED33)"+"\n"+
			      "\t\t"+""  +""+"--size"       +"\t\t\t"    + "[window size]" +"\t\t"+"Size of windows in bp  (default: "+stringify(sizeChunk)+")"+"\n"+	      
			      "\t\t"+""  +""+"--lambda"     +"\t\t"    + "[lambda]" +"\t\t"+"Skip coverage computation, specify lambda manually  (default: "+booleanAsString(lambdaCovSpecified)+")"+"\n"+	      
			      "\t\t"+""  +""+"--tstv"     +"\t\t\t"    + "[tstv]" +"\t\t\t"+"Ratio of transitions to transversions  (default: "+stringify(TStoTVratio)+")"+"\n"+	      



                              // "\n\tSample options:\n"+
                              // "\t\t"+""  +""+"--cont"  +"\t\t\t"    +  "[cont rate:0-1]" +"\t\t"+"Present-day human contamination rate (default: "+stringify(contrate)+")"+"\n"+
                              // // "\t\t"+"--phred64" +"\t\t\t\t"+"Use PHRED 64 as the offset for QC scores (default : PHRED33)"+"\n"+
			      
			      
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


	// if( string(argv[i]) == "--ingeno" ){
	//     genoFileAsInput     = string(argv[i+1]);
	//     genoFileAsInputFlag = true;
        //     i++;
        //     continue;
	// }
	
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

        if(string(argv[i]) == "--lambda"  ){
	    lambdaCovSpecified=true;
            lambdaCov         =destringify<double>(argv[i+1]);
	    i++;
            continue;
        }

        if(string(argv[i]) == "--tstv"  ){
            TStoTVratio         =destringify<long double>(argv[i+1]);
	    i++;
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

    cerr<<"Parsing arguments ...";

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
    for(int b1=0;b1<4;b1++){
	for(int b2=0;b2<4;b2++){
	    int dinucIndex = b1 *4+ b2;
	    illuminaErrorsProbMatrix.p[b1][b2] = illuminaErrorsProb.s[dinucIndex];
	}
    }

#ifdef DEBUGILLUMINAFREQ
    for(unsigned int i=0;i<16;i++){
    	cout<<i<<"\t"<<illuminaErrorsProb.s[i]<<endl;	
    }

    for(int b1=0;b1<4;b1++)
	cerr<<"ACGT"[b1]<<"\t";
    cerr<<endl;
    for(int b1=0;b1<4;b1++){
	cerr<<"ACGT"[b1]<<"\t";
	for(int b2=0;b2<4;b2++){
	    cerr<<illuminaErrorsProbMatrix.p[b1][b2]<<"\t";
	}
	cerr<<endl;
    }
    cerr<<endl;
    


    return 1;
#endif

    cerr<<"...done"<<endl;
    ///////////////////////////////////////
    // END read Illumina error profile //
    ///////////////////////////////////////    


    cerr<<"Begin pre-computation ...";
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

    cerr<<"..done"<<endl;
    //    return 1;

    cerr<<"Computing average coverage.."<<endl;

    int                   rc=0;
    pthread_t             threadCov[numberOfThreads];

    int    bpToExtract       = sizeChunk;
    



    GenomicWindows     rw  (fastaIndex,false);


    vector<GenomicRange> v = rw.getGenomicWindows(bpToExtract,0);
    if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
    
    unsigned int      rank  = 0;
    int          lastRank   =-1;
    unsigned int sizeGenome = 0;

    for(unsigned int i=0;i<v.size();i++){
	//cout<<"genomic region #"<<i<<" "<<v[i]<<endl;
	DataChunk * currentChunk = new DataChunk();
	
	currentChunk->rangeGen   = v[i];
	currentChunk->rank       = rank;
	lastRank                 = rank;
	//sizeGenome             += v[i].getLength();
 
	queueDataToprocess.push(currentChunk);
	rank++;
    }
    readDataDone=true;
    //    return 1;

    vector<GenoResults *> vectorGenoResults;


    if( !lambdaCovSpecified ){

	// if(!genoFileAsInputFlag){

	/////////////////////////////
	// BEGIN  Compute coverage //
	/////////////////////////////
	int bpToComputeCoverage = 10000000;
	int genomicRegionsToUse = bpToComputeCoverage/bpToExtract;
	if( genomicRegionsToUse > int(queueDataToprocess.size())){
	    genomicRegionsToUse = int(queueDataToprocess.size());
	}
    

	//TODO to renable
	//queueDataForCoverage = randomSubQueue( queueDataToprocess,genomicRegionsToUse);
	queueDataForCoverage = subFirstElemsQueue( queueDataToprocess,genomicRegionsToUse);


	pthread_mutex_init(&mutexQueue,   NULL);
	pthread_mutex_init(&mutexCounter, NULL);
	pthread_mutex_init(&mutexRank ,   NULL);

	for(int i=0;i<numberOfThreads;i++){
	    rc = pthread_create(&threadCov[i], NULL, mainCoverageComputationThread, NULL);
	    checkResults("pthread_create()\n", rc);
	}

	cerr<<"Creating threads for coverage calculation, need to process="<<queueDataForCoverage.size()<<" out of a total of "<<queueDataToprocess.size()<<endl;


	while(queueDataForCoverage.size()!=0){
	    cerr<<getDateString()<<" "<<getTimeString()<<" # of slices left to process: "<<queueDataForCoverage.size()<<"/"<<queueDataToprocess.size()<<endl;
	    sleep(timeThreadSleep);
	}
    
	//waiting for threads to finish
	for (int i=0; i <numberOfThreads; ++i) {	
	    rc = pthread_join(threadCov[i], NULL);
	    checkResults("pthread_join()\n", rc);
	}
	cerr<<"Coverage computations are done"<<endl;
	pthread_mutex_destroy(&mutexRank);
	pthread_mutex_destroy(&mutexQueue);
	pthread_mutex_destroy(&mutexCounter);

	//    cout<<"Final" <<" "<<totalBasesSum<<"\t"<<totalSitesSum<<"\t"<<double(totalBasesSum)/double(totalSitesSum)<<endl;
	pthread_exit(NULL);
    
	rateForPoissonCov    = ((long double)totalBasesSum)/((long double)totalSitesSum);
    }else{ //not lambdaCovSpecified
	//use the rate specified via the command line
	rateForPoissonCov = lambdaCov;
    }

    long double rateForPoissonCovFloor = floorl(rateForPoissonCov);
    long double rateForPoissonCovCeil  =  ceill(rateForPoissonCov);
    
#ifdef DEBUGCOV
    cout<<"rateForPoissonCovFloor "<<rateForPoissonCovFloor<<endl;
    cout<<"rateForPoissonCovCeil "<<rateForPoissonCovCeil<<endl;
#endif

    cov2probPoisson = new vector<long double> (MAXCOV+1,0);

    for(int i=0;i<=MAXCOV;i++){
	long double florPPMF=poisson_pmfl( (long double)(i), rateForPoissonCovFloor)/poisson_pmfl( rateForPoissonCovFloor, rateForPoissonCovFloor);
	long double ceilPPMF=poisson_pmfl( (long double)(i), rateForPoissonCovCeil) /poisson_pmfl( rateForPoissonCovCeil,  rateForPoissonCovCeil);
		
	cov2probPoisson->at(i) = (1.0-(rateForPoissonCov-rateForPoissonCovFloor))*florPPMF  + (1.0-(rateForPoissonCovCeil-rateForPoissonCov))*ceilPPMF;

#ifdef DEBUGCOV
	cerr<<"florPPMF "<<i<<" = "<<poisson_pmfl( (long double)(i), rateForPoissonCovFloor)<<"/"<<poisson_pmfl( rateForPoissonCovFloor, rateForPoissonCovFloor)<<" "<<florPPMF<<" w= "<<(1.0-(rateForPoissonCov-rateForPoissonCovFloor))<<endl;
	cerr<<"ceilPPMF "<<i<<" = "<<poisson_pmfl( (long double)(i), rateForPoissonCovCeil) <<"/"<<poisson_pmfl( rateForPoissonCovCeil,  rateForPoissonCovCeil)<<" "<<ceilPPMF<<" w= "<<(1.0-(rateForPoissonCovCeil-rateForPoissonCov))<<endl;
	cerr<<"cov2probPoisson["<<i<<"] = "<<cov2probPoisson->at(i)<<endl;
#endif

    }


    
    cerr<<"Results\tbp="<<totalBasesSum<<"\tsites="<<totalSitesSum<<"\tlambda="<<rateForPoissonCov<<endl;
    // for(int i=0;i<20;i++){
    // 	cout<<i<<"\t"<<pdfPoisson( (long double)i, rateForPoissonCov)/pdfPoisson( rateForPoissonCov, rateForPoissonCov)<<endl;
    // }


    // for(int i=0;i<100;i++){
    // 	cout<<i<<"\t"<<pdfPoisson( (long double)i, 20)/pdfPoisson( 20, 20)<<endl;
    // }

    cerr<<"..done"<<endl;

    // return 1;

    // pthread_exit(NULL);

    ////////////////////////////
    // END   Compute coverage //
    ////////////////////////////


















    // doneReading=true;    

    ///////////////////////
    //  Compute hetero   //
    ///////////////////////
    pthread_t             threadHet[numberOfThreads];
    threadID2Rank=map<unsigned int, int> ();


    cerr<<"Creating threads for heterozygosity calculation"<<endl;
    pthread_mutex_init(&mutexQueue,   NULL);
    pthread_mutex_init(&mutexCounter, NULL);
    pthread_mutex_init(&mutexRank ,   NULL);
    

    for(int i=0;i<numberOfThreads;i++){
	rc = pthread_create(&threadHet[i], NULL, mainHeteroComputationThread, NULL);
	checkResults("pthread_create()\n", rc);
    }
    //    return 1;
    // 	//threads are running here

    // unsigned int originalSize = queueDataToprocess.size();
    // while(queueDataToprocess.size()!=0){
    // 	cout<<"# of slices left to process: "<<queueDataToprocess.size()<<"/"<<originalSize<<endl;
    // 	sleep(timeThreadSleep);
    // }
    cerr<<"done creating threads "<<endl;
    cerr<<"outFileSiteLLFlag ="<<outFileSiteLLFlag<<endl;
    ///////////////////
    //Writing data out/
    ///////////////////


    Internal::BgzfStream  bgzipWriter;

    // if(outFileSiteLLFlag){
    // 	bgzipWriter.Open(outFileSiteLL, IBamIODevice::WriteOnly);
    // 	if(!bgzipWriter.IsOpen()){
    // 	    cerr<<"Cannot open file "<<outFileSiteLL<<" in bgzip writer"<<endl;
    // 	    return 1;
    // 	}
    

    // 	string headerOutFile;
    // 	if(useVCFoutput){
    // 	    headerOutFile="#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+sampleName+"\n";	
    // 	}else{
    // 	    headerOutFile="#CHROM\tPOS\tA\tC\tG\tT\tGENO\tGENOS\tQualL\tCovL\tAA\tAC\tAG\tAT\tCC\tCG\tCT\tGG\tGT\tTT\n";	
    // 	}

    // 	bgzipWriter.Write(headerOutFile.c_str(),headerOutFile.size());
    // }

    bool wroteEverything=false;
    int lastWrittenChunk=-1;   
    cerr<<"test wroteEverything="<<wroteEverything<<endl;
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
	rc = pthread_join(threadHet[i], NULL);
	checkResults("pthread_join()\n", rc);
    }

    cerr<<"Heterozygosity computations are done"<<endl;
    pthread_mutex_destroy(&mutexRank);
    pthread_mutex_destroy(&mutexQueue);
    pthread_mutex_destroy(&mutexCounter);


    if(outFileSiteLLFlag){
	bgzipWriter.Close();
    }

    // //}else{//if we suply the geno file as input

    // string    lineGENOL;
    // igzstream myFileGENOL;
    // 	cerr<<"Reading genotypes: "<<genoFileAsInput<<" ...";
    // 	myFileGENOL.open(genoFileAsInput.c_str(), ios::in);

    // 	if (myFileGENOL.good()){
    // 	    getline (myFileGENOL,lineGENOL);//header
    // 	    while ( getline (myFileGENOL,lineGENOL)){
    // 		//cout<<lineGENOL<<endl;
    // 		GenoResults * toadd =  new GenoResults( lineGENOL );
    // 		vectorGenoResults.push_back(toadd);
    // 		sizeGenome++;
    // 	    }
    // 	    myFileGENOL.close();
    // 	}else{
    // 	    cerr << "Error: unable to open file with genotype likelihoods: "<<genoFileAsInput<<endl;
    // 	    return 1;
    // 	}
    // 	cerr<<"..done"<<endl;

    // }

    // 


    //    return 1;
    //////////////////////////////////
    //                              //
    // BEGIN COMPUTE HETERO RATE    //
    //                              //
    //////////////////////////////////


#ifdef NODEF

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
        



    for(unsigned int i=0;i<vectorGenoResults.size();i++){
	delete( vectorGenoResults[i] );
    }

    pthread_exit(NULL);

#endif    


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

    delete(cov2probPoisson);
    
    return 0;
}

