#include <iostream>
#include <fstream>
#include <queue>
#include <algorithm>   
//#include <random>



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

#include "utils.h"

using namespace std;
using namespace BamTools;

//#define COVERAGETVERBOSE

#define MAXMAPPINGQUAL 257     // maximal mapping quality, should be sufficient as mapping qualities are encoded using 8 bits
#define MAXCOV          50     // maximal coverage

char offset=33;
long double likeMatch        [MAXMAPPINGQUAL];
long double likeMismatch     [MAXMAPPINGQUAL];

long double likeMatchProb    [MAXMAPPINGQUAL];
long double likeMismatchProb [MAXMAPPINGQUAL];

vector< vector<long double> > binomVec (MAXCOV,vector<long double>(MAXCOV,0)) ;
unsigned int totalBasesSum;
unsigned int totalSitesSum;


// // Returns logl( expl(x)+expl(y) )
// inline long double oplusl(long double x, long double y ){
//     return x > y 
//         ? x + log1pl( expl( y-x ) )
//         : y + log1pl( expl( x-y ) )  ;
// }




long double pdfPoisson(const long double l,const long double k ) {
    return expl(k*logl(l)-lgammal(k+1.0)-l);
}


//! A method to initialize various probability scores to avoid recomputation
/*!
  This method is called by the main after capturing the arguments
*/
void initScores(){
    totalBasesSum=0;
    totalSitesSum=0;

    //Computing for quality scores 2 and up
    for(int i=0;i<MAXMAPPINGQUAL;i++){
        likeMatch[i]        = log1pl(    -pow(10.0,i/-10.0) );          
        likeMismatch[i]     = logl  (     pow(10.0,i/-10.0) );

        likeMatchProb[i]           = 1.0-pow(10.0,i/-10.0);
        likeMismatchProb[i]        =     pow(10.0,i/-10.0);
    }



    for(int i=1;i<MAXCOV;i++){
	//cout<<i<<endl;

	for(int j=0;j<=i;j++){	    
	    binomVec[i][j] = ( logl(nChoosek(i,j))+logl(powl(0.5,i)) );	    
	    //cout<<j<<"\t"<<(logl(nChoosek(i,j))+logl(powl(0.5,i)))<<endl;//todo: precompute that line
	}
    }

}//end initScores


long double computeLL(const char                 al1Current,
		      const char                 al2Current,		      
		      const vector<int>          obsBase   ,
		      const vector<long double>  probDeam  ,
		      const vector<int>          obsQual   ,
		      const long double          contRate  ,
		      const char                 alContCurrent ,
		      const vector<long double>  mismappingProb
		      ){

    long double llik=0;
    long double llik1=0;
    long double llik2=0;
    long double llikC=0;
    

    for(int i=0;i<int(obsBase.size());i++){

	//contaminant
	if(obsBase[i] == alContCurrent){
	    llikC    =      (1.0-mismappingProb[i])*(likeMatchProb[obsQual[i]]*(1.0) + likeMismatchProb[obsQual[i]]*(0.5))+mismappingProb[i]*0.5;
	}else{
	    llikC    =      (1.0-mismappingProb[i])*(likeMatchProb[obsQual[i]]*(0.0) + likeMismatchProb[obsQual[i]]*(0.5))+mismappingProb[i]*0.5;
	}

	long double llikAl1t=0;
	long double llikAl2t=0;

	//replace prob deam with proper value from matrix
	if(obsBase[i] == al1Current){
	    llikAl1t =     (1.0-mismappingProb[i])*(likeMatchProb[obsQual[i]]*(1.0-probDeam[i])  + likeMismatchProb[obsQual[i]]*(0.5))+mismappingProb[i]*0.5;;
	}else{
	    llikAl1t =     (1.0-mismappingProb[i])*(likeMatchProb[obsQual[i]]*(    probDeam[i])  + likeMismatchProb[obsQual[i]]*(0.5))+mismappingProb[i]*0.5;
	}


	if(obsBase[i] == al2Current){
	    llikAl2t =     (1.0-mismappingProb[i])*(likeMatchProb[obsQual[i]]*(1.0-probDeam[i])  + likeMismatchProb[obsQual[i]]*(0.5))+mismappingProb[i]*0.5;;
	}else{
	    llikAl2t =     (1.0-mismappingProb[i])*(likeMatchProb[obsQual[i]]*(    probDeam[i])  + likeMismatchProb[obsQual[i]]*(0.5))+mismappingProb[i]*0.5;
	}

	// cout<<i<<"\t"<<likeMatchProb[obsQual[i]]<<"\t"<<likeMismatchProb[obsQual[i]]<<endl;
	// cout<<i<<"\t"<<llikAl1t <<"\t"<<llikAl2t<<endl;
	// exit(1);
	long double llikT  = (1.0-contRate)*( 0.5*llikAl1t + 0.5*llikAl2t ) + (contRate)*llikC  ;	
	llik              += logl(llikT);

 	llik1=oplusInitnatl( llik1, logl(llikAl1t) );
	llik2=oplusInitnatl( llik2, logl(llikAl2t) );

	// cout<<i<<"\t"<<llikAl1t <<"\t"<<llikAl2t<<"\t"<<llikC<<"\t"<<llikT<<endl;
	// llik1+=logl( llikAl1t );
	// llik2+=logl( llikAl2t );
    }
    //cout<<"CC\t"<<llikCC<<"\t"<<llikCC1<<"\t"<<llikCC2<<endl;    
    //cout<<llik1<<"\t"<<llik2<<"\t"<<expl(llik1)<<"\t"<<expl(llik2)<<endl;
    long double expal1  = roundl(  int(obsBase.size()) * ( expl(llik1) / expl(oplusnatl(llik1,llik2))) );
    //long double expal1=roundl(  sizeAr * ( expl(llik1) / llik1+llik2)) );
    //long double expal2  = sizeAr-expal1;
    long double binomE2 = binomVec[int(obsBase.size())][expal1]; //logl(nChoosek(sizeAr,expal1)*powl(0.5,expal1+expal2));
    //cout<<al1Current<<""<<al2Current<<"\t"<<llik<<"\t"<<expal1<<"\t"<<expal2<<"\t"<<binomE2<<"\t"<<(binomE2+llik)<<endl;

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
	bool foundOneFragment=false;
	if(pileupData.Position < int(m_leftCoord)   || 
	   pileupData.Position > int(m_rightCoord) ){
	    return ;
	}
	//cout<<m_leftCoord<<"\t"<<m_rightCoord<<"\t"<<pileupData.Position<<endl;

	for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){
	    // if( pileupData.PileupAlignments[i].IsCurrentDeletion &&
	    // 	pileupData.PileupAlignments[i].IsNextInsertion ){
	    // 	continue;
	    // }
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
    heteroComputerVisitor(const RefVector& references,unsigned int leftCoord, unsigned int rightCoord)
	: PileupVisitor()
	, m_references(references)
	, m_leftCoord(leftCoord)
	, m_rightCoord(rightCoord)
    { 
    }
    ~heteroComputerVisitor(void) { }
  
    // PileupVisitor interface implementation

    
    void Visit(const PileupPosition& pileupData) {   
	
	int                 counterB  [4];
	vector<int>         obsBase   ;
	vector<int>         obsQual   ;
	vector<long double> probDeam  ;
	vector<long double> mmProb    ;
	unsigned int posAlign = pileupData.Position+1;
    
	for(unsigned int i=0;i<pileupData.PileupAlignments.size();i++){
	    if(i>=MAXCOV){
		break;
	    }

	    char  b   =     pileupData.PileupAlignments[i].Alignment.QueryBases[ pileupData.PileupAlignments[i].PositionInAlignment ];
	    if(!isResolvedDNA(b)){ continue; }//avoid Ns
	    int bIndex = baseResolved2int(b);
	    int   q   = int(pileupData.PileupAlignments[i].Alignment.Qualities[  pileupData.PileupAlignments[i].PositionInAlignment ]-offset); 
	    int   m   = int(pileupData.PileupAlignments[i].Alignment.MapQuality);
	    
	    counterB[ bIndex ]++;
	    
	    obsBase.push_back( bIndex );
	    obsQual.push_back( q      );
	    mmProb.push_back(  m      );
	    
	    //Fill deamination vector, do not forget fragment orientation
	    //put proper probabilities
	    probDeam.push_back( 0.0   );
	}
	
	int counterUnique=0;
	for(int i=0;i<4;i++){
	    if(counterB[i]!=0) counterUnique++;
	}
	
	//skip sites with no defined bases and tri/tetra allelic sites
	if(counterUnique==0 ||
	   counterUnique>=3){
	    return;
	}

	int ref=-1;
	int alt=-1;

	if(counterUnique==1){
	    for(int i=0;i<4;i++){
		if(counterB[i]!=0) ref=i;
	    }
	    alt=randomBPExceptInt(ref);	    //maybe put a better dna sub model here?
	}


	if(counterUnique==2){
	    for(int i=0;i<4;i++){
		if(counterB[i]!=0){
		    if(ref == -1) 
			ref=i;
		    else         
			alt=i;
		}
	    }
	}
	


	
	cout<<posAlign<<"\t"<<"ACGT"[ref]<<","<<"ACGT"[alt]<<"\t"
	    <<"\t"<<computeLL(ref,
			      ref,		      		  
			      obsBase   ,
			      probDeam  ,
			      obsQual   ,
			      0.0       ,//todo set contamination rate
			      ref       ,
			      mmProb)
	    <<"\t"<<computeLL(ref,
			      alt,		      		  
			      obsBase   ,
			      probDeam  ,
			      obsQual   ,
			      0.0       ,//todo set contamination rate
			      ref       ,
			      mmProb)
	    <<"\t"<<computeLL(alt,
			      alt,		      		  
			      obsBase   ,
			      probDeam  ,
			      obsQual   ,
			      0.0       ,//todo set contamination rate
			      ref       ,
			      mmProb)<<endl;
	
    }
    

private:
    RefVector m_references;
    //Fasta * m_fastaReference;
    // unsigned int totalBases;
    // unsigned int totalSites;
    unsigned int m_leftCoord;
    unsigned int m_rightCoord;

};//heteroComputerVisitor




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

DataChunk::DataChunk(){
    //cerr<<"Constructor addr: "<<this<<endl;
}

DataChunk::~DataChunk(){
    //cerr<<"Destructor  addr: "<<this<<endl;
}

class CompareDataChunk {
public:
    bool operator() ( DataChunk * cd1, DataChunk * cd2)  {
        //comparison code here
	return ( cd1->rank > cd2->rank );
    }
};


int    timeThreadSleep =    10;
bool      readDataDone = false;
unsigned int sizeChunk =  5000;

string                                                             bamFileToOpen;
queue< DataChunk * >                                               queueDataToprocess;
queue< DataChunk * >                                               queueDataForCoverage;

priority_queue<DataChunk *, vector<DataChunk *>, CompareDataChunk> queueDataTowrite;

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
 



//! Method called for each thread
/*!
  

*/				
void *mainHeteroComputationThread(void * argc){

    int   rc;
    int rankThread=0;

    rc = pthread_mutex_lock(&mutexRank);
    checkResults("pthread_mutex_lock()\n", rc);

    threadID2Rank[*(int *)pthread_self()]  = threadID2Rank.size()+1;
    rankThread = threadID2Rank[*(int *)pthread_self()];

    
    rc = pthread_mutex_unlock(&mutexRank);
    checkResults("pthread_mutex_unlock()\n", rc);

 checkqueue:    
    // stackIndex=-1;
    //check stack

    
    rc = pthread_mutex_lock(&mutexQueue);
    checkResults("pthread_mutex_lock()\n", rc);


    bool foundData=false;
    
    cerr<<"Thread #"<<rankThread <<" started and is requesting data"<<endl;


    DataChunk * currentChunk;


    if(!queueDataToprocess.empty()){    
 	foundData=true;
 	currentChunk = queueDataToprocess.front();
 	queueDataToprocess.pop();
 	cerr<<"Thread #"<<rankThread<<" is reading "<<currentChunk->rank<<endl;
	//cout<<"rank "<< &(currentChunk->dataToProcess) <<endl;
    }

    
  

    if(!foundData){
	rc = pthread_mutex_unlock(&mutexQueue);
	checkResults("pthread_mutex_unlock()\n", rc);


	if(readDataDone){
	    cerr<<"Thread #"<<rankThread<<" is done"<<endl;
	    return NULL;	
	}else{
	    cerr<<"Thread #"<<rankThread<<" sleeping for "<<timeThreadSleep<<endl;
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
    cerr<<"Thread #"<<rankThread<<" is reading "<<currentChunk->rangeGen<<endl;
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
    	cerr << "Could not set region "<<currentChunk->rangeGen<<" for BAM file:" << bamFileToOpen <<" "<<currentChunk->rangeGen.getStartCoord()<<" "<< currentChunk->rangeGen.getEndCoord()<< endl;
    	exit(1);
    }

   
    heteroComputerVisitor* cv = new heteroComputerVisitor(references,
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
    
    //cerr<<"Thread #"<<rankThread <<" "<<cv->getTotalBases()<<"\t"<<cv->getTotalSites()<<"\t"<<double(cv->getTotalBases())/double(cv->getTotalSites())<<endl;

    delete cv;

	

    cerr<<"Thread #"<<rankThread <<" is done with computations"<<endl;

    //////////////////////////////////////////////////////////////
    //                END   COMPUTATION                         //
    //////////////////////////////////////////////////////////////
	

    //COUNTERS
    rc = pthread_mutex_lock(&mutexCounter);
    checkResults("pthread_mutex_lock()\n", rc);
    
    //queueDataTowrite.push(currentChunk);

    rc = pthread_mutex_unlock(&mutexCounter);
    checkResults("pthread_mutex_unlock()\n", rc);

    cerr<<"Thread #"<<rankThread <<" is re-starting"<<endl;

    goto checkqueue;	   


    

    
    cerr<<"Thread "<<rankThread<<" ended "<<endl;
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



void *mainCoverageComputationThread(void * argc){
    initScores();
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
		       currentChunk->rangeGen.getStartCoord(), 
		       refID, 
		       currentChunk->rangeGen.getEndCoord()   );

    bool setRegionRes=reader.SetRegion( bregion   );

    // cerr<<"Thread #"<<rankThread<<" "<<references[0].RefName<<"\t"<<setRegionRes<<endl;

    if( refID!=-1 &&
       !setRegionRes){
    	cerr << "Could not set region "<<currentChunk->rangeGen<<" for BAM file:" << bamFileToOpen <<" "<<currentChunk->rangeGen.getStartCoord()<<" "<< currentChunk->rangeGen.getEndCoord()<< endl;
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





int main (int argc, char *argv[]) {

    if(argc!=4){
	cerr<<"usage:"<<endl<<"\t"<<argv[0]<<" [bamfile in] [faidx] [# of threads]"<<endl<<endl;
	return 1;
    }
    initScores();
    bamFileToOpen            = string(argv[1]);
    string fastaIndex        = string(argv[2]);
    int    numberOfThreads   = destringify<int>( string(argv[3]) );
    int    bpToExtract       = 500;
    
    pthread_t             thread[numberOfThreads];
    int                   rc=0;



    GenomicWindows     rw  (fastaIndex,false);
    //TODO add genomic ranges in queue
    vector<GenomicRange> v = rw.getGenomicWindows(bpToExtract,0);
    if( v.size() == 0 ){    cerr<<"No range found using these chr/loci settings"<<endl; return 1;}
    unsigned int rank=0;

    for(unsigned int i=0;i<v.size();i++){
	//cout<<v[i]<<endl;
	DataChunk * currentChunk = new DataChunk();

	currentChunk->rangeGen = v[i];
	currentChunk->rank     = rank;

	queueDataToprocess.push(currentChunk);
	rank++;
    }
    readDataDone=true;
    //return 1;

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
    cout<<"creating threads"<<endl;

    //waiting for threads to finish
    for (int i=0; i <numberOfThreads; ++i) {
	rc = pthread_join(thread[i], NULL);
	checkResults("pthread_join()\n", rc);
    }
    cout<<"done threads"<<endl;
    pthread_mutex_destroy(&mutexRank);
    pthread_mutex_destroy(&mutexQueue);
    pthread_mutex_destroy(&mutexCounter);

    //    cout<<"Final" <<" "<<totalBasesSum<<"\t"<<totalSitesSum<<"\t"<<double(totalBasesSum)/double(totalSitesSum)<<endl;
    //    pthread_exit(NULL);
    
    long double rateForPoissonCov = ((long double)totalBasesSum)/((long double)totalSitesSum);

    cout<<"Final" <<" "<<totalBasesSum<<"\t"<<totalSitesSum<<"\t"<<double(totalBasesSum)/double(totalSitesSum)<<endl;
    // for(int i=0;i<20;i++){
    // 	cout<<i<<"\t"<<pdfPoisson( (long double)i, rateForPoissonCov)/pdfPoisson( rateForPoissonCov, rateForPoissonCov)<<endl;
    // }


    for(int i=0;i<100;i++){
	cout<<i<<"\t"<<pdfPoisson( (long double)i, 20)/pdfPoisson( 20, 20)<<endl;
    }

    return 1;


    ////////////////////////////
    // END   Compute coverage //
    ////////////////////////////
    // doneReading=true;    

    ///////////////////////
    //  Compute hetero   //
    ///////////////////////

    pthread_mutex_init(&mutexQueue,   NULL);
    pthread_mutex_init(&mutexCounter, NULL);
    pthread_mutex_init(&mutexRank ,   NULL);

    for(int i=0;i<numberOfThreads;i++){
	rc = pthread_create(&thread[i], NULL, mainHeteroComputationThread, NULL);
	checkResults("pthread_create()\n", rc);
    }


    //waiting for threads to finish
    for (int i=0; i <numberOfThreads; ++i) {
	rc = pthread_join(thread[i], NULL);
	checkResults("pthread_join()\n", rc);
    }

    pthread_mutex_destroy(&mutexRank);
    pthread_mutex_destroy(&mutexQueue);
    pthread_mutex_destroy(&mutexCounter);

    //writer.Close();
    //reader.Close();    
    pthread_exit(NULL);

    return 0;
}

