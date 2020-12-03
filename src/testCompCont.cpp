/*
 * testComp
 * Date: Mar-06-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "libgab.h"



using namespace std;

#define MAXMAPPINGQUAL 257     // maximal mapping quality, should be sufficient as mapping qualities are encoded using 8 bits
#define MAXCOV          50     // maximal coverage

long double likeMatch        [MAXMAPPINGQUAL];
long double likeMismatch     [MAXMAPPINGQUAL];

long double likeMatchProb    [MAXMAPPINGQUAL];
long double likeMismatchProb [MAXMAPPINGQUAL];

vector< vector<long double> > binomVec;
// // Returns logl( expl(x)+expl(y) )
// inline long double oplusl(long double x, long double y ){
//     return x > y 
//         ? x + log1pl( expl( y-x ) )
//         : y + log1pl( expl( x-y ) )  ;
// }


//! A method to initialize various probability scores to avoid recomputation
/*!
  This method is called by the main after capturing the arguments
*/
void initScores(){

    // for(int i=0;i<2;i++){
    //     likeMatch[i]        = log1pl(    -pow(10.0,2.0/-10.0) )   ;         
    //     likeMismatch[i]     = logl  (     pow(10.0,2.0/-10.0)/1.0 );

    // 	likeMatchProb[i]           = 1.0-pow(10.0,2.0/-10.0) ;
    //     likeMismatchProb[i]        =     pow(10.0,2.0/-10.0)/1.0 ;
    // }


    //Computing for quality scores 2 and up
    for(int i=0;i<MAXMAPPINGQUAL;i++){
        likeMatch[i]        = log1pl(    -pow(10.0,i/-10.0) );          
        likeMismatch[i]     = logl  (     pow(10.0,i/-10.0) );

        likeMatchProb[i]           = 1.0-pow(10.0,i/-10.0);
        likeMismatchProb[i]        =     pow(10.0,i/-10.0);
    }


    binomVec.resize(MAXCOV);
    for(int i=1;i<MAXCOV;i++){
	binomVec[i].resize(MAXCOV);
    }

    for(int i=1;i<MAXCOV;i++){
	//cout<<i<<endl;

	for(int j=0;j<=i;j++){
	    //todo avoid precision loss
	    binomVec[i][j] = ( logl(nChoosek(i,j))+logl(powl(0.5,i)) );	    
	    //cout<<j<<"\t"<<(logl(nChoosek(i,j))+logl(powl(0.5,i)))<<endl;//todo: precompute that line
	}
    }


}

long double computeLL(const char         al1Current,
		      const char         al2Current,		      
		      const int          sizeAr,
		      const char         obsBase [],
		      const long double  probDeam [],
		      const int          obsQual [],
		      const long double  contRate,
		      const char         alContCurrent ,
		      const long double  mismappingProb []
		      ){

    long double llik=0;
    long double llik1=0;
    long double llik2=0;
    long double llikC=0;
    
    // for(int i=0;i<sizeAr;i++){
    // 	cout<<i<<"\t"<<obsBase[i]<<"\t"<<obsQual[i]<<"\t"<<probDeam[i]<<endl;
    // }

    for(int i=0;i<sizeAr;i++){

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
    long double expal1  = roundl(  sizeAr * ( expl(llik1) / expl(oplusnatl(llik1,llik2))) );
    //long double expal1=roundl(  sizeAr * ( expl(llik1) / llik1+llik2)) );
    //long double expal2  = sizeAr-expal1;
    long double binomE2 = binomVec[sizeAr][expal1]; //logl(nChoosek(sizeAr,expal1)*powl(0.5,expal1+expal2));
    //cout<<al1Current<<""<<al2Current<<"\t"<<llik<<"\t"<<expal1<<"\t"<<expal2<<"\t"<<binomE2<<"\t"<<(binomE2+llik)<<endl;

    return (binomE2+llik);
}

int main (int argc, char *argv[]) {
    initScores();


    // cout<<expl(oplusl( logl(0.4),logl(0.5)))<<endl;
    // for(int i=0;i<40;i++){
    // 	cout<<i<<"\t"<<likeMatchProb[i]<<"\t"<<likeMismatchProb[i]<<endl;
    // }
    // return 1;
    //cout<<nChoosek(20,10)<<endl;

    char alleles [] = {'C','T'};

    char al1Current   = ' ';
    char al2Current   = ' ';



    // for(int i=0;i<=sizeAr;i++){
    // 	cout<<sizeAr<<"\t"<<i<<"\t"<<nChoosek(sizeAr,i)<<"\t"<<(nChoosek(sizeAr,i)*powl(0.5,sizeAr-i)*powl(0.5,i))<<endl;
    // }
    // int sizeAr =  10;
    // char obsBase          [ ] =  { 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C' }; 
    // int  obsQual          [ ] =  {  40,  40,  40,  40,  40,  40,  40,  40,  40,  40 }; 
    // //int  obsQual          [ ] =  {  0,   0,   0,   0,   0,   0,   0,   0,   0,   0  }; 
    // long double probDeam  [ ] =  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; 
    // long double mmProb    [ ] =  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; //mismapping prob

    int sizeAr =  6;
    //char obsBase [ ] =  { 'C', 'C', 'C', 'C', 'C', 'C' }; 
    char obsBase          [ ] =  { 'C', 'C', 'C', 'C', 'T', 'T' }; 
    int  obsQual          [ ] =  {  40,  40,  40,  40,  40,  40 }; 
    long double probDeam  [ ] =  { 0.0, 0.0, 0.0, 0.0, 0.2, 0.2 }; 
    long double mmProb    [ ] =  { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; //mismapping prob

    // char obsBase          [ ] =  { 'C', 'C', 'C', 'C', 'T' }; 
    // int  obsQual          [ ] =  {  40,  40,  40,  40,  40 }; 
    // //int  obsQual          [ ] =  {   0,   0,   0,   0,   0 }; 
    // long double probDeam  [ ] =  { 0.0, 0.0, 0.0, 0.0, 0.0 }; 
    //int  obsQual          [ ] =  {  2,   2,   2,   2,   2 }; 
    long double contRate=1.0;
    //int  obsQual [ ] =  {  2,   2,   2,   2,   2,   2 }; 

    // long double expal1;
    // long double expal2;
    // long double binomE2;


    // for(int i=0;i<sizeAr;i++){
    // 	cout<<i<<"\t"<<obsBase[i]<<"\t"<<obsQual[i]<<"\t"<<probDeam[i]<<endl;
    // }

    long double llikCCtotal=0.0;
    long double llikCTtotal=0.0;
    long double llikTTtotal=0.0;

    for(int i=0;i<2;i++){
	char al1   = alleles[0];
	char al2   = alleles[1];
	char alC   = alleles[i];


	//CC
	al1Current   = al1;
	al2Current   = al1;
	// long double llikCC=0;
	// long double llikCC1=0;
	// long double llikCC2=0;

	long double llikCC=computeLL(al1Current,
				     al2Current,		      
				     sizeAr    ,
				     obsBase   ,
				     probDeam  ,
				     obsQual   ,
				     contRate  ,
				     alC       ,
				     mmProb);

    

	//CT
	al1Current   = al1;
	al2Current   = al2;

	long double llikCT=computeLL(al1Current,
				     al2Current,		      
				     sizeAr    ,
				     obsBase   ,
				     probDeam  ,
				     obsQual   ,
				     contRate  ,
				     alC       ,
				     mmProb);


	//TT
	al1Current   = al2;
	al2Current   = al2;

	long double llikTT=computeLL(al1Current,
				     al2Current,		      
				     sizeAr    ,
				     obsBase   ,
				     probDeam  ,
				     obsQual   ,
				     contRate  ,
				     alC       ,
				     mmProb);

	cout<<alC<<"\t"<<llikCC<<","<<llikCT<<","<<llikTT<<endl;
	llikCCtotal = oplusInitnatl(llikCCtotal,llikCC)+logl(0.5); //for *0.5
	llikCTtotal = oplusInitnatl(llikCTtotal,llikCT)+logl(0.5); //for *0.5
	llikTTtotal = oplusInitnatl(llikTTtotal,llikTT)+logl(0.5); //for *0.5

	//
	// long double probHeto  = expl( llikCT-sumllik );
	// long double probHomo  = 1-probHeto;

	// // cout<<expl( llikCT-sumllik )<<endl;
	// // cout<<expl( llikTT-sumllik )<<endl;
	// cout<<probHomo<<"\t"<<probHeto<<endl;
    }
    cout<<llikCCtotal<<","<<llikCTtotal<<","<<llikTTtotal<<endl;
    long double bestGTllk     = llikCCtotal>llikTTtotal?llikCCtotal:llikTTtotal;

    long double sumHetPlusHom = oplusnatl(bestGTllk,llikCTtotal);
    long double ratioHet      = expl(llikCTtotal)/expl(sumHetPlusHom);
    
    cout<<bestGTllk<<"\t"<<sumHetPlusHom<<"\tp(het)=\t"<<ratioHet<<endl;
    cout<<abs(bestGTllk-llikCTtotal)<<"\t"<<1.0-expl(-1.0*abs(bestGTllk-llikCTtotal))<<"\t"<<1-expl(-50)<<"\t"<<1-expl(-3)<<endl;
    long double probAcc = 1.0-expl(-1.0*abs(bestGTllk-llikCTtotal));
    cout<<"qual\t"<<probAcc<<endl;

    //estimate h
    for(long double h=0;h<1;h+=0.01){
	long double ll=0.0;

	//for(int i=0;i<size;i++){
	long double llT = logl( 
			       (1-h)*(probAcc*(1-ratioHet) + (1-probAcc)*0.5)
			       + 
			       (  h)*(probAcc*(  ratioHet) + (1-probAcc)*0.5)
				);
	ll+=llT;
	//}
	
	cout<<h<<"\t"<<(1-h)<<"\t"<<ll<<endl;	
    }



    //idea: if qual bad, revert to 0.5. if qual good, revert to ratioHet
    return 0;
}

