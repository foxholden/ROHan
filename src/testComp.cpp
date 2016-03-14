/*
 * testComp
 * Date: Mar-06-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>

#include "utils.h"

using namespace std;

#define MAXMAPPINGQUAL 257     // maximal mapping quality, should be sufficient as mapping qualities are encoded using 8 bits

long double likeMatch        [MAXMAPPINGQUAL];
long double likeMismatch     [MAXMAPPINGQUAL];

long double likeMatchProb    [MAXMAPPINGQUAL];
long double likeMismatchProb [MAXMAPPINGQUAL];

// Returns logl( expl(x)+expl(y) )
inline long double oplusl(long double x, long double y ){
    return x > y 
        ? x + log1pl( expl( y-x ) )
        : y + log1pl( expl( x-y ) )  ;
}


//! A method to initialize various probability scores to avoid recomputation
/*!
  This method is called by the main after capturing the arguments
*/
void initScores(){

    for(int i=0;i<2;i++){
        likeMatch[i]        = log1pl(    -pow(10.0,2.0/-10.0) )   ;         
        likeMismatch[i]     = logl  (     pow(10.0,2.0/-10.0)/3.0 );

	likeMatchProb[i]           = 1.0-pow(10.0,2.0/-10.0) ;
        likeMismatchProb[i]        =     pow(10.0,2.0/-10.0)/3.0 ;
    }


    //Computing for quality scores 2 and up
    for(int i=2;i<MAXMAPPINGQUAL;i++){
        likeMatch[i]        = log1pl(    -pow(10.0,i/-10.0) )     ;          
        likeMismatch[i]     = logl  (     pow(10.0,i/-10.0)/3.0  );

        likeMatchProb[i]           = 1.0-pow(10.0,i/-10.0);
        likeMismatchProb[i]        =     pow(10.0,i/-10.0)/3.0;
    }



}

int main (int argc, char *argv[]) {
    initScores();
    // cout<<expl(oplusl( logl(0.4),logl(0.5)))<<endl;

    // return 1;
    //cout<<nChoosek(20,10)<<endl;
    char al1   = 'C';
    char al2   = 'T';
    char al1Current   = ' ';
    char al2Current   = ' ';

    int sizeAr =  5;
    for(int i=0;i<=sizeAr;i++){
	cout<<sizeAr<<"\t"<<i<<"\t"<<nChoosek(sizeAr,i)<<"\t"<<(nChoosek(sizeAr,i)*powl(0.5,sizeAr-i)*powl(0.5,i))<<endl;
    }

    //char obsBase [ ] =  { 'C', 'C', 'C', 'C', 'C', 'C' }; 
    // char obsBase          [ ] =  { 'C', 'C', 'T', 'T', 'T', 'T' }; 
    // int  obsQual          [ ] =  {  40,  40,  40,  40,  40,  40 }; 
    // long double probDeam  [ ] =  { 0.0, 0.0, 0.0, 0.0, 0.2, 0.2 }; 
    char obsBase          [ ] =  { 'C', 'C', 'C', 'C', 'T' }; 
    //int  obsQual          [ ] =  {  40,  40,  40,  40,  40 }; 
    long double probDeam  [ ] =  { 0.0, 0.0, 0.0, 0.0, 0.0 }; 
    int  obsQual          [ ] =  {  2,   2,   2,   2,   2 }; 

    //int  obsQual [ ] =  {  2,   2,   2,   2,   2,   2 }; 

    long double expal1;
    long double expal2;
    long double binomE2;
    
    //CC
    al1Current   = al1;
    al2Current   = al1;
    long double llikCC=0;
    long double llikCC1=0;
    long double llikCC2=0;

    for(int i=0;i<sizeAr;i++){
	cout<<i<<"\t"<<obsBase[i]<<"\t"<<obsQual[i]<<"\t"<<probDeam[i]<<endl;
    }

    for(int i=0;i<sizeAr;i++){
	//cout<<i<<"\t"<<obsBase[i]<<"\t"<<obsQual[i]<<endl;
	long double llikAl1t=0;
	long double llikAl2t=0;

	if(obsBase[i] == al1Current){
	    llikAl1t = (1.0-probDeam[i])*likeMatchProb[obsQual[i]] + 0*likeMismatchProb[obsQual[i]];
	}else{
	    llikAl1t =       probDeam[i]*likeMatchProb[obsQual[i]] + 1*likeMismatchProb[obsQual[i]];
	}


	if(obsBase[i] == al2Current){
	    llikAl2t = (1.0-probDeam[i])*likeMatchProb[obsQual[i]] + 0*likeMismatchProb[obsQual[i]];
	}else{
	    llikAl2t =      probDeam[i]*likeMatchProb[obsQual[i]] + 1*likeMismatchProb[obsQual[i]];
	}


	llikCC +=logl( 0.5*llikAl1t + 0.5*llikAl2t );
	llikCC1+=logl( llikAl1t );
	llikCC2+=logl( llikAl2t );
    }
    //cout<<"CC\t"<<llikCC<<"\t"<<llikCC1<<"\t"<<llikCC2<<endl;    
    expal1=roundl(sizeAr*(llikCC1 / (llikCC1+llikCC2)));
    expal2=sizeAr-expal1;
    binomE2= logl(nChoosek(sizeAr,expal1)*powl(0.5,expal1)*powl(0.5,expal2));
    cout<<"CC\t"<<llikCC<<"\t"<<expal1<<"\t"<<expal2<<"\t"<<binomE2<<"\t"<<(binomE2+llikCC)<<endl;
    

    //CT
    al1Current   = al1;
    al2Current   = al2;

    long double llikCT=0;
    long double llikCT1=0;
    long double llikCT2=0;

    for(int i=0;i<sizeAr;i++){
	///cout<<i<<"\t"<<obsBase[i]<<"\t"<<obsQual[i]<<endl;

	long double llikAl1t=0;
	long double llikAl2t=0;


	if(obsBase[i] == al1Current){
	    llikAl1t = (1.0-probDeam[i])*likeMatchProb[obsQual[i]] + 0*likeMismatchProb[obsQual[i]];
	}else{
	    llikAl1t =       probDeam[i]*likeMatchProb[obsQual[i]] + 1*likeMismatchProb[obsQual[i]];
	}

	if(obsBase[i] == al2Current){
	    llikAl2t = (1.0-probDeam[i])*likeMatchProb[obsQual[i]] + 0*likeMismatchProb[obsQual[i]];
	}else{
	    llikAl2t =       probDeam[i]*likeMatchProb[obsQual[i]] + 1*likeMismatchProb[obsQual[i]];
	}

	llikCT+=logl( 0.5*llikAl1t + 0.5*llikAl2t );
	llikCT1+=logl( llikAl1t );
	llikCT2+=logl( llikAl2t );
    }
    //cout<<"CT\t"<<llikCT<<"\t"<<llikCT1<<"\t"<<llikCT2<<endl;

    // cout<<oplusl(llikCT1,llikCT2)<<endl;
    // cout<<llikCT1 / (llikCT1+llikCT2)<<endl;
    // cout<<llikCT2 / (llikCT1+llikCT2)<<endl;
    //cout<<sizeAr*(llikCT1 / (llikCT1+llikCT2))<<endl;
    // cout<<sizeAr*(llikCT2 / (llikCT1+llikCT2))<<endl;

    // cout<<llikCT1 - oplusl(llikCT1,llikCT2)<<endl;
    // cout<<expl( llikCT1 - oplusl(llikCT1,llikCT2))<<endl;
    // cout<<expl( llikCT1 - oplusl(llikCT1,llikCT2))*sizeAr<<endl;

    expal1=roundl(sizeAr*(llikCT1 / (llikCT1+llikCT2)));
    expal2=sizeAr-expal1;
    binomE2= logl(nChoosek(sizeAr,expal1)*powl(0.5,expal1)*powl(0.5,expal2));
    cout<<"CT\t"<<llikCT<<"\t"<<expal1<<"\t"<<expal2<<"\t"<<binomE2<<"\t"<<(binomE2+llikCT)<<endl;


    //    cout<<expal1<<"\t"<<expal2<<endl;

    //TT
    al1Current   = al2;
    al2Current   = al2;

    long double llikTT=0;
    long double llikTT1=0;
    long double llikTT2=0;

    for(int i=0;i<sizeAr;i++){
	//cout<<i<<"\t"<<obsBase[i]<<"\t"<<obsQual[i]<<endl;

	long double llikAl1t=0;
	long double llikAl2t=0;


	if(obsBase[i] == al1Current){
	    llikAl1t = (1.0-probDeam[i])*likeMatchProb[obsQual[i]] + 0*likeMismatchProb[obsQual[i]];
	}else{
	    llikAl1t =      probDeam[i]*likeMatchProb[obsQual[i]] + 1*likeMismatchProb[obsQual[i]];
	}

	if(obsBase[i] == al2Current){
	    llikAl2t = (1.0-probDeam[i])*likeMatchProb[obsQual[i]] + 0*likeMismatchProb[obsQual[i]];
	}else{
	    llikAl2t =       probDeam[i]*likeMatchProb[obsQual[i]] + 1*likeMismatchProb[obsQual[i]];
	}

	llikTT+=logl( 0.5*llikAl1t + 0.5*llikAl2t );
	llikTT1+=logl( llikAl1t );
	llikTT2+=logl( llikAl2t );
    }
    //cout<<"TT\t"<<llikTT<<"\t"<<llikTT1<<"\t"<<llikTT2<<endl;
    expal1=roundl(sizeAr*(llikTT1 / (llikTT1+llikTT2)));
    expal2=sizeAr-expal1;
    //    cout<<expal1<<"\t"<<expal2<<endl;
    binomE2= logl(nChoosek(sizeAr,expal1)*powl(0.5,expal1)*powl(0.5,expal2));
    cout<<"TT\t"<<llikTT<<"\t"<<expal1<<"\t"<<expal2<<"\t"<<binomE2<<"\t"<<(binomE2+llikTT)<<endl;





    return 0;
}

