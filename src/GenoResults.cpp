/*
 * GenoResults
 * Date: Mar-21-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "GenoResults.h"

GenoResults::GenoResults(){
    //cerr<<"Constructor addr: "<<this<<endl;
}

GenoResults::GenoResults(const PositionResult * pr){
    rrll   = pr->rrll ;
    rall   = pr->rall ;
    aall   = pr->aall ;
    lqual  = pr->lqual;
    llCov  = pr->llCov;
    geno   = pr->geno ;    
}

GenoResults::GenoResults(const string lineToAdd){
    vector<string> tokensS=allTokens(lineToAdd,'\t');

    rrll   = destringify<long double>(tokensS[6+1]);
    rall   = destringify<long double>(tokensS[6+2]);
    aall   = destringify<long double>(tokensS[6+3]);
    lqual  = destringify<long double>(tokensS[6+4]);
    llCov  = destringify<long double>(tokensS[6+5]);


    if(tokensS[6+0] == "0/0"){
	geno=0;
    }else{
	if(tokensS[6+0] == "0/1"){
	    geno=1;
	}else{
	    if(tokensS[6+0] == "1/1"){
		geno=2;
	    }else{
		cerr<<"Cannot convert line "<<lineToAdd<<" into a GenoResults object, wrong genotype:"<<tokensS[6+0]<<endl;
		exit(1);
	    }	    
	}
    }

}

GenoResults::GenoResults(long double  rrll ,
			 long double  rall ,
			 long double  aall ,
			 long double  lqual,
			 long double  llCov,
			 int          geno): 
    rrll(  rrll  ),
    rall(  rall  ),
    aall(  aall  ),
    lqual( lqual ),
    llCov( llCov ),
    geno(  geno  ){
    
}

GenoResults::~GenoResults(){
    //cerr<<"Destructor  addr: "<<this<<endl;
}


