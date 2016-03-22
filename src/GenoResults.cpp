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


