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
    cov=0;
    for(int n=0;n<4;n++){
	cov+=pr->baseC[n];
    }


    for(int g=0;g<10;g++)
	ll[g] = pr->ll[g];
    // rrll   = pr->rrll ;
    // rall   = pr->rall ;
    // aall   = pr->aall ;
    lqual  = pr->lqual;
    llCov  = pr->llCov;
    geno   = pr->geno ;    
}

GenoResults::GenoResults(const string lineToAdd){
    vector<string> tokensS=allTokens(lineToAdd,'\t');

    cov=0;
    for(int n=0;n<4;n++){
	cov+=destringify<int>(tokensS[2+n]);
    }
   
    if(         tokensS[6+0] == "0/0"){
	geno=0;
    }else{
	if(     tokensS[6+0] == "0/1"){
	    geno=1;
	}else{
	    if( tokensS[6+0] == "1/1"){
		geno=2;
	    }else{
		cerr<<"Cannot convert line "<<lineToAdd<<" into a GenoResults object, wrong genotype:"<<tokensS[6+0]<<endl;
		exit(1);
	    }	    
	}
    }

    genoS[0]= tokensS[6+1][0];
    genoS[1]= tokensS[6+1][1];

    lqual  = destringify<long double>(tokensS[6+2]);
    llCov  = destringify<long double>(tokensS[6+3]);
    
    for(int g=0;g<10;g++){
	ll[g] = destringify<long double>(tokensS[6+4+g]);
    }
    // rrll   = destringify<long double>(tokensS[6+1]);
    // rall   = destringify<long double>(tokensS[6+2]);
    // aall   = destringify<long double>(tokensS[6+3]);
    





}

// GenoResults::GenoResults(long double  rrll ,
// 			 long double  rall ,
// 			 long double  aall ,
// 			 long double  lqual,
// 			 long double  llCov,
// 			 int          geno): 
//     rrll(  rrll  ),
//     rall(  rall  ),
//     aall(  aall  ),
//     lqual( lqual ),
//     llCov( llCov ),
//     geno(  geno  ){
    
// }

GenoResults::~GenoResults(){
    //cerr<<"Destructor  addr: "<<this<<endl;
}




ostream & operator << (ostream & os, const GenoResults & gr){
    os<<"gen= "<<gr.geno<<"\t"
	"ges= "<<gr.genoS<<"\t"

	"q=   "<<gr.lqual<<"\t"
	"cov= "<<gr.llCov<<"\t"


	"AA= "<<gr.ll[0]<<"\t"
 	"AC= "<<gr.ll[1]<<"\t"
	"AG= "<<gr.ll[2]<<"\t"
	"AT= "<<gr.ll[3]<<"\t"

 	"CC= "<<gr.ll[4]<<"\t"
	"CG= "<<gr.ll[5]<<"\t"
	"CT= "<<gr.ll[6]<<"\t"

	"GG= "<<gr.ll[7]<<"\t"
	"GT= "<<gr.ll[8]<<"\t"

	"TT= "<<gr.ll[9];
	

	// "pacc= "<<gr.probAccurate<<"\t"
	// "exp= "<<gr.expectedH;

    return os;
}
