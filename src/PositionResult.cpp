/*
 * PositionResult
 * Date: Mar-21-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "PositionResult.h"

PositionResult::PositionResult(){
    //cerr<<"Constructor addr: "<<this<<endl;
}

PositionResult::~PositionResult(){
    //cerr<<"Destructor  addr: "<<this<<endl;
}


string PositionResult::toString(const RefVector  references) const{
    stringstream s;
    s.precision(13);
    //cerr<<"Constructor addr: "<<this<<endl;
    //string toReturn="";
    s<<references[refID].RefName<<"\t";
    s<<pos<<"\t";

    for(int n=0;n<4;n++)
	s<<baseC[n]<<"\t";
    // s<<refB<<"\t";
    // s<<altB<<"\t";

    // s<<refC<<"\t";
    // s<<altC<<"\t";

    if(geno==0){
	s<<"0/0\t";
    }else{
	if(geno==1){
	    s<<"0/1\t";
	}else{
	    if(geno==2){
		s<<"1/1\t";
	    }else{
		cerr<<"Internal error for genotype ="<<geno<<endl;
		exit(1);
	    }	    
	}
    }

    s<<genoS[0]<<genoS[1]<<"\t";

    s<<lqual<<"\t";
    s<<llCov<<"\t";

    for(int g=0;g<10;g++)
	s<<ll[g]<<"\t";	
    // s<<rrll<<"\t";
    // s<<rall<<"\t";
    // s<<aall<<"\t";


    s<<"\n";

    //cout<<"toString "<<toReturn<<endl;
    return s.str();
}

//TODO code to VCF
// string PositionResult::toString(const RefVector  references) const{
//     //cerr<<"Constructor addr: "<<this<<endl;
//     string toReturn="";
//     toReturn += ""+references[refID].RefName+"\t";
//     toReturn += ""+stringify(pos)+"\t";
//     toReturn += ""+stringify(refB)+"\t";
//     toReturn += ""+stringify(altB)+"\t";
//     toReturn += ".\t"; //ID
//     toReturn += "0\t"; //QUAL   
//     toReturn += "\t";

//     toReturn += "\n";

//     //cout<<"toString "<<toReturn<<endl;
//     return toReturn;
// }

// ostream & operator << (ostream & os, const PositionResult & pr){


//     os<<pr.toString();


//     return os;
// }


