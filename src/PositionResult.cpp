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
    //cerr<<"Constructor addr: "<<this<<endl;
    string toReturn="";
    toReturn += ""+references[refID].RefName+"\t";
    toReturn += ""+stringify(pos)+"\t";

    toReturn += ""+stringify(refB)+"\t";
    toReturn += ""+stringify(altB)+"\t";

    toReturn += ""+stringify(refC)+"\t";
    toReturn += ""+stringify(altC)+"\t";

    if(geno==0){
	toReturn += "0/0\t";
    }else{
	if(geno==1){
	    toReturn += "0/1\t";
	}else{
	    if(geno==2){
		toReturn += "1/1\t";
	    }else{
		cerr<<"Internal error for genotype"<<endl;
		exit(1);
	    }	    
	}
    }


    toReturn += ""+stringify(rrll)+"\t";
    toReturn += ""+stringify(rall)+"\t";
    toReturn += ""+stringify(aall)+"\t";

    toReturn += ""+stringify(lqual)+"\t";
    toReturn += ""+stringify(llCov)+"\t";

    toReturn += "\n";

    //cout<<"toString "<<toReturn<<endl;
    return toReturn;
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


