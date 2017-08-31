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

pair<char,char> PositionResult::hetIndex2Bases(){
    if(mostLikelyGenoHetIdx==1)
	return make_pair<char,char>('A','C');
    if(mostLikelyGenoHetIdx==2)
	return make_pair<char,char>('A','G');
    if(mostLikelyGenoHetIdx==3)
	return make_pair<char,char>('A','T');
    if(mostLikelyGenoHetIdx==5)
	return make_pair<char,char>('C','G');
    if(mostLikelyGenoHetIdx==6)
	return make_pair<char,char>('C','T');
    if(mostLikelyGenoHetIdx==8)
	return make_pair<char,char>('G','T');
    cerr<<"PositionResult: wrong state"<<endl;
    exit(1);
    
}


char PositionResult::homoIndex2Base(){
    if(mostLikelyGenoIdx==0)
	return 'A';
    if(mostLikelyGenoIdx==4)
	return 'C';
    if(mostLikelyGenoIdx==7)
	return 'G';
    if(mostLikelyGenoIdx==9)
	return 'T';
    cerr<<"PositionResult: wrong state"<<endl;
    exit(1);
    
}

string PositionResult::toString(const RefVector  * references, const int & refID) const{
    //cerr<<"toString: "<<endl;
    stringstream s;
    s.precision(13);
    // cerr<<refID<<endl;
    //string toReturn="";
    s<<references->at(refID).RefName<<"\t";
    s<<pos<<"\t";

    // for(int n=0;n<4;n++)
    // 	s<<baseC[n]<<"\t";
    //REF
    s<<refB<<"\t";

    string altB;
    
    if(mostLikelyGenoIdx ==  mostLikelyGenoHetIdx){//heterozygous
	bool hasRefAsAlt=false;
	pair<char,char> hetBase = hetIndex2Bases();
	hasRefAsAlt = (refB == hetBase.first) || (refB == hetBase.second) ;
	if(hasRefAsAlt){
	    if( refB == hetBase.first )
		altB = string(hetBase.second);
	    else
		altB = string(hetBase.first);
	    
	}else{
	    altB = string(hetBase.first)+","+string(hetBase.second);
	}
    }else{//homozygous
	char bhomo=homoIndex2Base();
	if(bhomo!=refB){
	    altB=string(bhomo);
	}else{
	    pair<char,char> hetBase = hetIndex2Bases();
	    hasRefAsAlt = (refB == hetBase.first) || (refB == hetBase.second) ;
	    if(hasRefAsAlt){
		if( refB == hetBase.first )
		    altB = string(hetBase.second);
		else
		    altB = string(hetBase.first);
	    }else{
		cerr<<"ERROR at site "<<references->at(refID).RefName<<":"<<pos<<endl;

	    }
	}
    }

    //find if has ref
    s<<altB<<"\t";

    s<<"GT:AD:DP:GQ:PL\t";

    s<<gq<<"\t.\t";//the . is the filter
    

    // s<<refC<<"\t";
    // s<<altC<<"\t";

    // if(geno==0){
    // 	s<<"0/0\t";
    // }else{
    // 	if(geno==1){
    // 	    s<<"0/1\t";
    // 	}else{
    // 	    if(geno==2){
    // 		s<<"1/1\t";
    // 	    }else{
    // 		cerr<<"Internal error for genotype ="<<geno<<endl;
    // 		exit(1);
    // 	    }	    
    // 	}
    // }

    // s<<genoS[0]<<genoS[1]<<"\t";

    // s<<lqual<<"\t";
    // s<<llCov<<"\t";

     for(int g=0;g<10;g++)
     	s<<ll[g]<<"\t";	
    // // s<<rrll<<"\t";
    // // s<<rall<<"\t";
    // // s<<aall<<"\t";


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



    // if(refB=='A'){
    // 	hasRefAsAlt = hasRefAsAlt || (mostLikelyGenoHetIdx==1) || (mostLikelyGenoHetIdx==2) || (mostLikelyGenoHetIdx==3) ;
    // }else{
    // if(refB=='C'){
    // 	hasRefAsAlt = hasRefAsAlt || (mostLikelyGenoHetIdx==1) || (mostLikelyGenoHetIdx==5) || (mostLikelyGenoHetIdx==6) ;
    // }else{
    // if(refB=='G'){
    // 	hasRefAsAlt = hasRefAsAlt || (mostLikelyGenoHetIdx==2) || (mostLikelyGenoHetIdx==5) || (mostLikelyGenoHetIdx==8) ;
    // }else{
    // if(refB=='T'){
    // 	hasRefAsAlt = hasRefAsAlt || (mostLikelyGenoHetIdx==3) || (mostLikelyGenoHetIdx==6) || (mostLikelyGenoHetIdx==8) ;
    // }else{
    // }}}}
