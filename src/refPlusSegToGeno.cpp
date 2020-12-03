/*
 * refPlusSegToGeno
 * Date: Apr-13-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>

#include "libgab.h"

using namespace std;

typedef struct{
    string ref;
    unsigned int coordinate;
    char al1;
    char al2;
} segsite;

int main (int argc, char *argv[]) {


   igzstream myFileSeg;
   igzstream myFileRef;

   string refFile  = string(argv[1]);
   string segSites = string(argv[2]);

   string lineRef;
   string lineSeg;
   
   myFileRef.open(refFile.c_str(), ios::in);
   myFileSeg.open(segSites.c_str(), ios::in);

   segsite nextSS;

   if (myFileSeg.good()){
       if(getline (myFileSeg,lineSeg)){
	   vector<string> token=allTokens(lineSeg,'\t');
			 
	   nextSS.ref        =                      token[0];
	   nextSS.coordinate = destringify<unsigned int>(token[1]);
	   nextSS.al1        = destringify<char>(        token[2]);
	   nextSS.al2        = destringify<char>(        token[3]);
       }else{
	   nextSS.ref        =                      "-1";
       }
       
   }else{
       cerr << "Unable to open file "<<segSites<<endl;
       return 1;
    }


   //   return 1;

   string currentRef="-1";
   unsigned int coordinate=0;

   if (myFileRef.good()){
     while ( getline (myFileRef,lineRef)){
	 //cout<<lineRef<<endl;
	 if(lineRef[0] == '>'){	     
	     currentRef=lineRef;
	     coordinate=0;
	 }else{
	     for(unsigned int i=0;i<lineRef.size();i++){
		 coordinate++;

		 if(currentRef == nextSS.ref 
		    &&
		    coordinate == nextSS.coordinate ){
		     if(lineRef[i] != nextSS.al1 &&
			lineRef[i] != nextSS.al2 ){
			 cerr<<"internal error"<<endl;
			 return 1;
		     }	
		     
		     if(nextSS.al1<nextSS.al2)
			 cout<<""<<currentRef<<"\t"<<coordinate<<"\t"<<nextSS.al1<<nextSS.al2<<endl;
		     else
			 cout<<""<<currentRef<<"\t"<<coordinate<<"\t"<<nextSS.al2<<nextSS.al1<<endl;
		     
		     if(getline (myFileSeg,lineSeg)){
			 vector<string> token=allTokens(lineSeg,'\t');
			 
			 nextSS.ref        =                      token[0];
			 nextSS.coordinate = destringify<unsigned int>(token[1]);
			 nextSS.al1        = destringify<char>(        token[2]);
			 nextSS.al2        = destringify<char>(        token[3]);
		     }else{
			 nextSS.ref        =                      "-1";
		     }

		 }else{
		     cout<<currentRef<<"\t"<<coordinate<<"\t"<<lineRef[i]<<lineRef[i]<<endl;
		 }

	     }
	     
	 }
     }
     myFileRef.close();
   }else{
       cerr << "Unable to open file "<<refFile<<endl;
       return 1;
    }
    return 0;
}

