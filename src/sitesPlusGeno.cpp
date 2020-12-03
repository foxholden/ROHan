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

typedef struct{
    string ref;
    unsigned int coordinate;
    char al1;
    char al2;
    double qual;
} geno;

int main (int argc, char *argv[]) {


    igzstream myFileGen;
    igzstream myFileSeg;
    

    string genFile   = string(argv[1]);
    string segSites  = string(argv[2]);

    string lineGen;
    string lineSeg;
   
    myFileGen.open(genFile.c_str(),  ios::in);
    myFileSeg.open(segSites.c_str(), ios::in);





    if (myFileGen.good()){
	
    }else{
    	cerr << "Unable to open file "<<genFile<<endl;
    	return 1;
    }

    if (myFileSeg.good()){

    }else{
	cerr << "Unable to open file "<<segSites<<endl;
	return 1;
    }



    //   return 1;
    getline (myFileGen,lineGen);//header 

    string       currentRef = "-1";
    //    unsigned int coordinate =    0;


    while(getline (myFileGen,lineGen)){
	vector<string> tokenG=allTokens(lineGen,'\t');
	//cout<<"G"<<lineGen<<endl;
	geno g;
	g.ref        =                           tokenG[0];
	g.coordinate = destringify<unsigned int>(tokenG[1]);
	g.al1        =                           tokenG[7][0];
	g.al2        =                           tokenG[7][1];
	g.qual       = destringify<double>(      tokenG[8] );


	while ( getline (myFileSeg,lineSeg)){
	    //cout<<"S"<<lineSeg<<endl;
	    vector<string> token=allTokens(lineSeg,'\t');

	    segsite s;//truth
	    s.ref        =                           token[0];
	    s.coordinate = destringify<unsigned int>(token[1]);
	    s.al1        =                           token[2][0];
	    s.al2        =                           token[2][1];
	    
	    if( (">"+g.ref)        == s.ref && 
		g.coordinate == s.coordinate){
		if(s.al1 == g.al1 //ok
		   &&
		   s.al2 == g.al2 ){
		    
		    if(s.al1 == s.al2)//truth is homo
			cout<<"CHOM2HOM\t"<<g.ref<<"\t"<<g.coordinate<<"\t"<<g.al1<<g.al2<<"\t"<<s.al1<<s.al2<<"\t"<<g.qual<<"\t"<<lineGen<<endl;
		    else
			cout<<"CHET2HET\t"<<g.ref<<"\t"<<g.coordinate<<"\t"<<g.al1<<g.al2<<"\t"<<s.al1<<s.al2<<"\t"<<g.qual<<"\t"<<lineGen<<endl;
		    
		}else{//error
		    
		    if(s.al1 == s.al2){//truth is homo
			if(g.al1 == g.al2){//predicted is homo
			    cout<<"WHOM2HOM\t"<<g.ref<<"\t"<<g.coordinate<<"\t"<<g.al1<<g.al2<<"\t"<<s.al1<<s.al2<<"\t"<<g.qual<<"\t"<<lineGen<<endl;
			}else{
			    cout<<"WHOM2HET\t"<<g.ref<<"\t"<<g.coordinate<<"\t"<<g.al1<<g.al2<<"\t"<<s.al1<<s.al2<<"\t"<<g.qual<<"\t"<<lineGen<<endl;
			}
		    }else{//truth is hetero

			if(g.al1 == g.al2){//predicted is homo
			    cout<<"WHET2HOM\t"<<g.ref<<"\t"<<g.coordinate<<"\t"<<g.al1<<g.al2<<"\t"<<s.al1<<s.al2<<"\t"<<g.qual<<"\t"<<lineGen<<endl;
			}else{
			    cout<<"WHET2HET\t"<<g.ref<<"\t"<<g.coordinate<<"\t"<<g.al1<<g.al2<<"\t"<<s.al1<<s.al2<<"\t"<<g.qual<<"\t"<<lineGen<<endl;
			}
			
		    }
		}
		break;
	    }
		  
	}

    }
	    //cout<<lineSeg<<endl;

    if (myFileSeg.good()){
	myFileSeg.close();
    }else{
	cerr << "Unable to open file "<<segSites<<endl;
	return 1;
    }

    return 0;
}

