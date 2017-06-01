/*
 * ref2DNAprof
 * Date: Jun-01-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include <iostream>
#include <fstream>
#include <gzstream.h>

#include "utils.h"

using namespace std;

int main (int argc, char *argv[]) {

    const string usage=string("\nThis program computes the natural occurrence of the 4 DNA bases in a reference file\n\n\t"+
                              string(argv[0])+                        
                              " [options] [fasta file]  "+"\n\n"+
                              "\twhere:\n"+
                              "\t\t[fasta file]\t\tThe fasta file used for alignement\n"+
                              "\n\n"+
                              
                              
                              // "\n\tDeamination and error options:\n"+                                   
                              // "\t\t"+""  +""+"--deam5p\t\t"+"[.prof file]" +"\t\t"+"5p deamination frequency for the endogenous\n\t\t\t\t\t\t\t\t(default: "+deam5pfreqE+")"+"\n"+
                              // "\t\t"+""  +""+"--deam3p\t\t"+"[.prof file]" +"\t\t"+"3p deamination frequency for the endogenous\n\t\t\t\t\t\t\t\t(default: "+deam3pfreqE+")"+"\n"+
                              // // "\t\t"+"-deam5pc [.prof file]" +"\t\t"+"5p deamination frequency for the contaminant (default: "+deam5pfreqC+")"+"\n"+
                              // // "\t\t"+"-deam3pc [.prof file]" +"\t\t"+"3p deamination frequency for the contaminant (default: "+deam3pfreqC+")"+"\n"+                       
                              // "\t\t"+""  +""+"--err\t\t\t"    +"[.prof file]"+"\t\t"    +" Illumina error profile (default: "+illuminafreq+")"+"\n"+
                              "");


    if( (argc== 1) ||
        (argc== 2 && string(argv[1]) == "-h") ||
        (argc== 2 && string(argv[1]) == "-help") ||
        (argc== 2 && string(argv[1]) == "--help") ){
        cout<<usage<<endl;
        return 1;
    }

    string line;
    igzstream myFile;
    string filename = string(argv[1]);
    myFile.open(filename.c_str(), ios::in);
    uint64_t totalBases=0;
    uint64_t totalA    =0;
    uint64_t totalC    =0;
    uint64_t totalG    =0;
    uint64_t totalT    =0;

    if (myFile.good()){
	while ( getline (myFile,line)){
	    if( (line.size()>1) &&
		(line[0] == '>')){
		continue;
	    }
	    for(unsigned int i=0;i<line.size();i++){
		char c = toupper(line[i]);
		if( c == 'A'){
		    totalA++;
		    totalBases++;		 
		}
		if( c == 'C'){
		    totalC++;
		    totalBases++;		 
		}
		if( c == 'G'){
		    totalG++;
		    totalBases++;		 
		}
		if( c == 'T'){
		    totalT++;
		    totalBases++;		 
		}

	    }
	}
	myFile.close();
    }else{
	cerr << "Unable to open file "<<filename<<endl;
	return 1;
    }

    cout<<double(totalA)/double(totalBases)<<endl;
    cout<<double(totalC)/double(totalBases)<<endl;
    cout<<double(totalG)/double(totalBases)<<endl;
    cout<<double(totalT)/double(totalBases)<<endl;

    return 0;
}

