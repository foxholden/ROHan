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


int main (int argc, char *argv[]) {
    int size           = 6;
    //long double het [] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
    long double het [] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    long double pEr [] = {1.0, 0.0, 0.0, 0.0, 1.0, 1.0};
    
    for(long double h=0;h<1;h+=0.01){
	long double ll=0.0;

	for(int i=0;i<size;i++){
	    long double llT = logl( 
				   (1-pEr[i])*((1-h)*(1-het[i]) + (h)*(het[i]) )
				   +
				   (pEr[i])*( 0.5 )
				    );
	    ll+=llT;
	}
	
	cout<<h<<"\t"<<(1-h)<<"\t"<<ll<<endl;	
    }

    return 0;
}

