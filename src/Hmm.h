/*
 * Hmm
 * Date: Dec-05-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef Hmm_h
#define Hmm_h

#include <iostream>
#include <vector>
#include "HmmState.h"

#include "utils.h"

#define NBRSTATES 2

using namespace std;

typedef struct { 
    double p;
    int idx;
} emission;


class Hmm{
private:
    double **trans;

public:

    long double startingState[NBRSTATES];
    long double probTrans[NBRSTATES];
    long double probStay[NBRSTATES];
    HmmState * hmmstates[NBRSTATES];
	
    Hmm();
    Hmm(const Hmm & other);
    ~Hmm();
    Hmm & operator= (const Hmm & other);

    vector<emission> generateStates(unsigned int N,unsigned int total) const;
    int getNumberStates() const;
    double getTrans(int i,int j) const;
    
};
#endif
