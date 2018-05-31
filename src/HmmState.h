/*
 * HmmState
 * Date: Dec-04-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef HmmState_h
#define HmmState_h
#include <iostream>
#include <vector>
#include <limits>

#include <ctime>//for time cmd
#include <math.h>

#include <gsl/gsl_randist.h>

#include "miscfunc.h"

using namespace std;

class HmmState{
private:
    long double h;
    long double theta;
    long double rateGeom;
    unsigned int nrwPerSizeChunk;
    
    int minSegSitesPerChunk;
    int maxSegSitesPerChunk;
    int sizeChunk;

    vector<long double> * probabilitiesForEmission;
    HmmState * second;
    gsl_rng * rng;
    int idx;
    
public:
    HmmState(int idx_,long double p_,int minSegSitesPerChunk_,int maxSegSitesPerChunk_,int sizeChunk_,unsigned int nrwPerSizeChunk_);
    HmmState(const HmmState & other);
    ~HmmState();
    HmmState & operator= (const HmmState & other);
    
    /* HmmState & operator= (const HmmState & other){ */
    /* 	p      = other.p; */
    /* 	second = other.second; */
	
    /* } */

    long double probEmission(     const int mutations,                          const int total) const;
    long double probEmissionRange(const int mutationsMin,const int mutationsMax,const int total) const;
    
    unsigned int randomEmission(int total) const;
    int getIdx() const;
    long double getH() const;
    long double getTheta() const;
    unsigned int getNrwPerSizeChunk() const;
    long double getRateGeom() const;
    
    void setH(long double newH);
    void setNrwPerSizeChunk(unsigned int nrwPerSizeChunk_); //nrw = non-recombining window
    void recomputeProbs(bool verbose=false);
    void setSecond(HmmState * second);
};
#endif
