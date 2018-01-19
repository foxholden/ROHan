/*
 * HmmState
 * Date: Dec-04-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef HmmState_h
#define HmmState_h
#include <iostream>
#include <ctime>//for time cmd
#include <math.h>

#include <gsl/gsl_randist.h>

using namespace std;

class HmmState{
private:
    long double h;
    long double theta;
    long double rateGeom;

    HmmState * second;
    gsl_rng * rng;
    int idx;
    
public:
    HmmState(int idx_,long double p_);
    HmmState(const HmmState & other);
    ~HmmState();
    HmmState & operator= (const HmmState & other);
    
    /* HmmState & operator= (const HmmState & other){ */
    /* 	p      = other.p; */
    /* 	second = other.second; */
	
    /* } */

    double probEmission(unsigned int mutations,unsigned int total) const;
    unsigned int randomEmission(int total) const;
    int getIdx() const;
    long double getH();
    long double getTheta();
    long double getRateGeom();

    void setH(long double newH);

    void setSecond(HmmState * second);
};
#endif
