#include "HmmState.h"


HmmState::HmmState(int idx_,long double h_){
    h        = h_;
    theta    = h/(1-h); //h = theta/theta+1
    rateGeom = 1/(theta+1); //rate of geom

    rng  = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(rng, time(NULL));
    idx = idx_;
}

HmmState::HmmState(const HmmState & other){
    h        = other.h;
    theta    = other.theta;
    rateGeom = other.rateGeom;
    idx      = other.idx;

    rng  = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(rng, time(NULL));
}

HmmState::~HmmState(){
    gsl_rng_free(rng);

}

void HmmState::setSecond(HmmState * second_){
    second = second_;
}

long double HmmState::probEmission(unsigned int mutations,unsigned int total) const{

    // double noMut=gsl_ran_geometric_pdf  (1, rateGeom);
    // cout<<"no mut="<<noMut<<" "<<(1-noMut)<<" "<<rateGeom<<endl;
    // if(mutations == 0){
    // 	return ( pow(noMut,(total)) );
    // }else{
    // 	return ( pow((1-noMut),mutations) + pow(noMut,(total-mutations)) );
    // }
    return gsl_ran_binomial_pdf(mutations,h,total);
}

unsigned int HmmState::randomEmission(int total) const{


    //geom
    unsigned int toreturn;
    unsigned int mut=0;
    for(int i=0;i<total;i++){
	mut+=gsl_ran_geometric(rng,rateGeom);
    }
    //cout<<(mut-total)<<endl;
    toreturn = (mut-total);
    //unsigned int toreturn =  gsl_ran_binomial(rng,p,total);
    //unsigned int toreturn =  gsl_ran_negative_binomial(rng,1-p,total);
    //cout<<"randomEmission() "<<total<<"\tp="<<p<<"\t"<<toreturn<<endl;
    return toreturn;
}

int HmmState::getIdx() const{
    return idx;
}

// HmmState::HmmState(const HmmState & other){
//     p      = other.p;
//     second = other.second; 
// }


long double HmmState::getH(){
    return h;
}

long double HmmState::getTheta(){
    return theta;
}

long double HmmState::getRateGeom(){
    return rateGeom;
}

void HmmState::setH(long double newH){
    h        = newH;
    theta    = h/(1-h); //h = theta/theta+1
    rateGeom = 1/(theta+1); //rate of geom
}
