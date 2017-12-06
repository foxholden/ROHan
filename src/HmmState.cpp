#include "HmmState.h"


HmmState::HmmState(int idx_,long double p_){
    p = p_;
    rng  = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(rng, time(NULL));
    idx = idx_;
}


HmmState::~HmmState(){
    gsl_rng_free(rng);

}

void HmmState::setSecond(HmmState * second_){
    second = second_;
}

double HmmState::probEmission(unsigned int mutations,unsigned int total){
    return gsl_ran_binomial_pdf(mutations,p,total);
}

unsigned int HmmState::randomEmission(int total){


    //geom
    double rate = 1/(p+1);
    //cout<<"geom "<<p<<" "<<rate<<endl;
    unsigned int toreturn;
    unsigned int mut=0;
    for(int i=0;i<total;i++){
	mut+=gsl_ran_geometric(rng,rate);
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
