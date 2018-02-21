#include "HmmState.h"


HmmState::HmmState(int idx_,long double h_,int minSegSitesPerChunk_,int maxSegSitesPerChunk_,int sizeChunk_){
    h        = h_;
    theta    = h/(1-h); //h = theta/theta+1
    rateGeom = 1/(theta+1); //rate of geom

    minSegSitesPerChunk = minSegSitesPerChunk_;
    maxSegSitesPerChunk = maxSegSitesPerChunk_;
    sizeChunk           = sizeChunk_;

    rng  = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(rng, time(NULL));
    idx = idx_;
    probabilitiesForEmission = new vector<long double> ();

    for(unsigned int mutations=0;mutations<=( (unsigned int)maxSegSitesPerChunk);mutations++){ //unsigned int total) const{
	long double probMut = gsl_ran_binomial_pdf(mutations,h,sizeChunk);
	//cout<<mutations<<"\t"<<probMut<<endl;
	probabilitiesForEmission->push_back( probMut );
    }
    // double noMut=gsl_ran_geometric_pdf  (1, rateGeom);
    // cout<<"no mut="<<noMut<<" "<<(1-noMut)<<" "<<rateGeom<<endl;
    // if(mutations == 0){
    // 	return ( pow(noMut,(total)) );
    // }else{
    // 	return ( pow((1-noMut),mutations) + pow(noMut,(total-mutations)) );
    // }

}
    


HmmState::HmmState(const HmmState & other){
    h                   = other.h;
    theta               = other.theta;
    rateGeom            = other.rateGeom;
    idx                 = other.idx;

    minSegSitesPerChunk = other.minSegSitesPerChunk;
    maxSegSitesPerChunk = other.maxSegSitesPerChunk;
    sizeChunk           = other.sizeChunk;

    delete probabilitiesForEmission;
    probabilitiesForEmission = new vector<long double> (maxSegSitesPerChunk+1,0);
    
    for(unsigned int mutations=0;mutations<=( (unsigned int)maxSegSitesPerChunk);mutations++){ //unsigned int total) const{
	long double probMut = gsl_ran_binomial_pdf(mutations,h,sizeChunk);
	probabilitiesForEmission->at(mutations) = probMut ;
    }

    rng  = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set(rng, time(NULL));
}

HmmState::~HmmState(){
    gsl_rng_free(rng);
    delete probabilitiesForEmission;
}

void HmmState::setSecond(HmmState * second_){
    second = second_;
}

long double HmmState::probEmission(     const int mutations,                          const int total) const{//probEmission(unsigned int mutations,unsigned int total) const{

    // double noMut=gsl_ran_geometric_pdf  (1, rateGeom);
    // cout<<"no mut="<<noMut<<" "<<(1-noMut)<<" "<<rateGeom<<endl;
    // if(mutations == 0){
    // 	return ( pow(noMut,(total)) );
    // }else{
    // 	return ( pow((1-noMut),mutations) + pow(noMut,(total-mutations)) );
    // }
    //long double normalVal =   gsl_ran_binomial_pdf(mutations,h,total);

    if(int(total) != sizeChunk){
	cerr<<"HmmState::probEmission("<<mutations<<","<<total<<") is different than the pre-computed size of chunk "<<sizeChunk<<endl;
	exit(1);
    }
    

    if(mutations<0){
	cerr<<"HmmState::probEmission("<<mutations<<","<<total<<") is lower than the defined minimum "<<minSegSitesPerChunk<<endl;
	exit(1);	
    }

    if(int(mutations)>maxSegSitesPerChunk){
	cerr<<"HmmState::probEmission("<<mutations<<","<<total<<") is higher than the defined minimum "<<maxSegSitesPerChunk<<endl;
	exit(1);	
    }

    
    long double toReturn = probabilitiesForEmission->at( mutations );
    // cerr<<h<<"\t"<<mutations<<"\t"<<normalVal<<"\t"<<toReturn<<endl;
    // if(normalVal != toReturn){
    // 	exit(1);
    // }
    return toReturn;
}


//if we have a range to marginalize over using a uniform prior
long double HmmState::probEmissionRange(const int mutationsMin,const int mutationsMax,const int total) const{

    int mutationsMax_ = mutationsMax;
    //cerr<<"probEmissionRange "<<mutationsMin<<" "<<mutationsMax<<endl;
    // exit(1);
    if(int(total) != sizeChunk){
	cerr<<"HmmState::probEmissionRange("<<mutationsMin<<","<<mutationsMax_<<","<<total<<") is different than the pre-computed size of chunk "<<sizeChunk<<endl;
	exit(1);
    }
    
    if(mutationsMin>mutationsMax_){
	cerr<<"HmmState::probEmissionRange("<<mutationsMin<<","<<mutationsMax_<<","<<total<<") is different than the pre-computed size of chunk "<<sizeChunk<<endl;
	exit(1);
    }

    
    // if(mutationsMin<0){
    // 	cerr<<"HmmState::probEmissionRange("<<mutationsMin<<","<<mutationsMax_<<","<<total<<") minimum is lower than the defined minimum "<<minSegSitesPerChunk<<endl;
    // 	exit(1);
    // }

    if(mutationsMin == mutationsMax_){
	mutationsMax_ = mutationsMin+1;
    }

    // if(int(mutationsMax_) > maxSegSitesPerChunk ){
    // 	cerr<<"HmmState::probEmissionRange("<<mutationsMin<<","<<mutationsMax_<<","<<total<<") is higher than the defined minimum "<<maxSegSitesPerChunk<<endl;
    // 	exit(1);	
    // }
    long double sumPFE=0.0;
    for(int m=mutationsMin;m<=mutationsMax_;m++){
        //int m_  = m;
        if(m < 0){
            sumPFE += probabilitiesForEmission->at( 0 )*(-1.0*m);
            m  = 0;
            continue;
        }

        if(m > maxSegSitesPerChunk){
            //m = maxSegSitesPerChunk;//went over the limit, set to max
            sumPFE += probabilitiesForEmission->at( maxSegSitesPerChunk )*(m-maxSegSitesPerChunk);
            break;
        }
        sumPFE += probabilitiesForEmission->at( m );
    }
    
    //if( (sumPFE / (mutationsMax-mutationsMin) ) == 1){
    //cerr<<endl<<"probabilitiesForEmission("<<h<<","<<idx<<") sumPFE = "<<sumPFE<<" "<<mutationsMin<<" "<<mutationsMax<<" "<<(mutationsMax-mutationsMin)<<" return="<<(sumPFE / (mutationsMax-mutationsMin) )<<endl;
    // 	exit(1);
    // }
    //to test
    return (sumPFE / (mutationsMax_-mutationsMin) );
    //return probabilitiesForEmission->at( (mutationsMax-mutationsMin)/2);
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

    //need to recompute the probabilitiesForEmission
    delete probabilitiesForEmission;
    probabilitiesForEmission = new vector<long double> (maxSegSitesPerChunk+1,0.0);
    
    for(unsigned int mutations=0;
	mutations<=( (unsigned int)maxSegSitesPerChunk);
	mutations++){ //unsigned int total) const{
	long double probMut = gsl_ran_binomial_pdf(mutations,h,sizeChunk);
	probabilitiesForEmission->at(mutations) = probMut;
    }

}
