/*
 * Hmm
 * Date: Dec-05-2017 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "Hmm.h"


Hmm::Hmm(){
    //cerr<<"Hmm constr"<<endl;
    for(int n=0;n<NBRSTATES;n++){
	if(n==0)
	    hmmstates[n] = new  HmmState (n,0.000000012);
	if(n==1)
	    hmmstates[n] = new  HmmState (n,0.00072);	
    }

    for(int n=0;n<NBRSTATES;n++){
	if(n==0)
	    hmmstates[n]->setSecond(hmmstates[n+1]);
	if(n==1)
	    hmmstates[n]->setSecond(hmmstates[n-1]);
    }
    // for(int m=0;m<1000;m++)
    // 	cout<<m<<"\tp="<<hmmstates[1]->probEmission(m,1000000)<<endl;
    double pTrans = 0.1;
    trans = new double*[NBRSTATES];
    
    for(int i = 0; i < NBRSTATES; ++i)
	trans[i] = new double[NBRSTATES];
    
    for(int n=0;n<NBRSTATES;n++){
	probTrans[n] =   pTrans;
	probStay[n]  = 1-pTrans;

	trans[n][0]  =   pTrans;
	trans[n][1]  = 1-pTrans;	
    }
    
    for(int n=0;n<NBRSTATES;n++){
	startingState[n] = 1/double(NBRSTATES);
    }
    // HmmState roh    (0.0007);
    // HmmState normal (0.000000012);
    // roh.setSecond(&normal);
    // normal.setSecond(&roh);

    // hmmstates[0]=roh;
    // hmmstates[1]=normal;

    // double pTrans = 0.1;
    // probTrans[0]=pTrans;
    // probTrans[1]=pTrans;

    // probStay[0]=1-pTrans;
    // probStay[1]=1-pTrans;
    
    //    cout<<"test"<<endl;
    
}

Hmm::~Hmm(){

    for(int n=0;n<NBRSTATES;n++){
	delete hmmstates[n];
    }

}


vector<emission> Hmm::generateStates(unsigned int N,unsigned int total) const{
    vector<emission> vecToReturn;
    //int initState = randomInt(0,NBRSTATES-1);
    int initState = 1;
    HmmState * currrentState = hmmstates[ initState ];
    
    for(unsigned int i=0;i<N;i++){
	//cerr<<"generateStates "<<i<<endl;
	unsigned int d_ = currrentState->randomEmission(total);
	emission eToAdd;
	eToAdd.p     = double(d_)/double(total);
	eToAdd.idx   = currrentState->getIdx();
	
	vecToReturn.push_back( eToAdd );
	double pt = randomProb();
	//cout<<"state#"<<currrentState->getIdx()<<"\t"<<d<<"\t"<<pt<<"\t"<<probTrans[currrentState->getIdx()]<<endl;	

	//if(0)
	if(pt < probTrans[currrentState->getIdx()]){//transit
	    if(currrentState->getIdx() == 0)
		currrentState = hmmstates[1];
	    else
		if(currrentState->getIdx()     == 1)
		    currrentState  = hmmstates[0];
	}else{
	    //do nothing
	}
    }
    return vecToReturn;
}

int  Hmm::getNumberStates() const{
    return NBRSTATES;
}


double Hmm::getTrans(int i,int j) const{
    return trans[i][j];
}
