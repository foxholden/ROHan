/*
 * DataToWrite
 * Date: Mar-22-2016 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "DataToWrite.h"



DataToWrite::DataToWrite(){
    vecPositionResults =  new vector<PositionResult *>();
    //cerr<<"Constructor addr: "<<this<<endl;
}

DataToWrite::~DataToWrite(){
    for (unsigned int i =0; i< vecPositionResults->size();i++){
	delete (vecPositionResults->at(i));
    } 
    vecPositionResults->clear();
    delete vecPositionResults;
    //cerr<<"Destructor  addr: "<<this<<endl;
}


