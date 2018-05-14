/*
 * miscfunc
 * Date: Jun-08-2014 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */


#include "miscfunc.h"





void readIlluminaError(const string errFile,probSubstition & illuminaErrorsProb){

    igzstream errFileSt;

    errFileSt.open(errFile.c_str(), ios::in);

    //    unsigned int counterCont=0;
    if (errFileSt.good()){
	vector<string> fields;
	string line;
	//header
	if ( !getline (errFileSt,line)){
	    cerr << "Unable to open file "<<errFile<<endl;
	    exit(1);
	}
	fields = allTokens(line,'\t');
	
	if(fields.size() != 12){
	    cerr << "line from error profile does not have 12 fields "<<line<<endl;
	    exit(1);
	}

	//raw sums
	if ( !getline (errFileSt,line)){
	    cerr << "Unable to open file "<<errFile<<endl;
	    exit(1);
	}
	
	fields = allTokens(line,'\t');

	if(fields.size() != 12){
	    cerr << "line from error profile does not have 12 fields "<<line<<endl;
	    exit(1);
	}

	//probs
	if ( !getline (errFileSt,line)){
	    cerr << "Unable to open file "<<errFile<<endl;
	    exit(1);
	}
	
	fields = allTokens(line,'\t');

	if(fields.size() != 12){
	    cerr << "line from error profile does not have 12 fields "<<line<<endl;
	    exit(1);
	}
	substitutionRates tempFreq;	
	    
	
	for(unsigned int k=0;k<=9;k+=3){	

	    for(unsigned int t=0;t<=2;t++){	
		tempFreq.s[k+t]=destringify<long double>(fields[k+t]);
		//cerr<<freqIlluminaError.s[k+t]<<endl;
	    }

	}


	int indexFirstArray =0;
	int indexSecondArray=0;

	for(int nuc1=0;nuc1<4;nuc1++){
	    for(int nuc2=0;nuc2<4;nuc2++){
		if(nuc1==nuc2) // prob of error is 0 if both nucleotides are identical
		    illuminaErrorsProb.s[indexFirstArray++]=0.0;
		else //           rely on the substitution frequency
		    illuminaErrorsProb.s[indexFirstArray++]=tempFreq.s[indexSecondArray++];
	    }
	}
	
	             	              
	errFileSt.close();
    }else{
	cerr << "Unable to open file "<<errFile<<endl;
	exit(1);
    }



}

void readNucSubstitionFreq(const string filename,vector<probSubstition> & subVec){
    igzstream subFP;

    subFP.open(filename.c_str(), ios::in);

    //    unsigned int counterCont=0;
    if (subFP.good()){
	vector<string> fields;
	string line;

	//header
	if ( !getline (subFP,line)){
	    cerr << "Unable to open file "<<filename<<endl;
	    exit(1);
	}
	fields = allTokens(line,'\t');
	
	if(fields.size() != 12){
	    cerr << "line from error profile does not have 12 fields "<<line<<endl;
	    exit(1);
	}


	//probs
	while ( getline (subFP,line)){
	    
	    fields = allTokens(line,'\t');

	    
	    if(fields.size() != 12){
		cerr << "line from error profile does not have 12 fields "<<line<<endl;
		exit(1);
	    }

	    substitutionRates tempFreq;	
	    probSubstition toaddSub;


	    for(unsigned int k=0;k<=9;k+=3){	

		for(unsigned int t=0;t<=2;t++){	
		    tempFreq.s[k+t]=destringify<long double>(fields[k+t]);
		}

	    }


	    int indexFirstArray =0;
	    int indexSecondArray=0;

	    for(int nuc1=0;nuc1<4;nuc1++){
		long double sumMismatchProb=0.0;
		int indexInArrayMatch=1;
		for(int nuc2=0;nuc2<4;nuc2++){
		    if(nuc1==nuc2){ // prob of error is 0 if both nucleotides are identical
			indexInArrayMatch                       = indexFirstArray;
			toaddSub.s[indexFirstArray++]           = 1.0;		    
		    }else{ //           rely on the substitution frequency
			sumMismatchProb                         += tempFreq.s[indexSecondArray];
			toaddSub.s[indexFirstArray++]            = tempFreq.s[indexSecondArray++];
		    }
		}

		toaddSub.s[indexInArrayMatch] = 1.0 - sumMismatchProb;
	    }

	    // for(int nuc1=0;nuc1<4;nuc1++){
	    // 	for(int nuc2=0;nuc2<4;nuc2++){
	    // 	    cout<<(nuc1*4+nuc2)<<"\t"<<toaddSub.s[nuc1*4+nuc2]<<endl;
	    // 	}	       
	    // }
	    
	    // exit(1);

	    subVec.push_back( toaddSub );
	}	             	              
	subFP.close();
    }else{
	cerr << "Unable to open file "<<filename<<endl;
	exit(1);
    }



}



void readNucSubstitionRatesFreq(const string filename,vector<substitutionRates> & subVec){
    igzstream subFP;

    subFP.open(filename.c_str(), ios::in);

    //    unsigned int counterCont=0;
    if (subFP.good()){
	vector<string> fields;
	string line;

	//header
	if ( !getline (subFP,line)){
	    cerr << "Unable to open file "<<filename<<endl;
	    exit(1);
	}
	fields = allTokens(line,'\t');
	
	if(fields.size() != 12){
	    cerr << "line from error profile does not have 12 fields "<<line<<endl;
	    exit(1);
	}


	//probs
	while ( getline (subFP,line)){
	    
	    fields = allTokens(line,'\t');

	    if(fields.size() != 12){
		cerr << "line from error profile does not have 12 fields "<<line<<endl;
		exit(1);
	    }

	    substitutionRates tempFreq;	


	    for(unsigned int k=0;k<12;k++){	
		//for(unsigned int t=0;t<=2;t++){	
		tempFreq.s[k]=destringify<long double>(fields[k]);
		//}
	    }



	    subVec.push_back( tempFreq );
	}	             	              
	subFP.close();
    }else{
	cerr << "Unable to open file "<<filename<<endl;
	exit(1);
    }



}

void readDNABaseFreq(const string filename, alleleFrequency & dnaDefaultFreq){
    
    igzstream           dnaProfFP;
    vector<long double> fields;
    
    dnaProfFP.open(filename.c_str(), ios::in);
    string line;
    if (dnaProfFP.good()){
	//probs

	while ( getline (dnaProfFP,line)){	   
	    long double f = destringify<long double>(line);	    
	    fields.push_back(f);
	}	             	              
	dnaProfFP.close();
    }else{
	cerr << "Unable to open file "<<filename<<endl;
	exit(1);
    }

    if(fields.size() != 4){
	cerr<<"File prof "<<filename<<" contains "<<fields.size()<<"lines, should be 4"<<endl;
	exit(1);
    }

    for(unsigned int i=0;i<4;i++)
	dnaDefaultFreq.f[i] = fields[i];


}


void populatedCoverateVector(const string  programName, vector<long double> * cov2ProbSite, long double rateForPoissonCov,int maxcov){
    //cerr<<"Computing coverage priors.";
    long double rateForPoissonCovFloor = floorl(rateForPoissonCov);
    long double rateForPoissonCovCeil  =  ceill(rateForPoissonCov);

    vector<long double>  cov2ProbSiteFloor;
    vector<long double>  cov2ProbSiteCeil;

    rateForPoissonCovFloor = MAX2(rateForPoissonCovFloor, (long double)1);
    rateForPoissonCovFloor = MIN2(rateForPoissonCovFloor, (long double)maxcov);
    rateForPoissonCovCeil  = MIN2(rateForPoissonCovCeil,  (long double)maxcov);
    string directoryProgram;
    string commandPath=string(programName);
    size_t posSlash=commandPath.find_last_of("/");
    if(posSlash == string::npos){
        directoryProgram="";
    }else{
        directoryProgram=commandPath.substr(0,posSlash);
    }
    directoryProgram = directoryProgram +"/";
    populatedCoverateVectorSingle(directoryProgram, &cov2ProbSiteFloor, rateForPoissonCovFloor ,  maxcov);
    populatedCoverateVectorSingle(directoryProgram, &cov2ProbSiteCeil , rateForPoissonCovCeil  ,  maxcov);

    //cov2ProbSite->push_back(0.0);//cov =0

    for(int cov=1;cov<=maxcov;cov++){
	
	long double florPPMF = cov2ProbSiteFloor[cov-1];
	long double ceilPPMF = cov2ProbSiteCeil[cov-1];
	if(rateForPoissonCov<=1){
	    rateForPoissonCovFloor = 0.0;
	    rateForPoissonCovCeil  = 1.0;
	    rateForPoissonCov      = 1.0;
	}
#ifdef DEBUGCOV
	cerr<<cov<<"\t"<<rateForPoissonCov<<"\t"<<(1.0-(rateForPoissonCov-rateForPoissonCovFloor))<<"\t"<<florPPMF<<"\t"<<(1.0-(rateForPoissonCovCeil-rateForPoissonCov))<<"\t"<<ceilPPMF<<"\t"<<(	    (1.0-(rateForPoissonCov-rateForPoissonCovFloor))*florPPMF
	    +
	    (1.0-(rateForPoissonCovCeil-rateForPoissonCov))*ceilPPMF
	)<<endl;
#endif
	cov2ProbSite->at(cov) =
	    (
		(1.0-(rateForPoissonCov-rateForPoissonCovFloor))*florPPMF
		+
		(1.0-(rateForPoissonCovCeil-rateForPoissonCov))*ceilPPMF
	    );

    }
    
#ifdef DEBUGCOV    
    cerr<<"..done"<<endl;
#endif

}

void populatedCoverateVectorSingle(const string  directoryProgram,vector<long double> * cov2ProbSite, long double lambda,int maxcov){
   
    
    string filein  = string(directoryProgram+"../preComputated/coverageprior/correctionCov_"+stringify(lambda)+".bin");

#ifdef DEBUGCOV
    cerr<<"file in "<<filein<<endl;
#endif
    if(!isFile(filein)){
	cerr<<"ERROR: Cannot find file "<<filein<<"  containing pre-computed coverage weights"<<endl;
	exit(1);
    }
    
    ifstream in(filein.c_str(), ios::in | ios::binary);
    if(!in.good()){
	cerr<<"ERROR: Opening file "<<filein<<" failed containing pre-computed coverage weights"<<endl;
	exit(1);
    }

    //TODO what to do about cross platform?
    // uint8_t sizeld;
    // in.read((char *) &sizeld, sizeof(sizeld) );
    //cerr<<"size long double "<<int(sizeld)<<endl;
    for(int cov=1;cov<=maxcov;cov++){
        long double fnum;
	//cerr<<sizeof(fnum)<<endl;
        in.read((char *) &fnum, sizeof(fnum) );

#ifdef DEBUGCOV
	cerr<<cov<<" cov = "<<fnum<<endl;
#endif

	cov2ProbSite->push_back(fnum);
        //cout<<setprecision(5)<<fnum<<endl;
    }
    in.close();

    //cout<<"populatedCoverateVectorSingle "<<lambda<<" "<<maxcov<<endl;
    //     cov2ProbSite->push_back(0.0);//cov =0
    //     long double   p3_ [3] = {double(1.0)/double(3.0),double(1.0)/double(3.0),double(1.0)/double(3.0)}; 
    //     long double   p4_ [4] = {double(1.0)/double(4.0),double(1.0)/double(4.0),double(1.0)/double(4.0),double(1.0)/double(4.0)}; 
    //     long double   p5_ [5] = {double(1.0)/double(5.0),double(1.0)/double(5.0),double(1.0)/double(5.0),double(1.0)/double(5.0),double(1.0)/double(5.0)}; 
    //     // long double   p6_ [6] = {double(1.0)/double(6.0),double(1.0)/double(6.0),double(1.0)/double(6.0),double(1.0)/double(6.0),double(1.0)/double(6.0),double(1.0)/double(6.0)}; 
    //     // long double   p7_ [7] = {double(1.0)/double(7.0),double(1.0)/double(7.0),double(1.0)/double(7.0),double(1.0)/double(7.0),double(1.0)/double(7.0),double(1.0)/double(7.0),double(1.0)/double(7.0)}; 
    //     // long double   p8_ [8] = {double(1.0)/double(8.0),double(1.0)/double(8.0),double(1.0)/double(8.0),double(1.0)/double(8.0),double(1.0)/double(8.0),double(1.0)/double(8.0),double(1.0)/double(8.0),double(1.0)/double(8.0)}; 
    
    //     for(int cov=1;cov<=maxcov;cov++){
    // 	cerr<<".";
    // #ifdef DEBUGCOV
    // 	cout<<"populatedCoverateVectorSingle cov = "<<cov<<endl;
    // #endif
    // 	//P[cov|no dup]
    // 	double pNoDup=poisson_pmfl(double(cov),double(lambda));
    
    // #ifdef DEBUG2
    // 	    cout<<"P[cov|no dup] "<<pNoDup<<endl;
    // 	    cout<<"----------------------------------"<<endl;
    // #endif
    // 	    //}
    // 	//P[cov|with collapse in reference]
    // 	//This means that there are 2 locations in the original genome
    // 	//that contribute to the same genomic location
    // 	double pCollapseInReference=0;
    // 	for(int cov1=0;cov1<=cov;cov1++){//coverage of 1st location
    
    // 	    int cov2=(cov-cov1); //coverage of second location
    // 	    double p1 = poisson_pmf(double(cov1),double(lambda));
    // 	    double p2 = poisson_pmf(double(cov2),double(lambda));
    // 	    double p1_2=  p1*p2;
    // 	    pCollapseInReference += p1_2;
    // #ifdef DEBUG2
    // 	    cout<<cov1<<"\t"<<cov2<<"\t"<<p1<<"\t"<<p2<<"\t"<<p1_2<<"\t"<<pCollapseInReference<<endl;
    // #endif	    
    // 	}
    // #ifdef DEBUG2
    // 	cout<<"P[cov|with collapse in reference] "<<pCollapseInReference<<endl;
    // 	cout<<"----------------------------------"<<endl;
    // #endif
    // 	//P[cov|with dup in reference]
    // 	//there are 2 genomic copies in the reference thus leading to "half" the coverage from the original location
    // 	//
    // 	double pDupInReference2=0;
    // 	for(int covOrig=0;covOrig<=100;covOrig++){//coverage in the original genome
    // 	    //2 copies
    // 	    double pOrig = poisson_pmf(double(covOrig),double(lambda));
    
    // #ifdef DEBUG2
    // 	    cout<<covOrig<<"\t"<<pOrig<<endl;
    // #endif

    // 	    for(int cov1=0;cov1<=covOrig;cov1++){//coverage of 1st location	
    // 		int cov2=(covOrig-cov1); //coverage of second location
    // 		if(cov1!=cov && cov2!=cov)
    // 		    continue;
    
    // 		double probBinom = 	nChoosek(covOrig,cov1) * pow(0.5,covOrig);
    // #ifdef DEBUG2
    // 		cout<<cov1<<"\t"<<cov2<<"\t"<<probBinom<<endl;
    // #endif
    // 		pDupInReference2+=probBinom*pOrig;
    // 	    }
    // 	}
    
    
    // 	//3 copies
    // 	long double pDupInReference3=0.0;
    // 	for(int covOrig=0;covOrig<=100;covOrig++){//coverage in the original genome
    // 	    //2 copies
    // 	    long double pOrig = poisson_pmfl((long double)covOrig,(long double)lambda);
    
    // #ifdef  DEBUG3CP
    // 	    cout<<"covOrig\t"<<covOrig<<"\tpOrig\t"<<pOrig<<"\t"<<pDupInReference3<<endl;
    // #endif
    // 	    long double pDupInReference3_=0.0;
    // 	    for(int cov1=0;cov1<=covOrig;cov1++){//coverage of 1st location	
    
    // 		for(int cov2=0;cov2<=(covOrig-cov1);cov2++){//coverage of 2nd location	
    // 		    int cov3=(covOrig-cov1-cov2); //coverage of 3rd location
    
    // 		    // 		int cov2=(covOrig-cov1); //coverage of second location
    // 		    if(cov1!=cov && cov2!=cov && cov3!=cov)
    // 			continue;
    
    
    // 		    unsigned int  n_ [3] = {(unsigned int)cov1,(unsigned int)cov2,(unsigned int)cov3};   
    // 		    double probMnom = 	 gsl_ran_multinomial_pdf(3,p3_,n_);
    
    // 		    if(isinf(probMnom))
    // 			continue;
    // #ifdef DEBUG3CP
    // 		    cout<<"covOrig\t"<<covOrig<<"\tc1\t"<<cov1<<"\tc2\t"<<cov2<<"\tc3\t"<<cov3<<"\t"<<probMnom<<"\t"<<probMnom*pOrig<<endl;
    // #endif
    
    // 		    pDupInReference3_ += probMnom*pOrig;
    
    // 		    // #ifdef DEBUG2
    // 		    // 		cout<<cov1<<"\t"<<cov2<<"\t"<<probBinom<<endl;
    // 		    // #endif
    // 		    // 		pDupInReference+=probBinom*pOrig;
    // 		}
    // 	    }
    // 	    pDupInReference3 += pDupInReference3_;
    // 	    if(pDupInReference3_>0 && pDupInReference3_< 1.0e-10){
    // 		break;
    // 	    }
    // 	}//end cov orig
    

    
    // 	//4 copies
    // 	long double pDupInReference4=0.0;
    // 	for(int covOrig=0;covOrig<=100;covOrig++){//coverage in the original genome
    // 	    //2 copies
    // 	    long double pOrig = poisson_pmfl((long double)covOrig,(long double)lambda);
    
    // #ifdef  DEBUG4CP
    // 	    cout<<"covOrig\t"<<covOrig<<"\tpOrig\t"<<pOrig<<"\t"<<pDupInReference4<<endl;
    // #endif
    // 	    long double pDupInReference4_=0.0;
    // 	    for(int cov1=0;cov1<=covOrig;cov1++){//coverage of 1st location	
    
    // 		for(int cov2=0;cov2<=(covOrig-cov1);cov2++){//coverage of 2nd location	
    // 		    for(int cov3=0;cov3<=(covOrig-cov1-cov2);cov3++){//coverage of 2nd location	
    // 			int cov4=(covOrig-cov1-cov2-cov3); //coverage of 3rd location
    // 			//cout<<"covOrig\t"<<covOrig<<"\tc1\t"<<cov1<<"\tc2\t"<<cov2<<"\tc3\t"<<cov3<<"\tc4\t"<<cov4<<endl;
    // 			// 		int cov2=(covOrig-cov1); //coverage of second location
    // 			if(cov1!=cov && cov2!=cov && cov3!=cov && cov4!=cov)
    // 			    continue;
    

    // 			unsigned int  n_ [4] = {(unsigned int)cov1,(unsigned int)cov2,(unsigned int)cov3,(unsigned int)cov4};   
    // 			double probMnom = 	 gsl_ran_multinomial_pdf(4,p4_,n_);
    		   
    // 			if(isinf(probMnom))
    // 			    continue;
    // #ifdef DEBUG4CP
    // 			cout<<"covOrig\t"<<covOrig<<"\tc1\t"<<cov1<<"\tc2\t"<<cov2<<"\tc3\t"<<cov3<<"\tc4\t"<<cov4<<"\t"<<probMnom<<"\t"<<probMnom*pOrig<<endl;
    // #endif
		    
    // 			pDupInReference4_ += probMnom*pOrig;
		
    // 			// #ifdef DEBUG2
    // 			// 		cout<<cov1<<"\t"<<cov2<<"\t"<<probBinom<<endl;
    // 			// #endif
    // 			// 		pDupInReference+=probBinom*pOrig;
    // 		    }
    // 		}
    // 	    }
    // 	    pDupInReference4 += pDupInReference4_;
    // 	    if(pDupInReference4_>0 && pDupInReference4_< 1.0e-10){
    // 		break;
    // 	    }

    // 	}//end cov orig

    // 	//5 copies
    // 	long double pDupInReference5=0.0;
    // 	for(int covOrig=0;covOrig<=100;covOrig++){//coverage in the original genome
    // 	    //2 copies
    // 	    long double pOrig = poisson_pmfl((long double)covOrig,(long double)lambda);

    // #ifdef  DEBUG5CP
    // 	    cout<<"covOrig\t"<<covOrig<<"\tpOrig\t"<<pOrig<<"\t"<<pDupInReference5<<endl;
    // #endif
    // 	    long double pDupInReference5_=0.0;
    // 	    for(int cov1=0;cov1<=covOrig;cov1++){//coverage of 1st location	

    // 		for(int cov2=0;cov2<=(covOrig-cov1);cov2++){//coverage of 2nd location	
    // 		    for(int cov3=0;cov3<=(covOrig-cov1-cov2);cov3++){//coverage of 3rd location	
    // 		    for(int cov4=0;cov4<=(covOrig-cov1-cov2-cov3);cov4++){//coverage of 4th location	
    // 			int cov5=(covOrig-cov1-cov2-cov3-cov4); //coverage of 5th location
    // 			//cout<<"covOrig\t"<<covOrig<<"\tc1\t"<<cov1<<"\tc2\t"<<cov2<<"\tc3\t"<<cov3<<"\tc4\t"<<cov4<<endl;
    // 			// 		int cov2=(covOrig-cov1); //coverage of second location
    // 			if(cov1!=cov && cov2!=cov && cov3!=cov && cov4!=cov && cov5!=cov)
    // 			    continue;
		

    // 			unsigned int  n_ [5] = {(unsigned int)cov1,(unsigned int)cov2,(unsigned int)cov3,(unsigned int)cov4,(unsigned int)cov5};   
    // 			double probMnom = 	 gsl_ran_multinomial_pdf(5,p5_,n_);
    		   
    // 			if(isinf(probMnom))
    // 			    continue;
    // #ifdef DEBUG5CP
    // 			cout<<"covOrig\t"<<covOrig<<"\tc1\t"<<cov1<<"\tc2\t"<<cov2<<"\tc3\t"<<cov3<<"\tc4\t"<<cov4<<"\tc5\t"<<cov5<<"\t"<<probMnom<<"\t"<<probMnom*pOrig<<endl;
    // #endif
		    
    // 			pDupInReference5_ += probMnom*pOrig;
		
    // 		    }
    // 		}
    // 	    }
    // 	    }
    // 	    pDupInReference5 += pDupInReference5_;
    // 	    if(pDupInReference5_>0 && pDupInReference5_< 1.0e-10){
    // 		break;
    // 	    }

    // 	}//end cov orig



    // 	//6 copies
    // 	long double pDupInReference6=0.0;
    // 	for(int covOrig=0;covOrig<=100;covOrig++){//coverage in the original genome
    // 	    //2 copies
    // 	    long double pOrig = poisson_pmfl((long double)covOrig,(long double)lambda);

	    
    // #ifdef  DEBUG6CP
    // 	    cout<<"6cp covOrig\t"<<covOrig<<"\tpOrig\t"<<pOrig<<"\t"<<pDupInReference6<<endl;
    // #endif

    // 	    long double pDupInReference6_=0.0;

    // 	    string filenameIn = "../preComputated/combinatorialCov/comb6_"+stringify(cov)+"_"+stringify(covOrig)+".bin";
    // 	    ifstream in(filenameIn.c_str(), ios::in | ios::binary);	    
    // 	    if(!in.good()){
    // 		cerr<<"Cannot open "<<filenameIn<<endl;
    // 		exit(1);
    // 	    }
    // 	    in.read((char *) &pDupInReference6_, sizeof(pDupInReference6_) );
    // 	    in.close();
    // 	    pDupInReference6 += pDupInReference6_*pOrig;

    // #ifdef  DEBUG6CP
    // 	    cout<<"6cp covOrig\t"<<covOrig<<"\tpOrig\t"<<pOrig<<"\t"<<pDupInReference6_<<"\t"<<pDupInReference6<<endl;
    // #endif


    // 	}//end cov orig

    // 	//return 1;


    // 	//7 copies
    // 	long double pDupInReference7=0.0;
    // 	for(int covOrig=0;covOrig<=100;covOrig++){//coverage in the original genome
    // 	    //2 copies
    // 	    long double pOrig = poisson_pmfl((long double)covOrig,(long double)lambda);

    // #ifdef  DEBUG7CP
    // 	    cout<<"7cp covOrig\t"<<covOrig<<"\tpOrig\t"<<pOrig<<"\t"<<pDupInReference7<<endl;
    // #endif
	    
    // 	    long double pDupInReference7_=0.0;
    // 	    string filenameIn = "../preComputated/combinatorialCov/comb7_"+stringify(cov)+"_"+stringify(covOrig)+".bin";
    // 	    ifstream in(filenameIn.c_str(), ios::in | ios::binary);	    
    // 	    if(!in.good()){
    // 		cerr<<"Cannot open "<<filenameIn<<endl;
    // 		exit(1);
    // 	    }
    // 	    in.read((char *) &pDupInReference7_, sizeof(pDupInReference7_) );
    // 	    in.close();
	    

    // 	    pDupInReference7 += pDupInReference7_*pOrig;

    // 	}//end cov orig

    // 	//8 copies 
    // 	long double pDupInReference8=0.0;
    // 	for(int covOrig=0;covOrig<=100;covOrig++){//coverage in the original genome
    // 	    //2 copies
    // 	    long double pOrig = poisson_pmfl((long double)covOrig,(long double)lambda);

    // #ifdef  DEBUG8CP
    // 	    cout<<"8cp covOrig\t"<<covOrig<<"\tpOrig\t"<<pOrig<<"\t"<<pDupInReference8<<endl;
    // #endif

    // 	    long double pDupInReference8_=0.0;
    // 	    string filenameIn = "../preComputated/combinatorialCov/comb8_"+stringify(cov)+"_"+stringify(covOrig)+".bin";
    // 	    ifstream in(filenameIn.c_str(), ios::in | ios::binary);	    
    // 	    if(!in.good()){
    // 		cerr<<"Cannot open "<<filenameIn<<endl;
    // 		exit(1);
    // 	    }
    // 	    in.read((char *) &pDupInReference8_, sizeof(pDupInReference8_) );
    // 	    in.close();


    // 	    pDupInReference8 += pDupInReference8_*pOrig;

    // 	}//end cov orig
	
    // 	double 	p2to2copies = 0.0;
    // 	for(int covOrig1=0;covOrig1<=100;covOrig1++){//coverage in the original genome
    // 	    double p2to2copies_=0.0;
    // 	    for(int covOrig2=0;covOrig2<=100;covOrig2++){//coverage in the original genome

    // 		double pOrig1 = poisson_pmf(double(covOrig1),double(lambda));
    // 		double pOrig2 = poisson_pmf(double(covOrig2),double(lambda));
    // 		double p1_2=  pOrig1*pOrig2;
    // 		//cout<<"c1\t"<<covOrig1<<"\tc2\t"<<covOrig2<<"\t"<<p1_2<<endl;
    // 		double p2to2copies__=0.0;
    // 		for(int cov1=0;cov1<=(covOrig1+covOrig2);cov1++){//coverage of 1st location
		    
    // 		    int cov2=(covOrig1+covOrig2-cov1); //coverage of second location
    // 		    if(cov1!=cov && cov2!=cov)
    // 			continue;
		    
    // 		    double probBinom = 	nChoosek((covOrig1+covOrig2),cov1) * pow(0.5,(covOrig1+covOrig2) );
    // 		    p2to2copies__ += probBinom*p1_2;
		    
    // 		    //cout<<cov1<<"\t"<<cov2<<"\t"<<probBinom<<"\t"<<probBinom*p1_2<<endl;
    // 		}
    // 		p2to2copies_ += p2to2copies__;

    // 		if(p2to2copies__>0 && p2to2copies__<1.0e-10)
    // 		    break;
    // 	    }

    // 	    p2to2copies += p2to2copies_;
	    
    // 	    if(p2to2copies_>0 && p2to2copies_<1.0e-10)
    // 		break;
    // 	}
	
    // #ifdef DEBUG
    // 	cout<<"P[cov|no dup] "<<pNoDup<<endl;
    // 	cout<<"P[cov|with collapse in reference] "<<pCollapseInReference<<endl;
    // 	cout<<"P[cov|with dup in reference] "<<pDupInReference<<endl;
    // 	cout<<"P[cov|3 with dup in reference] "<<pDupInReference3<<endl;
    // #endif
	
    // 	//double pFraction  = (pNoDup/(pNoDup+pCollapseInReference+pDupInReference));
    // 	//double pFraction3 = (pNoDup/(pNoDup+pCollapseInReference+pDupInReference+pDupInReference3));
    // 	double prior          = poisson_pmf(double(cov),double(lambda))/ poisson_pmf(double(lambda),double(lambda));
    // 	double pFraction82to2 = (pNoDup/(pNoDup+pCollapseInReference+pDupInReference2+pDupInReference3+pDupInReference4+pDupInReference5+pDupInReference6+pDupInReference7+pDupInReference8+p2to2copies));

    // 	cov2ProbSite->push_back(pFraction82to2);
	
    // 	//cout<<lambda<<"\t"<<cov<<"\t"<<pFraction3 <<endl;
    // 	//cout<<cov<<"\t"<<prior<<"\t"<< pFraction<<"\t"<<(pFraction*prior) <<"\t"<<pFraction3<<"\t"<<(pFraction3*prior) <<endl;

    // #ifdef DEBUG
    // 	cout<<"----------------------------------"<<endl;
    // 	cout<<"----------------------------------"<<endl;
    // #endif
    //     }//end cov
    //     //cout<<lambda<<"\t"<<vectorToString(pFinalVec,"\n")<<endl;

}

// void readMTConsensus(const string consensusFile,
// 		     map<int, PHREDgeno> & pos2phredgeno,
// 		     int & sizeGenome,
// 		     vector<int> & posOfIndels){

//     string line;
//     igzstream consensusFD;
//     consensusFD.open(consensusFile.c_str());
//     if (consensusFD.good()){
// 	getline (consensusFD,line);

// 	while ( getline (consensusFD,line)){
// 	    if (line.empty())
// 		continue;

// 	    vector<string> fields = allTokens(line,'\t');
// 	    PHREDgeno toadd;
// 	    // cerr<<line<<endl;


// 	    if(fields.size() != 11){
// 		cerr << "line "<<line<<"  in file  "<<consensusFile<<" does not have 11 fields"<<endl;
// 		exit(1);
// 	    }
	    

// 	    if(fields[0][fields[0].size()-1] == 'i'){ //skip insertion
// 		posOfIndels.push_back( destringify<int>( fields[0]) );
// 		continue;
// 	    }

// 	    if(fields[2] == "D"){ //skip deletions
// 		posOfIndels.push_back( destringify<int>( fields[0]) );
// 		continue;
// 	    }	    

// 	    toadd.consensus = fields[2][0];
// 	    for(int nuc=0;nuc<4;nuc++){		
// 		toadd.phred[nuc]  = destringify<double>(fields[nuc+7]);		
// 		toadd.perror[nuc] = pow(10.0,toadd.phred[nuc]/(-10.0));		
// 	    }

// 	    pos2phredgeno[     destringify<int>( fields[0])   ] = toadd;
// 	    sizeGenome =  max( destringify<int>( fields[0]), sizeGenome);
// 	    // cout<<destringify<int>( fields[0])<<endl;
	    
// 	}
// 	consensusFD.close();

//     }else{
// 	cerr << "Cannot open consensus file  "<<consensusFile<<""<<endl;
// 	exit(1);
//     }


// }

// void readMTAlleleFreq(const string freqFile,	map<int, alleleFrequency> & pos2allelefreq){
//     // map<int, alleleFrequency> pos2allelefreq;

//     string line;
//     igzstream freqAlleleFile;
//     freqAlleleFile.open(freqFile.c_str());
//     if (freqAlleleFile.good()){

// 	while ( getline (freqAlleleFile,line)){

// 	    vector<string> fields = allTokens(line,'\t');
// 	    alleleFrequency freqToadd;
	    
// 	    if(fields.size() != 5){
// 		cerr << "line "<<line<<"  in file  "<<freqFile<<" does not have 5 fields"<<endl;
// 		exit(1);
// 	    }
	   

// 	    for(int nuc=0;nuc<4;nuc++){
// 		freqToadd.f[nuc]=destringify<double>(fields[nuc+1]);
// 	    }

// 	    pos2allelefreq[ destringify<int>( fields[0])  ] = freqToadd;
	    	    
// 	}
// 	freqAlleleFile.close();

//     }else{
// 	cerr << "Cannot open allele frequency file  "<<freqFile<<""<<endl;
// 	exit(1);
//     }

// }
