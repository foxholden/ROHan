/*
 * PdfWriter
 * Date: Feb-01-2018 
 * Author : Gabriel Renaud gabriel.reno [at sign here ] gmail.com
 *
 */

#include "PdfWriter.h"
jmp_buf env;

#ifdef HPDF_DLL
void  __stdcall
#else
void
#endif
error_handler  (HPDF_STATUS   error_no,
                HPDF_STATUS   detail_no,
                void         *user_data)
{
    printf ("ERROR: error_no=%04X, detail_no=%u\n", (HPDF_UINT)error_no,
                (HPDF_UINT)detail_no);
    longjmp(env, 1);
}

PdfWriter::PdfWriter(const string fname_,
		     double heightChr_ //height of the rectangle set by the user
){
     fname     = fname_;
     heightChr = heightChr_;
     totalNumChrToDraw = 0;
    //string fname      = string(argv[indexOflastOpt]);
     string page_title = "";

     //gray = {0.6,0.6,0.6};
     // gray[0]  = 0.6;
     // gray[1]  = 0.6;
     // gray[2]  = 0.6;
     // cerr<<"test gray"<<endl;
     
	 
     pdf = HPDF_New (error_handler, NULL);
     if (!pdf) {
	 cerr<<"error: cannot create PdfDoc object for file="<<fname<<endl;
	 exit(1);
     }

     if(setjmp(env)){
	 HPDF_Free (pdf);
	 cerr<<"error: jump failed file="<<fname<<endl;
	 exit(1);
     }
     font = HPDF_GetFont (pdf, "Helvetica", NULL);

     /* add a new page object. */
     page = HPDF_AddPage (pdf);





     // /* print the lines of the page. */
     //HPDF_Page_SetLineWidth (page, 1);
     // HPDF_Page_Rectangle (page, 50, 50, HPDF_Page_GetWidth(page) - 100,
     //             HPDF_Page_GetHeight (page) - 110);
     // HPDF_Page_Stroke (page);

     /* print the title of the page (with positioning center). */
     HPDF_Page_SetFontAndSize (page, font, 24);
     tw = HPDF_Page_TextWidth (page, page_title.c_str());
     HPDF_Page_BeginText (page);
     HPDF_Page_MoveTextPos (page, (HPDF_Page_GetWidth(page) - tw) / 2,
			    HPDF_Page_GetHeight (page) - 50);
     //HPDF_Page_ShowText (page, page_title.c_str());
     HPDF_Page_EndText (page);
     HPDF_Page_SetFontAndSize (page, font, 10);
     //cerr<<"Creating a page in file  "<<fname_<<"w="<<HPDF_Page_GetWidth(page)<<" h="<<HPDF_Page_GetHeight (page)<<endl;
     //hmmPrevious = false;

     hmmprevdata  hmmpmin;
     hmmprevdata  hmmpmid;
     hmmprevdata  hmmpmax;

     hmmpmin.hmmPrevious =false;
     hmmpmid.hmmPrevious =false;
     hmmpmax.hmmPrevious =false;
     
     hmmpv.push_back( hmmpmin ); // [HMMCODEMIN]
     hmmpv.push_back( hmmpmid ); // [HMMCODEMID] 
     hmmpv.push_back( hmmpmax ); // [HMMCODEMAX] 
     
  
}

PdfWriter::~PdfWriter(){
    
     /* save the document to a file */
    cerr<<"Saving PDF file "<<fname<<endl;
    HPDF_SaveToFile (pdf, fname.c_str());

     /* clean up */
     HPDF_Free (pdf);

}

void PdfWriter::draw_rect (HPDF_Page     page,
			   double        x,
			   double        y,
			   double length,
			   const char   *label){    
    
    //cerr<<"draw_rect "<<x<<","<<y<<" l = "<<length<<" label="<<label<<" "<<GRAY<<endl;
    //text pos
    HPDF_Page_BeginText (page);
    HPDF_Page_MoveTextPos (page, x, y - heightLabel+2);
    HPDF_Page_ShowText (page, label);
    HPDF_Page_EndText (page);
    
    //rectangle
    HPDF_Page_SetRGBStroke  (page,
    			     GRAY,// GRAY,// gray[0],
    			     GRAY,// GRAY,// gray[1],
    			     GRAY// GRAY// gray[2]
    );
    //cerr<<"set1 "<<endl;
    HPDF_Page_Rectangle(page,
			x,             //x
			y - (heightChr+heightLabel),  //y  offset minus the heightLabel and height
			length,  //width
			heightChr);     //height
    ///cerr<<"set2 "<<endl;
     // HPDF_Page_SetRGBStroke  (page,
     // 			      0,// 0.0,
     // 			      0,// 0.0,
     // 			      0);// 0.0);
    //cerr<<"done"<<endl;
    // HPDF_Page_Fill (page);
}


void PdfWriter::draw_Simplerect (HPDF_Page     page,
				 double        x,
				 double        y,
				 double length){    
    // cout<<"x "<<x<<endl;
    // cout<<"y "<<y<<endl;

    length=max(length,0.4);
    // cout<<"l "<<length<<endl;
    // HPDF_Page_BeginText (page);
    // //text pos
    // HPDF_Page_MoveTextPos (page, x, y - 10);
    // HPDF_Page_ShowText (page, label);
    // HPDF_Page_EndText (page);
    //rectangle

    HPDF_Page_Rectangle(page, x, y - 30, length, 15);
}

inline void PdfWriter::addRange(HPDF_Page & page,double begin,double end, const chrScreenInfo & chrInfToUse ){
    draw_Simplerect (page, 
		     10+ chrInfToUse.lengthScreen*(begin/chrInfToUse.length ),
		     chrInfToUse.y,
		     chrInfToUse.lengthScreen * ((end-begin)/chrInfToUse.length) );			     
    // cout<<line<<endl;
    HPDF_Page_Fill (page);
}

inline void PdfWriter::addRangeCov(HPDF_Page & page,double begin,double end, const chrScreenInfo & chrInfToUse, double covFrac , int indexofinputF){
    // cout<<"1 "<<(covFrac)<<endl;
    // cout<<"2 " <<(1.0-covFrac)<<endl;
    //    covFrac=double(rand())/double(RAND_MAX);



    if(indexofinputF==1){
	HPDF_Page_SetRGBStroke (page, 1.00,        1.0-covFrac,  1.0-covFrac);
	HPDF_Page_SetRGBFill   (page, 1.00,        1.0-covFrac,  1.0-covFrac);
    }else if(indexofinputF==2){
	HPDF_Page_SetRGBStroke (page, 1.0-covFrac, 1.00,         1.0-covFrac);
	HPDF_Page_SetRGBFill   (page, 1.0-covFrac, 1.00,         1.0-covFrac);
    }else if(indexofinputF==3){
	HPDF_Page_SetRGBStroke (page, 1.0-covFrac, 1.0-covFrac,  1.00);
	HPDF_Page_SetRGBFill   (page, 1.0-covFrac, 1.0-covFrac,  1.00);
    }else{
	cerr<<"Color not defined for this file"<<endl;
	exit(1);
    }

    draw_Simplerect (page, 
		     10+ chrInfToUse.lengthScreen*(begin/chrInfToUse.length ),
		     chrInfToUse.y,
		     chrInfToUse.lengthScreen * ((end-begin)/chrInfToUse.length) );			     
    // cout<<line<<endl;
    HPDF_Page_Fill (page);
}


int PdfWriter::drawFrame(const string & fastaIndex,
			 const double   windowSizeForHest,
			 int indexFirstChr ,
			 const  set<string> * listAutosomes,
			 bool specifiedAutosomes){
    
     string line;
     igzstream myFaidxFile;
     // bool oneChr=false;
     // string oneChrName="";
     int indexLastChr=indexFirstChr;
     myFaidxFile.open(fastaIndex.c_str(), ios::in);


     vector<chrinfo> chrFound;
     unsigned int  genomeLength=0;
     unsigned int  maxLengthFound=0;
     if (myFaidxFile.good()){
	 while ( getline (myFaidxFile,line)){
	     chrinfo toadd;
	     vector<string> fields = allTokens(line,'\t');

	     toadd.name         =fields[0];
	     toadd.startIndexChr=genomeLength+1;
	     toadd.length       =destringify<unsigned int>(fields[1]);
	     if(specifiedAutosomes){//we have specified the autosomes		 
		 if(listAutosomes->find(toadd.name) == listAutosomes->end()){//not found in autosomes
		     continue;
		 }
		 totalNumChrToDraw++;
	     }else{
		 totalNumChrToDraw++;
	     }

	     if(toadd.length> maxLengthFound){
		 maxLengthFound = toadd.length;
	     }
	     toadd.endIndexChr  =genomeLength+toadd.length;

	     chrFound.push_back(toadd);

	 }
	 myFaidxFile.close();
     }else{
	 cerr << "Unable to open fasta index file "<<fastaIndex<<endl;
	 return 1;
     }

     
     //bool found=false;
     // if(oneChr){
     // 	 for(unsigned int i=0;
     // 	     i<chrFound.size();
     // 	     i++){
     // 	     if(chrFound[i].name == oneChrName){
     // 		 found=true;
     // 		 maxLengthFound=chrFound[i].length;
     // 	     }
     // 	 }

     // 	 if(!found){
     // 	     cerr<<"Chromosome you entered "<<oneChrName<<" was not found"<<endl;
     // 	     return 1;
     // 	 }
     // }


     //double sizeToUse=HPDF_Page_GetHeight(page)/double(2.0*chrFound.size());
     double sizeToUse = heightLabel + heightChr;
     //cerr<<HPDF_Page_GetHeight(page)

     double widthScreen= (HPDF_Page_GetWidth(page)-10.0);    
     for(int i=indexFirstChr;
	 i<int(chrFound.size());
	 i++){
	 double yOffset=sizeToUse*double(i-indexFirstChr);

	 //cerr<<"i1 ="<<i<<"\t"<<yOffset<<"\t"<<chrFound[i].name<<" "<<(yOffset+getTotalHeightChr())<<" "<<HPDF_Page_GetHeight(page) <<endl;
	 //reached the end
	 if( (yOffset+getTotalHeightChr()) > HPDF_Page_GetHeight(page)){
	     break;
	 }
	 //cerr<<"accepted"<<endl;
	 indexLastChr = i;

	 // if(oneChr){
	 //     if(chrFound[i].name != oneChrName)
	 // 	 continue;
	 //     yOffset=0;
	 // }
	
	 draw_rect (page, 
		    xmargin ,                            //x
		    HPDF_Page_GetHeight(page) - yOffset, //y
		    (widthScreen-xmargin)  * (double(chrFound[i].length)/double(maxLengthFound)),     //length
		    //heightChr,//height
		    chrFound[i].name.c_str());

	 //cerr<<indexFirstChr<<"\tchrFound["<<i<<"].name="<< chrFound[i].name<<"\t"<<yOffset<<endl;	 
	 name2chrScreenInfo[ chrFound[i].name.c_str() ].y            = HPDF_Page_GetHeight(page) -yOffset; //y offset
	 name2chrScreenInfo[ chrFound[i].name.c_str() ].length       = double(chrFound[i].length);   //actual length in bases 
	 name2chrScreenInfo[ chrFound[i].name.c_str() ].lengthScreen = ( (HPDF_Page_GetWidth(page)- 2*xmargin)  * (double(chrFound[i].length)/double(maxLengthFound)) );     //length on screen in pixel

	 //cerr<<"i3 ="<<i<<"\t"<<yOffset<<endl;



	 //top corner minus height label + height chr drawing

	 HPDF_Page_Stroke (page);
	 //cerr<<
     }
     //cerr<< HPDF_Page_GetHeight(page) <<endl;
     //exit(1);
     //cerr<<"test "<<endl;
     

     HPDF_Page_GSave (page);
     gstate = HPDF_CreateExtGState (pdf);
     HPDF_ExtGState_SetAlphaFill (gstate, alpha);
     HPDF_Page_SetExtGState (page, gstate);


     //write out x labels
     for(int i=indexFirstChr;
	 i<=indexLastChr;
	 i++){
    // for(int i=indexFirstChr;
    // 	 i<int(chrFound.size());
    // 	 i++){

	 
	 unsigned int mark = 100000000;
	 while(mark>=1){
	     //cerr<<"mark "<<mark<<endl;
	     if( name2chrScreenInfo[ chrFound[i].name.c_str() ].length > mark){
		 name2chrScreenInfo[ chrFound[i].name.c_str() ].markings = mark;
		 break;
	     }else{
		 mark           = mark / 10;
	     }
	 }
	 //cerr<<"mark "<<mark<<endl;
	 //DRAW MARKINGS
	 //cerr<<chrFound[i].name<<endl;
	 for(unsigned int m=0;m<name2chrScreenInfo[ chrFound[i].name.c_str() ].length;m+=mark ){
	     //name2chrScreenInfo[ chrFound[i].name.c_str() ].y - (heightChr+heightLabel)
	     double xfraction = ( m  /(name2chrScreenInfo[ chrFound[i].name.c_str() ].length+windowSizeForHest));
	     double x  = xmargin +   xfraction*name2chrScreenInfo[ chrFound[i].name.c_str() ].lengthScreen ;
	     double y1 = name2chrScreenInfo[ chrFound[i].name.c_str() ].y - (heightChr+heightLabel);
	     double y2 = name2chrScreenInfo[ chrFound[i].name.c_str() ].y - (heightChr+heightLabel)-heightLabel*fracHeightLabelMark;//1/5 of the height label
	     
	     //cerr<<"x "<<x<<" m "<<m<<" "<<(name2chrScreenInfo[ chrFound[i].name.c_str() ].y - (heightChr+heightLabel))<<endl;

	     // HPDF_Page_SetRGBStroke  (page,
	     // 			      gray[0],
	     // 			      gray[1],
	     // 			      gray[2]);
	     
	     drawVerticalLine( x,                                                                                                             // x offset
	      		       y1,  // baseline
	     		       y2,
	     		       GRAY,// gray[0],
	     		       GRAY,// gray[1],
	     		       GRAY,// gray[2],
	     		       1.0//widthtouse/3.0
	     ); // baseline

	     //text pos
	     string labl ;
	     if(m==0){
		 labl = "0";
	     }else{
		 if(m > 1000000){
		     labl     = stringify(m/1000000)+"Mb";
		 }else{
		     if(m > 1000){
			 labl = stringify(m/1000)   +"kb";
		     }else{
			 labl = stringify(m)        + "b";
		     }
		 }
		 //labl = "0";
	     }
	     if(labl == "0") 	 continue;
	     //cerr<<"labl "<<labl<<endl;
	     double fontsize =5;
	     HPDF_Page_SetFontAndSize (page, font, 5);

	     HPDF_Page_BeginText   (page);
	     HPDF_Page_MoveTextPos (page, x - double(labl.size())/2.0, y2-fontsize);
	     HPDF_Page_ShowText    (page, labl.c_str() );
	     HPDF_Page_EndText     (page);
	     
	 }
     }
     
     // for (map<string,chrScreenInfo >::iterator it=name2chrScreenInfo.begin(); it!=name2chrScreenInfo.end(); ++it)
     // 	 cerr << it->first  << '\n';

     return 0;
}

void PdfWriter::setFname(const string newfname){
    fname=newfname;
}

string PdfWriter::getFname() const{
    return fname;
}

int PdfWriter::getTotalNumChrToDraw() const{
    return totalNumChrToDraw;
}

int PdfWriter::getPageHeight() const{
    return int(HPDF_Page_GetHeight(page));
}



int PdfWriter::getPageWidth() const{
    return int(HPDF_Page_GetWidth(page));
}

int PdfWriter::getHeightLabel() const{//height of label
    return int(heightLabel);
}

int PdfWriter::getHeightChr()   const{//height of label
    return int(heightChr);
}

int PdfWriter::getTotalHeightChr()   const{//height of label
    return int(heightChr+heightLabel+heightLabel*fracHeightLabelMark);
}


int PdfWriter::drawYLabels(const long double minHFoundPlotting,
			   const long double maxHFoundPlotting,
			   const bool scientific){
    string formattouse;
    if(scientific){
	formattouse="%.0e\n";
    }else{
	formattouse="%.1f\n";
    }
    
    //write out y labels
    for(map<string , chrScreenInfo>::iterator it = name2chrScreenInfo.begin(); it != name2chrScreenInfo.end(); ++it) {

	
	// for(unsigned int i=0;
	//  	 i<chrFound.size();
	//  	 i++){
	
	// double xfraction = ( m  /(name2chrScreenInfo[ chrFound[i].name.c_str() ].length+windowSizeForHest));
	// double x = xmargin +   xfraction*name2chrScreenInfo[ chrFound[i].name.c_str() ].lengthScreen ;
	double x  = xmargin;
	double y1 = name2chrScreenInfo[ it->first ].y - (          heightLabel); //top
	double y2 = name2chrScreenInfo[ it->first ].y - (heightChr+heightLabel); //bottom

	//cerr<<"drawYLabels "<<x<<" "<<y1<<" "<<y2<<endl;
	double fontsize =2;
	stringstream stream;
	string labl;
	
	//labl=stringify(minHFoundPlotting);
	stream<<minHFoundPlotting;
	labl=stream.str();
	//cout<<"labl min "<<labl<<endl;
	HPDF_Page_SetFontAndSize (page, font, fontsize);

	
	HPDF_Page_BeginText   (page);
	HPDF_Page_MoveTextPos (page, x - double(labl.size())*fontsize*2, y2-fontsize/2.0); // x - double(labl.size())/2.
	HPDF_Page_ShowText    (page, labl.c_str() );
	HPDF_Page_EndText     (page);

	drawHorizontalLine((x - x/4.0), //0,//x/2.0  ,
			   x,
			   y1,
			   GRAY,
			   GRAY,
			   GRAY,
			   1.0       );

	//labl=stringify(maxHFoundPlotting);
	char buf[64];
	
	sprintf(buf,formattouse.c_str(), double(maxHFoundPlotting));
	
	labl = string(buf);
	//cout<<"labl max "<<labl<<endl;
	HPDF_Page_SetFontAndSize (page, font, fontsize);
	
	HPDF_Page_BeginText   (page);
	HPDF_Page_MoveTextPos (page, x - double(labl.size())*1.5, y1-fontsize/4.0);
	HPDF_Page_ShowText    (page, labl.c_str() );
	HPDF_Page_EndText     (page);
	
	
	drawHorizontalLine((x - x/4.0), //0,//x/2.0  ,
			   x,
			   y2,
			   GRAY,
			   GRAY,
			   GRAY,
			   1.0       );



	double factor=double(1)/double(3.0);
	double y3 = name2chrScreenInfo[ it->first ].y - (heightChr*(1-factor)+heightLabel); //bottom

	
	sprintf(buf,formattouse.c_str(), double(maxHFoundPlotting*factor));

	labl = string(buf);
	//cout<<"labl max "<<labl<<endl;
	HPDF_Page_SetFontAndSize (page, font, fontsize);
	
	HPDF_Page_BeginText   (page);
	HPDF_Page_MoveTextPos (page, x - double(labl.size())*1.5, y3-fontsize/4.0);
	HPDF_Page_ShowText    (page, labl.c_str() );
	HPDF_Page_EndText     (page);
	
	
	drawHorizontalLine((x - x/4.0), //0,//x/2.0  ,
			   x,
			   y3,
			   GRAY,
			   GRAY,
			   GRAY,
			   1.0       );


	factor=double(2)/double(3.0);
	double y4 = name2chrScreenInfo[ it->first ].y - (heightChr*(1-factor)+heightLabel); //bottom

	sprintf(buf,formattouse.c_str(), double(maxHFoundPlotting*factor));
	
	labl = string(buf);
	//cout<<"labl max "<<labl<<endl;
	HPDF_Page_SetFontAndSize (page, font, fontsize);
	
	HPDF_Page_BeginText   (page);
	HPDF_Page_MoveTextPos (page, x - double(labl.size())*1.5, y4-fontsize/4.0);
	HPDF_Page_ShowText    (page, labl.c_str() );
	HPDF_Page_EndText     (page);
	
	
	drawHorizontalLine((x - x/4.0), //0,//x/2.0  ,
			   x,
			   y4,
			   GRAY,
			   GRAY,
			   GRAY,
			   1.0       );


	

	
	// cerr<<"x "<<x<<" m "<<m<<" "<<(name2chrScreenInfo[ chrFound[i].name.c_str() ].y - (heightChr+heightLabel))<<endl;	 
	// // HPDF_Page_SetRGBStroke  (page,
	// // 			      gray[0],
	// // 			      gray[1],
	// // 			      gray[2]);	
	// drawVerticalLine( x,                                                                                                             // x offset
	//  		       y1,  // baseline
	// 		       y2,
	// 		       GRAY,// gray[0],
	// 		       GRAY,// gray[1],
	// 		       GRAY,// gray[2],
	// 		       1.0//widthtouse/3.0
	
    }
    return 0;
}
   


int PdfWriter::drawVerticalLine(const double x,
				const double y1,
				const double y2,
				double r,
				double g,
				double b,
				double w, //width
				bool   dash
){
    // HPDF_Page_BeginText (page);
    // HPDF_Page_MoveTextPos (page, x, y - 10);
    // HPDF_Page_ShowText (page, label);
    // HPDF_Page_EndText (page);
    //cerr<<"drawHorizontalLine "<<x<<" "<<y1<<" "<<y2<<endl;
    if(dash){
	HPDF_Page_SetDash (page, DASH_MODE1, 1, 1);//set line dash
    }

    HPDF_Page_SetRGBStroke (page, r, g, b);//set color
    HPDF_Page_SetLineWidth (page, w); //set line width
    
    HPDF_Page_MoveTo (page, x,   y1);
    HPDF_Page_LineTo (page, x,   y2);
    HPDF_Page_Stroke (page);
    
    HPDF_Page_SetRGBStroke (page, 0, 0, 0); //unset color
    HPDF_Page_SetLineWidth (page, 1.0); //unset line width

    HPDF_Page_SetDash (page, NULL, 0, 0); //unset line dash
    
    return 0;
}
#define NEWPLOT

bool  PdfWriter::chrIspresent(const string chrname) const{
    return ( name2chrScreenInfo.find( chrname ) != name2chrScreenInfo.end() );
}

int PdfWriter::drawHEst(const GenomicRange  cr,         // genomic range to plot
			const long double   h_,          // value of h
			const long double   hlow_,       // lower  conf. int. for h
			const long double   hhigh_,      // higher conf. int. for h
			const double        hLimLow,    // lower  limit for the h plot 
			const double        hLimHigh,
			const double        windowSizeForHest){  // higher limit for the h plot 
    // HPDF_Page_BeginText (page);
    // HPDF_Page_MoveTextPos (page, x, y - 10);
    // HPDF_Page_ShowText (page, label);
    // HPDF_Page_EndText (page);
    long double   h     = h_;
    long double   hlow  = hlow_;
    long double   hhigh = hhigh_;

    if(h < hLimLow){	 h  = hLimLow;    }
    if(h > hLimHigh){	 h  = hLimHigh;   }

    if(hlow < hLimLow){	 hlow  = hLimLow;    }
    if(hlow > hLimHigh){ hlow  = hLimHigh;   }
    
    //if the hLimHigh was bound to avoid large confidence interval
    if(hhigh > hLimHigh){hhigh  = hLimHigh;    }
    if(hhigh > hLimHigh){hhigh  = hLimHigh;    }
    
    //cerr<<"drawHEst "<<cr<<" "<<h<<" "<<hlow<<" "<<hhigh<<" hlims "<<hLimLow<<" "<<hLimHigh<<	endl;

    if( name2chrScreenInfo.find( cr.getChrName() ) == name2chrScreenInfo.end() ){
	cerr<<"PdfWriter chromosome not found: "<<cr.getChrName()<<endl;
	return 1;	
    }

    double xfraction = ( ( (cr.getEndCoord()+cr.getStartCoord())/2.0 )/(name2chrScreenInfo[ cr.getChrName() ].length+windowSizeForHest));
    
    //cerr<<"drawHEst2 "<<cr.getEndCoord()<<" "<<cr.getStartCoord()<<" "<<name2chrScreenInfo[ cr.getChrName() ].y<<" l="<<name2chrScreenInfo[ cr.getChrName() ].length<<" "<<name2chrScreenInfo[ cr.getChrName() ].lengthScreen<<" "<<(( (h-hLimLow)/hLimHigh) ) <<" xf="<<xfraction<<endl;

    double x = xmargin +   xfraction*name2chrScreenInfo[ cr.getChrName() ].lengthScreen ;
	
    // drawHorizontalLine( x,//x offset
    // 			name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel), //baseline
    // 			name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h-hLimLow)/hLimHigh) );//baseline
    //cerr<<"PX "<<(name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (hhigh-hLimLow)/hLimHigh))<<endl;
    // upper estimate

    // drawHorizontalLine( x,                                                                                                             // x offset
    //  			name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh) ,  // baseline
    //  			name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (hhigh-hLimLow)/hLimHigh),
    // 			0,
    // 			0,
    // 			0.5,
    // 			0.0
    // ); // baseline
    //cerr<<name2chrScreenInfo[ cr.getChrName() ].le
    double widthtouse = name2chrScreenInfo[ cr.getChrName() ].lengthScreen/(name2chrScreenInfo[ cr.getChrName() ].length/windowSizeForHest);
    //cerr<<"PX "<<cr<<" "<<name2chrScreenInfo[ cr.getChrName() ].y <<" "<<(name2chrScreenInfo[ cr.getChrName() ].y-(heightChr))<<" "<<(name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h-hLimLow)/hLimHigh))<<" "<<h<<" "<<hLimLow<<" "<<hLimHigh<<endl;
    
    // #ifdef NEWPLOT 
    drawVerticalLine( x,                                                                                                             // x offset
		      name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh) ,  // baseline
		      name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (hhigh-hLimLow)/hLimHigh),
		      0.0,
		      0.0,
		      0.0,
		      0.0//widthtouse/3.0
    ); // baseline

    

    drawCircle(x-widthtouse/2.0,                                                                                                             // x offset
    	       x+widthtouse/2.0,                                                                                                             // x offset			,
    	       name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh),
    	       0.9,
    	       0.0,
    	       0.0,
    	       widthtouse/10.0);
	
    // drawHorizontalLine( x-widthtouse/2.0,                                                                                                             // x offset
    // 			x+widthtouse/2.0,                                                                                                             // x offset			
    // 			name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh),                   // baseline
    // 			0,
    // 			0,
    // 			0,
    // 			0.2);
			


    // lower estimate
    drawVerticalLine( x,                                                                                                             // x offset
		      name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (hlow -hLimLow)/hLimHigh) ,  // baseline
		      name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh) ,
		      0.0,
		      0.0,
		      0.0,
		      0.0//widthtouse/3.0
    ); // baseline
    // #else
    //     drawVerticalLine( x,                                                                                                             // x offset
    // 		      name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh) ,  // baseline
    // 		      name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (hhigh-hLimLow)/hLimHigh),
    // 		      0.9,
    // 		      0,
    // 		      0,
    // 		      0.0//widthtouse/3.0
    //     ); // baseline
    
    

    //     drawHorizontalLine( x-widthtouse/2.0,                                                                                                             // x offset
    // 			x+widthtouse/2.0,                                                                                                             // x offset			
    // 			name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh),                   // baseline
    // 			0,
    // 			0,
    // 			0,
    // 			0.2);
    


    //     // lower estimate
    //     drawVerticalLine( x,                                                                                                             // x offset
    // 		      name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (hlow -hLimLow)/hLimHigh) ,  // baseline
    // 		      name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh) ,
    // 		      0.9,
    // 		      0.0,
    // 		      0.0,
    // 		      0.0//widthtouse/3.0
    //     ); // baseline
    // #endif


    
    // drawHorizontalLine( x,                                                                                                             // x offset
    //  			name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh)-1 ,  // baseline
    //  			name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh)+1,
    // 			0,
    // 			0,
    // 			0.9,
    // 			2.0
    // ); // baseline

    //  			name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh)-1 ,  // baseline
    //  			name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh)+1,
    // 			0,
    // 			0,
    // 			0.9,
    // 			2.0
    // ); 
    // drawVerticalLine( x,                                                                                                             // x offset
    // 		      name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh),  // baseline
    // 		      name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh),
    // 		      0.9,
    // 		      0,
    // 		      0,
    // 		      widthtouse
    // ); // baseline

    //name2chrScreenInfo[ chrFound[i].name.c_str() ].length       = double(chrFound[i].length);   //length
    // HPDF_Page_MoveTo (page, x,   y1);
    // HPDF_Page_LineTo (page, x,   y2);
    // HPDF_Page_Stroke (page);
    
    return 0;
}




int PdfWriter::drawHMM(const GenomicRange  cr,         // genomic range to plot
		       const long double   h_,          // value of h
		       const long double   hlow_,       // lower  conf. int. for h
		       const long double   hhigh_,      // higher conf. int. for h
		       const double        hLimLow,    // lower  limit for the h plot 
		       const double        hLimHigh,
		       const double        windowSizeForHest,
		       const bool          defined,
		       const bool          chrbreak,		       
		       const unsigned char useminmidmax){  // higher limit for the h plot 
    long double   h     = h_;
    //long double   hlow  = hlow_;
    long double   hhigh = hhigh_;

    //if the hLimHigh was bound to avoid large confidence interval
    if(hhigh > hLimHigh){
	hhigh  = hLimHigh;
    }
    //cerr<<"drawHMM "<<cr<<" "<<h<<" "<<hlow<<" "<<hhigh<<" hlims "<<hLimLow<<" "<<hLimHigh<<	endl;

    if( name2chrScreenInfo.find( cr.getChrName() ) == name2chrScreenInfo.end() ){
	cerr<<"PdfWriter chromosome not found: "<<cr.getChrName()<<endl;
	return 1;	
    }

    double xfraction  = ( ( (cr.getEndCoord()+cr.getStartCoord())/2.0 )/(name2chrScreenInfo[ cr.getChrName() ].length+windowSizeForHest));
    
    // cerr<<"drawHEst2 "<<cr.getEndCoord()<<" "<<cr.getStartCoord()<<" "<<name2chrScreenInfo[ cr.getChrName() ].y<<" l="<<name2chrScreenInfo[ cr.getChrName() ].length<<" "<<name2chrScreenInfo[ cr.getChrName() ].lengthScreen<<" "<<(( (h-hLimLow)/hLimHigh) ) <<" xf="<<xfraction<<endl;

    double x          = xmargin +   xfraction*name2chrScreenInfo[ cr.getChrName() ].lengthScreen ;
	
    double widthtouse = name2chrScreenInfo[ cr.getChrName() ].lengthScreen/(name2chrScreenInfo[ cr.getChrName() ].length/windowSizeForHest);

    
     // cerr<<"drawHMM "<<cr.getEndCoord()<<" "<<cr.getStartCoord()<<" "<<xfraction<<" "<<x<<" "<<widthtouse<<" "<<name2chrScreenInfo[ cr.getChrName() ].y <<" "<< (heightChr+heightLabel) <<" "<<heightChr <<" "<<  ( (h    -hLimLow)/hLimHigh)<<" "<<((name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh)))<<endl;

    double r,g,b,w;
    double xoffsetVisible=0;
    
    if(defined){
	if( useminmidmax == HMMCODEMIN){//red
	    r =   1.0; g = 0.0; b =   0.0; w = 1.0; xoffsetVisible=-1;
	}
	if( useminmidmax == HMMCODEMID){//green
	    r =   0.0; g = 1.0; b =   0.0; w = 1.0; xoffsetVisible= 0;
	}
	if( useminmidmax == HMMCODEMAX){//magenta
	    r =   1.0; g = 0.0; b =   1.0; w = 1.0; xoffsetVisible= 1;
	}
    }else{
	r = 1.0; g = 1.0; b =   0; w = 0.5;  //yellow
    }

    widthtouse = 0.9*widthtouse;

    double xlhmm = x-widthtouse/2.0;// x offset
    double xrhmm = x+widthtouse/2.0;// x offset			
    double yhmm  = (name2chrScreenInfo[ cr.getChrName() ].y - (heightChr+heightLabel) + heightChr*  ( (h    -hLimLow)/hLimHigh))  + xoffsetVisible; // baseline
    
    
    drawHorizontalLine( xlhmm,
			xrhmm,
			yhmm ,
			r,g,b,w);

    // double xlhmmPrevious;
    // double xrhmmPrevious;
    // double yhmmPrevious;
    // bool   hmmPrevious;



    // if( useminmidmax == HMMCODEMIN){ 	}
    // 	if( useminmidmax == HMMCODEMID){
    // 	}
    // 	if( useminmidmax == HMMCODEMAX){	}
    
    if(chrbreak){ hmmpv[useminmidmax].hmmPrevious =false; }
    // cout<<"draw hmm previous "<<hmmpv[useminmidmax].xlhmmPrevious<<" "<<hmmpv[useminmidmax].xrhmmPrevious<<" "<<hmmpv[useminmidmax].yhmmPrevious<<" "<<xlhmm<<" "<<xrhmm<<" "<<yhmm<<endl; 



    
    if(hmmpv[useminmidmax].hmmPrevious){
	// cout<<"draw hmm previous "<<hmmpv[useminmidmax].xlhmmPrevious<<" "<<hmmpv[useminmidmax].xrhmmPrevious<<" "<<hmmpv[useminmidmax].yhmmPrevious<<" "<<xlhmm<<" "<<xrhmm<<" "<<yhmm<<endl; 

	//cout<<"1"<<endl;
    	//draw connector to previous
    	HPDF_Page_SetLineWidth(page, w/2);
	
	//cout<<"2"<<endl;
	HPDF_Page_SetRGBStroke (page, r, g, b);//set color
	
	//cout<<"3"<<endl;
	//current 
	HPDF_Page_MoveTo(page,hmmpv[useminmidmax].xrhmmPrevious,hmmpv[useminmidmax].yhmmPrevious);
	//cout<<"4"<<endl;

	double curvature=1;
	//if(0) //todo remove
    	HPDF_Page_CurveTo( page,

			   //first point for derivative
    			   hmmpv[useminmidmax].xrhmmPrevious+curvature,//xlhmm,        //x1
    			   hmmpv[useminmidmax].yhmmPrevious, //y1

			   //second point for derivative
			   xlhmm-curvature, //xrhmmPrevious,//x2
    			   yhmm,         //y2

			   //next connecting point
    			   xlhmm,        //x3
    			   yhmm);        //y3

       	//cout<<"4"<<endl;
    	HPDF_Page_Stroke (page);
	//cout<<"5"<<endl;
    }
    
    hmmpv[useminmidmax].xlhmmPrevious = xlhmm;
    hmmpv[useminmidmax].xrhmmPrevious = xrhmm;
    hmmpv[useminmidmax].yhmmPrevious  = yhmm;   
    hmmpv[useminmidmax].hmmPrevious   = true;
    
    return 0;
}





int PdfWriter::drawGlobalHEst(//string chrname,
			      const long double   h_,
			      const long double   hlow_,
			      const long double   hhigh_,
			      const double hLimLow,
			      const double hLimHigh){
    //cerr<<"drawGlobalHEst "<<h<<" "<<hlow<<" "<<hhigh<<endl;
    // if( name2chrScreenInfo.find( chrname ) == name2chrScreenInfo.end() ){
    // 	cerr<<"PdfWriter chromosome not found: "<<chrname<<endl;
    // 	return 1;	
    // }
    long double   h     = h_;
    long double   hlow  = hlow_;
    long double   hhigh = hhigh_;

    if(h < hLimLow){	 h  = hLimLow;    }
    if(h > hLimHigh){	 h  = hLimHigh;   }

    if(hlow < hLimLow){	 hlow  = hLimLow;    }
    if(hlow > hLimHigh){ hlow  = hLimHigh;   }
    
    //if the hLimHigh was bound to avoid large confidence interval
    if(hhigh > hLimHigh){hhigh  = hLimHigh;    }
    if(hhigh > hLimHigh){hhigh  = hLimHigh;    }
    
    vector<string> allnames;
    for(map<string , chrScreenInfo>::iterator it = name2chrScreenInfo.begin(); it != name2chrScreenInfo.end(); ++it) {
	// cerr<<fname<<endl;
	// cerr<<"drawGlobalHEst "<<h<<" "<<hlow<<" "<<hhigh<<" "<<it->first<<endl;
	// cerr<<name2chrScreenInfo[ it->first ].y  <<" "<< heightChr<<" "<<heightLabel <<" "<<  ( (h     -hLimLow)/hLimHigh)<<" "<<(name2chrScreenInfo[ it->first ].y - (heightChr+heightLabel) + heightChr*  ( (h     -hLimLow)/hLimHigh))<<endl;
	// cerr<<xmargin<<" "<<(xmargin+name2chrScreenInfo[ it->first ].lengthScreen)<<" "<<HPDF_Page_GetWidth(page)<<endl;
	// //low
	// drawVerticalLine(xmargin,
	// 		 name2chrScreenInfo[ it->first ].lengthScreen,
	// 		 name2chrScreenInfo[ it->first ].y - (heightChr+heightLabel) + heightChr*  ( (hlow  -hLimLow)/hLimHigh) );//baseline
	//middle

	
	// drawHorizontalLine(xmargin,
	// 		   xmargin+100,
	// 		   name2chrScreenInfo[ it->first ].y - (heightChr+heightLabel) + heightChr*  ( (h     -hLimLow)/hLimHigh),
	// 		   0.0,
	// 		   0.0,
	// 		   0.7,
	// 		   0.3, //w
	// 		   true
	// );//baseline

	//	if(0){

	drawHorizontalLine(xmargin,
			   xmargin+name2chrScreenInfo[ it->first ].lengthScreen,
			   name2chrScreenInfo[ it->first ].y - (heightChr+heightLabel) + heightChr*  ( (h     -hLimLow)/hLimHigh),
			   0.0,
			   0.0,
			   0.7,
			   0.3, //w
			   true
	);//baseline

	drawHorizontalLine(xmargin,
			   xmargin+name2chrScreenInfo[ it->first ].lengthScreen,
			   name2chrScreenInfo[ it->first ].y - (heightChr+heightLabel) + heightChr*  ( (hlow  -hLimLow)/hLimHigh),
			   0.0,
			   0.0,
			   0.2,
			   0.1, //w
			   true
	);//baseline
	
	drawHorizontalLine(xmargin,
			   xmargin+name2chrScreenInfo[ it->first ].lengthScreen,
			   name2chrScreenInfo[ it->first ].y - (heightChr+heightLabel) + heightChr*  ( (hhigh -hLimLow)/hLimHigh),
			   0.0,
			   0.0,
			   0.2,
			   0.1, //w
			   true
	);//baseline
	// }
	// break;
	// //high
	// drawVerticalLine(xmargin,
	// 		 name2chrScreenInfo[ it->first ].lengthScreen,
	// 		 name2chrScreenInfo[ it->first ].y - (heightChr+heightLabel) + heightChr*  ( (hhigh -hLimLow)/hLimHigh) );//baseline
	
    }
    
    return 0;
}

int PdfWriter::drawHorizontalLine(const double x1,
				  const double x2,
				  const double y,
				  double r,
				  double g,
				  double b,
				  double w,
				  bool   dash
){

    if(dash){
	HPDF_Page_SetDash (page, DASH_MODE1, 1, 1);//set line dash
    }
    HPDF_Page_SetRGBStroke (page, r, g, b);//set color
    HPDF_Page_SetLineWidth (page, w); //set line width

    //cerr<<"drawHorizontalLine "<<x1<<" "<<x2<<" "<<y<<" "<<r<<" "<<g<<" "<<b<<" "<<w<<" "<<dash<<endl;

    HPDF_Page_MoveTo (page, x1,   y);
    HPDF_Page_LineTo (page, x2,   y);
    HPDF_Page_Stroke (page);

    //if(dash){
    HPDF_Page_SetDash (page, NULL, 0, 0); //unset line dash
    //}
    HPDF_Page_SetRGBStroke (page, 0.0, 0.0, 0.0 );//unset color
    HPDF_Page_SetLineWidth (page, 1.0); //unset line width

    return 0;
}




int PdfWriter::drawCircle(const double x1,
			  const double x2,
			  const double y,
			  double r,
			  double g,
			  double b,
			  double radius){
    HPDF_Page_SetLineWidth (page, 0.0);

    HPDF_Page_SetRGBStroke (page, r, g, b);//set color
    //    HPDF_Page_SetLineWidth (page, ra); //set line width
    HPDF_Page_SetRGBFill (page, r, g, b);

    //cerr<<"drawHorizontalLine "<<x1<<" "<<x2<<" "<<y<<" "<<r<<" "<<g<<" "<<b<<" "<<w<<" "<<dash<<endl;
    HPDF_Page_Circle ( page, (x1+x2)/2.0 , y, radius);
    // HPDF_Page_MoveTo (page, x1,   y);
    // HPDF_Page_LineTo (page, x2,   y);
    // HPDF_Page_Stroke (page);
    HPDF_Page_ClosePathFillStroke (page);

    HPDF_Page_SetRGBStroke (page, 0.0, 0.0, 0.0 );//unset color
    HPDF_Page_SetRGBFill (page, 0.0, 0.0, 0.0);
    //HPDF_Page_SetLineWidth (page, 1.0); //unset line width

    return 0;
}
