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

PdfWriter::PdfWriter(const string fname_){
     fname = fname_;

    //string fname      = string(argv[indexOflastOpt]);
     string page_title = "";

 

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
    

  
}

PdfWriter::~PdfWriter(){
    
     /* save the document to a file */
    cerr<<"Saving file "<<fname<<endl;
    HPDF_SaveToFile (pdf, fname.c_str());

     /* clean up */
     HPDF_Free (pdf);

}

void PdfWriter::draw_rect (HPDF_Page     page,
			   double        x,
			   double        y,
			   double length,
			   const char   *label){    
    HPDF_Page_BeginText (page);
    //text pos
    HPDF_Page_MoveTextPos (page, x, y - 10);
    HPDF_Page_ShowText (page, label);
    HPDF_Page_EndText (page);
    //rectangle
    HPDF_Page_Rectangle(page, x, y - 30, length, 15);
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


int PdfWriter::drawFrame(const string & fastaIndex){
    
     string line;
     igzstream myFaidxFile;
     // bool oneChr=false;
     // string oneChrName="";

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


     double sizeToUse=HPDF_Page_GetHeight(page)/double(2.0*chrFound.size());

     map<string, chrScreenInfo>  name2chrScreenInfo;

     double widthScreen= (HPDF_Page_GetWidth(page)-10.0);    
     for(unsigned int i=0;
	 i<chrFound.size();
	 i++){

	 double yOffset=sizeToUse*double(i*2);
	 // if(oneChr){
	 //     if(chrFound[i].name != oneChrName)
	 // 	 continue;
	 //     yOffset=0;
	 // }
	
	 draw_rect (page, 
		    10 , //x
		    HPDF_Page_GetHeight(page) - yOffset, //y
		    widthScreen  * (double(chrFound[i].length)/double(maxLengthFound)),     //length
		    chrFound[i].name.c_str());

	 name2chrScreenInfo[ chrFound[i].name.c_str() ].y            = HPDF_Page_GetHeight(page) -yOffset; //y offset
	 name2chrScreenInfo[ chrFound[i].name.c_str() ].length       = double(chrFound[i].length);   //length

	 name2chrScreenInfo[ chrFound[i].name.c_str() ].lengthScreen = ( (HPDF_Page_GetWidth(page)-10.0)  * (double(chrFound[i].length)/double(maxLengthFound)) );     //length on screen
	 HPDF_Page_Stroke (page);
     }




     HPDF_Page_GSave (page);
     gstate = HPDF_CreateExtGState (pdf);
     HPDF_ExtGState_SetAlphaFill (gstate, alpha);
     HPDF_Page_SetExtGState (page, gstate);

     return 0;
}
