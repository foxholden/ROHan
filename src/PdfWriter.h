/*
 * PdfWriter
 * Date: Feb-01-2018 
 * Author : Gabriel Renaud gabriel.reno [at sign here] gmail.com
 *
 */

#ifndef PdfWriter_h
#define PdfWriter_h

using namespace std;

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <gzstream.h>
#include <setjmp.h>
#include <vector>

#include "miscfunc.h"
#include "hpdf.h"
#include "GenomicWindows.h"

#define GRAY 0.6



/* typedef struct{ */
/*     string name; */
/*     unsigned int startIndexChr; */
/*     unsigned int endIndexChr; */
/*     unsigned int length; */
/* } chrinfo; */


//position on the screen
typedef struct{
    double  y;
    double  length;
    double  lengthScreen;
    unsigned int markings;
} chrScreenInfo;

const HPDF_UINT16 DASH_MODE1[] = {3};

class PdfWriter{
 private:
    double xlhmmPrevious;
    double xrhmmPrevious;
    double yhmmPrevious;
    
    bool   hmmPrevious;
    string fname;
    HPDF_Doc  pdf;
    HPDF_Font font;
    HPDF_Page page;
    HPDF_ExtGState gstate;
    float tw;
    double alpha=0.8;
    double xmargin=10;
    
    double fracHeightLabelMark=0.2;
    
    double heightLabel=15;//height of label
    double heightChr  =15;//height of label
    //double gray [3] = {0.6,0.6,0.6};
    int totalNumChrToDraw;
    map<string, chrScreenInfo>  name2chrScreenInfo;

    void draw_rect (HPDF_Page     page,
		    double        x,
		    double        y,		    
		    double length  ,
		    //double height  ,
		    const char   *label);


    void draw_Simplerect (HPDF_Page     page,
			  double        x,
			  double        y,
			  double length);
    inline void addRange(HPDF_Page & page,double begin,double end, const chrScreenInfo & chrInfToUse );
    inline void addRangeCov(HPDF_Page & page,double begin,double end, const chrScreenInfo & chrInfToUse, double covFrac , int indexofinputF);    
 public:
    PdfWriter(const string fname,const double heightChr);
    PdfWriter(const PdfWriter & other);
    ~PdfWriter();
    PdfWriter & operator= (const PdfWriter & other);
    int drawFrame(const string & fastaIndex,const double   windowSizeForHest,int indexFirstChr=0,const  set<string> * listAutosomes=NULL, bool specifiedAutosomes=false );

    int drawYLabels(const long double minHFoundPlotting,
		    const long double maxHFoundPlotting,
		    const bool scientific);

    int drawVerticalLine(const double x,
			 const double y1,
			 const double y2,
			 double r=0.0,
			 double g=0.0,
			 double b=0.0,
			 double w=1.0,
			 bool   dash=false
    );

    int drawHorizontalLine(const double x1,
			   const double x2,
			   const double y,
			   double r=0.0,
			   double g=0.0,
			   double b=0.0,
			   double w=1.0,
			   bool   dash=false
    );

    int drawHEst(const GenomicRange  cr,
		 const long double   h,
		 const long double   hlow,
		 const long double   hhigh,
		 const double hLimLow,
		 const double hLimHigh,
		 const double windowSizeForHest);

    int drawHMM(const GenomicRange  cr,         // genomic range to plot
		const long double   h_,          // value of h
		const long double   hlow_,       // lower  conf. int. for h
		const long double   hhigh_,      // higher conf. int. for h
		const double        hLimLow,    // lower  limit for the h plot 
		const double        hLimHigh,
		const double        windowSizeForHest,
		const bool          defined,
		const unsigned char useminmidmax);  // higher 
    
    int drawGlobalHEst(const long double   h,
		       const long double   hlow,
		       const long double   hhigh,
		       const double hLimLow,
		       const double hLimHigh);

    int drawCircle(const double x1,
		   const double x2,
		   const double y,
		   double r,
		   double g,
		   double b,
		   double radius);

    int getPageHeight()  const;
    int getPageWidth()   const;
    int getHeightLabel() const;//height of label
    int getHeightChr()   const;//height of label
    int getTotalHeightChr()   const;//height of label

    
    int getTotalNumChrToDraw() const;
    void setFname(const string newfname);
    string getFname() const;
    bool  chrIspresent(const string chrname) const;



};


#endif
