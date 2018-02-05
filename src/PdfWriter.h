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

#include "hpdf.h"
#include "GenomicWindows.h"





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
} chrScreenInfo;

const HPDF_UINT16 DASH_MODE1[] = {3};
class PdfWriter{
 private:
    string fname;
    HPDF_Doc  pdf;
    HPDF_Font font;
    HPDF_Page page;
    HPDF_ExtGState gstate;
    float tw;
    double alpha=0.8;
    double xmargin=10;
    double heightLabel=15;//height of label
    double heightChr  =15;//height of label


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
    int drawFrame(const string & fastaIndex);
    int drawHorizontalLine(const double x,const double y1,const double y2);
    int drawVerticalLine(const double x1,const double x2,const double y);

    int drawHEst(const GenomicRange  cr,
		 const long double   h,
		 const long double   hlow,
		 const long double   hhigh,
		 const double hLimLow,
		 const double hLimHigh    );

    int drawGlobalHEst(const long double   h,
		       const long double   hlow,
		       const long double   hhigh,
		       const double hLimLow,
		       const double hLimHigh);
};


#endif
