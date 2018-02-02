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

#include "hpdf.h"






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


class PdfWriter{
 private:
    string fname;
    HPDF_Doc  pdf;

    void draw_rect (HPDF_Page     page,
		    double        x,
		    double        y,
		    double length,
		    const char   *label);


    void draw_Simplerect (HPDF_Page     page,
			  double        x,
			  double        y,
			  double length);
    inline void addRange(HPDF_Page & page,double begin,double end, const chrScreenInfo & chrInfToUse );
    inline void addRangeCov(HPDF_Page & page,double begin,double end, const chrScreenInfo & chrInfToUse, double covFrac , int indexofinputF);    
 public:
    PdfWriter(const string fname);
    PdfWriter(const PdfWriter & other);
    ~PdfWriter();
    PdfWriter & operator= (const PdfWriter & other);

};


#endif
