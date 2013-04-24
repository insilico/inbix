#ifndef __DCDFLIB_H__
#define __DCDFLIB_H__

/*
 * DCDFLIB is a FORTRAN90 library which evaluates the cumulative density 
 * function (CDF) associated with common probability distributions, by 
 * Barry Brown, James Lovato, Kathy Russell.
 * 
 * DCDFLIB includes routines for evaluating the cumulative density functions 
 * of a variety of standard probability distributions. An unusual feature of 
 * this library is its ability to easily compute any one parameter of the CDF 
 * given the others. This means that a single routine can evaluate the CDF 
 * given the usual parameters, or determine the value of a parameter that 
 * produced a given CDF value.
 * 
 * http://people.sc.fsu.edu/~jburkardt/f_src/dcdflib/dcdflib.html
 */

extern void cdfbet(int*,double*,double*,double*,double*,double*,double*,
	    int*,double*);
extern void cdfbin(int*,double*,double*,double*,double*,double*,double*,
	    int*,double*);
extern void cdfchi(int*,double*,double*,double*,double*,int*,double*);
extern void cdfchn(int*,double*,double*,double*,double*,double*,int*,double*);
extern void cdff(int*,double*,double*,double*,double*,double*,int*,double*);
extern void cdffnc(int*,double*,double*,double*,double*,double*,double*,
	    int*s,double*);
extern void cdfgam(int*,double*,double*,double*,double*,double*,int*,double*);
extern void cdfnbn(int*,double*,double*,double*,double*,double*,double*,
	    int*,double*);
extern void cdfnor(int*,double*,double*,double*,double*,double*,int*,double*);
extern void cdfpoi(int*,double*,double*,double*,double*,int*,double*);
extern void cdft(int*,double*,double*,double*,double*,int*,double*);

#endif

