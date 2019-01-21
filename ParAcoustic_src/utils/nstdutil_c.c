/* Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
 * Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
 * certain rights in this software.
 *
 * NOTICE:
 * For five (5) years from  the United States Government is granted for itself and others
 *  acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this
 *  data to reproduce, prepare derivative works, and perform publicly and display
 *  publicly, by or on behalf of the Government. There is provision for the possible
 *  extension of the term of this license. Subsequent to that period or any extension
 *  granted, the United States Government is granted for itself and others acting on its
 *  behalf a paid-up, nonexclusive, irrevocable worldwide license in this data to reproduce,
 *  prepare derivative works, distribute copies to the public, perform publicly and display
 *  publicly, and to permit others to do so. The specific term of the license can be
 *  identified by inquiry made to National Technology and Engineering Solutions of Sandia,
 *  LLC or DOE.
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR
 *  NATIONAL TECHNOLOGY AND ENGINEERING SOLUTIONS OF SANDIA, LLC, NOR ANY OF THEIR EMPLOYEES,
 *  MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL RESPONSIBILITY FOR THE
 *  ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS
 *  DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
 * Any licensee of this software has the obligation and responsibility to abide by the
 *  applicable export control laws, regulations, and general prohibitions relating to the
 *  export of technical data. Failure to obtain an export control license or other authority
 *  from the Government may result in criminal liability under U.S. laws.*/
/*
 *  nsrdutil_c.c
 *
 *
 *  Defines some functions that have general usability.
 *
 *  Defines the following functions:
 *  setErrorFunc
 *  assert
 *  deassert
 *  tEprintf
 *  lower2
 *  NUMREC_dfour1
 *  NUMREC_four1
 *
 *  Defines the following macros (only used internally):
 *  SWAP
 *
 */
#include <signal.h>
#include "nstdutil.h"
static FILE* _printFile=NULL;
static void (*_errorFunc_c_)(int)=NULL;

void setErrorFunc(void (*errorFunc)(int)){
  _errorFunc_c_=errorFunc;
}

int assert(int test,const char* msg,...){
  va_list args;
  if(test) return test;
  va_start(args,msg);
  printf("Failed Assertion \"");
  vprintf(msg,args);
  printf("\"\n***EXITING***\n");
  if(_printFile){
    fprintf(_printFile,"Failed Assertion \"");
    vfprintf(_printFile,msg,args);
    fprintf(_printFile,"\"\n***EXITING***\n");
  }
  va_end(args);
  if(_errorFunc_c_) _errorFunc_c_(test);
#ifdef INTERUPT_ASSERTION
  raise(SIGINT);
#endif
  exit(1);
  return 0;
}
  
int deassert(int test,const char* msg,...){
  va_list args;
  if(!test) return !test;
  va_start(args,msg);
  printf("Failed de-Assertion \"");
  vprintf(msg,args);
  printf("\"\n***EXITING***\n");
  if(_printFile){
    fprintf(_printFile,"Failed de-Assertion \"");
    vfprintf(_printFile,msg,args);
    fprintf(_printFile,"\"\n***EXITING***\n");
  }
  va_end(args);
  if(_errorFunc_c_) _errorFunc_c_(test);
#ifdef INTERUPT_ASSERTION
  raise(SIGINT);
#endif
  exit(1);
  return 0;
}

void tEprintf(int test,const char* printString,...){
  if(test){
    va_list args;
    va_start(args,printString);
    vfprintf(stderr,printString,args);
    fflush(stderr);
  }
}


///Here is a function to convert a string to all lower case in a work array. 
/// The string is unchanged. Note that work must be at least as long as str.
char* lower2(char* str,char* work){
  int i;
  strcpy(work,str);
  for(i=0;i<strlen(work);i++)
    work[i]=tolower(work[i]);
  return work;
}

//Here are some minimal numerical receipies subroutines that I actually use.
#ifndef SWAP
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#endif

/* (C) Copr. 1986-92 Numerical Recipes Software 21A24#21r1.. */
void NUMREC_dfour1(double data[], unsigned long nn, int isign)
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  double tempr,tempi;

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software 21A24#21r1.. */

/* (C) Copr. 1986-92 Numerical Recipes Software 21A24#21r1.. */
void NUMREC_four1(float data[], unsigned long nn, int isign)
{
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;

  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software 21A24#21r1.. */
#undef SWAP
