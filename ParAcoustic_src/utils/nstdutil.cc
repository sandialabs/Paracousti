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
 *  nstdutil.cc
 *
 *
 *  Defines functions that have general functionality for many applications.  These are C++
 *  versions.
 *
 *  Defines the following functions:
 *  setVectorValues
 *  setVectorValues
 *  allocateEnough
 *  assert
 *  deassert
 *  tEprintf
 *
 */
/*Utility subroutines for C++ code.
*/
#include "nstdutil.hh"
#include <signal.h>

//
//Utility functions
//

//Here is a subroutine to parse a vector input from the command line.
// The input can be either 3 numbers ix dx nx or a matlab style spec
// ix:dx:fx.
void setVectorValues(int argc,char* argv[],int& i,
		     float& min,float& dx,int& NX){
  assert(argc>i+1,"setVectorValues--no arguments after flag");
    
  char *colon=strchr(argv[++i],':');
  if(colon){
    //Matlab style definitions start:inc:end or start:end
    char* rest;
    min=strtod(argv[i],&rest);
    assert(*rest==':',
	   "setVectorValues--malformed vector \"%s\"",argv[i]);
    dx=strtod(rest+1,&rest);
    double max;
    if(*rest!=':'){
      max=dx;
      dx=1.0;
    }else{
      max=atof(rest+1);
    }

    if(max>min){
      NX=(int)rint((max-min)/dx)+1;
    }else{
      NX=1;
    }
  }else{
    assert(argc>i+2,"setVectorValues--not enough arguments for 3 number spec.");
    min=atof(argv[i]);
    dx=atof(argv[++i]);
    NX=atoi(argv[++i]);
  }
}
void setVectorValues(int argc,char* argv[],int& i,
		     int &NX,intPtr& values){
  assert(argc>i+1,"setVectorValues--no arguments after flag");
    
  char *colon=strchr(argv[++i],':');
  int min,dx,max;
  if(colon){
    char* rest;
    min=strtol(argv[i],&rest,10);
    assert(*rest==':',
	   "setVectorValues--malformed vector \"%s\"",argv[i]);
    dx=strtol(rest+1,&rest,10);
    if(*rest!=':'){
      max=dx;
      dx=1;
    }else{
      max=atoi(rest+1);
    }

    if(max>min){
      NX=(int)rint((max-min)/dx)+1;
    }else{
      NX=1;
    }
  }else{
    assert(argc>i+2,"setVectorValues--not enough arguments for 3 number spec.");
    min=atoi(argv[i]);
    dx=atoi(argv[++i]);
    NX=atoi(argv[++i]);
  }

  assert((values=(int*)malloc(NX*sizeof(int)))!=NULL,
	 "setVectorValues--unable to allocate %i ints for vector",
	 NX);
  for(int i=0;i<NX;i++)
    values[i]=min+i*dx;
}
  

void* allocateEnough(void* ptr,int& currSize,int requiredSize,size_t objSize){
  void* returnVal=ptr;
  if(!ptr){
    assert((returnVal=malloc(requiredSize*objSize))!=NULL,
	   "allocateEnough--unable to allocate %i objects of size %i",
	   requiredSize,objSize);
    currSize=requiredSize;
  }else if(currSize<requiredSize){
    assert((returnVal=realloc(ptr,requiredSize*objSize))!=NULL,
	   "allocateEnough--unable to reallocate ptr from %i objects to %i (size %i)",
	   currSize,requiredSize,objSize);
    currSize=requiredSize;
  }
  return returnVal;
}

static FILE* _printFile=NULL;
static void (*_errorFunc_)(int)=NULL;

int assert(int test,const char* msg,...){
  if(test) return test;
  va_list args;
  va_start(args,msg);
  fprintf(stderr,"Failed Assertion \"");
  vfprintf(stderr,msg,args);
  fprintf(stderr,"\"\n***EXITING***\n");
  if(_printFile){
    fprintf(_printFile,"Failed Assertion \"");
    vfprintf(_printFile,msg,args);
    fprintf(_printFile,"\"\n***EXITING***\n");
  }
  va_end(args);
  if(_errorFunc_) _errorFunc_(test);

#ifdef INTERUPT_ASSERTION
  raise(SIGINT);
#endif
  exit(1);
  return 0;
}
 
int deassert(int test,const char* msg,...){
  if(!test) return !test;
  va_list args;
  va_start(args,msg);
  fprintf(stderr,"Failed de-Assertion \"");
  vfprintf(stderr,msg,args);
  fprintf(stderr,"\"\n***EXITING***\n");
  if(_printFile){
    fprintf(_printFile,"Failed de-Assertion \"");
    vfprintf(_printFile,msg,args);
    fprintf(_printFile,"\"\n***EXITING***\n");
  }
  va_end(args);
  if(_errorFunc_) _errorFunc_(test);

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
