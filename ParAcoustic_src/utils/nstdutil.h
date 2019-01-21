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
 *  nstdutil.h
 *
 *
 *  Declares functions that have general usability.  These are C versions.
 *
 *  Declares the following functions:
 *  setErrorFunc
 *  assert
 *  deassert
 *  tEprintf
 *  lower2
 *  NUMREC_dfour1
 *  NUMREC_four1
 *
 */
/*Utility subroutine headers for C code.
*/
#ifndef _nstdutic
#define _nstdutic

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "constants.h"

void setErrorFunc(void (*errorFunc)(int));

///Test an assertion exit on failure.
int assert(int test,const char* msg,...);
///Test a de-assertion exit on failure.
int deassert(int test,const char* msg,...);
///Print to standard error if test passes.
void tEprintf(int test,const char* printString,...);

//
//These are functions that are part of the c (not c++) interface that might be
// called from C++ code. The previous subroutines are only expected to be called
// from C code.
#if defined(__cplusplus)
extern "C" {
#endif

  ///Here is a function to convert a string to all lower case in a work array. 
  /// The string is unchanged. Note that work must be at least as long as str.
  char* lower2(char* str,char* work);

  //Here are some minimal numerical receipies subroutines that I actually use.
  /* (C) Copr. 1986-92 Numerical Recipes Software 21A24#21r1.. */
  void NUMREC_dfour1(double data[], unsigned long nn, int isign);
  void NUMREC_four1(float data[], unsigned long nn, int isign);

#if defined(__cplusplus)
}
#endif

#endif
