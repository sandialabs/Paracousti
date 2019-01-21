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
 *  nstdutil.hh
 *
 *
 *  Declares several functions that have general applicability.  These are C++ versions.
 *
 *  Declares the following classes:
 *  setVectorValues
 *  setVectorValues
 *  allocateEnough
 *  assert
 *  deassert
 *  tEprintf
 *
 */
/*Utility subroutine headers for C++ code.
*/

#ifndef _nstdutil
#define _nstdutil

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "constants.h"

//
//Utility functions
//
/// \brief Smart allocation.
/// Allocate ptr to required size and set currSize if it is not
/// set; do nothing if currSize >= required size; reallocate if required.
void* allocateEnough(void* ptr,int& currSize,int requiredSize,size_t objSize);

///Test an assertion exit on failure.
int assert(int test,const char* msg,...);

///Test a de-assertion exit on failure.
int deassert(int test,const char* msg,...);
///Print to standard error if test passes.
void tEprintf(int test,const char* printString,...);

//Here is a subroutine to parse a vector input from the command line.
// The input can be either 3 numbers ix dx nx or a matlab style spec
// ix:dx:fx.
void setVectorValues(int argc,char* argv[],int& i,
		     float& min,float& dx,int& NX);
void setVectorValues(int argc,char* argv[],int& i,
		     int &NX,intPtr& values);
#endif

