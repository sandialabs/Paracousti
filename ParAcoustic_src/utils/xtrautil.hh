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
 *  xtrautil.hh
 *
 *
 *  Declares functions used to read NetCDF files.
 *
 *  Declares the following functions:
 *  readModelVariable
 *
 */
#ifndef _xtrautil
#define _xtrautil

#include "netcdf.h"
#include "io_procs.h"
#include "nstdutil.hh"

//
//CDF output
//

//
//Here are some new subroutines to read a variable from a model. The variable
// name is checked for upper and lower case versions. A possible 1D version of the
// variable is checked. The file is checked for versions that have been written
// with the variable name appended. And then we look for files that have been
// broken into subdomains.
float* readModelVariable(floatPtr &var,
			 char* modelName,const char* varName,
			 int *procLim,int globalNX,int globalNY,int globalNZ,
			 int modelFile=-1,int noFindError=TRUE);
int* readModelVariable(intPtr &var,
                       char* modelName,const char* varName,
                       int *procLim,int globalNX,int globalNY,int globalNZ,
                       int modelFile=-1,int noFindError=TRUE);

#endif //_xtrautil
