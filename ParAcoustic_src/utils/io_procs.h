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
 *  io_procs.h
 *
 *
 *  Declares some functions used in reading and writing to NetCDF files.
 *
 *  Declares the following functions:
 *  openCDFFile
 *  writeCDFHeader
 *  writeCDFVarFloat
 *  readCDFVarBlock
 *
 */
/*
//These are io procedures.
*/
#ifndef _io_procs_h_
#define _io_procs_h_

#if defined(__cplusplus)
extern "C" {
#endif
  /*open a file, this first checks if the file name ends in cdf
  // then opens the file and makes sure it is really open*/
  int openCDFFile(const char* fileName,int create,int mode);

  /*Write a generic header. This includes all four dimensions, limits and increments
  // Variables for x,y,z and t 
  // if numVars is non-zero then x-y-z variables are defined, names are expected
  // in varNames to correspond with numVars*/
  int writeCDFHeader(const char* fileName,
		     int nx,int ny,int nz,int nt,
		     float* origin,
		     float dx,float minX,
		     float dy,float minY,
		     float dz,float minZ,
		     float dt,float minT,
		     int numVars,char** varNames);

  /*Read and write basic elements, each of these will open and close
  // the file, so it is not the best way to make alot of reads or writes
  // but it is simple if only 1 or 2 values are needed*/
  float writeCDFVarFloat(const char* fileName,
			 const char* varName,float data);

  /*These are slightly more complicated read and write functions for a 3-D
  // block of data*/
  int readCDFVarBlock(const char* fileName,
		      const char* varName,float** varData,
		      int* procLim);
#if defined(__cplusplus)
}
#endif



#endif /*_io_procs_h_*/
