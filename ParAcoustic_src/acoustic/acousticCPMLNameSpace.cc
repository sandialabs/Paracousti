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
 *  acousticCPMLNameSpace.h
 *  
 *
 *  This is the declaration the special namespace I made for the CPML
 *  boundaries.  This keeps things better separated.
 *
 */
#include "acousticCPMLNameSpace.h"

namespace acoustic_CPML {
  float *__restrict__ kxTaper, *__restrict__ kyTaper, *__restrict__ kzTaper;
  float *__restrict__ alphaXTaper, *__restrict__ alphaYTaper, *__restrict__ alphaZTaper;
  float *__restrict__ xTaper, *__restrict__ yTaper, *__restrict__ zTaper;
  float *__restrict__ axTaper, *__restrict__ ayTaper, *__restrict__ azTaper;
  float *__restrict__ bxTaper, *__restrict__ byTaper, *__restrict__ bzTaper;
  float *__restrict__ kxTaperPH, *__restrict__ kyTaperPH, *__restrict__ kzTaperPH;
  float *__restrict__ alphaXTaperPH, *__restrict__ alphaYTaperPH, *__restrict__ alphaZTaperPH;
  float *__restrict__ xTaperPH, *__restrict__ yTaperPH, *__restrict__ zTaperPH;
  float *__restrict__ axTaperPH, *__restrict__ ayTaperPH, *__restrict__ azTaperPH;
  float *__restrict__ bxTaperPH, *__restrict__ byTaperPH, *__restrict__ bzTaperPH;
  float *__restrict__ xTaperRel, *__restrict__ yTaperRel;
  float *__restrict__ xTaperRelPH, *__restrict__ yTaperRelPH;
  float *__restrict__ del;
  
  float *__restrict__ vxxXmin, *__restrict__ vxxXmax, *__restrict__ vxxYmin, *__restrict__ vxxYmax, *__restrict__ vxxZmin, *__restrict__ vxxZmax;
  float *__restrict__ vyyXmin, *__restrict__ vyyXmax, *__restrict__ vyyYmin, *__restrict__ vyyYmax, *__restrict__ vyyZmin, *__restrict__ vyyZmax;
  float *__restrict__ vzzXmin, *__restrict__ vzzXmax, *__restrict__ vzzYmin, *__restrict__ vzzYmax, *__restrict__ vzzZmin, *__restrict__ vzzZmax;
  
  float *__restrict__ pxXmin, *__restrict__ pxXmax, *__restrict__ pxYmin, *__restrict__ pxYmax, *__restrict__ pxZmin, *__restrict__ pxZmax;
  float *__restrict__ pyXmin, *__restrict__ pyXmax, *__restrict__ pyYmin, *__restrict__ pyYmax, *__restrict__ pyZmin, *__restrict__ pyZmax;
  float *__restrict__ pzXmin, *__restrict__ pzXmax, *__restrict__ pzYmin, *__restrict__ pzYmax, *__restrict__ pzZmin, *__restrict__ pzZmax;
  
  int nyInt, nzInt;
  
  int xMinSize, xMaxSize, yMinSize, yMaxSize, zMinSize, zMaxSize;
  int _nXmin, _nXmax, _nYmin, _nYmax, _nZmin, _nZmax;
  int _xMaxStart, _yMaxStart, _zMaxStart;
  int xStart1, xStart2, yStart1, yStart2, zStart, zStart_z;
  int xEnd, yEnd, zEnd, xEnd1, yEnd1, zEnd1;
  int _surfaceBCMode;
  int _kstart;
  
}
