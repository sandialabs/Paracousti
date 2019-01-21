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
 *  updateAcousticMPMLFull.cpp
 *  
 *
 *  This file contains the setup and updating functions for all the
 *  dependent variables for the case where MPMLs serve as the absorbing
 *  BCs.  It also includes the specialized functions it uses if a pressure
 *  free surface is used.  For the updating functions, updating proceeds
 *  in exactly the order that the dependent variables are stored in
 *  memory.  It gives a big speed boost compared to updating interior
 *  points then positive x flanks, then negative x flanks, etc.
 *
 *  Defines the following functions:
 *  setupAcousticMPMLBoundsFull
 *  updateVxMPMLBounds
 *  updateVyMPMLBounds
 *  updateVzMPMLBounds
 *  updateVzPressFreeMPMLBounds
 *  updateAcousticPressureMPML
 *  updateAcousticPressurePressFreeMPML
 *
 */

#include "updateAcousticAttenMPMLFull.h"
#include "acousticAttenMPMLNameSpace.h"
#include "sgfd.hh"
#include "constants.h"

#define ELASTIC_FS_2_0 3
#define SURFACE_PRESS_FREE 2

//  This function sets up the extra variables and arrays needed for the
//  MPML BCs and fills in the values of the MPML sigmas, kappas and
//  alphas as a function of distance into the MPML zone
void setupAcousticAttenMPMLBoundsFull(int nXmin, float xMinVal, float xMinAVal, float xMinKVal,int nXmax, float xMaxVal, float xMaxAval, float xMaxKVal, int nYmin, float yMinVal, float yMinAVal, float yMinKVal,
                               int nYmax, float yMaxVal, float yMaxAVal, float yMaxKVal, int nZmin, float zMinVal, float zMinAVal, float zMinKVal,int nZmax, float zMaxVal, float zMaxAVal, float zMaxKVal,int kstart, int kstart_z, sgfdModel* model, int sbcmode, float xtap) {
  using namespace acousticAtten_MPML;
  DEF_MODEL_SIZE(model->modelDef());
  DEF_MODEL_LIMITS(model->modelDef());
  DEF_PARALLEL(model->parallelDef());
  
  //variables used for atten
  tcr = new float[NX];
  tap = new float[NX];
  tdvx = new float[NX];
  tdvy = new float[NX];
  tdvz = new float[NX];

  //define the thickness of the MPML zones for this domain
  _nXmin = MAX(0,MIN(NX,nXmin-procLim[0]));
  _nXmax = MAX(0,MIN(NX,procLim[1]-globalNX+nXmax));
  _nYmin = MAX(0,MIN(NY,nYmin-procLim[2]));
  _nYmax = MAX(0,MIN(NY,procLim[3]-globalNY+nYmax));
  _nZmin = sbcmode==SURFACE_PRESS_FREE?0:(MAX(0,MIN(NZ,nZmin-procLim[4])));
  _nZmax = MAX(0,MIN(NZ,procLim[5]-globalNZ+nZmax));
  
  //Only if this is the top of the model do we pay attention to the free surface
  if(sbcmode==SURFACE_PRESS_FREE && kstart==3) {
    _surfaceBCMode = sbcmode;
  } else
    _surfaceBCMode = -1;
  
  //set up various variables needed to track starting and stopping indices for all the dependent variables.  Versions like xStart1 versus xStart2 are needed for those variables that live at half-integer or interger node points
  _xMaxStart = NX-_nXmax;
  _yMaxStart = NY-_nYmax;
  _zMaxStart = NZ-_nZmax;
  
  xStart1 = MIN(MAX(_nXmin,1),NX-2);
  xStart2 = MIN(MAX(_nXmin,2),NX-2);
  xEnd = MAX(MIN(_xMaxStart,NX-2),xStart2);
  xEnd1 = MAX(MIN(_xMaxStart-1,NX-2),xStart1);
  
  yStart1 = MIN(MAX(_nYmin,1),NY-2);
  yStart2 = MIN(MAX(_nYmin,2),NY-2);
  yEnd = MAX(MIN(_yMaxStart,NY-2),yStart2);
  yEnd1 = MAX(MIN(_yMaxStart-1,NY-2),yStart1);
  
  zStart = MIN(MAX(_nZmin,kstart),NZ-2);
  zStart_z = MIN(MAX(_nZmin,kstart_z),NZ-2);
  zEnd = MAX(MIN(_zMaxStart,NZ-2),zStart);
  zEnd1 = MAX(MIN(_zMaxStart-1,NZ-2),zStart_z);
  
  _kstart = kstart;

  //compute the flank values for sigma needed at each of the flanks
  float xfac = 0;
  
  xMinVal = reflToD(xMinVal,nXmin*dx,model->_vMax)/(1.+2.*xfac);
  xMaxVal = reflToD(xMaxVal,nXmax*dx,model->_vMax)/(1.+2.*xfac);
  yMinVal = reflToD(yMinVal,nYmin*dy,model->_vMax)/(1.+2.*xfac);
  yMaxVal = reflToD(yMaxVal,nYmax*dy,model->_vMax)/(1.+2.*xfac);
  zMinVal = reflToD(zMinVal,nZmin*dz,model->_vMax)/(1.+2.*xfac);
  zMaxVal = reflToD(zMaxVal,nZmax*dz,model->_vMax)/(1.+2.*xfac);
  
  //Allocate space for all the 1-D profile functions needed for the MPMLs
  //k tapers are for kappa and will actually contain 1/kappa.
  //alpha tapers are for alpha; plain Tapers are for sigma
  //a and b tapers are functions of alpha, kappa and sigma tapers and are what are actually used in computations
  //Those tapers ending in PH are for half-interger locations, while those without are for integer locations
  kxTaper = new float[NX];
  kyTaper = new float[NY];
  kzTaper = new float[NZ];
  alphaXTaper = new float[NX];
  alphaYTaper = new float[NY];
  alphaZTaper = new float[NZ];
  xTaper = new float[NX];
  yTaper = new float[NY];
  zTaper = new float[NZ];
  axTaper = new float[NX];
  ayTaper = new float[NY];
  azTaper = new float[NZ];
  bxTaper = new float[NX];
  byTaper = new float[NY];
  bzTaper = new float[NZ];
  //    xTaperFac = new float[NX];
  //    yTaperFac = new float[NY];
  //    zTaperFac = new float[NZ];
  kxTaperPH = new float[NX];
  kyTaperPH= new float[NY];
  kzTaperPH = new float[NZ];
  alphaXTaperPH = new float[NX];
  alphaYTaperPH= new float[NY];
  alphaZTaperPH = new float[NZ];
  xTaperPH = new float[NX];
  yTaperPH= new float[NY];
  zTaperPH = new float[NZ];
  axTaperPH = new float[NX];
  ayTaperPH= new float[NY];
  azTaperPH = new float[NZ];
  bxTaperPH = new float[NX];
  byTaperPH= new float[NY];
  bzTaperPH = new float[NZ];
  //    xTaperPHFac = new float[NX];
  //    yTaperPHFac= new float[NY];
  //    zTaperPHFac = new float[NZ];
  if(_surfaceBCMode==ELASTIC_FS_2_0) {
    xTaperRel = new float[NX];
    yTaperRel = new float[NY];
    xTaperRelPH = new float[NX];
    yTaperRelPH = new float[NY];
  }
  
  //Now define array sizes for the flanks.  Only allocate what is necessary for this domain, so some (or all) of these could have zero size
  //These parameters define the sizes for the MPML memory variables
  //Note that Z boundaries cover the entire portion of the top or bottom parts of the model.
  //Y boundaries cover the all of the north and south sides of a model that is not already covered by the Z boundaries
  //X boundaries cover the east and west sides of the model that are not already covered by the Y or Z boundaries
  nyInt = _yMaxStart-_nYmin;
  nzInt = _zMaxStart-_nZmin;
  
  xMinSize = _nXmin*nyInt*nzInt;
  xMaxSize = _nXmax*nyInt*nzInt;
  yMinSize = NX*_nYmin*nzInt;
  yMaxSize = NX*_nYmax*nzInt;
  zMinSize = NXY*_nZmin;
  zMaxSize = NXY*_nZmax;
  
  //allocate arrays for memory variables
  //The variable names end with which flank they belong to (e.g., Zmin)
  //The beginning of the variable name gives the type of memory variable
  //So, vxx is the memory variable belonging to the x derivative of vx
  //v variables are used in the pressure updating equations, while p variables are used in the velocity updating equations
  vxxZmin = new float[NXY*_nZmin];
  vxxZmax = new float[NXY*_nZmax];
  vxxYmin = new float[NX*_nYmin*nzInt];
  vxxYmax = new float[NX*_nYmax*nzInt];
  vxxXmin = new float[_nXmin*nyInt*nzInt];
  vxxXmax = new float[_nXmax*nyInt*nzInt];
  
  vyyZmin = new float[NXY*_nZmin];
  vyyZmax = new float[NXY*_nZmax];
  vyyYmin = new float[NX*_nYmin*nzInt];
  vyyYmax = new float[NX*_nYmax*nzInt];
  vyyXmin = new float[_nXmin*nyInt*nzInt];
  vyyXmax = new float[_nXmax*nyInt*nzInt];
  
  vzzZmin = new float[NXY*_nZmin];
  vzzZmax = new float[NXY*_nZmax];
  vzzYmin = new float[NX*_nYmin*nzInt];
  vzzYmax = new float[NX*_nYmax*nzInt];
  vzzXmin = new float[_nXmin*nyInt*nzInt];
  vzzXmax = new float[_nXmax*nyInt*nzInt];
  
  pxZmin = new float[NXY*_nZmin];
  pxZmax = new float[NXY*_nZmax];
  pxYmin = new float[NX*_nYmin*nzInt];
  pxYmax = new float[NX*_nYmax*nzInt];
  pxXmin = new float[_nXmin*nyInt*nzInt];
  pxXmax = new float[_nXmax*nyInt*nzInt];
  
  pyZmin = new float[NXY*_nZmin];
  pyZmax = new float[NXY*_nZmax];
  pyYmin = new float[NX*_nYmin*nzInt];
  pyYmax = new float[NX*_nYmax*nzInt];
  pyXmin = new float[_nXmin*nyInt*nzInt];
  pyXmax = new float[_nXmax*nyInt*nzInt];
  
  pzZmin = new float[NXY*_nZmin];
  pzZmax = new float[NXY*_nZmax];
  pzYmin = new float[NX*_nYmin*nzInt];
  pzYmax = new float[NX*_nYmax*nzInt];
  pzXmin = new float[_nXmin*nyInt*nzInt];
  pzXmax = new float[_nXmax*nyInt*nzInt];

ayX0TaperPH = new float [_nXmin];
byX0TaperPH = new float [_nXmin];
ayE210TaperxPH = new float*__restrict__ [_nYmin];
byE210TaperxPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  ayE210TaperxPH[j] = new float [NX-xEnd1];
  byE210TaperxPH[j] = new float [NX-xEnd1];
}
ayE021TaperyPH = new float*__restrict__ [_nZmin];
byE021TaperyPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE021TaperyPH[k] = new float [NY-yEnd1];
  byE021TaperyPH[k] = new float [NY-yEnd1];
}
azE022TaperyPH = new float*__restrict__ [NZ-zEnd1];
bzE022TaperyPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE022TaperyPH[k] = new float [NY-yEnd1];
  bzE022TaperyPH[k] = new float [NY-yEnd1];
}
azC010Taper = new float*__restrict__*__restrict__ [_nZmin];
bzC010Taper = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC010Taper[k] = new float*__restrict__ [NY-yEnd1];
  bzC010Taper[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC010Taper[k][j] = new float [_nXmin];
    bzC010Taper[k][j] = new float [_nXmin];
  }
}
ayE202TaperxPH = new float*__restrict__ [NZ-zEnd1];
byE202TaperxPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE202TaperxPH[k] = new float [NX-xEnd1];
  byE202TaperxPH[k] = new float [NX-xEnd1];
}
axC001TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC001TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC001TaperxPH[k] = new float*__restrict__ [_nYmin];
  bxC001TaperxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC001TaperxPH[k][j] = new float [_nXmin];
    bxC001TaperxPH[k][j] = new float [_nXmin];
  }
}
axC001TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC001TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC001TaperyPH[k] = new float*__restrict__ [_nYmin];
  bxC001TaperyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC001TaperyPH[k][j] = new float [_nXmin];
    bxC001TaperyPH[k][j] = new float [_nXmin];
  }
}
axC011Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC011Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC011Taper[k] = new float*__restrict__ [NY-yEnd1];
  bxC011Taper[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC011Taper[k][j] = new float [_nXmin];
    bxC011Taper[k][j] = new float [_nXmin];
  }
}
azE011TaperzPH = new float*__restrict__ [_nZmin];
bzE011TaperzPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE011TaperzPH[k] = new float [_nYmin];
  bzE011TaperzPH[k] = new float [_nYmin];
}
ayC010TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
byC010TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC010TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  byC010TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC010TaperyPHxPH[k][j] = new float [_nXmin];
    byC010TaperyPHxPH[k][j] = new float [_nXmin];
  }
}
azC000TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
bzC000TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC000TaperzPH[k] = new float*__restrict__ [_nYmin];
  bzC000TaperzPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC000TaperzPH[k][j] = new float [_nXmin];
    bzC000TaperzPH[k][j] = new float [_nXmin];
  }
}
azC110TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
bzC110TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC110TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC110TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC110TaperzPH[k][j] = new float [NX-xEnd1];
    bzC110TaperzPH[k][j] = new float [NX-xEnd1];
  }
}
ayC101TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC101TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC101TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  byC101TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC101TaperzPHyPH[k][j] = new float [NX-xEnd1];
    byC101TaperzPHyPH[k][j] = new float [NX-xEnd1];
  }
}
azC100TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
bzC100TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC100TaperzPH[k] = new float*__restrict__ [_nYmin];
  bzC100TaperzPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC100TaperzPH[k][j] = new float [NX-xEnd1];
    bzC100TaperzPH[k][j] = new float [NX-xEnd1];
  }
}
axE012TaperyPH = new float*__restrict__ [NZ-zEnd1];
bxE012TaperyPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE012TaperyPH[k] = new float [_nYmin];
  bxE012TaperyPH[k] = new float [_nYmin];
}
ayE012TaperzPHyPH = new float*__restrict__ [NZ-zEnd1];
byE012TaperzPHyPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE012TaperzPHyPH[k] = new float [_nYmin];
  byE012TaperzPHyPH[k] = new float [_nYmin];
}
axC110TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bxC110TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC110TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC110TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC110TaperzPHxPH[k][j] = new float [NX-xEnd1];
    bxC110TaperzPHxPH[k][j] = new float [NX-xEnd1];
  }
}
azC001TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC001TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC001TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  bzC001TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC001TaperyPHxPH[k][j] = new float [_nXmin];
    bzC001TaperyPHxPH[k][j] = new float [_nXmin];
  }
}
azY1TaperPH = new float [NY-yEnd1];
bzY1TaperPH = new float [NY-yEnd1];
azE220TaperyPHxPH = new float*__restrict__ [NY-yEnd1];
bzE220TaperyPHxPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  azE220TaperyPHxPH[j] = new float [NX-xEnd1];
  bzE220TaperyPHxPH[j] = new float [NX-xEnd1];
}
azE220Taper = new float*__restrict__ [NY-yEnd1];
bzE220Taper = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  azE220Taper[j] = new float [NX-xEnd1];
  bzE220Taper[j] = new float [NX-xEnd1];
}
axE201TaperzPHxPH = new float*__restrict__ [_nZmin];
bxE201TaperzPHxPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE201TaperzPHxPH[k] = new float [NX-xEnd1];
  bxE201TaperzPHxPH[k] = new float [NX-xEnd1];
}
ayC001TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC001TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC001TaperxPH[k] = new float*__restrict__ [_nYmin];
  byC001TaperxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC001TaperxPH[k][j] = new float [_nXmin];
    byC001TaperxPH[k][j] = new float [_nXmin];
  }
}
axE201TaperzPH = new float*__restrict__ [_nZmin];
bxE201TaperzPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE201TaperzPH[k] = new float [NX-xEnd1];
  bxE201TaperzPH[k] = new float [NX-xEnd1];
}
azE101TaperxPH = new float*__restrict__ [_nZmin];
bzE101TaperxPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE101TaperxPH[k] = new float [_nXmin];
  bzE101TaperxPH[k] = new float [_nXmin];
}
ayC001TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC001TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC001TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  byC001TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC001TaperzPHxPH[k][j] = new float [_nXmin];
    byC001TaperzPHxPH[k][j] = new float [_nXmin];
  }
}
ayE120TaperyPH = new float*__restrict__ [NY-yEnd1];
byE120TaperyPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  ayE120TaperyPH[j] = new float [_nXmin];
  byE120TaperyPH[j] = new float [_nXmin];
}
axE101TaperzPH = new float*__restrict__ [_nZmin];
bxE101TaperzPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE101TaperzPH[k] = new float [_nXmin];
  bxE101TaperzPH[k] = new float [_nXmin];
}
ayE201Taper = new float*__restrict__ [_nZmin];
byE201Taper = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE201Taper[k] = new float [NX-xEnd1];
  byE201Taper[k] = new float [NX-xEnd1];
}
azE021TaperzPHyPH = new float*__restrict__ [_nZmin];
bzE021TaperzPHyPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE021TaperzPHyPH[k] = new float [NY-yEnd1];
  bzE021TaperzPHyPH[k] = new float [NY-yEnd1];
}
azC000TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
bzC000TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC000TaperxPH[k] = new float*__restrict__ [_nYmin];
  bzC000TaperxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC000TaperxPH[k][j] = new float [_nXmin];
    bzC000TaperxPH[k][j] = new float [_nXmin];
  }
}
azC011TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC011TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC011TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC011TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC011TaperyPHxPH[k][j] = new float [_nXmin];
    bzC011TaperyPHxPH[k][j] = new float [_nXmin];
  }
}
azE120Taper = new float*__restrict__ [NY-yEnd1];
bzE120Taper = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  azE120Taper[j] = new float [_nXmin];
  bzE120Taper[j] = new float [_nXmin];
}
ayC100Taper = new float*__restrict__*__restrict__ [_nZmin];
byC100Taper = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC100Taper[k] = new float*__restrict__ [_nYmin];
  byC100Taper[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC100Taper[k][j] = new float [NX-xEnd1];
    byC100Taper[k][j] = new float [NX-xEnd1];
  }
}
axC100TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
bxC100TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC100TaperyPH[k] = new float*__restrict__ [_nYmin];
  bxC100TaperyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC100TaperyPH[k][j] = new float [NX-xEnd1];
    bxC100TaperyPH[k][j] = new float [NX-xEnd1];
  }
}
axE102Taper = new float*__restrict__ [NZ-zEnd1];
bxE102Taper = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE102Taper[k] = new float [_nXmin];
  bxE102Taper[k] = new float [_nXmin];
}
ayE101TaperzPHxPH = new float*__restrict__ [_nZmin];
byE101TaperzPHxPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE101TaperzPHxPH[k] = new float [_nXmin];
  byE101TaperzPHxPH[k] = new float [_nXmin];
}
ayE021Taper = new float*__restrict__ [_nZmin];
byE021Taper = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE021Taper[k] = new float [NY-yEnd1];
  byE021Taper[k] = new float [NY-yEnd1];
}
ayX1TaperPH = new float [NX-xEnd1];
byX1TaperPH = new float [NX-xEnd1];
azE110TaperxPH = new float*__restrict__ [_nYmin];
bzE110TaperxPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  azE110TaperxPH[j] = new float [_nXmin];
  bzE110TaperxPH[j] = new float [_nXmin];
}
ayC010TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
byC010TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC010TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  byC010TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC010TaperzPH[k][j] = new float [_nXmin];
    byC010TaperzPH[k][j] = new float [_nXmin];
  }
}
axE120TaperyPHxPH = new float*__restrict__ [NY-yEnd1];
bxE120TaperyPHxPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  axE120TaperyPHxPH[j] = new float [_nXmin];
  bxE120TaperyPHxPH[j] = new float [_nXmin];
}
ayE102Taper = new float*__restrict__ [NZ-zEnd1];
byE102Taper = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE102Taper[k] = new float [_nXmin];
  byE102Taper[k] = new float [_nXmin];
}
axC100Taper = new float*__restrict__*__restrict__ [_nZmin];
bxC100Taper = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC100Taper[k] = new float*__restrict__ [_nYmin];
  bxC100Taper[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC100Taper[k][j] = new float [NX-xEnd1];
    bxC100Taper[k][j] = new float [NX-xEnd1];
  }
}
axE022TaperzPH = new float*__restrict__ [NZ-zEnd1];
bxE022TaperzPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE022TaperzPH[k] = new float [NY-yEnd1];
  bxE022TaperzPH[k] = new float [NY-yEnd1];
}
axC010TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bxC010TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC010TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC010TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC010TaperzPHxPH[k][j] = new float [_nXmin];
    bxC010TaperzPHxPH[k][j] = new float [_nXmin];
  }
}
azC000TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
bzC000TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC000TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  bzC000TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC000TaperzPHyPH[k][j] = new float [_nXmin];
    bzC000TaperzPHyPH[k][j] = new float [_nXmin];
  }
}
axE012Taper = new float*__restrict__ [NZ-zEnd1];
bxE012Taper = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE012Taper[k] = new float [_nYmin];
  bxE012Taper[k] = new float [_nYmin];
}
ayC101TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC101TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC101TaperxPH[k] = new float*__restrict__ [_nYmin];
  byC101TaperxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC101TaperxPH[k][j] = new float [NX-xEnd1];
    byC101TaperxPH[k][j] = new float [NX-xEnd1];
  }
}
ayC010Taper = new float*__restrict__*__restrict__ [_nZmin];
byC010Taper = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC010Taper[k] = new float*__restrict__ [NY-yEnd1];
  byC010Taper[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC010Taper[k][j] = new float [_nXmin];
    byC010Taper[k][j] = new float [_nXmin];
  }
}
azC001TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC001TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC001TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  bzC001TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC001TaperzPHyPH[k][j] = new float [_nXmin];
    bzC001TaperzPHyPH[k][j] = new float [_nXmin];
  }
}
azE021TaperzPH = new float*__restrict__ [_nZmin];
bzE021TaperzPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE021TaperzPH[k] = new float [NY-yEnd1];
  bzE021TaperzPH[k] = new float [NY-yEnd1];
}
ayC100TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
byC100TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC100TaperzPH[k] = new float*__restrict__ [_nYmin];
  byC100TaperzPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC100TaperzPH[k][j] = new float [NX-xEnd1];
    byC100TaperzPH[k][j] = new float [NX-xEnd1];
  }
}
azE102TaperzPH = new float*__restrict__ [NZ-zEnd1];
bzE102TaperzPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE102TaperzPH[k] = new float [_nXmin];
  bzE102TaperzPH[k] = new float [_nXmin];
}
ayC111TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC111TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC111TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  byC111TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC111TaperzPH[k][j] = new float [NX-xEnd1];
    byC111TaperzPH[k][j] = new float [NX-xEnd1];
  }
}
ayE120TaperxPH = new float*__restrict__ [NY-yEnd1];
byE120TaperxPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  ayE120TaperxPH[j] = new float [_nXmin];
  byE120TaperxPH[j] = new float [_nXmin];
}
ayE210Taper = new float*__restrict__ [_nYmin];
byE210Taper = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  ayE210Taper[j] = new float [NX-xEnd1];
  byE210Taper[j] = new float [NX-xEnd1];
}
axC000TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bxC000TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC000TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  bxC000TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC000TaperyPHxPH[k][j] = new float [_nXmin];
    bxC000TaperyPHxPH[k][j] = new float [_nXmin];
  }
}
axC011TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC011TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC011TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC011TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC011TaperxPH[k][j] = new float [_nXmin];
    bxC011TaperxPH[k][j] = new float [_nXmin];
  }
}
ayC011TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC011TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC011TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  byC011TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC011TaperzPH[k][j] = new float [_nXmin];
    byC011TaperzPH[k][j] = new float [_nXmin];
  }
}
ayE022TaperyPH = new float*__restrict__ [NZ-zEnd1];
byE022TaperyPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE022TaperyPH[k] = new float [NY-yEnd1];
  byE022TaperyPH[k] = new float [NY-yEnd1];
}
ayC000TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
byC000TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC000TaperzPH[k] = new float*__restrict__ [_nYmin];
  byC000TaperzPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC000TaperzPH[k][j] = new float [_nXmin];
    byC000TaperzPH[k][j] = new float [_nXmin];
  }
}
azC100TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bzC100TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC100TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  bzC100TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC100TaperzPHxPH[k][j] = new float [NX-xEnd1];
    bzC100TaperzPHxPH[k][j] = new float [NX-xEnd1];
  }
}
axC110TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
bxC110TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC110TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC110TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC110TaperzPHyPH[k][j] = new float [NX-xEnd1];
    bxC110TaperzPHyPH[k][j] = new float [NX-xEnd1];
  }
}
ayE102TaperzPHxPH = new float*__restrict__ [NZ-zEnd1];
byE102TaperzPHxPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE102TaperzPHxPH[k] = new float [_nXmin];
  byE102TaperzPHxPH[k] = new float [_nXmin];
}
axC010TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
bxC010TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC010TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC010TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC010TaperzPHyPH[k][j] = new float [_nXmin];
    bxC010TaperzPHyPH[k][j] = new float [_nXmin];
  }
}
azC000TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bzC000TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC000TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  bzC000TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC000TaperyPHxPH[k][j] = new float [_nXmin];
    bzC000TaperyPHxPH[k][j] = new float [_nXmin];
  }
}
azE011TaperyPH = new float*__restrict__ [_nZmin];
bzE011TaperyPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE011TaperyPH[k] = new float [_nYmin];
  bzE011TaperyPH[k] = new float [_nYmin];
}
axE110TaperyPHxPH = new float*__restrict__ [_nYmin];
bxE110TaperyPHxPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  axE110TaperyPHxPH[j] = new float [_nXmin];
  bxE110TaperyPHxPH[j] = new float [_nXmin];
}
ayC001TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC001TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC001TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  byC001TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC001TaperzPHyPH[k][j] = new float [_nXmin];
    byC001TaperzPHyPH[k][j] = new float [_nXmin];
  }
}
axC100TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
bxC100TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC100TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  bxC100TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC100TaperzPHyPH[k][j] = new float [NX-xEnd1];
    bxC100TaperzPHyPH[k][j] = new float [NX-xEnd1];
  }
}
axE120TaperyPH = new float*__restrict__ [NY-yEnd1];
bxE120TaperyPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  axE120TaperyPH[j] = new float [_nXmin];
  bxE120TaperyPH[j] = new float [_nXmin];
}
azC001TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC001TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC001TaperzPH[k] = new float*__restrict__ [_nYmin];
  bzC001TaperzPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC001TaperzPH[k][j] = new float [_nXmin];
    bzC001TaperzPH[k][j] = new float [_nXmin];
  }
}
axE210Taper = new float*__restrict__ [_nYmin];
bxE210Taper = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  axE210Taper[j] = new float [NX-xEnd1];
  bxE210Taper[j] = new float [NX-xEnd1];
}
azE210TaperyPHxPH = new float*__restrict__ [_nYmin];
bzE210TaperyPHxPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  azE210TaperyPHxPH[j] = new float [NX-xEnd1];
  bzE210TaperyPHxPH[j] = new float [NX-xEnd1];
}
axE210TaperxPH = new float*__restrict__ [_nYmin];
bxE210TaperxPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  axE210TaperxPH[j] = new float [NX-xEnd1];
  bxE210TaperxPH[j] = new float [NX-xEnd1];
}
axE102TaperzPH = new float*__restrict__ [NZ-zEnd1];
bxE102TaperzPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE102TaperzPH[k] = new float [_nXmin];
  bxE102TaperzPH[k] = new float [_nXmin];
}
ayC111TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC111TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC111TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  byC111TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC111TaperzPHyPH[k][j] = new float [NX-xEnd1];
    byC111TaperzPHyPH[k][j] = new float [NX-xEnd1];
  }
}
ayE110TaperyPHxPH = new float*__restrict__ [_nYmin];
byE110TaperyPHxPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  ayE110TaperyPHxPH[j] = new float [_nXmin];
  byE110TaperyPHxPH[j] = new float [_nXmin];
}
axC110TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
bxC110TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC110TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC110TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC110TaperxPH[k][j] = new float [NX-xEnd1];
    bxC110TaperxPH[k][j] = new float [NX-xEnd1];
  }
}
azC111TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC111TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC111TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC111TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC111TaperzPHyPH[k][j] = new float [NX-xEnd1];
    bzC111TaperzPHyPH[k][j] = new float [NX-xEnd1];
  }
}
azE210TaperyPH = new float*__restrict__ [_nYmin];
bzE210TaperyPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  azE210TaperyPH[j] = new float [NX-xEnd1];
  bzE210TaperyPH[j] = new float [NX-xEnd1];
}
azC111TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC111TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC111TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC111TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC111TaperxPH[k][j] = new float [NX-xEnd1];
    bzC111TaperxPH[k][j] = new float [NX-xEnd1];
  }
}
axC101TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC101TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC101TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  bxC101TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC101TaperzPHxPH[k][j] = new float [NX-xEnd1];
    bxC101TaperzPHxPH[k][j] = new float [NX-xEnd1];
  }
}
axE220TaperyPHxPH = new float*__restrict__ [NY-yEnd1];
bxE220TaperyPHxPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  axE220TaperyPHxPH[j] = new float [NX-xEnd1];
  bxE220TaperyPHxPH[j] = new float [NX-xEnd1];
}
axC010TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
bxC010TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC010TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC010TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC010TaperyPH[k][j] = new float [_nXmin];
    bxC010TaperyPH[k][j] = new float [_nXmin];
  }
}
axC000TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bxC000TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC000TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  bxC000TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC000TaperzPHxPH[k][j] = new float [_nXmin];
    bxC000TaperzPHxPH[k][j] = new float [_nXmin];
  }
}
axC100TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
bxC100TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC100TaperzPH[k] = new float*__restrict__ [_nYmin];
  bxC100TaperzPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC100TaperzPH[k][j] = new float [NX-xEnd1];
    bxC100TaperzPH[k][j] = new float [NX-xEnd1];
  }
}
axC101TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC101TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC101TaperxPH[k] = new float*__restrict__ [_nYmin];
  bxC101TaperxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC101TaperxPH[k][j] = new float [NX-xEnd1];
    bxC101TaperxPH[k][j] = new float [NX-xEnd1];
  }
}
azE201TaperxPH = new float*__restrict__ [_nZmin];
bzE201TaperxPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE201TaperxPH[k] = new float [NX-xEnd1];
  bzE201TaperxPH[k] = new float [NX-xEnd1];
}
axE210TaperyPH = new float*__restrict__ [_nYmin];
bxE210TaperyPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  axE210TaperyPH[j] = new float [NX-xEnd1];
  bxE210TaperyPH[j] = new float [NX-xEnd1];
}
azE120TaperyPH = new float*__restrict__ [NY-yEnd1];
bzE120TaperyPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  azE120TaperyPH[j] = new float [_nXmin];
  bzE120TaperyPH[j] = new float [_nXmin];
}
axC111TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC111TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC111TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC111TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC111TaperxPH[k][j] = new float [NX-xEnd1];
    bxC111TaperxPH[k][j] = new float [NX-xEnd1];
  }
}
azC110Taper = new float*__restrict__*__restrict__ [_nZmin];
bzC110Taper = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC110Taper[k] = new float*__restrict__ [NY-yEnd1];
  bzC110Taper[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC110Taper[k][j] = new float [NX-xEnd1];
    bzC110Taper[k][j] = new float [NX-xEnd1];
  }
}
ayZ0TaperPH = new float [_nZmin];
byZ0TaperPH = new float [_nZmin];
ayC010TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
byC010TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC010TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  byC010TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC010TaperzPHyPH[k][j] = new float [_nXmin];
    byC010TaperzPHyPH[k][j] = new float [_nXmin];
  }
}
azE201Taper = new float*__restrict__ [_nZmin];
bzE201Taper = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE201Taper[k] = new float [NX-xEnd1];
  bzE201Taper[k] = new float [NX-xEnd1];
}
azC001Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC001Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC001Taper[k] = new float*__restrict__ [_nYmin];
  bzC001Taper[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC001Taper[k][j] = new float [_nXmin];
    bzC001Taper[k][j] = new float [_nXmin];
  }
}
azE012Taper = new float*__restrict__ [NZ-zEnd1];
bzE012Taper = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE012Taper[k] = new float [_nYmin];
  bzE012Taper[k] = new float [_nYmin];
}
azC101TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC101TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC101TaperzPH[k] = new float*__restrict__ [_nYmin];
  bzC101TaperzPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC101TaperzPH[k][j] = new float [NX-xEnd1];
    bzC101TaperzPH[k][j] = new float [NX-xEnd1];
  }
}
ayX0Taper = new float [_nXmin];
byX0Taper = new float [_nXmin];
axE011Taper = new float*__restrict__ [_nZmin];
bxE011Taper = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE011Taper[k] = new float [_nYmin];
  bxE011Taper[k] = new float [_nYmin];
}
axE021Taper = new float*__restrict__ [_nZmin];
bxE021Taper = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE021Taper[k] = new float [NY-yEnd1];
  bxE021Taper[k] = new float [NY-yEnd1];
}
azC100TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bzC100TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC100TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  bzC100TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC100TaperyPHxPH[k][j] = new float [NX-xEnd1];
    bzC100TaperyPHxPH[k][j] = new float [NX-xEnd1];
  }
}
azY0Taper = new float [_nYmin];
bzY0Taper = new float [_nYmin];
azC011TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC011TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC011TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC011TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC011TaperzPHyPH[k][j] = new float [_nXmin];
    bzC011TaperzPHyPH[k][j] = new float [_nXmin];
  }
}
ayC000TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
byC000TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC000TaperyPH[k] = new float*__restrict__ [_nYmin];
  byC000TaperyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC000TaperyPH[k][j] = new float [_nXmin];
    byC000TaperyPH[k][j] = new float [_nXmin];
  }
}
axC110TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
bxC110TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC110TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC110TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC110TaperzPH[k][j] = new float [NX-xEnd1];
    bxC110TaperzPH[k][j] = new float [NX-xEnd1];
  }
}
axE021TaperzPH = new float*__restrict__ [_nZmin];
bxE021TaperzPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE021TaperzPH[k] = new float [NY-yEnd1];
  bxE021TaperzPH[k] = new float [NY-yEnd1];
}
axY1TaperPH = new float [NY-yEnd1];
bxY1TaperPH = new float [NY-yEnd1];
ayC111TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC111TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC111TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  byC111TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC111TaperxPH[k][j] = new float [NX-xEnd1];
    byC111TaperxPH[k][j] = new float [NX-xEnd1];
  }
}
axE012TaperzPH = new float*__restrict__ [NZ-zEnd1];
bxE012TaperzPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE012TaperzPH[k] = new float [_nYmin];
  bxE012TaperzPH[k] = new float [_nYmin];
}
axE101TaperxPH = new float*__restrict__ [_nZmin];
bxE101TaperxPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE101TaperxPH[k] = new float [_nXmin];
  bxE101TaperxPH[k] = new float [_nXmin];
}
axE202TaperzPHxPH = new float*__restrict__ [NZ-zEnd1];
bxE202TaperzPHxPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE202TaperzPHxPH[k] = new float [NX-xEnd1];
  bxE202TaperzPHxPH[k] = new float [NX-xEnd1];
}
azC101TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC101TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC101TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  bzC101TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC101TaperyPHxPH[k][j] = new float [NX-xEnd1];
    bzC101TaperyPHxPH[k][j] = new float [NX-xEnd1];
  }
}
ayE021TaperzPH = new float*__restrict__ [_nZmin];
byE021TaperzPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE021TaperzPH[k] = new float [NY-yEnd1];
  byE021TaperzPH[k] = new float [NY-yEnd1];
}
axC101TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC101TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC101TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  bxC101TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC101TaperzPHyPH[k][j] = new float [NX-xEnd1];
    bxC101TaperzPHyPH[k][j] = new float [NX-xEnd1];
  }
}
ayE110TaperxPH = new float*__restrict__ [_nYmin];
byE110TaperxPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  ayE110TaperxPH[j] = new float [_nXmin];
  byE110TaperxPH[j] = new float [_nXmin];
}
azC101TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC101TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC101TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  bzC101TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC101TaperzPHxPH[k][j] = new float [NX-xEnd1];
    bzC101TaperzPHxPH[k][j] = new float [NX-xEnd1];
  }
}
axE102TaperzPHxPH = new float*__restrict__ [NZ-zEnd1];
bxE102TaperzPHxPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE102TaperzPHxPH[k] = new float [_nXmin];
  bxE102TaperzPHxPH[k] = new float [_nXmin];
}
azX1TaperPH = new float [NX-xEnd1];
bzX1TaperPH = new float [NX-xEnd1];
azE110TaperyPHxPH = new float*__restrict__ [_nYmin];
bzE110TaperyPHxPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  azE110TaperyPHxPH[j] = new float [_nXmin];
  bzE110TaperyPHxPH[j] = new float [_nXmin];
}
ayC101TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC101TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC101TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  byC101TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC101TaperzPHxPH[k][j] = new float [NX-xEnd1];
    byC101TaperzPHxPH[k][j] = new float [NX-xEnd1];
  }
}
ayC110TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
byC110TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC110TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  byC110TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC110TaperzPH[k][j] = new float [NX-xEnd1];
    byC110TaperzPH[k][j] = new float [NX-xEnd1];
  }
}
ayC010TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
byC010TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC010TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  byC010TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC010TaperzPHxPH[k][j] = new float [_nXmin];
    byC010TaperzPHxPH[k][j] = new float [_nXmin];
  }
}
axC100TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bxC100TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC100TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  bxC100TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC100TaperzPHxPH[k][j] = new float [NX-xEnd1];
    bxC100TaperzPHxPH[k][j] = new float [NX-xEnd1];
  }
}
axC111TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC111TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC111TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC111TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC111TaperzPH[k][j] = new float [NX-xEnd1];
    bxC111TaperzPH[k][j] = new float [NX-xEnd1];
  }
}
ayE201TaperxPH = new float*__restrict__ [_nZmin];
byE201TaperxPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE201TaperxPH[k] = new float [NX-xEnd1];
  byE201TaperxPH[k] = new float [NX-xEnd1];
}
axE220TaperxPH = new float*__restrict__ [NY-yEnd1];
bxE220TaperxPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  axE220TaperxPH[j] = new float [NX-xEnd1];
  bxE220TaperxPH[j] = new float [NX-xEnd1];
}
axC101TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC101TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC101TaperyPH[k] = new float*__restrict__ [_nYmin];
  bxC101TaperyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC101TaperyPH[k][j] = new float [NX-xEnd1];
    bxC101TaperyPH[k][j] = new float [NX-xEnd1];
  }
}
azC010TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bzC010TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC010TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC010TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC010TaperyPHxPH[k][j] = new float [_nXmin];
    bzC010TaperyPHxPH[k][j] = new float [_nXmin];
  }
}
azC011TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC011TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC011TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC011TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC011TaperyPH[k][j] = new float [_nXmin];
    bzC011TaperyPH[k][j] = new float [_nXmin];
  }
}
axC100TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bxC100TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC100TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  bxC100TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC100TaperyPHxPH[k][j] = new float [NX-xEnd1];
    bxC100TaperyPHxPH[k][j] = new float [NX-xEnd1];
  }
}
azX0Taper = new float [_nXmin];
bzX0Taper = new float [_nXmin];
azC010TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
bzC010TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC010TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC010TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC010TaperxPH[k][j] = new float [_nXmin];
    bzC010TaperxPH[k][j] = new float [_nXmin];
  }
}
axC001TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC001TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC001TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  bxC001TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC001TaperzPHxPH[k][j] = new float [_nXmin];
    bxC001TaperzPHxPH[k][j] = new float [_nXmin];
  }
}
ayC111TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC111TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC111TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  byC111TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC111TaperyPHxPH[k][j] = new float [NX-xEnd1];
    byC111TaperyPHxPH[k][j] = new float [NX-xEnd1];
  }
}
azE201TaperzPHxPH = new float*__restrict__ [_nZmin];
bzE201TaperzPHxPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE201TaperzPHxPH[k] = new float [NX-xEnd1];
  bzE201TaperzPHxPH[k] = new float [NX-xEnd1];
}
ayE210TaperyPHxPH = new float*__restrict__ [_nYmin];
byE210TaperyPHxPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  ayE210TaperyPHxPH[j] = new float [NX-xEnd1];
  byE210TaperyPHxPH[j] = new float [NX-xEnd1];
}
axC010TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
bxC010TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC010TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC010TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC010TaperxPH[k][j] = new float [_nXmin];
    bxC010TaperxPH[k][j] = new float [_nXmin];
  }
}
axC111TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC111TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC111TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC111TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC111TaperzPHxPH[k][j] = new float [NX-xEnd1];
    bxC111TaperzPHxPH[k][j] = new float [NX-xEnd1];
  }
}
axE201Taper = new float*__restrict__ [_nZmin];
bxE201Taper = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE201Taper[k] = new float [NX-xEnd1];
  bxE201Taper[k] = new float [NX-xEnd1];
}
azC110TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
bzC110TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC110TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC110TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC110TaperyPH[k][j] = new float [NX-xEnd1];
    bzC110TaperyPH[k][j] = new float [NX-xEnd1];
  }
}
azC111TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC111TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC111TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC111TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC111TaperyPHxPH[k][j] = new float [NX-xEnd1];
    bzC111TaperyPHxPH[k][j] = new float [NX-xEnd1];
  }
}
ayC100TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
byC100TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC100TaperxPH[k] = new float*__restrict__ [_nYmin];
  byC100TaperxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC100TaperxPH[k][j] = new float [NX-xEnd1];
    byC100TaperxPH[k][j] = new float [NX-xEnd1];
  }
}
ayE011TaperyPH = new float*__restrict__ [_nZmin];
byE011TaperyPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE011TaperyPH[k] = new float [_nYmin];
  byE011TaperyPH[k] = new float [_nYmin];
}
azE110Taper = new float*__restrict__ [_nYmin];
bzE110Taper = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  azE110Taper[j] = new float [_nXmin];
  bzE110Taper[j] = new float [_nXmin];
}
axC111TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC111TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC111TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC111TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC111TaperzPHyPH[k][j] = new float [NX-xEnd1];
    bxC111TaperzPHyPH[k][j] = new float [NX-xEnd1];
  }
}
azE101Taper = new float*__restrict__ [_nZmin];
bzE101Taper = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE101Taper[k] = new float [_nXmin];
  bzE101Taper[k] = new float [_nXmin];
}
axZ1TaperPH = new float [NZ-zEnd1];
bxZ1TaperPH = new float [NZ-zEnd1];
axE011TaperzPHyPH = new float*__restrict__ [_nZmin];
bxE011TaperzPHyPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE011TaperzPHyPH[k] = new float [_nYmin];
  bxE011TaperzPHyPH[k] = new float [_nYmin];
}
axE202TaperzPH = new float*__restrict__ [NZ-zEnd1];
bxE202TaperzPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE202TaperzPH[k] = new float [NX-xEnd1];
  bxE202TaperzPH[k] = new float [NX-xEnd1];
}
ayC100TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
byC100TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC100TaperyPH[k] = new float*__restrict__ [_nYmin];
  byC100TaperyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC100TaperyPH[k][j] = new float [NX-xEnd1];
    byC100TaperyPH[k][j] = new float [NX-xEnd1];
  }
}
azE101TaperzPH = new float*__restrict__ [_nZmin];
bzE101TaperzPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE101TaperzPH[k] = new float [_nXmin];
  bzE101TaperzPH[k] = new float [_nXmin];
}
ayE210TaperyPH = new float*__restrict__ [_nYmin];
byE210TaperyPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  ayE210TaperyPH[j] = new float [NX-xEnd1];
  byE210TaperyPH[j] = new float [NX-xEnd1];
}
axE120Taper = new float*__restrict__ [NY-yEnd1];
bxE120Taper = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  axE120Taper[j] = new float [_nXmin];
  bxE120Taper[j] = new float [_nXmin];
}
azE102Taper = new float*__restrict__ [NZ-zEnd1];
bzE102Taper = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE102Taper[k] = new float [_nXmin];
  bzE102Taper[k] = new float [_nXmin];
}
ayC011TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC011TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC011TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  byC011TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC011TaperzPHxPH[k][j] = new float [_nXmin];
    byC011TaperzPHxPH[k][j] = new float [_nXmin];
  }
}
azC001TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC001TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC001TaperyPH[k] = new float*__restrict__ [_nYmin];
  bzC001TaperyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC001TaperyPH[k][j] = new float [_nXmin];
    bzC001TaperyPH[k][j] = new float [_nXmin];
  }
}
axC000TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
bxC000TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC000TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  bxC000TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC000TaperzPHyPH[k][j] = new float [_nXmin];
    bxC000TaperzPHyPH[k][j] = new float [_nXmin];
  }
}
azE022Taper = new float*__restrict__ [NZ-zEnd1];
bzE022Taper = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE022Taper[k] = new float [NY-yEnd1];
  bzE022Taper[k] = new float [NY-yEnd1];
}
ayE012Taper = new float*__restrict__ [NZ-zEnd1];
byE012Taper = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE012Taper[k] = new float [_nYmin];
  byE012Taper[k] = new float [_nYmin];
}
axE220TaperyPH = new float*__restrict__ [NY-yEnd1];
bxE220TaperyPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  axE220TaperyPH[j] = new float [NX-xEnd1];
  bxE220TaperyPH[j] = new float [NX-xEnd1];
}
axC011TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC011TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC011TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC011TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC011TaperyPHxPH[k][j] = new float [_nXmin];
    bxC011TaperyPHxPH[k][j] = new float [_nXmin];
  }
}
ayE202Taper = new float*__restrict__ [NZ-zEnd1];
byE202Taper = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE202Taper[k] = new float [NX-xEnd1];
  byE202Taper[k] = new float [NX-xEnd1];
}
ayE220TaperyPH = new float*__restrict__ [NY-yEnd1];
byE220TaperyPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  ayE220TaperyPH[j] = new float [NX-xEnd1];
  byE220TaperyPH[j] = new float [NX-xEnd1];
}
axZ0Taper = new float [_nZmin];
bxZ0Taper = new float [_nZmin];
ayZ1Taper = new float [NZ-zEnd1];
byZ1Taper = new float [NZ-zEnd1];
ayE101TaperzPH = new float*__restrict__ [_nZmin];
byE101TaperzPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE101TaperzPH[k] = new float [_nXmin];
  byE101TaperzPH[k] = new float [_nXmin];
}
azE012TaperzPHyPH = new float*__restrict__ [NZ-zEnd1];
bzE012TaperzPHyPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE012TaperzPHyPH[k] = new float [_nYmin];
  bzE012TaperzPHyPH[k] = new float [_nYmin];
}
azC010TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
bzC010TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC010TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC010TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC010TaperzPHyPH[k][j] = new float [_nXmin];
    bzC010TaperzPHyPH[k][j] = new float [_nXmin];
  }
}
ayE021TaperzPHyPH = new float*__restrict__ [_nZmin];
byE021TaperzPHyPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE021TaperzPHyPH[k] = new float [NY-yEnd1];
  byE021TaperzPHyPH[k] = new float [NY-yEnd1];
}
azC000TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bzC000TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC000TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  bzC000TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC000TaperzPHxPH[k][j] = new float [_nXmin];
    bzC000TaperzPHxPH[k][j] = new float [_nXmin];
  }
}
azC111TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC111TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC111TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC111TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC111TaperyPH[k][j] = new float [NX-xEnd1];
    bzC111TaperyPH[k][j] = new float [NX-xEnd1];
  }
}
azE011TaperzPHyPH = new float*__restrict__ [_nZmin];
bzE011TaperzPHyPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE011TaperzPHyPH[k] = new float [_nYmin];
  bzE011TaperzPHyPH[k] = new float [_nYmin];
}
azC000TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
bzC000TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC000TaperyPH[k] = new float*__restrict__ [_nYmin];
  bzC000TaperyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC000TaperyPH[k][j] = new float [_nXmin];
    bzC000TaperyPH[k][j] = new float [_nXmin];
  }
}
axZ1Taper = new float [NZ-zEnd1];
bxZ1Taper = new float [NZ-zEnd1];
ayC011TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC011TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC011TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  byC011TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC011TaperxPH[k][j] = new float [_nXmin];
    byC011TaperxPH[k][j] = new float [_nXmin];
  }
}
ayE101Taper = new float*__restrict__ [_nZmin];
byE101Taper = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE101Taper[k] = new float [_nXmin];
  byE101Taper[k] = new float [_nXmin];
}
ayZ1TaperPH = new float [NZ-zEnd1];
byZ1TaperPH = new float [NZ-zEnd1];
ayC110TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
byC110TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC110TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  byC110TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC110TaperxPH[k][j] = new float [NX-xEnd1];
    byC110TaperxPH[k][j] = new float [NX-xEnd1];
  }
}
axC001TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC001TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC001TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  bxC001TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC001TaperyPHxPH[k][j] = new float [_nXmin];
    bxC001TaperyPHxPH[k][j] = new float [_nXmin];
  }
}
ayE110Taper = new float*__restrict__ [_nYmin];
byE110Taper = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  ayE110Taper[j] = new float [_nXmin];
  byE110Taper[j] = new float [_nXmin];
}
ayC110TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
byC110TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC110TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  byC110TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC110TaperyPH[k][j] = new float [NX-xEnd1];
    byC110TaperyPH[k][j] = new float [NX-xEnd1];
  }
}
azE120TaperxPH = new float*__restrict__ [NY-yEnd1];
bzE120TaperxPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  azE120TaperxPH[j] = new float [_nXmin];
  bzE120TaperxPH[j] = new float [_nXmin];
}
azC011TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC011TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC011TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC011TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC011TaperzPHxPH[k][j] = new float [_nXmin];
    bzC011TaperzPHxPH[k][j] = new float [_nXmin];
  }
}
ayC010TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
byC010TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC010TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  byC010TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC010TaperyPH[k][j] = new float [_nXmin];
    byC010TaperyPH[k][j] = new float [_nXmin];
  }
}
axC101TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC101TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC101TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  bxC101TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC101TaperyPHxPH[k][j] = new float [NX-xEnd1];
    bxC101TaperyPHxPH[k][j] = new float [NX-xEnd1];
  }
}
azC100TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
bzC100TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC100TaperxPH[k] = new float*__restrict__ [_nYmin];
  bzC100TaperxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC100TaperxPH[k][j] = new float [NX-xEnd1];
    bzC100TaperxPH[k][j] = new float [NX-xEnd1];
  }
}
ayC011TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC011TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC011TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  byC011TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC011TaperyPHxPH[k][j] = new float [_nXmin];
    byC011TaperyPHxPH[k][j] = new float [_nXmin];
  }
}
azE101TaperzPHxPH = new float*__restrict__ [_nZmin];
bzE101TaperzPHxPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE101TaperzPHxPH[k] = new float [_nXmin];
  bzE101TaperzPHxPH[k] = new float [_nXmin];
}
ayC011Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC011Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC011Taper[k] = new float*__restrict__ [NY-yEnd1];
  byC011Taper[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC011Taper[k][j] = new float [_nXmin];
    byC011Taper[k][j] = new float [_nXmin];
  }
}
axE110TaperyPH = new float*__restrict__ [_nYmin];
bxE110TaperyPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  axE110TaperyPH[j] = new float [_nXmin];
  bxE110TaperyPH[j] = new float [_nXmin];
}
ayE022Taper = new float*__restrict__ [NZ-zEnd1];
byE022Taper = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE022Taper[k] = new float [NY-yEnd1];
  byE022Taper[k] = new float [NY-yEnd1];
}
ayE012TaperyPH = new float*__restrict__ [NZ-zEnd1];
byE012TaperyPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE012TaperyPH[k] = new float [_nYmin];
  byE012TaperyPH[k] = new float [_nYmin];
}
axE022Taper = new float*__restrict__ [NZ-zEnd1];
bxE022Taper = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE022Taper[k] = new float [NY-yEnd1];
  bxE022Taper[k] = new float [NY-yEnd1];
}
ayC111TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC111TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC111TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  byC111TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC111TaperyPH[k][j] = new float [NX-xEnd1];
    byC111TaperyPH[k][j] = new float [NX-xEnd1];
  }
}
ayC100TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
byC100TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC100TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  byC100TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC100TaperyPHxPH[k][j] = new float [NX-xEnd1];
    byC100TaperyPHxPH[k][j] = new float [NX-xEnd1];
  }
}
ayE202TaperzPH = new float*__restrict__ [NZ-zEnd1];
byE202TaperzPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE202TaperzPH[k] = new float [NX-xEnd1];
  byE202TaperzPH[k] = new float [NX-xEnd1];
}
axE202TaperxPH = new float*__restrict__ [NZ-zEnd1];
bxE202TaperxPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE202TaperxPH[k] = new float [NX-xEnd1];
  bxE202TaperxPH[k] = new float [NX-xEnd1];
}
axE202Taper = new float*__restrict__ [NZ-zEnd1];
bxE202Taper = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE202Taper[k] = new float [NX-xEnd1];
  bxE202Taper[k] = new float [NX-xEnd1];
}
azC100TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
bzC100TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC100TaperyPH[k] = new float*__restrict__ [_nYmin];
  bzC100TaperyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC100TaperyPH[k][j] = new float [NX-xEnd1];
    bzC100TaperyPH[k][j] = new float [NX-xEnd1];
  }
}
ayC110TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
byC110TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC110TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  byC110TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC110TaperzPHxPH[k][j] = new float [NX-xEnd1];
    byC110TaperzPHxPH[k][j] = new float [NX-xEnd1];
  }
}
azE022TaperzPHyPH = new float*__restrict__ [NZ-zEnd1];
bzE022TaperzPHyPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE022TaperzPHyPH[k] = new float [NY-yEnd1];
  bzE022TaperzPHyPH[k] = new float [NY-yEnd1];
}
azC010TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bzC010TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC010TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC010TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC010TaperzPHxPH[k][j] = new float [_nXmin];
    bzC010TaperzPHxPH[k][j] = new float [_nXmin];
  }
}
axC110Taper = new float*__restrict__*__restrict__ [_nZmin];
bxC110Taper = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC110Taper[k] = new float*__restrict__ [NY-yEnd1];
  bxC110Taper[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC110Taper[k][j] = new float [NX-xEnd1];
    bxC110Taper[k][j] = new float [NX-xEnd1];
  }
}
azE201TaperzPH = new float*__restrict__ [_nZmin];
bzE201TaperzPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE201TaperzPH[k] = new float [NX-xEnd1];
  bzE201TaperzPH[k] = new float [NX-xEnd1];
}
azE102TaperzPHxPH = new float*__restrict__ [NZ-zEnd1];
bzE102TaperzPHxPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE102TaperzPHxPH[k] = new float [_nXmin];
  bzE102TaperzPHxPH[k] = new float [_nXmin];
}
azC100TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
bzC100TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC100TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  bzC100TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC100TaperzPHyPH[k][j] = new float [NX-xEnd1];
    bzC100TaperzPHyPH[k][j] = new float [NX-xEnd1];
  }
}
axC000Taper = new float*__restrict__*__restrict__ [_nZmin];
bxC000Taper = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC000Taper[k] = new float*__restrict__ [_nYmin];
  bxC000Taper[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC000Taper[k][j] = new float [_nXmin];
    bxC000Taper[k][j] = new float [_nXmin];
  }
}
azE202Taper = new float*__restrict__ [NZ-zEnd1];
bzE202Taper = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE202Taper[k] = new float [NX-xEnd1];
  bzE202Taper[k] = new float [NX-xEnd1];
}
azC101TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC101TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC101TaperxPH[k] = new float*__restrict__ [_nYmin];
  bzC101TaperxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC101TaperxPH[k][j] = new float [NX-xEnd1];
    bzC101TaperxPH[k][j] = new float [NX-xEnd1];
  }
}
axC111TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC111TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC111TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC111TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC111TaperyPH[k][j] = new float [NX-xEnd1];
    bxC111TaperyPH[k][j] = new float [NX-xEnd1];
  }
}
axC000TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
bxC000TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC000TaperzPH[k] = new float*__restrict__ [_nYmin];
  bxC000TaperzPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC000TaperzPH[k][j] = new float [_nXmin];
    bxC000TaperzPH[k][j] = new float [_nXmin];
  }
}
axE101Taper = new float*__restrict__ [_nZmin];
bxE101Taper = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE101Taper[k] = new float [_nXmin];
  bxE101Taper[k] = new float [_nXmin];
}
azC111Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC111Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC111Taper[k] = new float*__restrict__ [NY-yEnd1];
  bzC111Taper[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC111Taper[k][j] = new float [NX-xEnd1];
    bzC111Taper[k][j] = new float [NX-xEnd1];
  }
}
ayC001TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC001TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC001TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  byC001TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC001TaperyPHxPH[k][j] = new float [_nXmin];
    byC001TaperyPHxPH[k][j] = new float [_nXmin];
  }
}
axE102TaperxPH = new float*__restrict__ [NZ-zEnd1];
bxE102TaperxPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE102TaperxPH[k] = new float [_nXmin];
  bxE102TaperxPH[k] = new float [_nXmin];
}
ayE012TaperzPH = new float*__restrict__ [NZ-zEnd1];
byE012TaperzPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE012TaperzPH[k] = new float [_nYmin];
  byE012TaperzPH[k] = new float [_nYmin];
}
axE022TaperzPHyPH = new float*__restrict__ [NZ-zEnd1];
bxE022TaperzPHyPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE022TaperzPHyPH[k] = new float [NY-yEnd1];
  bxE022TaperzPHyPH[k] = new float [NY-yEnd1];
}
azE110TaperyPH = new float*__restrict__ [_nYmin];
bzE110TaperyPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  azE110TaperyPH[j] = new float [_nXmin];
  bzE110TaperyPH[j] = new float [_nXmin];
}
axC101Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC101Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC101Taper[k] = new float*__restrict__ [_nYmin];
  bxC101Taper[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC101Taper[k][j] = new float [NX-xEnd1];
    bxC101Taper[k][j] = new float [NX-xEnd1];
  }
}
azC001TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC001TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC001TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  bzC001TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC001TaperzPHxPH[k][j] = new float [_nXmin];
    bzC001TaperzPHxPH[k][j] = new float [_nXmin];
  }
}
axE011TaperyPH = new float*__restrict__ [_nZmin];
bxE011TaperyPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE011TaperyPH[k] = new float [_nYmin];
  bxE011TaperyPH[k] = new float [_nYmin];
}
azC011TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC011TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC011TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC011TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC011TaperxPH[k][j] = new float [_nXmin];
    bzC011TaperxPH[k][j] = new float [_nXmin];
  }
}
axE110TaperxPH = new float*__restrict__ [_nYmin];
bxE110TaperxPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  axE110TaperxPH[j] = new float [_nXmin];
  bxE110TaperxPH[j] = new float [_nXmin];
}
axY1Taper = new float [NY-yEnd1];
bxY1Taper = new float [NY-yEnd1];
ayC101TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC101TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC101TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  byC101TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC101TaperyPHxPH[k][j] = new float [NX-xEnd1];
    byC101TaperyPHxPH[k][j] = new float [NX-xEnd1];
  }
}
axC011TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC011TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC011TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC011TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC011TaperzPH[k][j] = new float [_nXmin];
    bxC011TaperzPH[k][j] = new float [_nXmin];
  }
}
ayC110TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
byC110TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC110TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  byC110TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC110TaperzPHyPH[k][j] = new float [NX-xEnd1];
    byC110TaperzPHyPH[k][j] = new float [NX-xEnd1];
  }
}
ayC011TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC011TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC011TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  byC011TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC011TaperzPHyPH[k][j] = new float [_nXmin];
    byC011TaperzPHyPH[k][j] = new float [_nXmin];
  }
}
ayE201TaperzPHxPH = new float*__restrict__ [_nZmin];
byE201TaperzPHxPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE201TaperzPHxPH[k] = new float [NX-xEnd1];
  byE201TaperzPHxPH[k] = new float [NX-xEnd1];
}
azC100Taper = new float*__restrict__*__restrict__ [_nZmin];
bzC100Taper = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC100Taper[k] = new float*__restrict__ [_nYmin];
  bzC100Taper[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC100Taper[k][j] = new float [NX-xEnd1];
    bzC100Taper[k][j] = new float [NX-xEnd1];
  }
}
axC001TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC001TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC001TaperzPH[k] = new float*__restrict__ [_nYmin];
  bxC001TaperzPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC001TaperzPH[k][j] = new float [_nXmin];
    bxC001TaperzPH[k][j] = new float [_nXmin];
  }
}
azC011Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC011Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC011Taper[k] = new float*__restrict__ [NY-yEnd1];
  bzC011Taper[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC011Taper[k][j] = new float [_nXmin];
    bzC011Taper[k][j] = new float [_nXmin];
  }
}
azC010TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
bzC010TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC010TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC010TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC010TaperzPH[k][j] = new float [_nXmin];
    bzC010TaperzPH[k][j] = new float [_nXmin];
  }
}
ayE120TaperyPHxPH = new float*__restrict__ [NY-yEnd1];
byE120TaperyPHxPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  ayE120TaperyPHxPH[j] = new float [_nXmin];
  byE120TaperyPHxPH[j] = new float [_nXmin];
}
ayC101TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC101TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC101TaperzPH[k] = new float*__restrict__ [_nYmin];
  byC101TaperzPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC101TaperzPH[k][j] = new float [NX-xEnd1];
    byC101TaperzPH[k][j] = new float [NX-xEnd1];
  }
}
azC000Taper = new float*__restrict__*__restrict__ [_nZmin];
bzC000Taper = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC000Taper[k] = new float*__restrict__ [_nYmin];
  bzC000Taper[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC000Taper[k][j] = new float [_nXmin];
    bzC000Taper[k][j] = new float [_nXmin];
  }
}
axC010TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
bxC010TaperzPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC010TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC010TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC010TaperzPH[k][j] = new float [_nXmin];
    bxC010TaperzPH[k][j] = new float [_nXmin];
  }
}
ayC111Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC111Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC111Taper[k] = new float*__restrict__ [NY-yEnd1];
  byC111Taper[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC111Taper[k][j] = new float [NX-xEnd1];
    byC111Taper[k][j] = new float [NX-xEnd1];
  }
}
axE220Taper = new float*__restrict__ [NY-yEnd1];
bxE220Taper = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  axE220Taper[j] = new float [NX-xEnd1];
  bxE220Taper[j] = new float [NX-xEnd1];
}
axC100TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
bxC100TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC100TaperxPH[k] = new float*__restrict__ [_nYmin];
  bxC100TaperxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC100TaperxPH[k][j] = new float [NX-xEnd1];
    bxC100TaperxPH[k][j] = new float [NX-xEnd1];
  }
}
azC001TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC001TaperxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC001TaperxPH[k] = new float*__restrict__ [_nYmin];
  bzC001TaperxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC001TaperxPH[k][j] = new float [_nXmin];
    bzC001TaperxPH[k][j] = new float [_nXmin];
  }
}
azE012TaperyPH = new float*__restrict__ [NZ-zEnd1];
bzE012TaperyPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE012TaperyPH[k] = new float [_nYmin];
  bzE012TaperyPH[k] = new float [_nYmin];
}
azY1Taper = new float [NY-yEnd1];
bzY1Taper = new float [NY-yEnd1];
azE210Taper = new float*__restrict__ [_nYmin];
bzE210Taper = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  azE210Taper[j] = new float [NX-xEnd1];
  bzE210Taper[j] = new float [NX-xEnd1];
}
axE210TaperyPHxPH = new float*__restrict__ [_nYmin];
bxE210TaperyPHxPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  axE210TaperyPHxPH[j] = new float [NX-xEnd1];
  bxE210TaperyPHxPH[j] = new float [NX-xEnd1];
}
ayC001TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC001TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC001TaperyPH[k] = new float*__restrict__ [_nYmin];
  byC001TaperyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC001TaperyPH[k][j] = new float [_nXmin];
    byC001TaperyPH[k][j] = new float [_nXmin];
  }
}
axY0TaperPH = new float [_nYmin];
bxY0TaperPH = new float [_nYmin];
ayE110TaperyPH = new float*__restrict__ [_nYmin];
byE110TaperyPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  ayE110TaperyPH[j] = new float [_nXmin];
  byE110TaperyPH[j] = new float [_nXmin];
}
axC010TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bxC010TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC010TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC010TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC010TaperyPHxPH[k][j] = new float [_nXmin];
    bxC010TaperyPHxPH[k][j] = new float [_nXmin];
  }
}
ayC110TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
byC110TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC110TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  byC110TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC110TaperyPHxPH[k][j] = new float [NX-xEnd1];
    byC110TaperyPHxPH[k][j] = new float [NX-xEnd1];
  }
}
ayE101TaperxPH = new float*__restrict__ [_nZmin];
byE101TaperxPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE101TaperxPH[k] = new float [_nXmin];
  byE101TaperxPH[k] = new float [_nXmin];
}
axC111Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC111Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC111Taper[k] = new float*__restrict__ [NY-yEnd1];
  bxC111Taper[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC111Taper[k][j] = new float [NX-xEnd1];
    bxC111Taper[k][j] = new float [NX-xEnd1];
  }
}
ayC011TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC011TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC011TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  byC011TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC011TaperyPH[k][j] = new float [_nXmin];
    byC011TaperyPH[k][j] = new float [_nXmin];
  }
}
axE201TaperxPH = new float*__restrict__ [_nZmin];
bxE201TaperxPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE201TaperxPH[k] = new float [NX-xEnd1];
  bxE201TaperxPH[k] = new float [NX-xEnd1];
}
axC001TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC001TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC001TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  bxC001TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC001TaperzPHyPH[k][j] = new float [_nXmin];
    bxC001TaperzPHyPH[k][j] = new float [_nXmin];
  }
}
axC010Taper = new float*__restrict__*__restrict__ [_nZmin];
bxC010Taper = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC010Taper[k] = new float*__restrict__ [NY-yEnd1];
  bxC010Taper[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC010Taper[k][j] = new float [_nXmin];
    bxC010Taper[k][j] = new float [_nXmin];
  }
}
azE220TaperxPH = new float*__restrict__ [NY-yEnd1];
bzE220TaperxPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  azE220TaperxPH[j] = new float [NX-xEnd1];
  bzE220TaperxPH[j] = new float [NX-xEnd1];
}
azC011TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC011TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC011TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC011TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC011TaperzPH[k][j] = new float [_nXmin];
    bzC011TaperzPH[k][j] = new float [_nXmin];
  }
}
azE022TaperzPH = new float*__restrict__ [NZ-zEnd1];
bzE022TaperzPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE022TaperzPH[k] = new float [NY-yEnd1];
  bzE022TaperzPH[k] = new float [NY-yEnd1];
}
ayC000Taper = new float*__restrict__*__restrict__ [_nZmin];
byC000Taper = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC000Taper[k] = new float*__restrict__ [_nYmin];
  byC000Taper[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC000Taper[k][j] = new float [_nXmin];
    byC000Taper[k][j] = new float [_nXmin];
  }
}
ayC100TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
byC100TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC100TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  byC100TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC100TaperzPHxPH[k][j] = new float [NX-xEnd1];
    byC100TaperzPHxPH[k][j] = new float [NX-xEnd1];
  }
}
axC001Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC001Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC001Taper[k] = new float*__restrict__ [_nYmin];
  bxC001Taper[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC001Taper[k][j] = new float [_nXmin];
    bxC001Taper[k][j] = new float [_nXmin];
  }
}
ayC111TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC111TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC111TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  byC111TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC111TaperzPHxPH[k][j] = new float [NX-xEnd1];
    byC111TaperzPHxPH[k][j] = new float [NX-xEnd1];
  }
}
ayC000TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
byC000TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC000TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  byC000TaperyPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC000TaperyPHxPH[k][j] = new float [_nXmin];
    byC000TaperyPHxPH[k][j] = new float [_nXmin];
  }
}
azE202TaperxPH = new float*__restrict__ [NZ-zEnd1];
bzE202TaperxPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE202TaperxPH[k] = new float [NX-xEnd1];
  bzE202TaperxPH[k] = new float [NX-xEnd1];
}
ayX1Taper = new float [NX-xEnd1];
byX1Taper = new float [NX-xEnd1];
azC110TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
bzC110TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC110TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC110TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC110TaperzPHyPH[k][j] = new float [NX-xEnd1];
    bzC110TaperzPHyPH[k][j] = new float [NX-xEnd1];
  }
}
azC111TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC111TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC111TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC111TaperzPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC111TaperzPH[k][j] = new float [NX-xEnd1];
    bzC111TaperzPH[k][j] = new float [NX-xEnd1];
  }
}
azC101Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC101Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC101Taper[k] = new float*__restrict__ [_nYmin];
  bzC101Taper[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC101Taper[k][j] = new float [NX-xEnd1];
    bzC101Taper[k][j] = new float [NX-xEnd1];
  }
}
ayC100TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
byC100TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC100TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  byC100TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC100TaperzPHyPH[k][j] = new float [NX-xEnd1];
    byC100TaperzPHyPH[k][j] = new float [NX-xEnd1];
  }
}
azX1Taper = new float [NX-xEnd1];
bzX1Taper = new float [NX-xEnd1];
axC110TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bxC110TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC110TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC110TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC110TaperyPHxPH[k][j] = new float [NX-xEnd1];
    bxC110TaperyPHxPH[k][j] = new float [NX-xEnd1];
  }
}
axC110TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
bxC110TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC110TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC110TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC110TaperyPH[k][j] = new float [NX-xEnd1];
    bxC110TaperyPH[k][j] = new float [NX-xEnd1];
  }
}
azE210TaperxPH = new float*__restrict__ [_nYmin];
bzE210TaperxPH = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  azE210TaperxPH[j] = new float [NX-xEnd1];
  bzE210TaperxPH[j] = new float [NX-xEnd1];
}
axC011TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC011TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC011TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC011TaperzPHyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC011TaperzPHyPH[k][j] = new float [_nXmin];
    bxC011TaperzPHyPH[k][j] = new float [_nXmin];
  }
}
ayC000TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
byC000TaperzPHyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC000TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  byC000TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC000TaperzPHyPH[k][j] = new float [_nXmin];
    byC000TaperzPHyPH[k][j] = new float [_nXmin];
  }
}
azE120TaperyPHxPH = new float*__restrict__ [NY-yEnd1];
bzE120TaperyPHxPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  azE120TaperyPHxPH[j] = new float [_nXmin];
  bzE120TaperyPHxPH[j] = new float [_nXmin];
}
azE202TaperzPH = new float*__restrict__ [NZ-zEnd1];
bzE202TaperzPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE202TaperzPH[k] = new float [NX-xEnd1];
  bzE202TaperzPH[k] = new float [NX-xEnd1];
}
axC101TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC101TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC101TaperzPH[k] = new float*__restrict__ [_nYmin];
  bxC101TaperzPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC101TaperzPH[k][j] = new float [NX-xEnd1];
    bxC101TaperzPH[k][j] = new float [NX-xEnd1];
  }
}
azE102TaperxPH = new float*__restrict__ [NZ-zEnd1];
bzE102TaperxPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE102TaperxPH[k] = new float [_nXmin];
  bzE102TaperxPH[k] = new float [_nXmin];
}
ayC001TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC001TaperzPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC001TaperzPH[k] = new float*__restrict__ [_nYmin];
  byC001TaperzPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC001TaperzPH[k][j] = new float [_nXmin];
    byC001TaperzPH[k][j] = new float [_nXmin];
  }
}
ayC000TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
byC000TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC000TaperxPH[k] = new float*__restrict__ [_nYmin];
  byC000TaperxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC000TaperxPH[k][j] = new float [_nXmin];
    byC000TaperxPH[k][j] = new float [_nXmin];
  }
}
azC110TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bzC110TaperyPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC110TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC110TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC110TaperyPHxPH[k][j] = new float [NX-xEnd1];
    bzC110TaperyPHxPH[k][j] = new float [NX-xEnd1];
  }
}
axE101TaperzPHxPH = new float*__restrict__ [_nZmin];
bxE101TaperzPHxPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE101TaperzPHxPH[k] = new float [_nXmin];
  bxE101TaperzPHxPH[k] = new float [_nXmin];
}
azE220TaperyPH = new float*__restrict__ [NY-yEnd1];
bzE220TaperyPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  azE220TaperyPH[j] = new float [NX-xEnd1];
  bzE220TaperyPH[j] = new float [NX-xEnd1];
}
ayE022TaperzPH = new float*__restrict__ [NZ-zEnd1];
byE022TaperzPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE022TaperzPH[k] = new float [NY-yEnd1];
  byE022TaperzPH[k] = new float [NY-yEnd1];
}
ayC110Taper = new float*__restrict__*__restrict__ [_nZmin];
byC110Taper = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC110Taper[k] = new float*__restrict__ [NY-yEnd1];
  byC110Taper[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC110Taper[k][j] = new float [NX-xEnd1];
    byC110Taper[k][j] = new float [NX-xEnd1];
  }
}
axE022TaperyPH = new float*__restrict__ [NZ-zEnd1];
bxE022TaperyPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE022TaperyPH[k] = new float [NY-yEnd1];
  bxE022TaperyPH[k] = new float [NY-yEnd1];
}
azE021Taper = new float*__restrict__ [_nZmin];
bzE021Taper = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE021Taper[k] = new float [NY-yEnd1];
  bzE021Taper[k] = new float [NY-yEnd1];
}
axC000TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
bxC000TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC000TaperxPH[k] = new float*__restrict__ [_nYmin];
  bxC000TaperxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC000TaperxPH[k][j] = new float [_nXmin];
    bxC000TaperxPH[k][j] = new float [_nXmin];
  }
}
axY0Taper = new float [_nYmin];
bxY0Taper = new float [_nYmin];
ayE022TaperzPHyPH = new float*__restrict__ [NZ-zEnd1];
byE022TaperzPHyPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE022TaperzPHyPH[k] = new float [NY-yEnd1];
  byE022TaperzPHyPH[k] = new float [NY-yEnd1];
}
ayE011Taper = new float*__restrict__ [_nZmin];
byE011Taper = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE011Taper[k] = new float [_nYmin];
  byE011Taper[k] = new float [_nYmin];
}
ayE102TaperzPH = new float*__restrict__ [NZ-zEnd1];
byE102TaperzPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE102TaperzPH[k] = new float [_nXmin];
  byE102TaperzPH[k] = new float [_nXmin];
}
axZ0TaperPH = new float [_nZmin];
bxZ0TaperPH = new float [_nZmin];
azC111TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC111TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC111TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC111TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC111TaperzPHxPH[k][j] = new float [NX-xEnd1];
    bzC111TaperzPHxPH[k][j] = new float [NX-xEnd1];
  }
}
axE120TaperxPH = new float*__restrict__ [NY-yEnd1];
bxE120TaperxPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  axE120TaperxPH[j] = new float [_nXmin];
  bxE120TaperxPH[j] = new float [_nXmin];
}
ayZ0Taper = new float [_nZmin];
byZ0Taper = new float [_nZmin];
axC000TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
bxC000TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axC000TaperyPH[k] = new float*__restrict__ [_nYmin];
  bxC000TaperyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    axC000TaperyPH[k][j] = new float [_nXmin];
    bxC000TaperyPH[k][j] = new float [_nXmin];
  }
}
azE012TaperzPH = new float*__restrict__ [NZ-zEnd1];
bzE012TaperzPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE012TaperzPH[k] = new float [_nYmin];
  bzE012TaperzPH[k] = new float [_nYmin];
}
azC101TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC101TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC101TaperyPH[k] = new float*__restrict__ [_nYmin];
  bzC101TaperyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC101TaperyPH[k][j] = new float [NX-xEnd1];
    bzC101TaperyPH[k][j] = new float [NX-xEnd1];
  }
}
azE202TaperzPHxPH = new float*__restrict__ [NZ-zEnd1];
bzE202TaperzPHxPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azE202TaperzPHxPH[k] = new float [NX-xEnd1];
  bzE202TaperzPHxPH[k] = new float [NX-xEnd1];
}
ayC001Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC001Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC001Taper[k] = new float*__restrict__ [_nYmin];
  byC001Taper[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC001Taper[k][j] = new float [_nXmin];
    byC001Taper[k][j] = new float [_nXmin];
  }
}
ayC101Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC101Taper = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC101Taper[k] = new float*__restrict__ [_nYmin];
  byC101Taper[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC101Taper[k][j] = new float [NX-xEnd1];
    byC101Taper[k][j] = new float [NX-xEnd1];
  }
}
azE011Taper = new float*__restrict__ [_nZmin];
bzE011Taper = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE011Taper[k] = new float [_nYmin];
  bzE011Taper[k] = new float [_nYmin];
}
ayC010TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
byC010TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC010TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  byC010TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    ayC010TaperxPH[k][j] = new float [_nXmin];
    byC010TaperxPH[k][j] = new float [_nXmin];
  }
}
azC110TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
bzC110TaperxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC110TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC110TaperxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC110TaperxPH[k][j] = new float [NX-xEnd1];
    bzC110TaperxPH[k][j] = new float [NX-xEnd1];
  }
}
axE012TaperzPHyPH = new float*__restrict__ [NZ-zEnd1];
bxE012TaperzPHyPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axE012TaperzPHyPH[k] = new float [_nYmin];
  bxE012TaperzPHyPH[k] = new float [_nYmin];
}
axE021TaperyPH = new float*__restrict__ [_nZmin];
bxE021TaperyPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE021TaperyPH[k] = new float [NY-yEnd1];
  bxE021TaperyPH[k] = new float [NY-yEnd1];
}
ayC000TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
byC000TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayC000TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  byC000TaperzPHxPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC000TaperzPHxPH[k][j] = new float [_nXmin];
    byC000TaperzPHxPH[k][j] = new float [_nXmin];
  }
}
axC011TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC011TaperzPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC011TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC011TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC011TaperzPHxPH[k][j] = new float [_nXmin];
    bxC011TaperzPHxPH[k][j] = new float [_nXmin];
  }
}
ayE011TaperzPH = new float*__restrict__ [_nZmin];
byE011TaperzPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE011TaperzPH[k] = new float [_nYmin];
  byE011TaperzPH[k] = new float [_nYmin];
}
ayE102TaperxPH = new float*__restrict__ [NZ-zEnd1];
byE102TaperxPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE102TaperxPH[k] = new float [_nXmin];
  byE102TaperxPH[k] = new float [_nXmin];
}
azX0TaperPH = new float [_nXmin];
bzX0TaperPH = new float [_nXmin];
ayE202TaperzPHxPH = new float*__restrict__ [NZ-zEnd1];
byE202TaperzPHxPH = new float*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayE202TaperzPHxPH[k] = new float [NX-xEnd1];
  byE202TaperzPHxPH[k] = new float [NX-xEnd1];
}
axC111TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC111TaperyPHxPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC111TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC111TaperyPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC111TaperyPHxPH[k][j] = new float [NX-xEnd1];
    bxC111TaperyPHxPH[k][j] = new float [NX-xEnd1];
  }
}
ayE220TaperxPH = new float*__restrict__ [NY-yEnd1];
byE220TaperxPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  ayE220TaperxPH[j] = new float [NX-xEnd1];
  byE220TaperxPH[j] = new float [NX-xEnd1];
}
azC110TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
bzC110TaperzPHxPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC110TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC110TaperzPHxPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC110TaperzPHxPH[k][j] = new float [NX-xEnd1];
    bzC110TaperzPHxPH[k][j] = new float [NX-xEnd1];
  }
}
ayE011TaperzPHyPH = new float*__restrict__ [_nZmin];
byE011TaperzPHyPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE011TaperzPHyPH[k] = new float [_nYmin];
  byE011TaperzPHyPH[k] = new float [_nYmin];
}
ayC101TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
byC101TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  ayC101TaperyPH[k] = new float*__restrict__ [_nYmin];
  byC101TaperyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    ayC101TaperyPH[k][j] = new float [NX-xEnd1];
    byC101TaperyPH[k][j] = new float [NX-xEnd1];
  }
}
ayE220TaperyPHxPH = new float*__restrict__ [NY-yEnd1];
byE220TaperyPHxPH = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  ayE220TaperyPHxPH[j] = new float [NX-xEnd1];
  byE220TaperyPHxPH[j] = new float [NX-xEnd1];
}
ayE220Taper = new float*__restrict__ [NY-yEnd1];
byE220Taper = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  ayE220Taper[j] = new float [NX-xEnd1];
  byE220Taper[j] = new float [NX-xEnd1];
}
azE021TaperyPH = new float*__restrict__ [_nZmin];
bzE021TaperyPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azE021TaperyPH[k] = new float [NY-yEnd1];
  bzE021TaperyPH[k] = new float [NY-yEnd1];
}
axE011TaperzPH = new float*__restrict__ [_nZmin];
bxE011TaperzPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE011TaperzPH[k] = new float [_nYmin];
  bxE011TaperzPH[k] = new float [_nYmin];
}
ayE201TaperzPH = new float*__restrict__ [_nZmin];
byE201TaperzPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  ayE201TaperzPH[k] = new float [NX-xEnd1];
  byE201TaperzPH[k] = new float [NX-xEnd1];
}
azC101TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bzC101TaperzPHyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  azC101TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  bzC101TaperzPHyPH[k] = new float*__restrict__ [_nYmin];
  for(int j=0;j<_nYmin;++j) {
    azC101TaperzPHyPH[k][j] = new float [NX-xEnd1];
    bzC101TaperzPHyPH[k][j] = new float [NX-xEnd1];
  }
}
axE110Taper = new float*__restrict__ [_nYmin];
bxE110Taper = new float*__restrict__ [_nYmin];
for(int j=0;j<_nYmin;++j) {
  axE110Taper[j] = new float [_nXmin];
  bxE110Taper[j] = new float [_nXmin];
}
ayE120Taper = new float*__restrict__ [NY-yEnd1];
byE120Taper = new float*__restrict__ [NY-yEnd1];
for(int j=0;j<NY-yEnd1;++j) {
  ayE120Taper[j] = new float [_nXmin];
  byE120Taper[j] = new float [_nXmin];
}
azC010TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
bzC010TaperyPH = new float*__restrict__*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  azC010TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  bzC010TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    azC010TaperyPH[k][j] = new float [_nXmin];
    bzC010TaperyPH[k][j] = new float [_nXmin];
  }
}
axE021TaperzPHyPH = new float*__restrict__ [_nZmin];
bxE021TaperzPHyPH = new float*__restrict__ [_nZmin];
for(int k=0;k<_nZmin;++k) {
  axE021TaperzPHyPH[k] = new float [NY-yEnd1];
  bxE021TaperzPHyPH[k] = new float [NY-yEnd1];
}
azY0TaperPH = new float [_nYmin];
bzY0TaperPH = new float [_nYmin];
axC011TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
bxC011TaperyPH = new float*__restrict__*__restrict__ [NZ-zEnd1];
for(int k=0;k<NZ-zEnd1;++k) {
  axC011TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  bxC011TaperyPH[k] = new float*__restrict__ [NY-yEnd1];
  for(int j=0;j<NY-yEnd1;++j) {
    axC011TaperyPH[k][j] = new float [_nXmin];
    bxC011TaperyPH[k][j] = new float [_nXmin];
  }
}


  
  //Initialize the memory variables to 0
  for(int i=0;i<NXY*_nZmin;i++) {vxxZmin[i]=vyyZmin[i]=vzzZmin[i]=pxZmin[i]=pyZmin[i]=pzZmin[i]=0.;}
  for(int i=0;i<NXY*_nZmax;i++) {vxxZmax[i]=vyyZmax[i]=vzzZmax[i]=pxZmax[i]=pyZmax[i]=pzZmax[i]=0.;}
  for(int i=0;i<NX*_nYmin*nzInt;i++) {vxxYmin[i]=vyyYmin[i]=vzzYmin[i]=pxYmin[i]=pyYmin[i]=pzYmin[i]=0.;}
  for(int i=0;i<NX*_nYmax*nzInt;i++) {vxxYmax[i]=vyyYmax[i]=vzzYmax[i]=pxYmax[i]=pyYmax[i]=pzYmax[i]=0.;}
  for(int i=0;i<_nXmin*nyInt*nzInt;i++) {vxxXmin[i]=vyyXmin[i]=vzzXmin[i]=pxXmin[i]=pyXmin[i]=pzXmin[i]=0.;}
  for(int i=0;i<_nXmax*nyInt*nzInt;i++) {vxxXmax[i]=vyyXmax[i]=vzzXmax[i]=pxXmax[i]=pyXmax[i]=pzXmax[i]=0.;}
  
  //Fill in the flank profile variables
  //fromStart and fromEnd are maximum at the min side or max side, respectively.  Appropriate values are filled in when outside the MPML layers; however, these values are not used.
  //sigma and kappa tapers vary quadratically.  sigma starts at zero at the interior edge of the MPML zone and increases to the flank sigma computed above at the flank of the model.  kappa starts at one at the interior and will decrease to the flank value at the edge.
  //alpha taper varies linearly from a maximum at the interior edge to zero at the flank
  
  //The X profile
  for(int i=0;i<NX;i++) {
    float fromStart = nXmin-i-procLim[0];
    float fromEnd = nXmax-globalNX+procLim[1]-(NX-1-i);
    if(fromStart>0.f) {
      xTaper[i] = xMinVal*sqr(fromStart/nXmin);
      xTaperPH[i] = xMinVal*sqr((fromStart-0.5)/nXmin);
      kxTaper[i] = 1.+(xMinKVal-1.)*sqr(fromStart/nXmin);
      kxTaperPH[i] = 1.+(xMinKVal-1.)*sqr((fromStart-0.5)/nXmin);
      alphaXTaper[i] = xMinAVal*(1.-fromStart/nXmin);
      alphaXTaperPH[i] = xMinAVal*(1.-(fromStart-0.5)/nXmin);
      bxTaper[i] = exp(-(xTaper[i]/kxTaper[i]+alphaXTaper[i])*dt);
      bxTaperPH[i] = exp(-(xTaperPH[i]/kxTaperPH[i]+alphaXTaperPH[i])*dt);
      axTaper[i] = xTaper[i]*(bxTaper[i]-1.)/(kxTaper[i]*(xTaper[i]+kxTaper[i]*alphaXTaper[i]));
      axTaperPH[i] = xTaperPH[i]*(bxTaperPH[i]-1.)/(kxTaperPH[i]*(xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
      //        bxTaper[i] *= dt;   //multiply by dt to save mults in loops
      //        bxTaperPH[i] *= dt;
      kxTaper[i] = 1./kxTaper[i]; //more easily handled
      kxTaperPH[i] = 1./kxTaperPH[i];
      //        xTaperFac[i] = xTaper[i]*xfac*sqr(fromStart/nXmin);
      //        xTaperPHFac[i] = xTaperPH[i]*xfac*sqr((fromStart-0.5)/nXmin);
      if(_surfaceBCMode==ELASTIC_FS_2_0) {
        xTaperRel[i] = sqr(fromStart/nXmin);
        xTaperRelPH[i] = sqr((fromStart-0.5)/nXmin);
      }
    } else if(fromEnd>-0.5f) {
      xTaper[i] = xMaxVal*sqr(fromEnd/nXmax);
      xTaperPH[i] = xMaxVal*sqr((fromEnd+0.5)/nXmax);
      kxTaper[i] = 1.+(xMinKVal-1.)*sqr(fromEnd/nXmax);
      kxTaperPH[i] = 1.+(xMinKVal-1.)*sqr((fromEnd+0.5)/nXmax);
      alphaXTaper[i] = xMinAVal*(1.-fromEnd/nXmax);
      alphaXTaperPH[i] = xMinAVal*(1.-(fromEnd+0.5)/nXmax);
      bxTaper[i] = exp(-(xTaper[i]/kxTaper[i]+alphaXTaper[i])*dt);
      bxTaperPH[i] = exp(-(xTaperPH[i]/kxTaperPH[i]+alphaXTaperPH[i])*dt);
      axTaper[i] = xTaper[i]==0.?0.:xTaper[i]*(bxTaper[i]-1.)/(kxTaper[i]*(xTaper[i]+kxTaper[i]*alphaXTaper[i]));
      axTaperPH[i] = xTaperPH[i]*(bxTaperPH[i]-1.)/(kxTaperPH[i]*(xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
      //        bxTaper[i] *= dt;   //multiply by dt to save mults in loops
      //        bxTaperPH[i] *= dt;
      kxTaper[i] = 1./kxTaper[i]; //more easily handled
      kxTaperPH[i] = 1./kxTaperPH[i];
      //        xTaperFac[i] = xTaper[i]*xfac*sqr(fromEnd/nXmax);
      //        xTaperPHFac[i] = xTaperPH[i]*xfac*sqr((fromEnd+0.5)/nXmax);
      if(_surfaceBCMode==ELASTIC_FS_2_0) {
        xTaperRel[i] = sqr(fromEnd/nXmax);
        xTaperRelPH[i] = sqr((fromEnd+0.5)/nXmax);
      }
    } else {
      xTaper[i] = xTaperPH[i] = alphaXTaper[i] = alphaXTaperPH[i] = axTaper[i] = axTaperPH[i] = 0.; //xTaperFac[i] = xTaperPHFac[i] = 0.f;
      kxTaper[i] = kxTaperPH[i] = 1.;
      bxTaper[i] = bxTaperPH[i] = 1.;  //dt;
      if(_surfaceBCMode==ELASTIC_FS_2_0) xTaperRel[i] = xTaperRelPH[i] = 0.f;
    }
    //      fprintf(stderr,"%f\t%f\t",xTaper[i]/xTaper[0],xTaperRel[i]);
  }
  //The Y profiles
  for(int i=0;i<NY;i++) {
    float fromStart = nYmin-i-procLim[2];
    float fromEnd = nYmax-globalNY+procLim[3]-(NY-1-i);
    if(fromStart>0.f) {
      yTaper[i] = yMinVal*sqr(fromStart/nYmin);
      yTaperPH[i] = yMinVal*sqr((fromStart-0.5)/nYmin);
      kyTaper[i] = 1.+(yMinKVal-1.)*sqr(fromStart/nYmin);
      kyTaperPH[i] = 1.+(yMinKVal-1.)*sqr((fromStart-0.5)/nYmin);
      alphaYTaper[i] = yMinAVal*(1.-fromStart/nYmin);
      alphaYTaperPH[i] = yMinAVal*(1.-(fromStart-0.5)/nYmin);
      byTaper[i] = exp(-(yTaper[i]/kyTaper[i]+alphaYTaper[i])*dt);
      byTaperPH[i] = exp(-(yTaperPH[i]/kyTaperPH[i]+alphaYTaperPH[i])*dt);
      ayTaper[i] = yTaper[i]*(byTaper[i]-1.)/(kyTaper[i]*(yTaper[i]+kyTaper[i]*alphaYTaper[i]));
      ayTaperPH[i] = yTaperPH[i]*(byTaperPH[i]-1.)/(kyTaperPH[i]*(yTaperPH[i]+kyTaperPH[i]*alphaYTaperPH[i]));
      //        byTaper[i] *= dt;   //multiply by dt to save mults in loops
      //        byTaperPH[i] *= dt;
      kyTaper[i] = 1./kyTaper[i]; //more easily handled
      kyTaperPH[i] = 1./kyTaperPH[i];
      //        yTaperFac[i] = yTaper[i]*xfac*sqr(fromStart/nYmin);
      //        yTaperPHFac[i] = yTaperPH[i]*xfac*sqr((fromStart-0.5)/nYmin);
      if(_surfaceBCMode==ELASTIC_FS_2_0) {
        yTaperRel[i] = sqr(fromStart/nYmin);
        yTaperRelPH[i] = sqr((fromStart-0.5)/nYmin);
      }
    } else if(fromEnd>-0.5f){
      yTaper[i] = yMaxVal*sqr(fromEnd/nYmax);
      yTaperPH[i] = sqr((fromEnd+0.5)/nYmax)*yMaxVal;
      kyTaper[i] = 1.+(yMinKVal-1.)*sqr(fromEnd/nYmax);
      kyTaperPH[i] = 1.+(yMinKVal-1.)*sqr((fromEnd+0.5)/nYmax);
      alphaYTaper[i] = yMinAVal*(1.-fromEnd/nYmax);
      alphaYTaperPH[i] = yMinAVal*(1.-(fromEnd+0.5)/nYmax);
      byTaper[i] = exp(-(yTaper[i]/kyTaper[i]+alphaYTaper[i])*dt);
      byTaperPH[i] = exp(-(yTaperPH[i]/kyTaperPH[i]+alphaYTaperPH[i])*dt);
      ayTaper[i] = yTaper[i]==0.?0.:yTaper[i]*(byTaper[i]-1.)/(kyTaper[i]*(yTaper[i]+kyTaper[i]*alphaYTaper[i]));
      ayTaperPH[i] = yTaperPH[i]*(byTaperPH[i]-1.)/(kyTaperPH[i]*(yTaperPH[i]+kyTaperPH[i]*alphaYTaperPH[i]));
      //        byTaper[i] *= dt;   //multiply by dt to save mults in loops
      //        byTaperPH[i] *= dt;
      kyTaper[i] = 1./kyTaper[i]; //more easily handled
      kyTaperPH[i] = 1./kyTaperPH[i];
      //        yTaperFac[i] = yTaper[i]*xfac*sqr(fromEnd/nYmax);
      //        yTaperPHFac[i] = yTaperPH[i]*xfac*sqr((fromEnd+0.5)/nYmax);
      if(_surfaceBCMode==ELASTIC_FS_2_0) {
        yTaperRel[i] = sqr(fromEnd/nYmax);
        yTaperRelPH[i] = sqr((fromEnd+0.5)/nYmax);
      }
    } else {
      yTaper[i] = yTaperPH[i] = alphaYTaper[i] = alphaYTaperPH[i] = ayTaper[i] = ayTaperPH[i] = 0.; //yTaperFac[i] = yTaperPHFac[i] = 0.f;
      kyTaper[i] = kyTaperPH[i] = 1.;
      byTaper[i] = byTaperPH[i] = 1.; //dt;
      if(_surfaceBCMode==ELASTIC_FS_2_0) yTaperRel[i] = yTaperRelPH[i] = 0.f;
    }
    //      fprintf(stderr,":%d %f %f %f %f %f %f %f %f %f %f\t",i,yTaper[i],yTaperPH[i],kyTaper[i],kyTaperPH[i],alphaYTaper[i],alphaYTaperPH[i],ayTaper[i],ayTaperPH[i],
    //              byTaper[i],byTaperPH[i]);
    //      fprintf(stderr,"%f\t%f\t",yTaper[i]/yTaper[0],yTaperRel[i]);
  }
  //The Z profiles
  for(int i=0;i<NZ;i++) {
    float fromStart = nZmin-i-procLim[4];
    float fromEnd = nZmax-globalNZ+procLim[5]-(NZ-1-i);
    if(fromStart>0.f) {
      zTaper[i] = sqr(fromStart/nZmin)*zMinVal;
      zTaperPH[i] = sqr((fromStart-0.5)/nZmin)*zMinVal;
      kzTaper[i] = 1.+(zMinKVal-1.)*sqr(fromStart/nZmin);
      kzTaperPH[i] = 1.+(zMinKVal-1.)*sqr((fromStart-0.5)/nZmin);
      alphaZTaper[i] = zMinAVal*(1.-fromStart/nZmin);
      alphaZTaperPH[i] = zMinAVal*(1.-(fromStart-0.5)/nZmin);
      bzTaper[i] = exp(-(zTaper[i]/kzTaper[i]+alphaZTaper[i])*dt);
      bzTaperPH[i] = exp(-(zTaperPH[i]/kzTaperPH[i]+alphaZTaperPH[i])*dt);
      azTaper[i] = zTaper[i]*(bzTaper[i]-1.)/(kzTaper[i]*(zTaper[i]+kzTaper[i]*alphaZTaper[i]));
      azTaperPH[i] = zTaperPH[i]*(bzTaperPH[i]-1.)/(kzTaperPH[i]*(zTaperPH[i]+kzTaperPH[i]*alphaZTaperPH[i]));
      //        bzTaper[i] *= dt;   //multiply by dt to save mults in loops
      //        bzTaperPH[i] *= dt;
      kzTaper[i] = 1./kzTaper[i]; //more easily handled
      kzTaperPH[i] = 1./kzTaperPH[i];
      //        zTaperFac[i] = zTaper[i]*xfac*sqr(fromStart/nZmin);
      //        zTaperPHFac[i] = zTaperPH[i]*xfac*sqr((fromStart-0.5)/nZmin);
    } else if(fromEnd>-0.5f) {
      zTaper[i] = sqr(fromEnd/nZmax)*zMaxVal;
      zTaperPH[i] = sqr((fromEnd+0.5)/nZmax)*zMaxVal;
      kzTaper[i] = 1.+(zMinKVal-1.)*sqr(fromEnd/nZmax);
      kzTaperPH[i] = 1.+(zMinKVal-1.)*sqr((fromEnd+0.5)/nZmax);
      alphaZTaper[i] = zMinAVal*(1.-fromEnd/nZmax);
      alphaZTaperPH[i] = zMinAVal*(1.-(fromEnd+0.5)/nZmax);
      bzTaper[i] = exp(-(zTaper[i]/kzTaper[i]+alphaZTaper[i])*dt);
      bzTaperPH[i] = exp(-(zTaperPH[i]/kzTaperPH[i]+alphaZTaperPH[i])*dt);
      azTaper[i] = zTaper[i]==0.?0.:zTaper[i]*(bzTaper[i]-1.)/(kzTaper[i]*(zTaper[i]+kzTaper[i]*alphaZTaper[i]));
      azTaperPH[i] = zTaperPH[i]*(bzTaperPH[i]-1.)/(kzTaperPH[i]*(zTaperPH[i]+kzTaperPH[i]*alphaZTaperPH[i]));
      //        bzTaper[i] *= dt;   //multiply by dt to save mults in loops
      //        bzTaperPH[i] *= dt;
      kzTaper[i] = 1./kzTaper[i]; //more easily handled
      kzTaperPH[i] = 1./kzTaperPH[i];
      //        zTaperFac[i] = zTaper[i]*xfac*sqr(fromEnd/nZmax);
      //        zTaperPHFac[i] = zTaperPH[i]*xfac*sqr((fromEnd+0.5)/nZmax);
    } else {
      zTaper[i] = zTaperPH[i] = alphaZTaper[i] = alphaZTaperPH[i] = azTaper[i] = azTaperPH[i] = 0.; //zTaperFac[i] = zTaperPHFac[i] = 0.f;
      kzTaper[i] = kzTaperPH[i] = 1.;
      bzTaper[i] = bzTaperPH[i] = 1.; //dt;
    }
    //      fprintf(stderr,"%f\t%f\t",zTaper[i],zTaperPH[i]);
  }

float amax = -1e10, amin = 1e10, bmax = -1e10, bmin = 1e10;
for(int i=0;i<_nXmin;++i) {
  byX0Taper[i] = exp(-xtap*xTaper[i]*dt);
  ayX0Taper[i] = byX0Taper[i]-1.;
amax = ayX0Taper[i]>amax?ayX0Taper[i]:amax;
amin = ayX0Taper[i]<amin?ayX0Taper[i]:amin;
bmax = byX0Taper[i]>bmax?byX0Taper[i]:bmax;
bmin = byX0Taper[i]<bmin?byX0Taper[i]:bmin;
  byX0TaperPH[i] = exp(-xtap*xTaperPH[i]*dt);
  ayX0TaperPH[i] = byX0TaperPH[i]-1.;
amax = ayX0TaperPH[i]>amax?ayX0TaperPH[i]:amax;
amin = ayX0TaperPH[i]<amin?ayX0TaperPH[i]:amin;
bmax = byX0TaperPH[i]>bmax?byX0TaperPH[i]:bmax;
bmin = byX0TaperPH[i]<bmin?byX0TaperPH[i]:bmin;
  bzX0Taper[i] = exp(-xtap*xTaper[i]*dt);
  azX0Taper[i] = bzX0Taper[i]-1.;
amax = azX0Taper[i]>amax?azX0Taper[i]:amax;
amin = azX0Taper[i]<amin?azX0Taper[i]:amin;
bmax = bzX0Taper[i]>bmax?bzX0Taper[i]:bmax;
bmin = bzX0Taper[i]<bmin?bzX0Taper[i]:bmin;
  bzX0TaperPH[i] = exp(-xtap*xTaperPH[i]*dt);
  azX0TaperPH[i] = bzX0TaperPH[i]-1.;
amax = azX0TaperPH[i]>amax?azX0TaperPH[i]:amax;
amin = azX0TaperPH[i]<amin?azX0TaperPH[i]:amin;
bmax = bzX0TaperPH[i]>bmax?bzX0TaperPH[i]:bmax;
bmin = bzX0TaperPH[i]<bmin?bzX0TaperPH[i]:bmin;
}
for(int i=xEnd1;i<NX;++i) {
  int ti = i-xEnd1;
  byX1Taper[ti] = exp(-xtap*xTaper[i]*dt);
  ayX1Taper[ti] = byX1Taper[ti]-1.;
amax = ayX1Taper[ti]>amax?ayX1Taper[ti]:amax;
amin = ayX1Taper[ti]<amin?ayX1Taper[ti]:amin;
bmax = byX1Taper[ti]>bmax?byX1Taper[ti]:bmax;
bmin = byX1Taper[ti]<bmin?byX1Taper[ti]:bmin;
  byX1TaperPH[ti] = exp(-xtap*xTaperPH[i]*dt);
  ayX1TaperPH[ti] = byX1TaperPH[ti]-1.;
amax = ayX1TaperPH[ti]>amax?ayX1TaperPH[ti]:amax;
amin = ayX1TaperPH[ti]<amin?ayX1TaperPH[ti]:amin;
bmax = byX1TaperPH[ti]>bmax?byX1TaperPH[ti]:bmax;
bmin = byX1TaperPH[ti]<bmin?byX1TaperPH[ti]:bmin;
  bzX1Taper[ti] = exp(-xtap*xTaper[i]*dt);
  azX1Taper[ti] = bzX1Taper[ti]-1.;
amax = azX1Taper[ti]>amax?azX1Taper[ti]:amax;
amin = azX1Taper[ti]<amin?azX1Taper[ti]:amin;
bmax = bzX1Taper[ti]>bmax?bzX1Taper[ti]:bmax;
bmin = bzX1Taper[ti]<bmin?bzX1Taper[ti]:bmin;
  bzX1TaperPH[ti] = exp(-xtap*xTaperPH[i]*dt);
  azX1TaperPH[ti] = bzX1TaperPH[ti]-1.;
amax = azX1TaperPH[ti]>amax?azX1TaperPH[ti]:amax;
amin = azX1TaperPH[ti]<amin?azX1TaperPH[ti]:amin;
bmax = bzX1TaperPH[ti]>bmax?bzX1TaperPH[ti]:bmax;
bmin = bzX1TaperPH[ti]<bmin?bzX1TaperPH[ti]:bmin;
}
for(int j=0;j<_nYmin;++j) {
  bxY0Taper[j] = exp(-xtap*yTaper[j]*dt);
  axY0Taper[j] = bxY0Taper[j]-1.;
amax = axY0Taper[j]>amax?axY0Taper[j]:amax;
amin = axY0Taper[j]<amin?axY0Taper[j]:amin;
bmax = bxY0Taper[j]>bmax?bxY0Taper[j]:bmax;
bmin = bxY0Taper[j]<bmin?bxY0Taper[j]:bmin;
  bxY0TaperPH[j] = exp(-xtap*yTaperPH[j]*dt);
  axY0TaperPH[j] = bxY0TaperPH[j]-1.;
amax = axY0TaperPH[j]>amax?axY0TaperPH[j]:amax;
amin = axY0TaperPH[j]<amin?axY0TaperPH[j]:amin;
bmax = bxY0TaperPH[j]>bmax?bxY0TaperPH[j]:bmax;
bmin = bxY0TaperPH[j]<bmin?bxY0TaperPH[j]:bmin;
  bzY0Taper[j] = exp(-xtap*yTaper[j]*dt);
  azY0Taper[j] = bzY0Taper[j]-1.;
amax = azY0Taper[j]>amax?azY0Taper[j]:amax;
amin = azY0Taper[j]<amin?azY0Taper[j]:amin;
bmax = bzY0Taper[j]>bmax?bzY0Taper[j]:bmax;
bmin = bzY0Taper[j]<bmin?bzY0Taper[j]:bmin;
  bzY0TaperPH[j] = exp(-xtap*yTaperPH[j]*dt);
  azY0TaperPH[j] = bzY0TaperPH[j]-1.;
amax = azY0TaperPH[j]>amax?azY0TaperPH[j]:amax;
amin = azY0TaperPH[j]<amin?azY0TaperPH[j]:amin;
bmax = bzY0TaperPH[j]>bmax?bzY0TaperPH[j]:bmax;
bmin = bzY0TaperPH[j]<bmin?bzY0TaperPH[j]:bmin;
}
for(int j=yEnd1;j<NY;++j) {
  int tj = j-yEnd1;
  bxY1Taper[tj] = exp(-xtap*yTaper[j]*dt);
  axY1Taper[tj] = bxY1Taper[tj]-1.;
amax = axY1Taper[tj]>amax?axY1Taper[tj]:amax;
amin = axY1Taper[tj]<amin?axY1Taper[tj]:amin;
bmax = bxY1Taper[tj]>bmax?bxY1Taper[tj]:bmax;
bmin = bxY1Taper[tj]<bmin?bxY1Taper[tj]:bmin;
  bxY1TaperPH[tj] = exp(-xtap*yTaperPH[j]*dt);
  axY1TaperPH[tj] = bxY1TaperPH[tj]-1.;
amax = axY1TaperPH[tj]>amax?axY1TaperPH[tj]:amax;
amin = axY1TaperPH[tj]<amin?axY1TaperPH[tj]:amin;
bmax = bxY1TaperPH[tj]>bmax?bxY1TaperPH[tj]:bmax;
bmin = bxY1TaperPH[tj]<bmin?bxY1TaperPH[tj]:bmin;
  bzY1Taper[tj] = exp(-xtap*yTaper[j]*dt);
  azY1Taper[tj] = bzY1Taper[tj]-1.;
amax = azY1Taper[tj]>amax?azY1Taper[tj]:amax;
amin = azY1Taper[tj]<amin?azY1Taper[tj]:amin;
bmax = bzY1Taper[tj]>bmax?bzY1Taper[tj]:bmax;
bmin = bzY1Taper[tj]<bmin?bzY1Taper[tj]:bmin;
  bzY1TaperPH[tj] = exp(-xtap*yTaperPH[j]*dt);
  azY1TaperPH[tj] = bzY1TaperPH[tj]-1.;
amax = azY1TaperPH[tj]>amax?azY1TaperPH[tj]:amax;
amin = azY1TaperPH[tj]<amin?azY1TaperPH[tj]:amin;
bmax = bzY1TaperPH[tj]>bmax?bzY1TaperPH[tj]:bmax;
bmin = bzY1TaperPH[tj]<bmin?bzY1TaperPH[tj]:bmin;
}
for(int k=0;k<_nZmin;++k) {
  bxZ0Taper[k] = exp(-xtap*zTaper[k]*dt);
  axZ0Taper[k] = bxZ0Taper[k]-1.;
amax = axZ0Taper[k]>amax?axZ0Taper[k]:amax;
amin = axZ0Taper[k]<amin?axZ0Taper[k]:amin;
bmax = bxZ0Taper[k]>bmax?bxZ0Taper[k]:bmax;
bmin = bxZ0Taper[k]<bmin?bxZ0Taper[k]:bmin;
  bxZ0TaperPH[k] = exp(-xtap*zTaperPH[k]*dt);
  axZ0TaperPH[k] = bxZ0TaperPH[k]-1.;
amax = axZ0TaperPH[k]>amax?axZ0TaperPH[k]:amax;
amin = axZ0TaperPH[k]<amin?axZ0TaperPH[k]:amin;
bmax = bxZ0TaperPH[k]>bmax?bxZ0TaperPH[k]:bmax;
bmin = bxZ0TaperPH[k]<bmin?bxZ0TaperPH[k]:bmin;
  byZ0Taper[k] = exp(-xtap*zTaper[k]*dt);
  ayZ0Taper[k] = byZ0Taper[k]-1.;
amax = ayZ0Taper[k]>amax?ayZ0Taper[k]:amax;
amin = ayZ0Taper[k]<amin?ayZ0Taper[k]:amin;
bmax = byZ0Taper[k]>bmax?byZ0Taper[k]:bmax;
bmin = byZ0Taper[k]<bmin?byZ0Taper[k]:bmin;
  byZ0TaperPH[k] = exp(-xtap*zTaperPH[k]*dt);
  ayZ0TaperPH[k] = byZ0TaperPH[k]-1.;
amax = ayZ0TaperPH[k]>amax?ayZ0TaperPH[k]:amax;
amin = ayZ0TaperPH[k]<amin?ayZ0TaperPH[k]:amin;
bmax = byZ0TaperPH[k]>bmax?byZ0TaperPH[k]:bmax;
bmin = byZ0TaperPH[k]<bmin?byZ0TaperPH[k]:bmin;
}
for(int k=zEnd1;k<NZ;++k) {
  int tk = k-zEnd1;
  bxZ1Taper[tk] = exp(-xtap*zTaper[k]*dt);
  axZ1Taper[tk] = bxZ1Taper[tk]-1.;
amax = axZ1Taper[tk]>amax?axZ1Taper[tk]:amax;
amin = axZ1Taper[tk]<amin?axZ1Taper[tk]:amin;
bmax = bxZ1Taper[tk]>bmax?bxZ1Taper[tk]:bmax;
bmin = bxZ1Taper[tk]<bmin?bxZ1Taper[tk]:bmin;
  bxZ1TaperPH[tk] = exp(-xtap*zTaperPH[k]*dt);
  axZ1TaperPH[tk] = bxZ1TaperPH[tk]-1.;
amax = axZ1TaperPH[tk]>amax?axZ1TaperPH[tk]:amax;
amin = axZ1TaperPH[tk]<amin?axZ1TaperPH[tk]:amin;
bmax = bxZ1TaperPH[tk]>bmax?bxZ1TaperPH[tk]:bmax;
bmin = bxZ1TaperPH[tk]<bmin?bxZ1TaperPH[tk]:bmin;
  byZ1Taper[tk] = exp(-xtap*zTaper[k]*dt);
  ayZ1Taper[tk] = byZ1Taper[tk]-1.;
amax = ayZ1Taper[tk]>amax?ayZ1Taper[tk]:amax;
amin = ayZ1Taper[tk]<amin?ayZ1Taper[tk]:amin;
bmax = byZ1Taper[tk]>bmax?byZ1Taper[tk]:bmax;
bmin = byZ1Taper[tk]<bmin?byZ1Taper[tk]:bmin;
  byZ1TaperPH[tk] = exp(-xtap*zTaperPH[k]*dt);
  ayZ1TaperPH[tk] = byZ1TaperPH[tk]-1.;
amax = ayZ1TaperPH[tk]>amax?ayZ1TaperPH[tk]:amax;
amin = ayZ1TaperPH[tk]<amin?ayZ1TaperPH[tk]:amin;
bmax = byZ1TaperPH[tk]>bmax?byZ1TaperPH[tk]:bmax;
bmin = byZ1TaperPH[tk]<bmin?byZ1TaperPH[tk]:bmin;
}
for(int k=0;k<_nZmin;++k) {
  for(int i=0;i<_nXmin;++i) {
    bxE101Taper[k][i] = bxTaper[i]*exp(-xtap*zTaper[k]*dt);
    axE101Taper[k][i] = (xfac*zTaper[k]+xTaper[i])*(bxE101Taper[k][i]-1.)/(kxTaper[i]*(xfac*zTaper[k]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE101Taper[k][i]>amax?axE101Taper[k][i]:amax;
amin = axE101Taper[k][i]<amin?axE101Taper[k][i]:amin;
bmax = bxE101Taper[k][i]>bmax?bxE101Taper[k][i]:bmax;
bmin = bxE101Taper[k][i]<bmin?bxE101Taper[k][i]:bmin;
    bxE101TaperzPH[k][i] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt);
    axE101TaperzPH[k][i] = (xfac*zTaperPH[k]+xTaper[i])*(bxE101TaperzPH[k][i]-1.)/(kxTaper[i]*(xfac*zTaperPH[k]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE101TaperzPH[k][i]>amax?axE101TaperzPH[k][i]:amax;
amin = axE101TaperzPH[k][i]<amin?axE101TaperzPH[k][i]:amin;
bmax = bxE101TaperzPH[k][i]>bmax?bxE101TaperzPH[k][i]:bmax;
bmin = bxE101TaperzPH[k][i]<bmin?bxE101TaperzPH[k][i]:bmin;
    bxE101TaperzPHxPH[k][i] = bxTaperPH[i]*exp(-xtap*zTaperPH[k]*dt);
    axE101TaperzPHxPH[k][i] = (xfac*zTaperPH[k]+xTaperPH[i])*(bxE101TaperzPHxPH[k][i]-1.)/(kxTaperPH[i]*(xfac*zTaperPH[k]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE101TaperzPHxPH[k][i]>amax?axE101TaperzPHxPH[k][i]:amax;
amin = axE101TaperzPHxPH[k][i]<amin?axE101TaperzPHxPH[k][i]:amin;
bmax = bxE101TaperzPHxPH[k][i]>bmax?bxE101TaperzPHxPH[k][i]:bmax;
bmin = bxE101TaperzPHxPH[k][i]<bmin?bxE101TaperzPHxPH[k][i]:bmin;
    bxE101TaperxPH[k][i] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt);
    axE101TaperxPH[k][i] = (xfac*zTaper[k]+xTaperPH[i])*(bxE101TaperxPH[k][i]-1.)/(kxTaperPH[i]*(xfac*zTaper[k]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE101TaperxPH[k][i]>amax?axE101TaperxPH[k][i]:amax;
amin = axE101TaperxPH[k][i]<amin?axE101TaperxPH[k][i]:amin;
bmax = bxE101TaperxPH[k][i]>bmax?bxE101TaperxPH[k][i]:bmax;
bmin = bxE101TaperxPH[k][i]<bmin?bxE101TaperxPH[k][i]:bmin;
    byE101Taper[k][i] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
    ayE101Taper[k][i] = byE101Taper[k][i]-1.;
amax = ayE101Taper[k][i]>amax?ayE101Taper[k][i]:amax;
amin = ayE101Taper[k][i]<amin?ayE101Taper[k][i]:amin;
bmax = byE101Taper[k][i]>bmax?byE101Taper[k][i]:bmax;
bmin = byE101Taper[k][i]<bmin?byE101Taper[k][i]:bmin;
    byE101TaperzPH[k][i] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
    ayE101TaperzPH[k][i] = byE101TaperzPH[k][i]-1.;
amax = ayE101TaperzPH[k][i]>amax?ayE101TaperzPH[k][i]:amax;
amin = ayE101TaperzPH[k][i]<amin?ayE101TaperzPH[k][i]:amin;
bmax = byE101TaperzPH[k][i]>bmax?byE101TaperzPH[k][i]:bmax;
bmin = byE101TaperzPH[k][i]<bmin?byE101TaperzPH[k][i]:bmin;
    byE101TaperzPHxPH[k][i] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
    ayE101TaperzPHxPH[k][i] = byE101TaperzPHxPH[k][i]-1.;
amax = ayE101TaperzPHxPH[k][i]>amax?ayE101TaperzPHxPH[k][i]:amax;
amin = ayE101TaperzPHxPH[k][i]<amin?ayE101TaperzPHxPH[k][i]:amin;
bmax = byE101TaperzPHxPH[k][i]>bmax?byE101TaperzPHxPH[k][i]:bmax;
bmin = byE101TaperzPHxPH[k][i]<bmin?byE101TaperzPHxPH[k][i]:bmin;
    byE101TaperxPH[k][i] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
    ayE101TaperxPH[k][i] = byE101TaperxPH[k][i]-1.;
amax = ayE101TaperxPH[k][i]>amax?ayE101TaperxPH[k][i]:amax;
amin = ayE101TaperxPH[k][i]<amin?ayE101TaperxPH[k][i]:amin;
bmax = byE101TaperxPH[k][i]>bmax?byE101TaperxPH[k][i]:bmax;
bmin = byE101TaperxPH[k][i]<bmin?byE101TaperxPH[k][i]:bmin;
    bzE101Taper[k][i] = bzTaper[k]*exp(-xtap*xTaper[i]*dt);
    azE101Taper[k][i] = (xfac*xTaper[i]+zTaper[k])*(bzE101Taper[k][i]-1.)/(kzTaper[k]*(xfac*xTaper[i]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE101Taper[k][i]>amax?azE101Taper[k][i]:amax;
amin = azE101Taper[k][i]<amin?azE101Taper[k][i]:amin;
bmax = bzE101Taper[k][i]>bmax?bzE101Taper[k][i]:bmax;
bmin = bzE101Taper[k][i]<bmin?bzE101Taper[k][i]:bmin;
    bzE101TaperzPH[k][i] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt);
    azE101TaperzPH[k][i] = (xfac*xTaper[i]+zTaperPH[k])*(bzE101TaperzPH[k][i]-1.)/(kzTaperPH[k]*(xfac*xTaper[i]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE101TaperzPH[k][i]>amax?azE101TaperzPH[k][i]:amax;
amin = azE101TaperzPH[k][i]<amin?azE101TaperzPH[k][i]:amin;
bmax = bzE101TaperzPH[k][i]>bmax?bzE101TaperzPH[k][i]:bmax;
bmin = bzE101TaperzPH[k][i]<bmin?bzE101TaperzPH[k][i]:bmin;
    bzE101TaperzPHxPH[k][i] = bzTaperPH[k]*exp(-xtap*xTaperPH[i]*dt);
    azE101TaperzPHxPH[k][i] = (xfac*xTaperPH[i]+zTaperPH[k])*(bzE101TaperzPHxPH[k][i]-1.)/(kzTaperPH[k]*(xfac*xTaperPH[i]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE101TaperzPHxPH[k][i]>amax?azE101TaperzPHxPH[k][i]:amax;
amin = azE101TaperzPHxPH[k][i]<amin?azE101TaperzPHxPH[k][i]:amin;
bmax = bzE101TaperzPHxPH[k][i]>bmax?bzE101TaperzPHxPH[k][i]:bmax;
bmin = bzE101TaperzPHxPH[k][i]<bmin?bzE101TaperzPHxPH[k][i]:bmin;
    bzE101TaperxPH[k][i] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt);
    azE101TaperxPH[k][i] = (xfac*xTaperPH[i]+zTaper[k])*(bzE101TaperxPH[k][i]-1.)/(kzTaper[k]*(xfac*xTaperPH[i]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE101TaperxPH[k][i]>amax?azE101TaperxPH[k][i]:amax;
amin = azE101TaperxPH[k][i]<amin?azE101TaperxPH[k][i]:amin;
bmax = bzE101TaperxPH[k][i]>bmax?bzE101TaperxPH[k][i]:bmax;
bmin = bzE101TaperxPH[k][i]<bmin?bzE101TaperxPH[k][i]:bmin;
  }
  for(int i=xEnd1;i<NX;++i) {
    int ti = i-xEnd1;
    bxE201Taper[k][ti] = bxTaper[i]*exp(-xtap*zTaper[k]*dt);
    axE201Taper[k][ti] = (xfac*zTaper[k]+xTaper[i])*(bxE201Taper[k][ti]-1.)/(kxTaper[i]*(xfac*zTaper[k]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE201Taper[k][ti]>amax?axE201Taper[k][ti]:amax;
amin = axE201Taper[k][ti]<amin?axE201Taper[k][ti]:amin;
bmax = bxE201Taper[k][ti]>bmax?bxE201Taper[k][ti]:bmax;
bmin = bxE201Taper[k][ti]<bmin?bxE201Taper[k][ti]:bmin;
    bxE201TaperzPH[k][ti] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt);
    axE201TaperzPH[k][ti] = (xfac*zTaperPH[k]+xTaper[i])*(bxE201TaperzPH[k][ti]-1.)/(kxTaper[i]*(xfac*zTaperPH[k]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE201TaperzPH[k][ti]>amax?axE201TaperzPH[k][ti]:amax;
amin = axE201TaperzPH[k][ti]<amin?axE201TaperzPH[k][ti]:amin;
bmax = bxE201TaperzPH[k][ti]>bmax?bxE201TaperzPH[k][ti]:bmax;
bmin = bxE201TaperzPH[k][ti]<bmin?bxE201TaperzPH[k][ti]:bmin;
    bxE201TaperzPHxPH[k][ti] = bxTaperPH[i]*exp(-xtap*zTaperPH[k]*dt);
    axE201TaperzPHxPH[k][ti] = (xfac*zTaperPH[k]+xTaperPH[i])*(bxE201TaperzPHxPH[k][ti]-1.)/(kxTaperPH[i]*(xfac*zTaperPH[k]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE201TaperzPHxPH[k][ti]>amax?axE201TaperzPHxPH[k][ti]:amax;
amin = axE201TaperzPHxPH[k][ti]<amin?axE201TaperzPHxPH[k][ti]:amin;
bmax = bxE201TaperzPHxPH[k][ti]>bmax?bxE201TaperzPHxPH[k][ti]:bmax;
bmin = bxE201TaperzPHxPH[k][ti]<bmin?bxE201TaperzPHxPH[k][ti]:bmin;
    bxE201TaperxPH[k][ti] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt);
    axE201TaperxPH[k][ti] = (xfac*zTaper[k]+xTaperPH[i])*(bxE201TaperxPH[k][ti]-1.)/(kxTaperPH[i]*(xfac*zTaper[k]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE201TaperxPH[k][ti]>amax?axE201TaperxPH[k][ti]:amax;
amin = axE201TaperxPH[k][ti]<amin?axE201TaperxPH[k][ti]:amin;
bmax = bxE201TaperxPH[k][ti]>bmax?bxE201TaperxPH[k][ti]:bmax;
bmin = bxE201TaperxPH[k][ti]<bmin?bxE201TaperxPH[k][ti]:bmin;
    byE201Taper[k][ti] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
    ayE201Taper[k][ti] = byE201Taper[k][ti]-1.;
amax = ayE201Taper[k][ti]>amax?ayE201Taper[k][ti]:amax;
amin = ayE201Taper[k][ti]<amin?ayE201Taper[k][ti]:amin;
bmax = byE201Taper[k][ti]>bmax?byE201Taper[k][ti]:bmax;
bmin = byE201Taper[k][ti]<bmin?byE201Taper[k][ti]:bmin;
    byE201TaperzPH[k][ti] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
    ayE201TaperzPH[k][ti] = byE201TaperzPH[k][ti]-1.;
amax = ayE201TaperzPH[k][ti]>amax?ayE201TaperzPH[k][ti]:amax;
amin = ayE201TaperzPH[k][ti]<amin?ayE201TaperzPH[k][ti]:amin;
bmax = byE201TaperzPH[k][ti]>bmax?byE201TaperzPH[k][ti]:bmax;
bmin = byE201TaperzPH[k][ti]<bmin?byE201TaperzPH[k][ti]:bmin;
    byE201TaperzPHxPH[k][ti] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
    ayE201TaperzPHxPH[k][ti] = byE201TaperzPHxPH[k][ti]-1.;
amax = ayE201TaperzPHxPH[k][ti]>amax?ayE201TaperzPHxPH[k][ti]:amax;
amin = ayE201TaperzPHxPH[k][ti]<amin?ayE201TaperzPHxPH[k][ti]:amin;
bmax = byE201TaperzPHxPH[k][ti]>bmax?byE201TaperzPHxPH[k][ti]:bmax;
bmin = byE201TaperzPHxPH[k][ti]<bmin?byE201TaperzPHxPH[k][ti]:bmin;
    byE201TaperxPH[k][ti] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
    ayE201TaperxPH[k][ti] = byE201TaperxPH[k][ti]-1.;
amax = ayE201TaperxPH[k][ti]>amax?ayE201TaperxPH[k][ti]:amax;
amin = ayE201TaperxPH[k][ti]<amin?ayE201TaperxPH[k][ti]:amin;
bmax = byE201TaperxPH[k][ti]>bmax?byE201TaperxPH[k][ti]:bmax;
bmin = byE201TaperxPH[k][ti]<bmin?byE201TaperxPH[k][ti]:bmin;
    bzE201Taper[k][ti] = bzTaper[k]*exp(-xtap*xTaper[i]*dt);
    azE201Taper[k][ti] = (xfac*xTaper[i]+zTaper[k])*(bzE201Taper[k][ti]-1.)/(kzTaper[k]*(xfac*xTaper[i]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE201Taper[k][ti]>amax?azE201Taper[k][ti]:amax;
amin = azE201Taper[k][ti]<amin?azE201Taper[k][ti]:amin;
bmax = bzE201Taper[k][ti]>bmax?bzE201Taper[k][ti]:bmax;
bmin = bzE201Taper[k][ti]<bmin?bzE201Taper[k][ti]:bmin;
    bzE201TaperzPH[k][ti] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt);
    azE201TaperzPH[k][ti] = (xfac*xTaper[i]+zTaperPH[k])*(bzE201TaperzPH[k][ti]-1.)/(kzTaperPH[k]*(xfac*xTaper[i]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE201TaperzPH[k][ti]>amax?azE201TaperzPH[k][ti]:amax;
amin = azE201TaperzPH[k][ti]<amin?azE201TaperzPH[k][ti]:amin;
bmax = bzE201TaperzPH[k][ti]>bmax?bzE201TaperzPH[k][ti]:bmax;
bmin = bzE201TaperzPH[k][ti]<bmin?bzE201TaperzPH[k][ti]:bmin;
    bzE201TaperzPHxPH[k][ti] = bzTaperPH[k]*exp(-xtap*xTaperPH[i]*dt);
    azE201TaperzPHxPH[k][ti] = (xfac*xTaperPH[i]+zTaperPH[k])*(bzE201TaperzPHxPH[k][ti]-1.)/(kzTaperPH[k]*(xfac*xTaperPH[i]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE201TaperzPHxPH[k][ti]>amax?azE201TaperzPHxPH[k][ti]:amax;
amin = azE201TaperzPHxPH[k][ti]<amin?azE201TaperzPHxPH[k][ti]:amin;
bmax = bzE201TaperzPHxPH[k][ti]>bmax?bzE201TaperzPHxPH[k][ti]:bmax;
bmin = bzE201TaperzPHxPH[k][ti]<bmin?bzE201TaperzPHxPH[k][ti]:bmin;
    bzE201TaperxPH[k][ti] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt);
    azE201TaperxPH[k][ti] = (xfac*xTaperPH[i]+zTaper[k])*(bzE201TaperxPH[k][ti]-1.)/(kzTaper[k]*(xfac*xTaperPH[i]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE201TaperxPH[k][ti]>amax?azE201TaperxPH[k][ti]:amax;
amin = azE201TaperxPH[k][ti]<amin?azE201TaperxPH[k][ti]:amin;
bmax = bzE201TaperxPH[k][ti]>bmax?bzE201TaperxPH[k][ti]:bmax;
bmin = bzE201TaperxPH[k][ti]<bmin?bzE201TaperxPH[k][ti]:bmin;
  }
}
for(int k=zEnd1;k<NZ;++k) {
  int tk = k-zEnd1;
  for(int i=0;i<_nXmin;++i) {
    bxE102Taper[tk][i] = bxTaper[i]*exp(-xtap*zTaper[k]*dt);
    axE102Taper[tk][i] = (xfac*zTaper[k]+xTaper[i])*(bxE102Taper[tk][i]-1.)/(kxTaper[i]*(xfac*zTaper[k]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE102Taper[tk][i]>amax?axE102Taper[tk][i]:amax;
amin = axE102Taper[tk][i]<amin?axE102Taper[tk][i]:amin;
bmax = bxE102Taper[tk][i]>bmax?bxE102Taper[tk][i]:bmax;
bmin = bxE102Taper[tk][i]<bmin?bxE102Taper[tk][i]:bmin;
    bxE102TaperzPH[tk][i] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt);
    axE102TaperzPH[tk][i] = (xfac*zTaperPH[k]+xTaper[i])*(bxE102TaperzPH[tk][i]-1.)/(kxTaper[i]*(xfac*zTaperPH[k]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE102TaperzPH[tk][i]>amax?axE102TaperzPH[tk][i]:amax;
amin = axE102TaperzPH[tk][i]<amin?axE102TaperzPH[tk][i]:amin;
bmax = bxE102TaperzPH[tk][i]>bmax?bxE102TaperzPH[tk][i]:bmax;
bmin = bxE102TaperzPH[tk][i]<bmin?bxE102TaperzPH[tk][i]:bmin;
    bxE102TaperzPHxPH[tk][i] = bxTaperPH[i]*exp(-xtap*zTaperPH[k]*dt);
    axE102TaperzPHxPH[tk][i] = (xfac*zTaperPH[k]+xTaperPH[i])*(bxE102TaperzPHxPH[tk][i]-1.)/(kxTaperPH[i]*(xfac*zTaperPH[k]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE102TaperzPHxPH[tk][i]>amax?axE102TaperzPHxPH[tk][i]:amax;
amin = axE102TaperzPHxPH[tk][i]<amin?axE102TaperzPHxPH[tk][i]:amin;
bmax = bxE102TaperzPHxPH[tk][i]>bmax?bxE102TaperzPHxPH[tk][i]:bmax;
bmin = bxE102TaperzPHxPH[tk][i]<bmin?bxE102TaperzPHxPH[tk][i]:bmin;
    bxE102TaperxPH[tk][i] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt);
    axE102TaperxPH[tk][i] = (xfac*zTaper[k]+xTaperPH[i])*(bxE102TaperxPH[tk][i]-1.)/(kxTaperPH[i]*(xfac*zTaper[k]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE102TaperxPH[tk][i]>amax?axE102TaperxPH[tk][i]:amax;
amin = axE102TaperxPH[tk][i]<amin?axE102TaperxPH[tk][i]:amin;
bmax = bxE102TaperxPH[tk][i]>bmax?bxE102TaperxPH[tk][i]:bmax;
bmin = bxE102TaperxPH[tk][i]<bmin?bxE102TaperxPH[tk][i]:bmin;
    byE102Taper[tk][i] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
    ayE102Taper[tk][i] = byE102Taper[tk][i]-1.;
amax = ayE102Taper[tk][i]>amax?ayE102Taper[tk][i]:amax;
amin = ayE102Taper[tk][i]<amin?ayE102Taper[tk][i]:amin;
bmax = byE102Taper[tk][i]>bmax?byE102Taper[tk][i]:bmax;
bmin = byE102Taper[tk][i]<bmin?byE102Taper[tk][i]:bmin;
    byE102TaperzPH[tk][i] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
    ayE102TaperzPH[tk][i] = byE102TaperzPH[tk][i]-1.;
amax = ayE102TaperzPH[tk][i]>amax?ayE102TaperzPH[tk][i]:amax;
amin = ayE102TaperzPH[tk][i]<amin?ayE102TaperzPH[tk][i]:amin;
bmax = byE102TaperzPH[tk][i]>bmax?byE102TaperzPH[tk][i]:bmax;
bmin = byE102TaperzPH[tk][i]<bmin?byE102TaperzPH[tk][i]:bmin;
    byE102TaperzPHxPH[tk][i] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
    ayE102TaperzPHxPH[tk][i] = byE102TaperzPHxPH[tk][i]-1.;
amax = ayE102TaperzPHxPH[tk][i]>amax?ayE102TaperzPHxPH[tk][i]:amax;
amin = ayE102TaperzPHxPH[tk][i]<amin?ayE102TaperzPHxPH[tk][i]:amin;
bmax = byE102TaperzPHxPH[tk][i]>bmax?byE102TaperzPHxPH[tk][i]:bmax;
bmin = byE102TaperzPHxPH[tk][i]<bmin?byE102TaperzPHxPH[tk][i]:bmin;
    byE102TaperxPH[tk][i] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
    ayE102TaperxPH[tk][i] = byE102TaperxPH[tk][i]-1.;
amax = ayE102TaperxPH[tk][i]>amax?ayE102TaperxPH[tk][i]:amax;
amin = ayE102TaperxPH[tk][i]<amin?ayE102TaperxPH[tk][i]:amin;
bmax = byE102TaperxPH[tk][i]>bmax?byE102TaperxPH[tk][i]:bmax;
bmin = byE102TaperxPH[tk][i]<bmin?byE102TaperxPH[tk][i]:bmin;
    bzE102Taper[tk][i] = bzTaper[k]*exp(-xtap*xTaper[i]*dt);
    azE102Taper[tk][i] = (xfac*xTaper[i]+zTaper[k])*(bzE102Taper[tk][i]-1.)/(kzTaper[k]*(xfac*xTaper[i]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE102Taper[tk][i]>amax?azE102Taper[tk][i]:amax;
amin = azE102Taper[tk][i]<amin?azE102Taper[tk][i]:amin;
bmax = bzE102Taper[tk][i]>bmax?bzE102Taper[tk][i]:bmax;
bmin = bzE102Taper[tk][i]<bmin?bzE102Taper[tk][i]:bmin;
    bzE102TaperzPH[tk][i] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt);
    azE102TaperzPH[tk][i] = (xfac*xTaper[i]+zTaperPH[k])*(bzE102TaperzPH[tk][i]-1.)/(kzTaperPH[k]*(xfac*xTaper[i]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE102TaperzPH[tk][i]>amax?azE102TaperzPH[tk][i]:amax;
amin = azE102TaperzPH[tk][i]<amin?azE102TaperzPH[tk][i]:amin;
bmax = bzE102TaperzPH[tk][i]>bmax?bzE102TaperzPH[tk][i]:bmax;
bmin = bzE102TaperzPH[tk][i]<bmin?bzE102TaperzPH[tk][i]:bmin;
    bzE102TaperzPHxPH[tk][i] = bzTaperPH[k]*exp(-xtap*xTaperPH[i]*dt);
    azE102TaperzPHxPH[tk][i] = (xfac*xTaperPH[i]+zTaperPH[k])*(bzE102TaperzPHxPH[tk][i]-1.)/(kzTaperPH[k]*(xfac*xTaperPH[i]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE102TaperzPHxPH[tk][i]>amax?azE102TaperzPHxPH[tk][i]:amax;
amin = azE102TaperzPHxPH[tk][i]<amin?azE102TaperzPHxPH[tk][i]:amin;
bmax = bzE102TaperzPHxPH[tk][i]>bmax?bzE102TaperzPHxPH[tk][i]:bmax;
bmin = bzE102TaperzPHxPH[tk][i]<bmin?bzE102TaperzPHxPH[tk][i]:bmin;
    bzE102TaperxPH[tk][i] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt);
    azE102TaperxPH[tk][i] = (xfac*xTaperPH[i]+zTaper[k])*(bzE102TaperxPH[tk][i]-1.)/(kzTaper[k]*(xfac*xTaperPH[i]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE102TaperxPH[tk][i]>amax?azE102TaperxPH[tk][i]:amax;
amin = azE102TaperxPH[tk][i]<amin?azE102TaperxPH[tk][i]:amin;
bmax = bzE102TaperxPH[tk][i]>bmax?bzE102TaperxPH[tk][i]:bmax;
bmin = bzE102TaperxPH[tk][i]<bmin?bzE102TaperxPH[tk][i]:bmin;
  }
  for(int i=xEnd1;i<NX;++i) {
    int ti = i-xEnd1;
    bxE202Taper[tk][ti] = bxTaper[i]*exp(-xtap*zTaper[k]*dt);
    axE202Taper[tk][ti] = (xfac*zTaper[k]+xTaper[i])*(bxE202Taper[tk][ti]-1.)/(kxTaper[i]*(xfac*zTaper[k]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE202Taper[tk][ti]>amax?axE202Taper[tk][ti]:amax;
amin = axE202Taper[tk][ti]<amin?axE202Taper[tk][ti]:amin;
bmax = bxE202Taper[tk][ti]>bmax?bxE202Taper[tk][ti]:bmax;
bmin = bxE202Taper[tk][ti]<bmin?bxE202Taper[tk][ti]:bmin;
    bxE202TaperzPH[tk][ti] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt);
    axE202TaperzPH[tk][ti] = (xfac*zTaperPH[k]+xTaper[i])*(bxE202TaperzPH[tk][ti]-1.)/(kxTaper[i]*(xfac*zTaperPH[k]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE202TaperzPH[tk][ti]>amax?axE202TaperzPH[tk][ti]:amax;
amin = axE202TaperzPH[tk][ti]<amin?axE202TaperzPH[tk][ti]:amin;
bmax = bxE202TaperzPH[tk][ti]>bmax?bxE202TaperzPH[tk][ti]:bmax;
bmin = bxE202TaperzPH[tk][ti]<bmin?bxE202TaperzPH[tk][ti]:bmin;
    bxE202TaperzPHxPH[tk][ti] = bxTaperPH[i]*exp(-xtap*zTaperPH[k]*dt);
    axE202TaperzPHxPH[tk][ti] = (xfac*zTaperPH[k]+xTaperPH[i])*(bxE202TaperzPHxPH[tk][ti]-1.)/(kxTaperPH[i]*(xfac*zTaperPH[k]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE202TaperzPHxPH[tk][ti]>amax?axE202TaperzPHxPH[tk][ti]:amax;
amin = axE202TaperzPHxPH[tk][ti]<amin?axE202TaperzPHxPH[tk][ti]:amin;
bmax = bxE202TaperzPHxPH[tk][ti]>bmax?bxE202TaperzPHxPH[tk][ti]:bmax;
bmin = bxE202TaperzPHxPH[tk][ti]<bmin?bxE202TaperzPHxPH[tk][ti]:bmin;
    bxE202TaperxPH[tk][ti] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt);
    axE202TaperxPH[tk][ti] = (xfac*zTaper[k]+xTaperPH[i])*(bxE202TaperxPH[tk][ti]-1.)/(kxTaperPH[i]*(xfac*zTaper[k]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE202TaperxPH[tk][ti]>amax?axE202TaperxPH[tk][ti]:amax;
amin = axE202TaperxPH[tk][ti]<amin?axE202TaperxPH[tk][ti]:amin;
bmax = bxE202TaperxPH[tk][ti]>bmax?bxE202TaperxPH[tk][ti]:bmax;
bmin = bxE202TaperxPH[tk][ti]<bmin?bxE202TaperxPH[tk][ti]:bmin;
    byE202Taper[tk][ti] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
    ayE202Taper[tk][ti] = byE202Taper[tk][ti]-1.;
amax = ayE202Taper[tk][ti]>amax?ayE202Taper[tk][ti]:amax;
amin = ayE202Taper[tk][ti]<amin?ayE202Taper[tk][ti]:amin;
bmax = byE202Taper[tk][ti]>bmax?byE202Taper[tk][ti]:bmax;
bmin = byE202Taper[tk][ti]<bmin?byE202Taper[tk][ti]:bmin;
    byE202TaperzPH[tk][ti] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
    ayE202TaperzPH[tk][ti] = byE202TaperzPH[tk][ti]-1.;
amax = ayE202TaperzPH[tk][ti]>amax?ayE202TaperzPH[tk][ti]:amax;
amin = ayE202TaperzPH[tk][ti]<amin?ayE202TaperzPH[tk][ti]:amin;
bmax = byE202TaperzPH[tk][ti]>bmax?byE202TaperzPH[tk][ti]:bmax;
bmin = byE202TaperzPH[tk][ti]<bmin?byE202TaperzPH[tk][ti]:bmin;
    byE202TaperzPHxPH[tk][ti] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
    ayE202TaperzPHxPH[tk][ti] = byE202TaperzPHxPH[tk][ti]-1.;
amax = ayE202TaperzPHxPH[tk][ti]>amax?ayE202TaperzPHxPH[tk][ti]:amax;
amin = ayE202TaperzPHxPH[tk][ti]<amin?ayE202TaperzPHxPH[tk][ti]:amin;
bmax = byE202TaperzPHxPH[tk][ti]>bmax?byE202TaperzPHxPH[tk][ti]:bmax;
bmin = byE202TaperzPHxPH[tk][ti]<bmin?byE202TaperzPHxPH[tk][ti]:bmin;
    byE202TaperxPH[tk][ti] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
    ayE202TaperxPH[tk][ti] = byE202TaperxPH[tk][ti]-1.;
amax = ayE202TaperxPH[tk][ti]>amax?ayE202TaperxPH[tk][ti]:amax;
amin = ayE202TaperxPH[tk][ti]<amin?ayE202TaperxPH[tk][ti]:amin;
bmax = byE202TaperxPH[tk][ti]>bmax?byE202TaperxPH[tk][ti]:bmax;
bmin = byE202TaperxPH[tk][ti]<bmin?byE202TaperxPH[tk][ti]:bmin;
    bzE202Taper[tk][ti] = bzTaper[k]*exp(-xtap*xTaper[i]*dt);
    azE202Taper[tk][ti] = (xfac*xTaper[i]+zTaper[k])*(bzE202Taper[tk][ti]-1.)/(kzTaper[k]*(xfac*xTaper[i]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE202Taper[tk][ti]>amax?azE202Taper[tk][ti]:amax;
amin = azE202Taper[tk][ti]<amin?azE202Taper[tk][ti]:amin;
bmax = bzE202Taper[tk][ti]>bmax?bzE202Taper[tk][ti]:bmax;
bmin = bzE202Taper[tk][ti]<bmin?bzE202Taper[tk][ti]:bmin;
    bzE202TaperzPH[tk][ti] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt);
    azE202TaperzPH[tk][ti] = (xfac*xTaper[i]+zTaperPH[k])*(bzE202TaperzPH[tk][ti]-1.)/(kzTaperPH[k]*(xfac*xTaper[i]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE202TaperzPH[tk][ti]>amax?azE202TaperzPH[tk][ti]:amax;
amin = azE202TaperzPH[tk][ti]<amin?azE202TaperzPH[tk][ti]:amin;
bmax = bzE202TaperzPH[tk][ti]>bmax?bzE202TaperzPH[tk][ti]:bmax;
bmin = bzE202TaperzPH[tk][ti]<bmin?bzE202TaperzPH[tk][ti]:bmin;
    bzE202TaperzPHxPH[tk][ti] = bzTaperPH[k]*exp(-xtap*xTaperPH[i]*dt);
    azE202TaperzPHxPH[tk][ti] = (xfac*xTaperPH[i]+zTaperPH[k])*(bzE202TaperzPHxPH[tk][ti]-1.)/(kzTaperPH[k]*(xfac*xTaperPH[i]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE202TaperzPHxPH[tk][ti]>amax?azE202TaperzPHxPH[tk][ti]:amax;
amin = azE202TaperzPHxPH[tk][ti]<amin?azE202TaperzPHxPH[tk][ti]:amin;
bmax = bzE202TaperzPHxPH[tk][ti]>bmax?bzE202TaperzPHxPH[tk][ti]:bmax;
bmin = bzE202TaperzPHxPH[tk][ti]<bmin?bzE202TaperzPHxPH[tk][ti]:bmin;
    bzE202TaperxPH[tk][ti] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt);
    azE202TaperxPH[tk][ti] = (xfac*xTaperPH[i]+zTaper[k])*(bzE202TaperxPH[tk][ti]-1.)/(kzTaper[k]*(xfac*xTaperPH[i]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE202TaperxPH[tk][ti]>amax?azE202TaperxPH[tk][ti]:amax;
amin = azE202TaperxPH[tk][ti]<amin?azE202TaperxPH[tk][ti]:amin;
bmax = bzE202TaperxPH[tk][ti]>bmax?bzE202TaperxPH[tk][ti]:bmax;
bmin = bzE202TaperxPH[tk][ti]<bmin?bzE202TaperxPH[tk][ti]:bmin;
  }
}
for(int k=0;k<_nZmin;++k) {
  for(int j=0;j<_nYmin;++j) {
    bxE011Taper[k][j] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
    axE011Taper[k][j] = bxE011Taper[k][j]-1.;
amax = axE011Taper[k][j]>amax?axE011Taper[k][j]:amax;
amin = axE011Taper[k][j]<amin?axE011Taper[k][j]:amin;
bmax = bxE011Taper[k][j]>bmax?bxE011Taper[k][j]:bmax;
bmin = bxE011Taper[k][j]<bmin?bxE011Taper[k][j]:bmin;
    bxE011TaperzPH[k][j] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
    axE011TaperzPH[k][j] = bxE011TaperzPH[k][j]-1.;
amax = axE011TaperzPH[k][j]>amax?axE011TaperzPH[k][j]:amax;
amin = axE011TaperzPH[k][j]<amin?axE011TaperzPH[k][j]:amin;
bmax = bxE011TaperzPH[k][j]>bmax?bxE011TaperzPH[k][j]:bmax;
bmin = bxE011TaperzPH[k][j]<bmin?bxE011TaperzPH[k][j]:bmin;
    bxE011TaperzPHyPH[k][j] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
    axE011TaperzPHyPH[k][j] = bxE011TaperzPHyPH[k][j]-1.;
amax = axE011TaperzPHyPH[k][j]>amax?axE011TaperzPHyPH[k][j]:amax;
amin = axE011TaperzPHyPH[k][j]<amin?axE011TaperzPHyPH[k][j]:amin;
bmax = bxE011TaperzPHyPH[k][j]>bmax?bxE011TaperzPHyPH[k][j]:bmax;
bmin = bxE011TaperzPHyPH[k][j]<bmin?bxE011TaperzPHyPH[k][j]:bmin;
    bxE011TaperyPH[k][j] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
    axE011TaperyPH[k][j] = bxE011TaperyPH[k][j]-1.;
amax = axE011TaperyPH[k][j]>amax?axE011TaperyPH[k][j]:amax;
amin = axE011TaperyPH[k][j]<amin?axE011TaperyPH[k][j]:amin;
bmax = bxE011TaperyPH[k][j]>bmax?bxE011TaperyPH[k][j]:bmax;
bmin = bxE011TaperyPH[k][j]<bmin?bxE011TaperyPH[k][j]:bmin;
    byE011Taper[k][j] = byTaper[j]*exp(-xtap*zTaper[k]*dt);
    ayE011Taper[k][j] = (xfac*zTaper[k]+yTaper[j])*(byE011Taper[k][j]-1.)/(kyTaper[j]*(xfac*zTaper[k]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE011Taper[k][j]>amax?ayE011Taper[k][j]:amax;
amin = ayE011Taper[k][j]<amin?ayE011Taper[k][j]:amin;
bmax = byE011Taper[k][j]>bmax?byE011Taper[k][j]:bmax;
bmin = byE011Taper[k][j]<bmin?byE011Taper[k][j]:bmin;
    byE011TaperzPH[k][j] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt);
    ayE011TaperzPH[k][j] = (xfac*zTaperPH[k]+yTaper[j])*(byE011TaperzPH[k][j]-1.)/(kyTaper[j]*(xfac*zTaperPH[k]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE011TaperzPH[k][j]>amax?ayE011TaperzPH[k][j]:amax;
amin = ayE011TaperzPH[k][j]<amin?ayE011TaperzPH[k][j]:amin;
bmax = byE011TaperzPH[k][j]>bmax?byE011TaperzPH[k][j]:bmax;
bmin = byE011TaperzPH[k][j]<bmin?byE011TaperzPH[k][j]:bmin;
    byE011TaperzPHyPH[k][j] = byTaperPH[j]*exp(-xtap*zTaperPH[k]*dt);
    ayE011TaperzPHyPH[k][j] = (xfac*zTaperPH[k]+yTaperPH[j])*(byE011TaperzPHyPH[k][j]-1.)/(kyTaperPH[j]*(xfac*zTaperPH[k]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE011TaperzPHyPH[k][j]>amax?ayE011TaperzPHyPH[k][j]:amax;
amin = ayE011TaperzPHyPH[k][j]<amin?ayE011TaperzPHyPH[k][j]:amin;
bmax = byE011TaperzPHyPH[k][j]>bmax?byE011TaperzPHyPH[k][j]:bmax;
bmin = byE011TaperzPHyPH[k][j]<bmin?byE011TaperzPHyPH[k][j]:bmin;
    byE011TaperyPH[k][j] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt);
    ayE011TaperyPH[k][j] = (xfac*zTaper[k]+yTaperPH[j])*(byE011TaperyPH[k][j]-1.)/(kyTaperPH[j]*(xfac*zTaper[k]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE011TaperyPH[k][j]>amax?ayE011TaperyPH[k][j]:amax;
amin = ayE011TaperyPH[k][j]<amin?ayE011TaperyPH[k][j]:amin;
bmax = byE011TaperyPH[k][j]>bmax?byE011TaperyPH[k][j]:bmax;
bmin = byE011TaperyPH[k][j]<bmin?byE011TaperyPH[k][j]:bmin;
    bzE011Taper[k][j] = bzTaper[k]*exp(-xtap*yTaper[j]*dt);
    azE011Taper[k][j] = (xfac*yTaper[j]+zTaper[k])*(bzE011Taper[k][j]-1.)/(kzTaper[k]*(xfac*yTaper[j]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE011Taper[k][j]>amax?azE011Taper[k][j]:amax;
amin = azE011Taper[k][j]<amin?azE011Taper[k][j]:amin;
bmax = bzE011Taper[k][j]>bmax?bzE011Taper[k][j]:bmax;
bmin = bzE011Taper[k][j]<bmin?bzE011Taper[k][j]:bmin;
    bzE011TaperzPH[k][j] = bzTaperPH[k]*exp(-xtap*yTaper[j]*dt);
    azE011TaperzPH[k][j] = (xfac*yTaper[j]+zTaperPH[k])*(bzE011TaperzPH[k][j]-1.)/(kzTaperPH[k]*(xfac*yTaper[j]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE011TaperzPH[k][j]>amax?azE011TaperzPH[k][j]:amax;
amin = azE011TaperzPH[k][j]<amin?azE011TaperzPH[k][j]:amin;
bmax = bzE011TaperzPH[k][j]>bmax?bzE011TaperzPH[k][j]:bmax;
bmin = bzE011TaperzPH[k][j]<bmin?bzE011TaperzPH[k][j]:bmin;
    bzE011TaperzPHyPH[k][j] = bzTaperPH[k]*exp(-xtap*yTaperPH[j]*dt);
    azE011TaperzPHyPH[k][j] = (xfac*yTaperPH[j]+zTaperPH[k])*(bzE011TaperzPHyPH[k][j]-1.)/(kzTaperPH[k]*(xfac*yTaperPH[j]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE011TaperzPHyPH[k][j]>amax?azE011TaperzPHyPH[k][j]:amax;
amin = azE011TaperzPHyPH[k][j]<amin?azE011TaperzPHyPH[k][j]:amin;
bmax = bzE011TaperzPHyPH[k][j]>bmax?bzE011TaperzPHyPH[k][j]:bmax;
bmin = bzE011TaperzPHyPH[k][j]<bmin?bzE011TaperzPHyPH[k][j]:bmin;
    bzE011TaperyPH[k][j] = bzTaper[k]*exp(-xtap*yTaperPH[j]*dt);
    azE011TaperyPH[k][j] = (xfac*yTaperPH[j]+zTaper[k])*(bzE011TaperyPH[k][j]-1.)/(kzTaper[k]*(xfac*yTaperPH[j]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE011TaperyPH[k][j]>amax?azE011TaperyPH[k][j]:amax;
amin = azE011TaperyPH[k][j]<amin?azE011TaperyPH[k][j]:amin;
bmax = bzE011TaperyPH[k][j]>bmax?bzE011TaperyPH[k][j]:bmax;
bmin = bzE011TaperyPH[k][j]<bmin?bzE011TaperyPH[k][j]:bmin;
  }
  for(int j=yEnd1;j<NY;++j) {
    int tj = j-yEnd1;
    bxE021Taper[k][tj] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
    axE021Taper[k][tj] = bxE021Taper[k][tj]-1.;
amax = axE021Taper[k][tj]>amax?axE021Taper[k][tj]:amax;
amin = axE021Taper[k][tj]<amin?axE021Taper[k][tj]:amin;
bmax = bxE021Taper[k][tj]>bmax?bxE021Taper[k][tj]:bmax;
bmin = bxE021Taper[k][tj]<bmin?bxE021Taper[k][tj]:bmin;
    bxE021TaperzPH[k][tj] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
    axE021TaperzPH[k][tj] = bxE021TaperzPH[k][tj]-1.;
amax = axE021TaperzPH[k][tj]>amax?axE021TaperzPH[k][tj]:amax;
amin = axE021TaperzPH[k][tj]<amin?axE021TaperzPH[k][tj]:amin;
bmax = bxE021TaperzPH[k][tj]>bmax?bxE021TaperzPH[k][tj]:bmax;
bmin = bxE021TaperzPH[k][tj]<bmin?bxE021TaperzPH[k][tj]:bmin;
    bxE021TaperzPHyPH[k][tj] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
    axE021TaperzPHyPH[k][tj] = bxE021TaperzPHyPH[k][tj]-1.;
amax = axE021TaperzPHyPH[k][tj]>amax?axE021TaperzPHyPH[k][tj]:amax;
amin = axE021TaperzPHyPH[k][tj]<amin?axE021TaperzPHyPH[k][tj]:amin;
bmax = bxE021TaperzPHyPH[k][tj]>bmax?bxE021TaperzPHyPH[k][tj]:bmax;
bmin = bxE021TaperzPHyPH[k][tj]<bmin?bxE021TaperzPHyPH[k][tj]:bmin;
    bxE021TaperyPH[k][tj] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
    axE021TaperyPH[k][tj] = bxE021TaperyPH[k][tj]-1.;
amax = axE021TaperyPH[k][tj]>amax?axE021TaperyPH[k][tj]:amax;
amin = axE021TaperyPH[k][tj]<amin?axE021TaperyPH[k][tj]:amin;
bmax = bxE021TaperyPH[k][tj]>bmax?bxE021TaperyPH[k][tj]:bmax;
bmin = bxE021TaperyPH[k][tj]<bmin?bxE021TaperyPH[k][tj]:bmin;
    byE021Taper[k][tj] = byTaper[j]*exp(-xtap*zTaper[k]*dt);
    ayE021Taper[k][tj] = (xfac*zTaper[k]+yTaper[j])*(byE021Taper[k][tj]-1.)/(kyTaper[j]*(xfac*zTaper[k]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE021Taper[k][tj]>amax?ayE021Taper[k][tj]:amax;
amin = ayE021Taper[k][tj]<amin?ayE021Taper[k][tj]:amin;
bmax = byE021Taper[k][tj]>bmax?byE021Taper[k][tj]:bmax;
bmin = byE021Taper[k][tj]<bmin?byE021Taper[k][tj]:bmin;
    byE021TaperzPH[k][tj] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt);
    ayE021TaperzPH[k][tj] = (xfac*zTaperPH[k]+yTaper[j])*(byE021TaperzPH[k][tj]-1.)/(kyTaper[j]*(xfac*zTaperPH[k]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE021TaperzPH[k][tj]>amax?ayE021TaperzPH[k][tj]:amax;
amin = ayE021TaperzPH[k][tj]<amin?ayE021TaperzPH[k][tj]:amin;
bmax = byE021TaperzPH[k][tj]>bmax?byE021TaperzPH[k][tj]:bmax;
bmin = byE021TaperzPH[k][tj]<bmin?byE021TaperzPH[k][tj]:bmin;
    byE021TaperzPHyPH[k][tj] = byTaperPH[j]*exp(-xtap*zTaperPH[k]*dt);
    ayE021TaperzPHyPH[k][tj] = (xfac*zTaperPH[k]+yTaperPH[j])*(byE021TaperzPHyPH[k][tj]-1.)/(kyTaperPH[j]*(xfac*zTaperPH[k]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE021TaperzPHyPH[k][tj]>amax?ayE021TaperzPHyPH[k][tj]:amax;
amin = ayE021TaperzPHyPH[k][tj]<amin?ayE021TaperzPHyPH[k][tj]:amin;
bmax = byE021TaperzPHyPH[k][tj]>bmax?byE021TaperzPHyPH[k][tj]:bmax;
bmin = byE021TaperzPHyPH[k][tj]<bmin?byE021TaperzPHyPH[k][tj]:bmin;
    byE021TaperyPH[k][tj] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt);
    ayE021TaperyPH[k][tj] = (xfac*zTaper[k]+yTaperPH[j])*(byE021TaperyPH[k][tj]-1.)/(kyTaperPH[j]*(xfac*zTaper[k]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE021TaperyPH[k][tj]>amax?ayE021TaperyPH[k][tj]:amax;
amin = ayE021TaperyPH[k][tj]<amin?ayE021TaperyPH[k][tj]:amin;
bmax = byE021TaperyPH[k][tj]>bmax?byE021TaperyPH[k][tj]:bmax;
bmin = byE021TaperyPH[k][tj]<bmin?byE021TaperyPH[k][tj]:bmin;
    bzE021Taper[k][tj] = bzTaper[k]*exp(-xtap*yTaper[j]*dt);
    azE021Taper[k][tj] = (xfac*yTaper[j]+zTaper[k])*(bzE021Taper[k][tj]-1.)/(kzTaper[k]*(xfac*yTaper[j]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE021Taper[k][tj]>amax?azE021Taper[k][tj]:amax;
amin = azE021Taper[k][tj]<amin?azE021Taper[k][tj]:amin;
bmax = bzE021Taper[k][tj]>bmax?bzE021Taper[k][tj]:bmax;
bmin = bzE021Taper[k][tj]<bmin?bzE021Taper[k][tj]:bmin;
    bzE021TaperzPH[k][tj] = bzTaperPH[k]*exp(-xtap*yTaper[j]*dt);
    azE021TaperzPH[k][tj] = (xfac*yTaper[j]+zTaperPH[k])*(bzE021TaperzPH[k][tj]-1.)/(kzTaperPH[k]*(xfac*yTaper[j]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE021TaperzPH[k][tj]>amax?azE021TaperzPH[k][tj]:amax;
amin = azE021TaperzPH[k][tj]<amin?azE021TaperzPH[k][tj]:amin;
bmax = bzE021TaperzPH[k][tj]>bmax?bzE021TaperzPH[k][tj]:bmax;
bmin = bzE021TaperzPH[k][tj]<bmin?bzE021TaperzPH[k][tj]:bmin;
    bzE021TaperzPHyPH[k][tj] = bzTaperPH[k]*exp(-xtap*yTaperPH[j]*dt);
    azE021TaperzPHyPH[k][tj] = (xfac*yTaperPH[j]+zTaperPH[k])*(bzE021TaperzPHyPH[k][tj]-1.)/(kzTaperPH[k]*(xfac*yTaperPH[j]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE021TaperzPHyPH[k][tj]>amax?azE021TaperzPHyPH[k][tj]:amax;
amin = azE021TaperzPHyPH[k][tj]<amin?azE021TaperzPHyPH[k][tj]:amin;
bmax = bzE021TaperzPHyPH[k][tj]>bmax?bzE021TaperzPHyPH[k][tj]:bmax;
bmin = bzE021TaperzPHyPH[k][tj]<bmin?bzE021TaperzPHyPH[k][tj]:bmin;
    bzE021TaperyPH[k][tj] = bzTaper[k]*exp(-xtap*yTaperPH[j]*dt);
    azE021TaperyPH[k][tj] = (xfac*yTaperPH[j]+zTaper[k])*(bzE021TaperyPH[k][tj]-1.)/(kzTaper[k]*(xfac*yTaperPH[j]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE021TaperyPH[k][tj]>amax?azE021TaperyPH[k][tj]:amax;
amin = azE021TaperyPH[k][tj]<amin?azE021TaperyPH[k][tj]:amin;
bmax = bzE021TaperyPH[k][tj]>bmax?bzE021TaperyPH[k][tj]:bmax;
bmin = bzE021TaperyPH[k][tj]<bmin?bzE021TaperyPH[k][tj]:bmin;
  }
}
for(int k=zEnd1;k<NZ;++k) {
  int tk = k-zEnd1;
  for(int j=0;j<_nYmin;++j) {
    bxE012Taper[tk][j] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
    axE012Taper[tk][j] = bxE012Taper[tk][j]-1.;
amax = axE012Taper[tk][j]>amax?axE012Taper[tk][j]:amax;
amin = axE012Taper[tk][j]<amin?axE012Taper[tk][j]:amin;
bmax = bxE012Taper[tk][j]>bmax?bxE012Taper[tk][j]:bmax;
bmin = bxE012Taper[tk][j]<bmin?bxE012Taper[tk][j]:bmin;
    bxE012TaperzPH[tk][j] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
    axE012TaperzPH[tk][j] = bxE012TaperzPH[tk][j]-1.;
amax = axE012TaperzPH[tk][j]>amax?axE012TaperzPH[tk][j]:amax;
amin = axE012TaperzPH[tk][j]<amin?axE012TaperzPH[tk][j]:amin;
bmax = bxE012TaperzPH[tk][j]>bmax?bxE012TaperzPH[tk][j]:bmax;
bmin = bxE012TaperzPH[tk][j]<bmin?bxE012TaperzPH[tk][j]:bmin;
    bxE012TaperzPHyPH[tk][j] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
    axE012TaperzPHyPH[tk][j] = bxE012TaperzPHyPH[tk][j]-1.;
amax = axE012TaperzPHyPH[tk][j]>amax?axE012TaperzPHyPH[tk][j]:amax;
amin = axE012TaperzPHyPH[tk][j]<amin?axE012TaperzPHyPH[tk][j]:amin;
bmax = bxE012TaperzPHyPH[tk][j]>bmax?bxE012TaperzPHyPH[tk][j]:bmax;
bmin = bxE012TaperzPHyPH[tk][j]<bmin?bxE012TaperzPHyPH[tk][j]:bmin;
    bxE012TaperyPH[tk][j] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
    axE012TaperyPH[tk][j] = bxE012TaperyPH[tk][j]-1.;
amax = axE012TaperyPH[tk][j]>amax?axE012TaperyPH[tk][j]:amax;
amin = axE012TaperyPH[tk][j]<amin?axE012TaperyPH[tk][j]:amin;
bmax = bxE012TaperyPH[tk][j]>bmax?bxE012TaperyPH[tk][j]:bmax;
bmin = bxE012TaperyPH[tk][j]<bmin?bxE012TaperyPH[tk][j]:bmin;
    byE012Taper[tk][j] = byTaper[j]*exp(-xtap*zTaper[k]*dt);
    ayE012Taper[tk][j] = (xfac*zTaper[k]+yTaper[j])*(byE012Taper[tk][j]-1.)/(kyTaper[j]*(xfac*zTaper[k]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE012Taper[tk][j]>amax?ayE012Taper[tk][j]:amax;
amin = ayE012Taper[tk][j]<amin?ayE012Taper[tk][j]:amin;
bmax = byE012Taper[tk][j]>bmax?byE012Taper[tk][j]:bmax;
bmin = byE012Taper[tk][j]<bmin?byE012Taper[tk][j]:bmin;
    byE012TaperzPH[tk][j] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt);
    ayE012TaperzPH[tk][j] = (xfac*zTaperPH[k]+yTaper[j])*(byE012TaperzPH[tk][j]-1.)/(kyTaper[j]*(xfac*zTaperPH[k]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE012TaperzPH[tk][j]>amax?ayE012TaperzPH[tk][j]:amax;
amin = ayE012TaperzPH[tk][j]<amin?ayE012TaperzPH[tk][j]:amin;
bmax = byE012TaperzPH[tk][j]>bmax?byE012TaperzPH[tk][j]:bmax;
bmin = byE012TaperzPH[tk][j]<bmin?byE012TaperzPH[tk][j]:bmin;
    byE012TaperzPHyPH[tk][j] = byTaperPH[j]*exp(-xtap*zTaperPH[k]*dt);
    ayE012TaperzPHyPH[tk][j] = (xfac*zTaperPH[k]+yTaperPH[j])*(byE012TaperzPHyPH[tk][j]-1.)/(kyTaperPH[j]*(xfac*zTaperPH[k]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE012TaperzPHyPH[tk][j]>amax?ayE012TaperzPHyPH[tk][j]:amax;
amin = ayE012TaperzPHyPH[tk][j]<amin?ayE012TaperzPHyPH[tk][j]:amin;
bmax = byE012TaperzPHyPH[tk][j]>bmax?byE012TaperzPHyPH[tk][j]:bmax;
bmin = byE012TaperzPHyPH[tk][j]<bmin?byE012TaperzPHyPH[tk][j]:bmin;
    byE012TaperyPH[tk][j] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt);
    ayE012TaperyPH[tk][j] = (xfac*zTaper[k]+yTaperPH[j])*(byE012TaperyPH[tk][j]-1.)/(kyTaperPH[j]*(xfac*zTaper[k]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE012TaperyPH[tk][j]>amax?ayE012TaperyPH[tk][j]:amax;
amin = ayE012TaperyPH[tk][j]<amin?ayE012TaperyPH[tk][j]:amin;
bmax = byE012TaperyPH[tk][j]>bmax?byE012TaperyPH[tk][j]:bmax;
bmin = byE012TaperyPH[tk][j]<bmin?byE012TaperyPH[tk][j]:bmin;
    bzE012Taper[tk][j] = bzTaper[k]*exp(-xtap*yTaper[j]*dt);
    azE012Taper[tk][j] = (xfac*yTaper[j]+zTaper[k])*(bzE012Taper[tk][j]-1.)/(kzTaper[k]*(xfac*yTaper[j]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE012Taper[tk][j]>amax?azE012Taper[tk][j]:amax;
amin = azE012Taper[tk][j]<amin?azE012Taper[tk][j]:amin;
bmax = bzE012Taper[tk][j]>bmax?bzE012Taper[tk][j]:bmax;
bmin = bzE012Taper[tk][j]<bmin?bzE012Taper[tk][j]:bmin;
    bzE012TaperzPH[tk][j] = bzTaperPH[k]*exp(-xtap*yTaper[j]*dt);
    azE012TaperzPH[tk][j] = (xfac*yTaper[j]+zTaperPH[k])*(bzE012TaperzPH[tk][j]-1.)/(kzTaperPH[k]*(xfac*yTaper[j]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE012TaperzPH[tk][j]>amax?azE012TaperzPH[tk][j]:amax;
amin = azE012TaperzPH[tk][j]<amin?azE012TaperzPH[tk][j]:amin;
bmax = bzE012TaperzPH[tk][j]>bmax?bzE012TaperzPH[tk][j]:bmax;
bmin = bzE012TaperzPH[tk][j]<bmin?bzE012TaperzPH[tk][j]:bmin;
    bzE012TaperzPHyPH[tk][j] = bzTaperPH[k]*exp(-xtap*yTaperPH[j]*dt);
    azE012TaperzPHyPH[tk][j] = (xfac*yTaperPH[j]+zTaperPH[k])*(bzE012TaperzPHyPH[tk][j]-1.)/(kzTaperPH[k]*(xfac*yTaperPH[j]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE012TaperzPHyPH[tk][j]>amax?azE012TaperzPHyPH[tk][j]:amax;
amin = azE012TaperzPHyPH[tk][j]<amin?azE012TaperzPHyPH[tk][j]:amin;
bmax = bzE012TaperzPHyPH[tk][j]>bmax?bzE012TaperzPHyPH[tk][j]:bmax;
bmin = bzE012TaperzPHyPH[tk][j]<bmin?bzE012TaperzPHyPH[tk][j]:bmin;
    bzE012TaperyPH[tk][j] = bzTaper[k]*exp(-xtap*yTaperPH[j]*dt);
    azE012TaperyPH[tk][j] = (xfac*yTaperPH[j]+zTaper[k])*(bzE012TaperyPH[tk][j]-1.)/(kzTaper[k]*(xfac*yTaperPH[j]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE012TaperyPH[tk][j]>amax?azE012TaperyPH[tk][j]:amax;
amin = azE012TaperyPH[tk][j]<amin?azE012TaperyPH[tk][j]:amin;
bmax = bzE012TaperyPH[tk][j]>bmax?bzE012TaperyPH[tk][j]:bmax;
bmin = bzE012TaperyPH[tk][j]<bmin?bzE012TaperyPH[tk][j]:bmin;
  }
  for(int j=yEnd1;j<NY;++j) {
    int tj = j-yEnd1;
    bxE022Taper[tk][tj] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
    axE022Taper[tk][tj] = bxE022Taper[tk][tj]-1.;
amax = axE022Taper[tk][tj]>amax?axE022Taper[tk][tj]:amax;
amin = axE022Taper[tk][tj]<amin?axE022Taper[tk][tj]:amin;
bmax = bxE022Taper[tk][tj]>bmax?bxE022Taper[tk][tj]:bmax;
bmin = bxE022Taper[tk][tj]<bmin?bxE022Taper[tk][tj]:bmin;
    bxE022TaperzPH[tk][tj] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
    axE022TaperzPH[tk][tj] = bxE022TaperzPH[tk][tj]-1.;
amax = axE022TaperzPH[tk][tj]>amax?axE022TaperzPH[tk][tj]:amax;
amin = axE022TaperzPH[tk][tj]<amin?axE022TaperzPH[tk][tj]:amin;
bmax = bxE022TaperzPH[tk][tj]>bmax?bxE022TaperzPH[tk][tj]:bmax;
bmin = bxE022TaperzPH[tk][tj]<bmin?bxE022TaperzPH[tk][tj]:bmin;
    bxE022TaperzPHyPH[tk][tj] = exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
    axE022TaperzPHyPH[tk][tj] = bxE022TaperzPHyPH[tk][tj]-1.;
amax = axE022TaperzPHyPH[tk][tj]>amax?axE022TaperzPHyPH[tk][tj]:amax;
amin = axE022TaperzPHyPH[tk][tj]<amin?axE022TaperzPHyPH[tk][tj]:amin;
bmax = bxE022TaperzPHyPH[tk][tj]>bmax?bxE022TaperzPHyPH[tk][tj]:bmax;
bmin = bxE022TaperzPHyPH[tk][tj]<bmin?bxE022TaperzPHyPH[tk][tj]:bmin;
    bxE022TaperyPH[tk][tj] = exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
    axE022TaperyPH[tk][tj] = bxE022TaperyPH[tk][tj]-1.;
amax = axE022TaperyPH[tk][tj]>amax?axE022TaperyPH[tk][tj]:amax;
amin = axE022TaperyPH[tk][tj]<amin?axE022TaperyPH[tk][tj]:amin;
bmax = bxE022TaperyPH[tk][tj]>bmax?bxE022TaperyPH[tk][tj]:bmax;
bmin = bxE022TaperyPH[tk][tj]<bmin?bxE022TaperyPH[tk][tj]:bmin;
    byE022Taper[tk][tj] = byTaper[j]*exp(-xtap*zTaper[k]*dt);
    ayE022Taper[tk][tj] = (xfac*zTaper[k]+yTaper[j])*(byE022Taper[tk][tj]-1.)/(kyTaper[j]*(xfac*zTaper[k]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE022Taper[tk][tj]>amax?ayE022Taper[tk][tj]:amax;
amin = ayE022Taper[tk][tj]<amin?ayE022Taper[tk][tj]:amin;
bmax = byE022Taper[tk][tj]>bmax?byE022Taper[tk][tj]:bmax;
bmin = byE022Taper[tk][tj]<bmin?byE022Taper[tk][tj]:bmin;
    byE022TaperzPH[tk][tj] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt);
    ayE022TaperzPH[tk][tj] = (xfac*zTaperPH[k]+yTaper[j])*(byE022TaperzPH[tk][tj]-1.)/(kyTaper[j]*(xfac*zTaperPH[k]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE022TaperzPH[tk][tj]>amax?ayE022TaperzPH[tk][tj]:amax;
amin = ayE022TaperzPH[tk][tj]<amin?ayE022TaperzPH[tk][tj]:amin;
bmax = byE022TaperzPH[tk][tj]>bmax?byE022TaperzPH[tk][tj]:bmax;
bmin = byE022TaperzPH[tk][tj]<bmin?byE022TaperzPH[tk][tj]:bmin;
    byE022TaperzPHyPH[tk][tj] = byTaperPH[j]*exp(-xtap*zTaperPH[k]*dt);
    ayE022TaperzPHyPH[tk][tj] = (xfac*zTaperPH[k]+yTaperPH[j])*(byE022TaperzPHyPH[tk][tj]-1.)/(kyTaperPH[j]*(xfac*zTaperPH[k]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE022TaperzPHyPH[tk][tj]>amax?ayE022TaperzPHyPH[tk][tj]:amax;
amin = ayE022TaperzPHyPH[tk][tj]<amin?ayE022TaperzPHyPH[tk][tj]:amin;
bmax = byE022TaperzPHyPH[tk][tj]>bmax?byE022TaperzPHyPH[tk][tj]:bmax;
bmin = byE022TaperzPHyPH[tk][tj]<bmin?byE022TaperzPHyPH[tk][tj]:bmin;
    byE022TaperyPH[tk][tj] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt);
    ayE022TaperyPH[tk][tj] = (xfac*zTaper[k]+yTaperPH[j])*(byE022TaperyPH[tk][tj]-1.)/(kyTaperPH[j]*(xfac*zTaper[k]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE022TaperyPH[tk][tj]>amax?ayE022TaperyPH[tk][tj]:amax;
amin = ayE022TaperyPH[tk][tj]<amin?ayE022TaperyPH[tk][tj]:amin;
bmax = byE022TaperyPH[tk][tj]>bmax?byE022TaperyPH[tk][tj]:bmax;
bmin = byE022TaperyPH[tk][tj]<bmin?byE022TaperyPH[tk][tj]:bmin;
    bzE022Taper[tk][tj] = bzTaper[k]*exp(-xtap*yTaper[j]*dt);
    azE022Taper[tk][tj] = (xfac*yTaper[j]+zTaper[k])*(bzE022Taper[tk][tj]-1.)/(kzTaper[k]*(xfac*yTaper[j]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE022Taper[tk][tj]>amax?azE022Taper[tk][tj]:amax;
amin = azE022Taper[tk][tj]<amin?azE022Taper[tk][tj]:amin;
bmax = bzE022Taper[tk][tj]>bmax?bzE022Taper[tk][tj]:bmax;
bmin = bzE022Taper[tk][tj]<bmin?bzE022Taper[tk][tj]:bmin;
    bzE022TaperzPH[tk][tj] = bzTaperPH[k]*exp(-xtap*yTaper[j]*dt);
    azE022TaperzPH[tk][tj] = (xfac*yTaper[j]+zTaperPH[k])*(bzE022TaperzPH[tk][tj]-1.)/(kzTaperPH[k]*(xfac*yTaper[j]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE022TaperzPH[tk][tj]>amax?azE022TaperzPH[tk][tj]:amax;
amin = azE022TaperzPH[tk][tj]<amin?azE022TaperzPH[tk][tj]:amin;
bmax = bzE022TaperzPH[tk][tj]>bmax?bzE022TaperzPH[tk][tj]:bmax;
bmin = bzE022TaperzPH[tk][tj]<bmin?bzE022TaperzPH[tk][tj]:bmin;
    bzE022TaperzPHyPH[tk][tj] = bzTaperPH[k]*exp(-xtap*yTaperPH[j]*dt);
    azE022TaperzPHyPH[tk][tj] = (xfac*yTaperPH[j]+zTaperPH[k])*(bzE022TaperzPHyPH[tk][tj]-1.)/(kzTaperPH[k]*(xfac*yTaperPH[j]+zTaperPH[k]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azE022TaperzPHyPH[tk][tj]>amax?azE022TaperzPHyPH[tk][tj]:amax;
amin = azE022TaperzPHyPH[tk][tj]<amin?azE022TaperzPHyPH[tk][tj]:amin;
bmax = bzE022TaperzPHyPH[tk][tj]>bmax?bzE022TaperzPHyPH[tk][tj]:bmax;
bmin = bzE022TaperzPHyPH[tk][tj]<bmin?bzE022TaperzPHyPH[tk][tj]:bmin;
    bzE022TaperyPH[tk][tj] = bzTaper[k]*exp(-xtap*yTaperPH[j]*dt);
    azE022TaperyPH[tk][tj] = (xfac*yTaperPH[j]+zTaper[k])*(bzE022TaperyPH[tk][tj]-1.)/(kzTaper[k]*(xfac*yTaperPH[j]+zTaper[k]+kzTaper[k]*alphaZTaper[k]));
amax = azE022TaperyPH[tk][tj]>amax?azE022TaperyPH[tk][tj]:amax;
amin = azE022TaperyPH[tk][tj]<amin?azE022TaperyPH[tk][tj]:amin;
bmax = bzE022TaperyPH[tk][tj]>bmax?bzE022TaperyPH[tk][tj]:bmax;
bmin = bzE022TaperyPH[tk][tj]<bmin?bzE022TaperyPH[tk][tj]:bmin;
  }
}
for(int j=0;j<_nYmin;++j) {
  for(int i=0;i<_nXmin;++i) {
    bxE110Taper[j][i] = bxTaper[i]*exp(-xtap*yTaper[j]*dt);
    axE110Taper[j][i] = (xfac*yTaper[j]+xTaper[i])*(bxE110Taper[j][i]-1.)/(kxTaper[i]*(xfac*yTaper[j]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE110Taper[j][i]>amax?axE110Taper[j][i]:amax;
amin = axE110Taper[j][i]<amin?axE110Taper[j][i]:amin;
bmax = bxE110Taper[j][i]>bmax?bxE110Taper[j][i]:bmax;
bmin = bxE110Taper[j][i]<bmin?bxE110Taper[j][i]:bmin;
    bxE110TaperyPH[j][i] = bxTaper[i]*exp(-xtap*yTaperPH[j]*dt);
    axE110TaperyPH[j][i] = (xfac*yTaperPH[j]+xTaper[i])*(bxE110TaperyPH[j][i]-1.)/(kxTaper[i]*(xfac*yTaperPH[j]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE110TaperyPH[j][i]>amax?axE110TaperyPH[j][i]:amax;
amin = axE110TaperyPH[j][i]<amin?axE110TaperyPH[j][i]:amin;
bmax = bxE110TaperyPH[j][i]>bmax?bxE110TaperyPH[j][i]:bmax;
bmin = bxE110TaperyPH[j][i]<bmin?bxE110TaperyPH[j][i]:bmin;
    bxE110TaperyPHxPH[j][i] = bxTaperPH[i]*exp(-xtap*yTaperPH[j]*dt);
    axE110TaperyPHxPH[j][i] = (xfac*yTaperPH[j]+xTaperPH[i])*(bxE110TaperyPHxPH[j][i]-1.)/(kxTaperPH[i]*(xfac*yTaperPH[j]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE110TaperyPHxPH[j][i]>amax?axE110TaperyPHxPH[j][i]:amax;
amin = axE110TaperyPHxPH[j][i]<amin?axE110TaperyPHxPH[j][i]:amin;
bmax = bxE110TaperyPHxPH[j][i]>bmax?bxE110TaperyPHxPH[j][i]:bmax;
bmin = bxE110TaperyPHxPH[j][i]<bmin?bxE110TaperyPHxPH[j][i]:bmin;
    bxE110TaperxPH[j][i] = bxTaperPH[i]*exp(-xtap*yTaper[j]*dt);
    axE110TaperxPH[j][i] = (xfac*yTaper[j]+xTaperPH[i])*(bxE110TaperxPH[j][i]-1.)/(kxTaperPH[i]*(xfac*yTaper[j]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE110TaperxPH[j][i]>amax?axE110TaperxPH[j][i]:amax;
amin = axE110TaperxPH[j][i]<amin?axE110TaperxPH[j][i]:amin;
bmax = bxE110TaperxPH[j][i]>bmax?bxE110TaperxPH[j][i]:bmax;
bmin = bxE110TaperxPH[j][i]<bmin?bxE110TaperxPH[j][i]:bmin;
    byE110Taper[j][i] = byTaper[j]*exp(-xtap*xTaper[i]*dt);
    ayE110Taper[j][i] = (xfac*xTaper[i]+yTaper[j])*(byE110Taper[j][i]-1.)/(kyTaper[j]*(xfac*xTaper[i]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE110Taper[j][i]>amax?ayE110Taper[j][i]:amax;
amin = ayE110Taper[j][i]<amin?ayE110Taper[j][i]:amin;
bmax = byE110Taper[j][i]>bmax?byE110Taper[j][i]:bmax;
bmin = byE110Taper[j][i]<bmin?byE110Taper[j][i]:bmin;
    byE110TaperyPH[j][i] = byTaperPH[j]*exp(-xtap*xTaper[i]*dt);
    ayE110TaperyPH[j][i] = (xfac*xTaper[i]+yTaperPH[j])*(byE110TaperyPH[j][i]-1.)/(kyTaperPH[j]*(xfac*xTaper[i]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE110TaperyPH[j][i]>amax?ayE110TaperyPH[j][i]:amax;
amin = ayE110TaperyPH[j][i]<amin?ayE110TaperyPH[j][i]:amin;
bmax = byE110TaperyPH[j][i]>bmax?byE110TaperyPH[j][i]:bmax;
bmin = byE110TaperyPH[j][i]<bmin?byE110TaperyPH[j][i]:bmin;
    byE110TaperyPHxPH[j][i] = byTaperPH[j]*exp(-xtap*xTaperPH[i]*dt);
    ayE110TaperyPHxPH[j][i] = (xfac*xTaperPH[i]+yTaperPH[j])*(byE110TaperyPHxPH[j][i]-1.)/(kyTaperPH[j]*(xfac*xTaperPH[i]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE110TaperyPHxPH[j][i]>amax?ayE110TaperyPHxPH[j][i]:amax;
amin = ayE110TaperyPHxPH[j][i]<amin?ayE110TaperyPHxPH[j][i]:amin;
bmax = byE110TaperyPHxPH[j][i]>bmax?byE110TaperyPHxPH[j][i]:bmax;
bmin = byE110TaperyPHxPH[j][i]<bmin?byE110TaperyPHxPH[j][i]:bmin;
    byE110TaperxPH[j][i] = byTaper[j]*exp(-xtap*xTaperPH[i]*dt);
    ayE110TaperxPH[j][i] = (xfac*xTaperPH[i]+yTaper[j])*(byE110TaperxPH[j][i]-1.)/(kyTaper[j]*(xfac*xTaperPH[i]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE110TaperxPH[j][i]>amax?ayE110TaperxPH[j][i]:amax;
amin = ayE110TaperxPH[j][i]<amin?ayE110TaperxPH[j][i]:amin;
bmax = byE110TaperxPH[j][i]>bmax?byE110TaperxPH[j][i]:bmax;
bmin = byE110TaperxPH[j][i]<bmin?byE110TaperxPH[j][i]:bmin;
    bzE110Taper[j][i] = exp(-xtap*yTaper[j]*dt)*exp(-xtap*xTaper[i]*dt);
    azE110Taper[j][i] = bzE110Taper[j][i]-1.;
amax = azE110Taper[j][i]>amax?azE110Taper[j][i]:amax;
amin = azE110Taper[j][i]<amin?azE110Taper[j][i]:amin;
bmax = bzE110Taper[j][i]>bmax?bzE110Taper[j][i]:bmax;
bmin = bzE110Taper[j][i]<bmin?bzE110Taper[j][i]:bmin;
    bzE110TaperyPH[j][i] = exp(-xtap*yTaperPH[j]*dt)*exp(-xtap*xTaper[i]*dt);
    azE110TaperyPH[j][i] = bzE110TaperyPH[j][i]-1.;
amax = azE110TaperyPH[j][i]>amax?azE110TaperyPH[j][i]:amax;
amin = azE110TaperyPH[j][i]<amin?azE110TaperyPH[j][i]:amin;
bmax = bzE110TaperyPH[j][i]>bmax?bzE110TaperyPH[j][i]:bmax;
bmin = bzE110TaperyPH[j][i]<bmin?bzE110TaperyPH[j][i]:bmin;
    bzE110TaperyPHxPH[j][i] = exp(-xtap*yTaperPH[j]*dt)*exp(-xtap*xTaperPH[i]*dt);
    azE110TaperyPHxPH[j][i] = bzE110TaperyPHxPH[j][i]-1.;
amax = azE110TaperyPHxPH[j][i]>amax?azE110TaperyPHxPH[j][i]:amax;
amin = azE110TaperyPHxPH[j][i]<amin?azE110TaperyPHxPH[j][i]:amin;
bmax = bzE110TaperyPHxPH[j][i]>bmax?bzE110TaperyPHxPH[j][i]:bmax;
bmin = bzE110TaperyPHxPH[j][i]<bmin?bzE110TaperyPHxPH[j][i]:bmin;
    bzE110TaperxPH[j][i] = exp(-xtap*yTaper[j]*dt)*exp(-xtap*xTaperPH[i]*dt);
    azE110TaperxPH[j][i] = bzE110TaperxPH[j][i]-1.;
amax = azE110TaperxPH[j][i]>amax?azE110TaperxPH[j][i]:amax;
amin = azE110TaperxPH[j][i]<amin?azE110TaperxPH[j][i]:amin;
bmax = bzE110TaperxPH[j][i]>bmax?bzE110TaperxPH[j][i]:bmax;
bmin = bzE110TaperxPH[j][i]<bmin?bzE110TaperxPH[j][i]:bmin;
  }
  for(int i=xEnd1;i<NX;++i) {
    int ti = i-xEnd1;
    bxE210Taper[j][ti] = bxTaper[i]*exp(-xtap*yTaper[j]*dt);
    axE210Taper[j][ti] = (xfac*yTaper[j]+xTaper[i])*(bxE210Taper[j][ti]-1.)/(kxTaper[i]*(xfac*yTaper[j]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE210Taper[j][ti]>amax?axE210Taper[j][ti]:amax;
amin = axE210Taper[j][ti]<amin?axE210Taper[j][ti]:amin;
bmax = bxE210Taper[j][ti]>bmax?bxE210Taper[j][ti]:bmax;
bmin = bxE210Taper[j][ti]<bmin?bxE210Taper[j][ti]:bmin;
    bxE210TaperyPH[j][ti] = bxTaper[i]*exp(-xtap*yTaperPH[j]*dt);
    axE210TaperyPH[j][ti] = (xfac*yTaperPH[j]+xTaper[i])*(bxE210TaperyPH[j][ti]-1.)/(kxTaper[i]*(xfac*yTaperPH[j]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE210TaperyPH[j][ti]>amax?axE210TaperyPH[j][ti]:amax;
amin = axE210TaperyPH[j][ti]<amin?axE210TaperyPH[j][ti]:amin;
bmax = bxE210TaperyPH[j][ti]>bmax?bxE210TaperyPH[j][ti]:bmax;
bmin = bxE210TaperyPH[j][ti]<bmin?bxE210TaperyPH[j][ti]:bmin;
    bxE210TaperyPHxPH[j][ti] = bxTaperPH[i]*exp(-xtap*yTaperPH[j]*dt);
    axE210TaperyPHxPH[j][ti] = (xfac*yTaperPH[j]+xTaperPH[i])*(bxE210TaperyPHxPH[j][ti]-1.)/(kxTaperPH[i]*(xfac*yTaperPH[j]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE210TaperyPHxPH[j][ti]>amax?axE210TaperyPHxPH[j][ti]:amax;
amin = axE210TaperyPHxPH[j][ti]<amin?axE210TaperyPHxPH[j][ti]:amin;
bmax = bxE210TaperyPHxPH[j][ti]>bmax?bxE210TaperyPHxPH[j][ti]:bmax;
bmin = bxE210TaperyPHxPH[j][ti]<bmin?bxE210TaperyPHxPH[j][ti]:bmin;
    bxE210TaperxPH[j][ti] = bxTaperPH[i]*exp(-xtap*yTaper[j]*dt);
    axE210TaperxPH[j][ti] = (xfac*yTaper[j]+xTaperPH[i])*(bxE210TaperxPH[j][ti]-1.)/(kxTaperPH[i]*(xfac*yTaper[j]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE210TaperxPH[j][ti]>amax?axE210TaperxPH[j][ti]:amax;
amin = axE210TaperxPH[j][ti]<amin?axE210TaperxPH[j][ti]:amin;
bmax = bxE210TaperxPH[j][ti]>bmax?bxE210TaperxPH[j][ti]:bmax;
bmin = bxE210TaperxPH[j][ti]<bmin?bxE210TaperxPH[j][ti]:bmin;
    byE210Taper[j][ti] = byTaper[j]*exp(-xtap*xTaper[i]*dt);
    ayE210Taper[j][ti] = (xfac*xTaper[i]+yTaper[j])*(byE210Taper[j][ti]-1.)/(kyTaper[j]*(xfac*xTaper[i]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE210Taper[j][ti]>amax?ayE210Taper[j][ti]:amax;
amin = ayE210Taper[j][ti]<amin?ayE210Taper[j][ti]:amin;
bmax = byE210Taper[j][ti]>bmax?byE210Taper[j][ti]:bmax;
bmin = byE210Taper[j][ti]<bmin?byE210Taper[j][ti]:bmin;
    byE210TaperyPH[j][ti] = byTaperPH[j]*exp(-xtap*xTaper[i]*dt);
    ayE210TaperyPH[j][ti] = (xfac*xTaper[i]+yTaperPH[j])*(byE210TaperyPH[j][ti]-1.)/(kyTaperPH[j]*(xfac*xTaper[i]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE210TaperyPH[j][ti]>amax?ayE210TaperyPH[j][ti]:amax;
amin = ayE210TaperyPH[j][ti]<amin?ayE210TaperyPH[j][ti]:amin;
bmax = byE210TaperyPH[j][ti]>bmax?byE210TaperyPH[j][ti]:bmax;
bmin = byE210TaperyPH[j][ti]<bmin?byE210TaperyPH[j][ti]:bmin;
    byE210TaperyPHxPH[j][ti] = byTaperPH[j]*exp(-xtap*xTaperPH[i]*dt);
    ayE210TaperyPHxPH[j][ti] = (xfac*xTaperPH[i]+yTaperPH[j])*(byE210TaperyPHxPH[j][ti]-1.)/(kyTaperPH[j]*(xfac*xTaperPH[i]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE210TaperyPHxPH[j][ti]>amax?ayE210TaperyPHxPH[j][ti]:amax;
amin = ayE210TaperyPHxPH[j][ti]<amin?ayE210TaperyPHxPH[j][ti]:amin;
bmax = byE210TaperyPHxPH[j][ti]>bmax?byE210TaperyPHxPH[j][ti]:bmax;
bmin = byE210TaperyPHxPH[j][ti]<bmin?byE210TaperyPHxPH[j][ti]:bmin;
    byE210TaperxPH[j][ti] = byTaper[j]*exp(-xtap*xTaperPH[i]*dt);
    ayE210TaperxPH[j][ti] = (xfac*xTaperPH[i]+yTaper[j])*(byE210TaperxPH[j][ti]-1.)/(kyTaper[j]*(xfac*xTaperPH[i]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE210TaperxPH[j][ti]>amax?ayE210TaperxPH[j][ti]:amax;
amin = ayE210TaperxPH[j][ti]<amin?ayE210TaperxPH[j][ti]:amin;
bmax = byE210TaperxPH[j][ti]>bmax?byE210TaperxPH[j][ti]:bmax;
bmin = byE210TaperxPH[j][ti]<bmin?byE210TaperxPH[j][ti]:bmin;
    bzE210Taper[j][ti] = exp(-xtap*yTaper[j]*dt)*exp(-xtap*xTaper[i]*dt);
    azE210Taper[j][ti] = bzE210Taper[j][ti]-1.;
amax = azE210Taper[j][ti]>amax?azE210Taper[j][ti]:amax;
amin = azE210Taper[j][ti]<amin?azE210Taper[j][ti]:amin;
bmax = bzE210Taper[j][ti]>bmax?bzE210Taper[j][ti]:bmax;
bmin = bzE210Taper[j][ti]<bmin?bzE210Taper[j][ti]:bmin;
    bzE210TaperyPH[j][ti] = exp(-xtap*yTaperPH[j]*dt)*exp(-xtap*xTaper[i]*dt);
    azE210TaperyPH[j][ti] = bzE210TaperyPH[j][ti]-1.;
amax = azE210TaperyPH[j][ti]>amax?azE210TaperyPH[j][ti]:amax;
amin = azE210TaperyPH[j][ti]<amin?azE210TaperyPH[j][ti]:amin;
bmax = bzE210TaperyPH[j][ti]>bmax?bzE210TaperyPH[j][ti]:bmax;
bmin = bzE210TaperyPH[j][ti]<bmin?bzE210TaperyPH[j][ti]:bmin;
    bzE210TaperyPHxPH[j][ti] = exp(-xtap*yTaperPH[j]*dt)*exp(-xtap*xTaperPH[i]*dt);
    azE210TaperyPHxPH[j][ti] = bzE210TaperyPHxPH[j][ti]-1.;
amax = azE210TaperyPHxPH[j][ti]>amax?azE210TaperyPHxPH[j][ti]:amax;
amin = azE210TaperyPHxPH[j][ti]<amin?azE210TaperyPHxPH[j][ti]:amin;
bmax = bzE210TaperyPHxPH[j][ti]>bmax?bzE210TaperyPHxPH[j][ti]:bmax;
bmin = bzE210TaperyPHxPH[j][ti]<bmin?bzE210TaperyPHxPH[j][ti]:bmin;
    bzE210TaperxPH[j][ti] = exp(-xtap*yTaper[j]*dt)*exp(-xtap*xTaperPH[i]*dt);
    azE210TaperxPH[j][ti] = bzE210TaperxPH[j][ti]-1.;
amax = azE210TaperxPH[j][ti]>amax?azE210TaperxPH[j][ti]:amax;
amin = azE210TaperxPH[j][ti]<amin?azE210TaperxPH[j][ti]:amin;
bmax = bzE210TaperxPH[j][ti]>bmax?bzE210TaperxPH[j][ti]:bmax;
bmin = bzE210TaperxPH[j][ti]<bmin?bzE210TaperxPH[j][ti]:bmin;
  }
}
for(int j=yEnd1;j<NY;++j) {
  int tj = j-yEnd1;
  for(int i=0;i<_nXmin;++i) {
    bxE120Taper[tj][i] = bxTaper[i]*exp(-xtap*yTaper[j]*dt);
    axE120Taper[tj][i] = (xfac*yTaper[j]+xTaper[i])*(bxE120Taper[tj][i]-1.)/(kxTaper[i]*(xfac*yTaper[j]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE120Taper[tj][i]>amax?axE120Taper[tj][i]:amax;
amin = axE120Taper[tj][i]<amin?axE120Taper[tj][i]:amin;
bmax = bxE120Taper[tj][i]>bmax?bxE120Taper[tj][i]:bmax;
bmin = bxE120Taper[tj][i]<bmin?bxE120Taper[tj][i]:bmin;
    bxE120TaperyPH[tj][i] = bxTaper[i]*exp(-xtap*yTaperPH[j]*dt);
    axE120TaperyPH[tj][i] = (xfac*yTaperPH[j]+xTaper[i])*(bxE120TaperyPH[tj][i]-1.)/(kxTaper[i]*(xfac*yTaperPH[j]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE120TaperyPH[tj][i]>amax?axE120TaperyPH[tj][i]:amax;
amin = axE120TaperyPH[tj][i]<amin?axE120TaperyPH[tj][i]:amin;
bmax = bxE120TaperyPH[tj][i]>bmax?bxE120TaperyPH[tj][i]:bmax;
bmin = bxE120TaperyPH[tj][i]<bmin?bxE120TaperyPH[tj][i]:bmin;
    bxE120TaperyPHxPH[tj][i] = bxTaperPH[i]*exp(-xtap*yTaperPH[j]*dt);
    axE120TaperyPHxPH[tj][i] = (xfac*yTaperPH[j]+xTaperPH[i])*(bxE120TaperyPHxPH[tj][i]-1.)/(kxTaperPH[i]*(xfac*yTaperPH[j]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE120TaperyPHxPH[tj][i]>amax?axE120TaperyPHxPH[tj][i]:amax;
amin = axE120TaperyPHxPH[tj][i]<amin?axE120TaperyPHxPH[tj][i]:amin;
bmax = bxE120TaperyPHxPH[tj][i]>bmax?bxE120TaperyPHxPH[tj][i]:bmax;
bmin = bxE120TaperyPHxPH[tj][i]<bmin?bxE120TaperyPHxPH[tj][i]:bmin;
    bxE120TaperxPH[tj][i] = bxTaperPH[i]*exp(-xtap*yTaper[j]*dt);
    axE120TaperxPH[tj][i] = (xfac*yTaper[j]+xTaperPH[i])*(bxE120TaperxPH[tj][i]-1.)/(kxTaperPH[i]*(xfac*yTaper[j]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE120TaperxPH[tj][i]>amax?axE120TaperxPH[tj][i]:amax;
amin = axE120TaperxPH[tj][i]<amin?axE120TaperxPH[tj][i]:amin;
bmax = bxE120TaperxPH[tj][i]>bmax?bxE120TaperxPH[tj][i]:bmax;
bmin = bxE120TaperxPH[tj][i]<bmin?bxE120TaperxPH[tj][i]:bmin;
    byE120Taper[tj][i] = byTaper[j]*exp(-xtap*xTaper[i]*dt);
    ayE120Taper[tj][i] = (xfac*xTaper[i]+yTaper[j])*(byE120Taper[tj][i]-1.)/(kyTaper[j]*(xfac*xTaper[i]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE120Taper[tj][i]>amax?ayE120Taper[tj][i]:amax;
amin = ayE120Taper[tj][i]<amin?ayE120Taper[tj][i]:amin;
bmax = byE120Taper[tj][i]>bmax?byE120Taper[tj][i]:bmax;
bmin = byE120Taper[tj][i]<bmin?byE120Taper[tj][i]:bmin;
    byE120TaperyPH[tj][i] = byTaperPH[j]*exp(-xtap*xTaper[i]*dt);
    ayE120TaperyPH[tj][i] = (xfac*xTaper[i]+yTaperPH[j])*(byE120TaperyPH[tj][i]-1.)/(kyTaperPH[j]*(xfac*xTaper[i]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE120TaperyPH[tj][i]>amax?ayE120TaperyPH[tj][i]:amax;
amin = ayE120TaperyPH[tj][i]<amin?ayE120TaperyPH[tj][i]:amin;
bmax = byE120TaperyPH[tj][i]>bmax?byE120TaperyPH[tj][i]:bmax;
bmin = byE120TaperyPH[tj][i]<bmin?byE120TaperyPH[tj][i]:bmin;
    byE120TaperyPHxPH[tj][i] = byTaperPH[j]*exp(-xtap*xTaperPH[i]*dt);
    ayE120TaperyPHxPH[tj][i] = (xfac*xTaperPH[i]+yTaperPH[j])*(byE120TaperyPHxPH[tj][i]-1.)/(kyTaperPH[j]*(xfac*xTaperPH[i]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE120TaperyPHxPH[tj][i]>amax?ayE120TaperyPHxPH[tj][i]:amax;
amin = ayE120TaperyPHxPH[tj][i]<amin?ayE120TaperyPHxPH[tj][i]:amin;
bmax = byE120TaperyPHxPH[tj][i]>bmax?byE120TaperyPHxPH[tj][i]:bmax;
bmin = byE120TaperyPHxPH[tj][i]<bmin?byE120TaperyPHxPH[tj][i]:bmin;
    byE120TaperxPH[tj][i] = byTaper[j]*exp(-xtap*xTaperPH[i]*dt);
    ayE120TaperxPH[tj][i] = (xfac*xTaperPH[i]+yTaper[j])*(byE120TaperxPH[tj][i]-1.)/(kyTaper[j]*(xfac*xTaperPH[i]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE120TaperxPH[tj][i]>amax?ayE120TaperxPH[tj][i]:amax;
amin = ayE120TaperxPH[tj][i]<amin?ayE120TaperxPH[tj][i]:amin;
bmax = byE120TaperxPH[tj][i]>bmax?byE120TaperxPH[tj][i]:bmax;
bmin = byE120TaperxPH[tj][i]<bmin?byE120TaperxPH[tj][i]:bmin;
    bzE120Taper[tj][i] = exp(-xtap*yTaper[j]*dt)*exp(-xtap*xTaper[i]*dt);
    azE120Taper[tj][i] = bzE120Taper[tj][i]-1.;
amax = azE120Taper[tj][i]>amax?azE120Taper[tj][i]:amax;
amin = azE120Taper[tj][i]<amin?azE120Taper[tj][i]:amin;
bmax = bzE120Taper[tj][i]>bmax?bzE120Taper[tj][i]:bmax;
bmin = bzE120Taper[tj][i]<bmin?bzE120Taper[tj][i]:bmin;
    bzE120TaperyPH[tj][i] = exp(-xtap*yTaperPH[j]*dt)*exp(-xtap*xTaper[i]*dt);
    azE120TaperyPH[tj][i] = bzE120TaperyPH[tj][i]-1.;
amax = azE120TaperyPH[tj][i]>amax?azE120TaperyPH[tj][i]:amax;
amin = azE120TaperyPH[tj][i]<amin?azE120TaperyPH[tj][i]:amin;
bmax = bzE120TaperyPH[tj][i]>bmax?bzE120TaperyPH[tj][i]:bmax;
bmin = bzE120TaperyPH[tj][i]<bmin?bzE120TaperyPH[tj][i]:bmin;
    bzE120TaperyPHxPH[tj][i] = exp(-xtap*yTaperPH[j]*dt)*exp(-xtap*xTaperPH[i]*dt);
    azE120TaperyPHxPH[tj][i] = bzE120TaperyPHxPH[tj][i]-1.;
amax = azE120TaperyPHxPH[tj][i]>amax?azE120TaperyPHxPH[tj][i]:amax;
amin = azE120TaperyPHxPH[tj][i]<amin?azE120TaperyPHxPH[tj][i]:amin;
bmax = bzE120TaperyPHxPH[tj][i]>bmax?bzE120TaperyPHxPH[tj][i]:bmax;
bmin = bzE120TaperyPHxPH[tj][i]<bmin?bzE120TaperyPHxPH[tj][i]:bmin;
    bzE120TaperxPH[tj][i] = exp(-xtap*yTaper[j]*dt)*exp(-xtap*xTaperPH[i]*dt);
    azE120TaperxPH[tj][i] = bzE120TaperxPH[tj][i]-1.;
amax = azE120TaperxPH[tj][i]>amax?azE120TaperxPH[tj][i]:amax;
amin = azE120TaperxPH[tj][i]<amin?azE120TaperxPH[tj][i]:amin;
bmax = bzE120TaperxPH[tj][i]>bmax?bzE120TaperxPH[tj][i]:bmax;
bmin = bzE120TaperxPH[tj][i]<bmin?bzE120TaperxPH[tj][i]:bmin;
  }
  for(int i=xEnd1;i<NX;++i) {
    int ti = i-xEnd1;
    bxE220Taper[tj][ti] = bxTaper[i]*exp(-xtap*yTaper[j]*dt);
    axE220Taper[tj][ti] = (xfac*yTaper[j]+xTaper[i])*(bxE220Taper[tj][ti]-1.)/(kxTaper[i]*(xfac*yTaper[j]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE220Taper[tj][ti]>amax?axE220Taper[tj][ti]:amax;
amin = axE220Taper[tj][ti]<amin?axE220Taper[tj][ti]:amin;
bmax = bxE220Taper[tj][ti]>bmax?bxE220Taper[tj][ti]:bmax;
bmin = bxE220Taper[tj][ti]<bmin?bxE220Taper[tj][ti]:bmin;
    bxE220TaperyPH[tj][ti] = bxTaper[i]*exp(-xtap*yTaperPH[j]*dt);
    axE220TaperyPH[tj][ti] = (xfac*yTaperPH[j]+xTaper[i])*(bxE220TaperyPH[tj][ti]-1.)/(kxTaper[i]*(xfac*yTaperPH[j]+xTaper[i]+kxTaper[i]*alphaXTaper[i]));
amax = axE220TaperyPH[tj][ti]>amax?axE220TaperyPH[tj][ti]:amax;
amin = axE220TaperyPH[tj][ti]<amin?axE220TaperyPH[tj][ti]:amin;
bmax = bxE220TaperyPH[tj][ti]>bmax?bxE220TaperyPH[tj][ti]:bmax;
bmin = bxE220TaperyPH[tj][ti]<bmin?bxE220TaperyPH[tj][ti]:bmin;
    bxE220TaperyPHxPH[tj][ti] = bxTaperPH[i]*exp(-xtap*yTaperPH[j]*dt);
    axE220TaperyPHxPH[tj][ti] = (xfac*yTaperPH[j]+xTaperPH[i])*(bxE220TaperyPHxPH[tj][ti]-1.)/(kxTaperPH[i]*(xfac*yTaperPH[j]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE220TaperyPHxPH[tj][ti]>amax?axE220TaperyPHxPH[tj][ti]:amax;
amin = axE220TaperyPHxPH[tj][ti]<amin?axE220TaperyPHxPH[tj][ti]:amin;
bmax = bxE220TaperyPHxPH[tj][ti]>bmax?bxE220TaperyPHxPH[tj][ti]:bmax;
bmin = bxE220TaperyPHxPH[tj][ti]<bmin?bxE220TaperyPHxPH[tj][ti]:bmin;
    bxE220TaperxPH[tj][ti] = bxTaperPH[i]*exp(-xtap*yTaper[j]*dt);
    axE220TaperxPH[tj][ti] = (xfac*yTaper[j]+xTaperPH[i])*(bxE220TaperxPH[tj][ti]-1.)/(kxTaperPH[i]*(xfac*yTaper[j]+xTaperPH[i]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axE220TaperxPH[tj][ti]>amax?axE220TaperxPH[tj][ti]:amax;
amin = axE220TaperxPH[tj][ti]<amin?axE220TaperxPH[tj][ti]:amin;
bmax = bxE220TaperxPH[tj][ti]>bmax?bxE220TaperxPH[tj][ti]:bmax;
bmin = bxE220TaperxPH[tj][ti]<bmin?bxE220TaperxPH[tj][ti]:bmin;
    byE220Taper[tj][ti] = byTaper[j]*exp(-xtap*xTaper[i]*dt);
    ayE220Taper[tj][ti] = (xfac*xTaper[i]+yTaper[j])*(byE220Taper[tj][ti]-1.)/(kyTaper[j]*(xfac*xTaper[i]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE220Taper[tj][ti]>amax?ayE220Taper[tj][ti]:amax;
amin = ayE220Taper[tj][ti]<amin?ayE220Taper[tj][ti]:amin;
bmax = byE220Taper[tj][ti]>bmax?byE220Taper[tj][ti]:bmax;
bmin = byE220Taper[tj][ti]<bmin?byE220Taper[tj][ti]:bmin;
    byE220TaperyPH[tj][ti] = byTaperPH[j]*exp(-xtap*xTaper[i]*dt);
    ayE220TaperyPH[tj][ti] = (xfac*xTaper[i]+yTaperPH[j])*(byE220TaperyPH[tj][ti]-1.)/(kyTaperPH[j]*(xfac*xTaper[i]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE220TaperyPH[tj][ti]>amax?ayE220TaperyPH[tj][ti]:amax;
amin = ayE220TaperyPH[tj][ti]<amin?ayE220TaperyPH[tj][ti]:amin;
bmax = byE220TaperyPH[tj][ti]>bmax?byE220TaperyPH[tj][ti]:bmax;
bmin = byE220TaperyPH[tj][ti]<bmin?byE220TaperyPH[tj][ti]:bmin;
    byE220TaperyPHxPH[tj][ti] = byTaperPH[j]*exp(-xtap*xTaperPH[i]*dt);
    ayE220TaperyPHxPH[tj][ti] = (xfac*xTaperPH[i]+yTaperPH[j])*(byE220TaperyPHxPH[tj][ti]-1.)/(kyTaperPH[j]*(xfac*xTaperPH[i]+yTaperPH[j]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayE220TaperyPHxPH[tj][ti]>amax?ayE220TaperyPHxPH[tj][ti]:amax;
amin = ayE220TaperyPHxPH[tj][ti]<amin?ayE220TaperyPHxPH[tj][ti]:amin;
bmax = byE220TaperyPHxPH[tj][ti]>bmax?byE220TaperyPHxPH[tj][ti]:bmax;
bmin = byE220TaperyPHxPH[tj][ti]<bmin?byE220TaperyPHxPH[tj][ti]:bmin;
    byE220TaperxPH[tj][ti] = byTaper[j]*exp(-xtap*xTaperPH[i]*dt);
    ayE220TaperxPH[tj][ti] = (xfac*xTaperPH[i]+yTaper[j])*(byE220TaperxPH[tj][ti]-1.)/(kyTaper[j]*(xfac*xTaperPH[i]+yTaper[j]+kyTaper[j]*alphaYTaper[j]));
amax = ayE220TaperxPH[tj][ti]>amax?ayE220TaperxPH[tj][ti]:amax;
amin = ayE220TaperxPH[tj][ti]<amin?ayE220TaperxPH[tj][ti]:amin;
bmax = byE220TaperxPH[tj][ti]>bmax?byE220TaperxPH[tj][ti]:bmax;
bmin = byE220TaperxPH[tj][ti]<bmin?byE220TaperxPH[tj][ti]:bmin;
    bzE220Taper[tj][ti] = exp(-xtap*yTaper[j]*dt)*exp(-xtap*xTaper[i]*dt);
    azE220Taper[tj][ti] = bzE220Taper[tj][ti]-1.;
amax = azE220Taper[tj][ti]>amax?azE220Taper[tj][ti]:amax;
amin = azE220Taper[tj][ti]<amin?azE220Taper[tj][ti]:amin;
bmax = bzE220Taper[tj][ti]>bmax?bzE220Taper[tj][ti]:bmax;
bmin = bzE220Taper[tj][ti]<bmin?bzE220Taper[tj][ti]:bmin;
    bzE220TaperyPH[tj][ti] = exp(-xtap*yTaperPH[j]*dt)*exp(-xtap*xTaper[i]*dt);
    azE220TaperyPH[tj][ti] = bzE220TaperyPH[tj][ti]-1.;
amax = azE220TaperyPH[tj][ti]>amax?azE220TaperyPH[tj][ti]:amax;
amin = azE220TaperyPH[tj][ti]<amin?azE220TaperyPH[tj][ti]:amin;
bmax = bzE220TaperyPH[tj][ti]>bmax?bzE220TaperyPH[tj][ti]:bmax;
bmin = bzE220TaperyPH[tj][ti]<bmin?bzE220TaperyPH[tj][ti]:bmin;
    bzE220TaperyPHxPH[tj][ti] = exp(-xtap*yTaperPH[j]*dt)*exp(-xtap*xTaperPH[i]*dt);
    azE220TaperyPHxPH[tj][ti] = bzE220TaperyPHxPH[tj][ti]-1.;
amax = azE220TaperyPHxPH[tj][ti]>amax?azE220TaperyPHxPH[tj][ti]:amax;
amin = azE220TaperyPHxPH[tj][ti]<amin?azE220TaperyPHxPH[tj][ti]:amin;
bmax = bzE220TaperyPHxPH[tj][ti]>bmax?bzE220TaperyPHxPH[tj][ti]:bmax;
bmin = bzE220TaperyPHxPH[tj][ti]<bmin?bzE220TaperyPHxPH[tj][ti]:bmin;
    bzE220TaperxPH[tj][ti] = exp(-xtap*yTaper[j]*dt)*exp(-xtap*xTaperPH[i]*dt);
    azE220TaperxPH[tj][ti] = bzE220TaperxPH[tj][ti]-1.;
amax = azE220TaperxPH[tj][ti]>amax?azE220TaperxPH[tj][ti]:amax;
amin = azE220TaperxPH[tj][ti]<amin?azE220TaperxPH[tj][ti]:amin;
bmax = bzE220TaperxPH[tj][ti]>bmax?bzE220TaperxPH[tj][ti]:bmax;
bmin = bzE220TaperxPH[tj][ti]<bmin?bzE220TaperxPH[tj][ti]:bmin;
  }
}
for(int k=0;k<_nZmin;++k) {
  for(int j=0;j<_nYmin;++j) {
    for(int i=0;i<_nXmin;++i) {
      bxC000Taper[k][j][i] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC000Taper[k][j][i] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC000Taper[k][j][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC000Taper[k][j][i]>amax?axC000Taper[k][j][i]:amax;
amin = axC000Taper[k][j][i]<amin?axC000Taper[k][j][i]:amin;
bmax = bxC000Taper[k][j][i]>bmax?bxC000Taper[k][j][i]:bmax;
bmin = bxC000Taper[k][j][i]<bmin?bxC000Taper[k][j][i]:bmin;
      bxC000TaperzPH[k][j][i] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC000TaperzPH[k][j][i] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC000TaperzPH[k][j][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC000TaperzPH[k][j][i]>amax?axC000TaperzPH[k][j][i]:amax;
amin = axC000TaperzPH[k][j][i]<amin?axC000TaperzPH[k][j][i]:amin;
bmax = bxC000TaperzPH[k][j][i]>bmax?bxC000TaperzPH[k][j][i]:bmax;
bmin = bxC000TaperzPH[k][j][i]<bmin?bxC000TaperzPH[k][j][i]:bmin;
      bxC000TaperzPHxPH[k][j][i] = bxTaperPH[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC000TaperzPHxPH[k][j][i] = (xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC000TaperzPHxPH[k][j][i]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC000TaperzPHxPH[k][j][i]>amax?axC000TaperzPHxPH[k][j][i]:amax;
amin = axC000TaperzPHxPH[k][j][i]<amin?axC000TaperzPHxPH[k][j][i]:amin;
bmax = bxC000TaperzPHxPH[k][j][i]>bmax?bxC000TaperzPHxPH[k][j][i]:bmax;
bmin = bxC000TaperzPHxPH[k][j][i]<bmin?bxC000TaperzPHxPH[k][j][i]:bmin;
      bxC000TaperzPHyPH[k][j][i] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC000TaperzPHyPH[k][j][i] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j])*(bxC000TaperzPHyPH[k][j][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC000TaperzPHyPH[k][j][i]>amax?axC000TaperzPHyPH[k][j][i]:amax;
amin = axC000TaperzPHyPH[k][j][i]<amin?axC000TaperzPHyPH[k][j][i]:amin;
bmax = bxC000TaperzPHyPH[k][j][i]>bmax?bxC000TaperzPHyPH[k][j][i]:bmax;
bmin = bxC000TaperzPHyPH[k][j][i]<bmin?bxC000TaperzPHyPH[k][j][i]:bmin;
      bxC000TaperyPH[k][j][i] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC000TaperyPH[k][j][i] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC000TaperyPH[k][j][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC000TaperyPH[k][j][i]>amax?axC000TaperyPH[k][j][i]:amax;
amin = axC000TaperyPH[k][j][i]<amin?axC000TaperyPH[k][j][i]:amin;
bmax = bxC000TaperyPH[k][j][i]>bmax?bxC000TaperyPH[k][j][i]:bmax;
bmin = bxC000TaperyPH[k][j][i]<bmin?bxC000TaperyPH[k][j][i]:bmin;
      bxC000TaperyPHxPH[k][j][i] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC000TaperyPHxPH[k][j][i] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC000TaperyPHxPH[k][j][i]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC000TaperyPHxPH[k][j][i]>amax?axC000TaperyPHxPH[k][j][i]:amax;
amin = axC000TaperyPHxPH[k][j][i]<amin?axC000TaperyPHxPH[k][j][i]:amin;
bmax = bxC000TaperyPHxPH[k][j][i]>bmax?bxC000TaperyPHxPH[k][j][i]:bmax;
bmin = bxC000TaperyPHxPH[k][j][i]<bmin?bxC000TaperyPHxPH[k][j][i]:bmin;
      bxC000TaperxPH[k][j][i] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC000TaperxPH[k][j][i] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC000TaperxPH[k][j][i]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC000TaperxPH[k][j][i]>amax?axC000TaperxPH[k][j][i]:amax;
amin = axC000TaperxPH[k][j][i]<amin?axC000TaperxPH[k][j][i]:amin;
bmax = bxC000TaperxPH[k][j][i]>bmax?bxC000TaperxPH[k][j][i]:bmax;
bmin = bxC000TaperxPH[k][j][i]<bmin?bxC000TaperxPH[k][j][i]:bmin;
      byC000Taper[k][j][i] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC000Taper[k][j][i] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC000Taper[k][j][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC000Taper[k][j][i]>amax?ayC000Taper[k][j][i]:amax;
amin = ayC000Taper[k][j][i]<amin?ayC000Taper[k][j][i]:amin;
bmax = byC000Taper[k][j][i]>bmax?byC000Taper[k][j][i]:bmax;
bmin = byC000Taper[k][j][i]<bmin?byC000Taper[k][j][i]:bmin;
      byC000TaperzPH[k][j][i] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC000TaperzPH[k][j][i] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC000TaperzPH[k][j][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC000TaperzPH[k][j][i]>amax?ayC000TaperzPH[k][j][i]:amax;
amin = ayC000TaperzPH[k][j][i]<amin?ayC000TaperzPH[k][j][i]:amin;
bmax = byC000TaperzPH[k][j][i]>bmax?byC000TaperzPH[k][j][i]:bmax;
bmin = byC000TaperzPH[k][j][i]<bmin?byC000TaperzPH[k][j][i]:bmin;
      byC000TaperzPHxPH[k][j][i] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC000TaperzPHxPH[k][j][i] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i])*(byC000TaperzPHxPH[k][j][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC000TaperzPHxPH[k][j][i]>amax?ayC000TaperzPHxPH[k][j][i]:amax;
amin = ayC000TaperzPHxPH[k][j][i]<amin?ayC000TaperzPHxPH[k][j][i]:amin;
bmax = byC000TaperzPHxPH[k][j][i]>bmax?byC000TaperzPHxPH[k][j][i]:bmax;
bmin = byC000TaperzPHxPH[k][j][i]<bmin?byC000TaperzPHxPH[k][j][i]:bmin;
      byC000TaperzPHyPH[k][j][i] = byTaperPH[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC000TaperzPHyPH[k][j][i] = (yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC000TaperzPHyPH[k][j][i]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC000TaperzPHyPH[k][j][i]>amax?ayC000TaperzPHyPH[k][j][i]:amax;
amin = ayC000TaperzPHyPH[k][j][i]<amin?ayC000TaperzPHyPH[k][j][i]:amin;
bmax = byC000TaperzPHyPH[k][j][i]>bmax?byC000TaperzPHyPH[k][j][i]:bmax;
bmin = byC000TaperzPHyPH[k][j][i]<bmin?byC000TaperzPHyPH[k][j][i]:bmin;
      byC000TaperyPH[k][j][i] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC000TaperyPH[k][j][i] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC000TaperyPH[k][j][i]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC000TaperyPH[k][j][i]>amax?ayC000TaperyPH[k][j][i]:amax;
amin = ayC000TaperyPH[k][j][i]<amin?ayC000TaperyPH[k][j][i]:amin;
bmax = byC000TaperyPH[k][j][i]>bmax?byC000TaperyPH[k][j][i]:bmax;
bmin = byC000TaperyPH[k][j][i]<bmin?byC000TaperyPH[k][j][i]:bmin;
      byC000TaperyPHxPH[k][j][i] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC000TaperyPHxPH[k][j][i] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC000TaperyPHxPH[k][j][i]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC000TaperyPHxPH[k][j][i]>amax?ayC000TaperyPHxPH[k][j][i]:amax;
amin = ayC000TaperyPHxPH[k][j][i]<amin?ayC000TaperyPHxPH[k][j][i]:amin;
bmax = byC000TaperyPHxPH[k][j][i]>bmax?byC000TaperyPHxPH[k][j][i]:bmax;
bmin = byC000TaperyPHxPH[k][j][i]<bmin?byC000TaperyPHxPH[k][j][i]:bmin;
      byC000TaperxPH[k][j][i] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC000TaperxPH[k][j][i] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC000TaperxPH[k][j][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC000TaperxPH[k][j][i]>amax?ayC000TaperxPH[k][j][i]:amax;
amin = ayC000TaperxPH[k][j][i]<amin?ayC000TaperxPH[k][j][i]:amin;
bmax = byC000TaperxPH[k][j][i]>bmax?byC000TaperxPH[k][j][i]:bmax;
bmin = byC000TaperxPH[k][j][i]<bmin?byC000TaperxPH[k][j][i]:bmin;
      bzC000Taper[k][j][i] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC000Taper[k][j][i] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC000Taper[k][j][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC000Taper[k][j][i]>amax?azC000Taper[k][j][i]:amax;
amin = azC000Taper[k][j][i]<amin?azC000Taper[k][j][i]:amin;
bmax = bzC000Taper[k][j][i]>bmax?bzC000Taper[k][j][i]:bmax;
bmin = bzC000Taper[k][j][i]<bmin?bzC000Taper[k][j][i]:bmin;
      bzC000TaperzPH[k][j][i] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC000TaperzPH[k][j][i] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC000TaperzPH[k][j][i]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC000TaperzPH[k][j][i]>amax?azC000TaperzPH[k][j][i]:amax;
amin = azC000TaperzPH[k][j][i]<amin?azC000TaperzPH[k][j][i]:amin;
bmax = bzC000TaperzPH[k][j][i]>bmax?bzC000TaperzPH[k][j][i]:bmax;
bmin = bzC000TaperzPH[k][j][i]<bmin?bzC000TaperzPH[k][j][i]:bmin;
      bzC000TaperzPHxPH[k][j][i] = bzTaperPH[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC000TaperzPHxPH[k][j][i] = (zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC000TaperzPHxPH[k][j][i]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC000TaperzPHxPH[k][j][i]>amax?azC000TaperzPHxPH[k][j][i]:amax;
amin = azC000TaperzPHxPH[k][j][i]<amin?azC000TaperzPHxPH[k][j][i]:amin;
bmax = bzC000TaperzPHxPH[k][j][i]>bmax?bzC000TaperzPHxPH[k][j][i]:bmax;
bmin = bzC000TaperzPHxPH[k][j][i]<bmin?bzC000TaperzPHxPH[k][j][i]:bmin;
      bzC000TaperzPHyPH[k][j][i] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC000TaperzPHyPH[k][j][i] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC000TaperzPHyPH[k][j][i]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC000TaperzPHyPH[k][j][i]>amax?azC000TaperzPHyPH[k][j][i]:amax;
amin = azC000TaperzPHyPH[k][j][i]<amin?azC000TaperzPHyPH[k][j][i]:amin;
bmax = bzC000TaperzPHyPH[k][j][i]>bmax?bzC000TaperzPHyPH[k][j][i]:bmax;
bmin = bzC000TaperzPHyPH[k][j][i]<bmin?bzC000TaperzPHyPH[k][j][i]:bmin;
      bzC000TaperyPH[k][j][i] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC000TaperyPH[k][j][i] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC000TaperyPH[k][j][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC000TaperyPH[k][j][i]>amax?azC000TaperyPH[k][j][i]:amax;
amin = azC000TaperyPH[k][j][i]<amin?azC000TaperyPH[k][j][i]:amin;
bmax = bzC000TaperyPH[k][j][i]>bmax?bzC000TaperyPH[k][j][i]:bmax;
bmin = bzC000TaperyPH[k][j][i]<bmin?bzC000TaperyPH[k][j][i]:bmin;
      bzC000TaperyPHxPH[k][j][i] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC000TaperyPHxPH[k][j][i] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j])*(bzC000TaperyPHxPH[k][j][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC000TaperyPHxPH[k][j][i]>amax?azC000TaperyPHxPH[k][j][i]:amax;
amin = azC000TaperyPHxPH[k][j][i]<amin?azC000TaperyPHxPH[k][j][i]:amin;
bmax = bzC000TaperyPHxPH[k][j][i]>bmax?bzC000TaperyPHxPH[k][j][i]:bmax;
bmin = bzC000TaperyPHxPH[k][j][i]<bmin?bzC000TaperyPHxPH[k][j][i]:bmin;
      bzC000TaperxPH[k][j][i] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC000TaperxPH[k][j][i] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC000TaperxPH[k][j][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC000TaperxPH[k][j][i]>amax?azC000TaperxPH[k][j][i]:amax;
amin = azC000TaperxPH[k][j][i]<amin?azC000TaperxPH[k][j][i]:amin;
bmax = bzC000TaperxPH[k][j][i]>bmax?bzC000TaperxPH[k][j][i]:bmax;
bmin = bzC000TaperxPH[k][j][i]<bmin?bzC000TaperxPH[k][j][i]:bmin;
    }
    for(int i=xEnd1;i<NX;++i) {
      int ti = i-xEnd1;
      bxC100Taper[k][j][ti] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC100Taper[k][j][ti] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC100Taper[k][j][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC100Taper[k][j][ti]>amax?axC100Taper[k][j][ti]:amax;
amin = axC100Taper[k][j][ti]<amin?axC100Taper[k][j][ti]:amin;
bmax = bxC100Taper[k][j][ti]>bmax?bxC100Taper[k][j][ti]:bmax;
bmin = bxC100Taper[k][j][ti]<bmin?bxC100Taper[k][j][ti]:bmin;
      bxC100TaperzPH[k][j][ti] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC100TaperzPH[k][j][ti] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC100TaperzPH[k][j][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC100TaperzPH[k][j][ti]>amax?axC100TaperzPH[k][j][ti]:amax;
amin = axC100TaperzPH[k][j][ti]<amin?axC100TaperzPH[k][j][ti]:amin;
bmax = bxC100TaperzPH[k][j][ti]>bmax?bxC100TaperzPH[k][j][ti]:bmax;
bmin = bxC100TaperzPH[k][j][ti]<bmin?bxC100TaperzPH[k][j][ti]:bmin;
      bxC100TaperzPHxPH[k][j][ti] = bxTaperPH[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC100TaperzPHxPH[k][j][ti] = (xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC100TaperzPHxPH[k][j][ti]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC100TaperzPHxPH[k][j][ti]>amax?axC100TaperzPHxPH[k][j][ti]:amax;
amin = axC100TaperzPHxPH[k][j][ti]<amin?axC100TaperzPHxPH[k][j][ti]:amin;
bmax = bxC100TaperzPHxPH[k][j][ti]>bmax?bxC100TaperzPHxPH[k][j][ti]:bmax;
bmin = bxC100TaperzPHxPH[k][j][ti]<bmin?bxC100TaperzPHxPH[k][j][ti]:bmin;
      bxC100TaperzPHyPH[k][j][ti] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC100TaperzPHyPH[k][j][ti] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j])*(bxC100TaperzPHyPH[k][j][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC100TaperzPHyPH[k][j][ti]>amax?axC100TaperzPHyPH[k][j][ti]:amax;
amin = axC100TaperzPHyPH[k][j][ti]<amin?axC100TaperzPHyPH[k][j][ti]:amin;
bmax = bxC100TaperzPHyPH[k][j][ti]>bmax?bxC100TaperzPHyPH[k][j][ti]:bmax;
bmin = bxC100TaperzPHyPH[k][j][ti]<bmin?bxC100TaperzPHyPH[k][j][ti]:bmin;
      bxC100TaperyPH[k][j][ti] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC100TaperyPH[k][j][ti] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC100TaperyPH[k][j][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC100TaperyPH[k][j][ti]>amax?axC100TaperyPH[k][j][ti]:amax;
amin = axC100TaperyPH[k][j][ti]<amin?axC100TaperyPH[k][j][ti]:amin;
bmax = bxC100TaperyPH[k][j][ti]>bmax?bxC100TaperyPH[k][j][ti]:bmax;
bmin = bxC100TaperyPH[k][j][ti]<bmin?bxC100TaperyPH[k][j][ti]:bmin;
      bxC100TaperyPHxPH[k][j][ti] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC100TaperyPHxPH[k][j][ti] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC100TaperyPHxPH[k][j][ti]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC100TaperyPHxPH[k][j][ti]>amax?axC100TaperyPHxPH[k][j][ti]:amax;
amin = axC100TaperyPHxPH[k][j][ti]<amin?axC100TaperyPHxPH[k][j][ti]:amin;
bmax = bxC100TaperyPHxPH[k][j][ti]>bmax?bxC100TaperyPHxPH[k][j][ti]:bmax;
bmin = bxC100TaperyPHxPH[k][j][ti]<bmin?bxC100TaperyPHxPH[k][j][ti]:bmin;
      bxC100TaperxPH[k][j][ti] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC100TaperxPH[k][j][ti] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC100TaperxPH[k][j][ti]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC100TaperxPH[k][j][ti]>amax?axC100TaperxPH[k][j][ti]:amax;
amin = axC100TaperxPH[k][j][ti]<amin?axC100TaperxPH[k][j][ti]:amin;
bmax = bxC100TaperxPH[k][j][ti]>bmax?bxC100TaperxPH[k][j][ti]:bmax;
bmin = bxC100TaperxPH[k][j][ti]<bmin?bxC100TaperxPH[k][j][ti]:bmin;
      byC100Taper[k][j][ti] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC100Taper[k][j][ti] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC100Taper[k][j][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC100Taper[k][j][ti]>amax?ayC100Taper[k][j][ti]:amax;
amin = ayC100Taper[k][j][ti]<amin?ayC100Taper[k][j][ti]:amin;
bmax = byC100Taper[k][j][ti]>bmax?byC100Taper[k][j][ti]:bmax;
bmin = byC100Taper[k][j][ti]<bmin?byC100Taper[k][j][ti]:bmin;
      byC100TaperzPH[k][j][ti] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC100TaperzPH[k][j][ti] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC100TaperzPH[k][j][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC100TaperzPH[k][j][ti]>amax?ayC100TaperzPH[k][j][ti]:amax;
amin = ayC100TaperzPH[k][j][ti]<amin?ayC100TaperzPH[k][j][ti]:amin;
bmax = byC100TaperzPH[k][j][ti]>bmax?byC100TaperzPH[k][j][ti]:bmax;
bmin = byC100TaperzPH[k][j][ti]<bmin?byC100TaperzPH[k][j][ti]:bmin;
      byC100TaperzPHxPH[k][j][ti] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC100TaperzPHxPH[k][j][ti] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i])*(byC100TaperzPHxPH[k][j][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC100TaperzPHxPH[k][j][ti]>amax?ayC100TaperzPHxPH[k][j][ti]:amax;
amin = ayC100TaperzPHxPH[k][j][ti]<amin?ayC100TaperzPHxPH[k][j][ti]:amin;
bmax = byC100TaperzPHxPH[k][j][ti]>bmax?byC100TaperzPHxPH[k][j][ti]:bmax;
bmin = byC100TaperzPHxPH[k][j][ti]<bmin?byC100TaperzPHxPH[k][j][ti]:bmin;
      byC100TaperzPHyPH[k][j][ti] = byTaperPH[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC100TaperzPHyPH[k][j][ti] = (yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC100TaperzPHyPH[k][j][ti]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC100TaperzPHyPH[k][j][ti]>amax?ayC100TaperzPHyPH[k][j][ti]:amax;
amin = ayC100TaperzPHyPH[k][j][ti]<amin?ayC100TaperzPHyPH[k][j][ti]:amin;
bmax = byC100TaperzPHyPH[k][j][ti]>bmax?byC100TaperzPHyPH[k][j][ti]:bmax;
bmin = byC100TaperzPHyPH[k][j][ti]<bmin?byC100TaperzPHyPH[k][j][ti]:bmin;
      byC100TaperyPH[k][j][ti] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC100TaperyPH[k][j][ti] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC100TaperyPH[k][j][ti]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC100TaperyPH[k][j][ti]>amax?ayC100TaperyPH[k][j][ti]:amax;
amin = ayC100TaperyPH[k][j][ti]<amin?ayC100TaperyPH[k][j][ti]:amin;
bmax = byC100TaperyPH[k][j][ti]>bmax?byC100TaperyPH[k][j][ti]:bmax;
bmin = byC100TaperyPH[k][j][ti]<bmin?byC100TaperyPH[k][j][ti]:bmin;
      byC100TaperyPHxPH[k][j][ti] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC100TaperyPHxPH[k][j][ti] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC100TaperyPHxPH[k][j][ti]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC100TaperyPHxPH[k][j][ti]>amax?ayC100TaperyPHxPH[k][j][ti]:amax;
amin = ayC100TaperyPHxPH[k][j][ti]<amin?ayC100TaperyPHxPH[k][j][ti]:amin;
bmax = byC100TaperyPHxPH[k][j][ti]>bmax?byC100TaperyPHxPH[k][j][ti]:bmax;
bmin = byC100TaperyPHxPH[k][j][ti]<bmin?byC100TaperyPHxPH[k][j][ti]:bmin;
      byC100TaperxPH[k][j][ti] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC100TaperxPH[k][j][ti] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC100TaperxPH[k][j][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC100TaperxPH[k][j][ti]>amax?ayC100TaperxPH[k][j][ti]:amax;
amin = ayC100TaperxPH[k][j][ti]<amin?ayC100TaperxPH[k][j][ti]:amin;
bmax = byC100TaperxPH[k][j][ti]>bmax?byC100TaperxPH[k][j][ti]:bmax;
bmin = byC100TaperxPH[k][j][ti]<bmin?byC100TaperxPH[k][j][ti]:bmin;
      bzC100Taper[k][j][ti] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC100Taper[k][j][ti] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC100Taper[k][j][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC100Taper[k][j][ti]>amax?azC100Taper[k][j][ti]:amax;
amin = azC100Taper[k][j][ti]<amin?azC100Taper[k][j][ti]:amin;
bmax = bzC100Taper[k][j][ti]>bmax?bzC100Taper[k][j][ti]:bmax;
bmin = bzC100Taper[k][j][ti]<bmin?bzC100Taper[k][j][ti]:bmin;
      bzC100TaperzPH[k][j][ti] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC100TaperzPH[k][j][ti] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC100TaperzPH[k][j][ti]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC100TaperzPH[k][j][ti]>amax?azC100TaperzPH[k][j][ti]:amax;
amin = azC100TaperzPH[k][j][ti]<amin?azC100TaperzPH[k][j][ti]:amin;
bmax = bzC100TaperzPH[k][j][ti]>bmax?bzC100TaperzPH[k][j][ti]:bmax;
bmin = bzC100TaperzPH[k][j][ti]<bmin?bzC100TaperzPH[k][j][ti]:bmin;
      bzC100TaperzPHxPH[k][j][ti] = bzTaperPH[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC100TaperzPHxPH[k][j][ti] = (zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC100TaperzPHxPH[k][j][ti]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC100TaperzPHxPH[k][j][ti]>amax?azC100TaperzPHxPH[k][j][ti]:amax;
amin = azC100TaperzPHxPH[k][j][ti]<amin?azC100TaperzPHxPH[k][j][ti]:amin;
bmax = bzC100TaperzPHxPH[k][j][ti]>bmax?bzC100TaperzPHxPH[k][j][ti]:bmax;
bmin = bzC100TaperzPHxPH[k][j][ti]<bmin?bzC100TaperzPHxPH[k][j][ti]:bmin;
      bzC100TaperzPHyPH[k][j][ti] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC100TaperzPHyPH[k][j][ti] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC100TaperzPHyPH[k][j][ti]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC100TaperzPHyPH[k][j][ti]>amax?azC100TaperzPHyPH[k][j][ti]:amax;
amin = azC100TaperzPHyPH[k][j][ti]<amin?azC100TaperzPHyPH[k][j][ti]:amin;
bmax = bzC100TaperzPHyPH[k][j][ti]>bmax?bzC100TaperzPHyPH[k][j][ti]:bmax;
bmin = bzC100TaperzPHyPH[k][j][ti]<bmin?bzC100TaperzPHyPH[k][j][ti]:bmin;
      bzC100TaperyPH[k][j][ti] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC100TaperyPH[k][j][ti] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC100TaperyPH[k][j][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC100TaperyPH[k][j][ti]>amax?azC100TaperyPH[k][j][ti]:amax;
amin = azC100TaperyPH[k][j][ti]<amin?azC100TaperyPH[k][j][ti]:amin;
bmax = bzC100TaperyPH[k][j][ti]>bmax?bzC100TaperyPH[k][j][ti]:bmax;
bmin = bzC100TaperyPH[k][j][ti]<bmin?bzC100TaperyPH[k][j][ti]:bmin;
      bzC100TaperyPHxPH[k][j][ti] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC100TaperyPHxPH[k][j][ti] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j])*(bzC100TaperyPHxPH[k][j][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC100TaperyPHxPH[k][j][ti]>amax?azC100TaperyPHxPH[k][j][ti]:amax;
amin = azC100TaperyPHxPH[k][j][ti]<amin?azC100TaperyPHxPH[k][j][ti]:amin;
bmax = bzC100TaperyPHxPH[k][j][ti]>bmax?bzC100TaperyPHxPH[k][j][ti]:bmax;
bmin = bzC100TaperyPHxPH[k][j][ti]<bmin?bzC100TaperyPHxPH[k][j][ti]:bmin;
      bzC100TaperxPH[k][j][ti] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC100TaperxPH[k][j][ti] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC100TaperxPH[k][j][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC100TaperxPH[k][j][ti]>amax?azC100TaperxPH[k][j][ti]:amax;
amin = azC100TaperxPH[k][j][ti]<amin?azC100TaperxPH[k][j][ti]:amin;
bmax = bzC100TaperxPH[k][j][ti]>bmax?bzC100TaperxPH[k][j][ti]:bmax;
bmin = bzC100TaperxPH[k][j][ti]<bmin?bzC100TaperxPH[k][j][ti]:bmin;
    }
  }
  for(int j=yEnd1;j<NY;++j) {
    int tj = j-yEnd1;
    for(int i=0;i<_nXmin;++i) {
      bxC010Taper[k][tj][i] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC010Taper[k][tj][i] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC010Taper[k][tj][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC010Taper[k][tj][i]>amax?axC010Taper[k][tj][i]:amax;
amin = axC010Taper[k][tj][i]<amin?axC010Taper[k][tj][i]:amin;
bmax = bxC010Taper[k][tj][i]>bmax?bxC010Taper[k][tj][i]:bmax;
bmin = bxC010Taper[k][tj][i]<bmin?bxC010Taper[k][tj][i]:bmin;
      bxC010TaperzPH[k][tj][i] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC010TaperzPH[k][tj][i] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC010TaperzPH[k][tj][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC010TaperzPH[k][tj][i]>amax?axC010TaperzPH[k][tj][i]:amax;
amin = axC010TaperzPH[k][tj][i]<amin?axC010TaperzPH[k][tj][i]:amin;
bmax = bxC010TaperzPH[k][tj][i]>bmax?bxC010TaperzPH[k][tj][i]:bmax;
bmin = bxC010TaperzPH[k][tj][i]<bmin?bxC010TaperzPH[k][tj][i]:bmin;
      bxC010TaperzPHxPH[k][tj][i] = bxTaperPH[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC010TaperzPHxPH[k][tj][i] = (xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC010TaperzPHxPH[k][tj][i]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC010TaperzPHxPH[k][tj][i]>amax?axC010TaperzPHxPH[k][tj][i]:amax;
amin = axC010TaperzPHxPH[k][tj][i]<amin?axC010TaperzPHxPH[k][tj][i]:amin;
bmax = bxC010TaperzPHxPH[k][tj][i]>bmax?bxC010TaperzPHxPH[k][tj][i]:bmax;
bmin = bxC010TaperzPHxPH[k][tj][i]<bmin?bxC010TaperzPHxPH[k][tj][i]:bmin;
      bxC010TaperzPHyPH[k][tj][i] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC010TaperzPHyPH[k][tj][i] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j])*(bxC010TaperzPHyPH[k][tj][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC010TaperzPHyPH[k][tj][i]>amax?axC010TaperzPHyPH[k][tj][i]:amax;
amin = axC010TaperzPHyPH[k][tj][i]<amin?axC010TaperzPHyPH[k][tj][i]:amin;
bmax = bxC010TaperzPHyPH[k][tj][i]>bmax?bxC010TaperzPHyPH[k][tj][i]:bmax;
bmin = bxC010TaperzPHyPH[k][tj][i]<bmin?bxC010TaperzPHyPH[k][tj][i]:bmin;
      bxC010TaperyPH[k][tj][i] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC010TaperyPH[k][tj][i] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC010TaperyPH[k][tj][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC010TaperyPH[k][tj][i]>amax?axC010TaperyPH[k][tj][i]:amax;
amin = axC010TaperyPH[k][tj][i]<amin?axC010TaperyPH[k][tj][i]:amin;
bmax = bxC010TaperyPH[k][tj][i]>bmax?bxC010TaperyPH[k][tj][i]:bmax;
bmin = bxC010TaperyPH[k][tj][i]<bmin?bxC010TaperyPH[k][tj][i]:bmin;
      bxC010TaperyPHxPH[k][tj][i] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC010TaperyPHxPH[k][tj][i] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC010TaperyPHxPH[k][tj][i]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC010TaperyPHxPH[k][tj][i]>amax?axC010TaperyPHxPH[k][tj][i]:amax;
amin = axC010TaperyPHxPH[k][tj][i]<amin?axC010TaperyPHxPH[k][tj][i]:amin;
bmax = bxC010TaperyPHxPH[k][tj][i]>bmax?bxC010TaperyPHxPH[k][tj][i]:bmax;
bmin = bxC010TaperyPHxPH[k][tj][i]<bmin?bxC010TaperyPHxPH[k][tj][i]:bmin;
      bxC010TaperxPH[k][tj][i] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC010TaperxPH[k][tj][i] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC010TaperxPH[k][tj][i]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC010TaperxPH[k][tj][i]>amax?axC010TaperxPH[k][tj][i]:amax;
amin = axC010TaperxPH[k][tj][i]<amin?axC010TaperxPH[k][tj][i]:amin;
bmax = bxC010TaperxPH[k][tj][i]>bmax?bxC010TaperxPH[k][tj][i]:bmax;
bmin = bxC010TaperxPH[k][tj][i]<bmin?bxC010TaperxPH[k][tj][i]:bmin;
      byC010Taper[k][tj][i] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC010Taper[k][tj][i] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC010Taper[k][tj][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC010Taper[k][tj][i]>amax?ayC010Taper[k][tj][i]:amax;
amin = ayC010Taper[k][tj][i]<amin?ayC010Taper[k][tj][i]:amin;
bmax = byC010Taper[k][tj][i]>bmax?byC010Taper[k][tj][i]:bmax;
bmin = byC010Taper[k][tj][i]<bmin?byC010Taper[k][tj][i]:bmin;
      byC010TaperzPH[k][tj][i] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC010TaperzPH[k][tj][i] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC010TaperzPH[k][tj][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC010TaperzPH[k][tj][i]>amax?ayC010TaperzPH[k][tj][i]:amax;
amin = ayC010TaperzPH[k][tj][i]<amin?ayC010TaperzPH[k][tj][i]:amin;
bmax = byC010TaperzPH[k][tj][i]>bmax?byC010TaperzPH[k][tj][i]:bmax;
bmin = byC010TaperzPH[k][tj][i]<bmin?byC010TaperzPH[k][tj][i]:bmin;
      byC010TaperzPHxPH[k][tj][i] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC010TaperzPHxPH[k][tj][i] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i])*(byC010TaperzPHxPH[k][tj][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC010TaperzPHxPH[k][tj][i]>amax?ayC010TaperzPHxPH[k][tj][i]:amax;
amin = ayC010TaperzPHxPH[k][tj][i]<amin?ayC010TaperzPHxPH[k][tj][i]:amin;
bmax = byC010TaperzPHxPH[k][tj][i]>bmax?byC010TaperzPHxPH[k][tj][i]:bmax;
bmin = byC010TaperzPHxPH[k][tj][i]<bmin?byC010TaperzPHxPH[k][tj][i]:bmin;
      byC010TaperzPHyPH[k][tj][i] = byTaperPH[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC010TaperzPHyPH[k][tj][i] = (yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC010TaperzPHyPH[k][tj][i]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC010TaperzPHyPH[k][tj][i]>amax?ayC010TaperzPHyPH[k][tj][i]:amax;
amin = ayC010TaperzPHyPH[k][tj][i]<amin?ayC010TaperzPHyPH[k][tj][i]:amin;
bmax = byC010TaperzPHyPH[k][tj][i]>bmax?byC010TaperzPHyPH[k][tj][i]:bmax;
bmin = byC010TaperzPHyPH[k][tj][i]<bmin?byC010TaperzPHyPH[k][tj][i]:bmin;
      byC010TaperyPH[k][tj][i] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC010TaperyPH[k][tj][i] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC010TaperyPH[k][tj][i]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC010TaperyPH[k][tj][i]>amax?ayC010TaperyPH[k][tj][i]:amax;
amin = ayC010TaperyPH[k][tj][i]<amin?ayC010TaperyPH[k][tj][i]:amin;
bmax = byC010TaperyPH[k][tj][i]>bmax?byC010TaperyPH[k][tj][i]:bmax;
bmin = byC010TaperyPH[k][tj][i]<bmin?byC010TaperyPH[k][tj][i]:bmin;
      byC010TaperyPHxPH[k][tj][i] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC010TaperyPHxPH[k][tj][i] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC010TaperyPHxPH[k][tj][i]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC010TaperyPHxPH[k][tj][i]>amax?ayC010TaperyPHxPH[k][tj][i]:amax;
amin = ayC010TaperyPHxPH[k][tj][i]<amin?ayC010TaperyPHxPH[k][tj][i]:amin;
bmax = byC010TaperyPHxPH[k][tj][i]>bmax?byC010TaperyPHxPH[k][tj][i]:bmax;
bmin = byC010TaperyPHxPH[k][tj][i]<bmin?byC010TaperyPHxPH[k][tj][i]:bmin;
      byC010TaperxPH[k][tj][i] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC010TaperxPH[k][tj][i] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC010TaperxPH[k][tj][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC010TaperxPH[k][tj][i]>amax?ayC010TaperxPH[k][tj][i]:amax;
amin = ayC010TaperxPH[k][tj][i]<amin?ayC010TaperxPH[k][tj][i]:amin;
bmax = byC010TaperxPH[k][tj][i]>bmax?byC010TaperxPH[k][tj][i]:bmax;
bmin = byC010TaperxPH[k][tj][i]<bmin?byC010TaperxPH[k][tj][i]:bmin;
      bzC010Taper[k][tj][i] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC010Taper[k][tj][i] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC010Taper[k][tj][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC010Taper[k][tj][i]>amax?azC010Taper[k][tj][i]:amax;
amin = azC010Taper[k][tj][i]<amin?azC010Taper[k][tj][i]:amin;
bmax = bzC010Taper[k][tj][i]>bmax?bzC010Taper[k][tj][i]:bmax;
bmin = bzC010Taper[k][tj][i]<bmin?bzC010Taper[k][tj][i]:bmin;
      bzC010TaperzPH[k][tj][i] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC010TaperzPH[k][tj][i] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC010TaperzPH[k][tj][i]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC010TaperzPH[k][tj][i]>amax?azC010TaperzPH[k][tj][i]:amax;
amin = azC010TaperzPH[k][tj][i]<amin?azC010TaperzPH[k][tj][i]:amin;
bmax = bzC010TaperzPH[k][tj][i]>bmax?bzC010TaperzPH[k][tj][i]:bmax;
bmin = bzC010TaperzPH[k][tj][i]<bmin?bzC010TaperzPH[k][tj][i]:bmin;
      bzC010TaperzPHxPH[k][tj][i] = bzTaperPH[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC010TaperzPHxPH[k][tj][i] = (zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC010TaperzPHxPH[k][tj][i]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC010TaperzPHxPH[k][tj][i]>amax?azC010TaperzPHxPH[k][tj][i]:amax;
amin = azC010TaperzPHxPH[k][tj][i]<amin?azC010TaperzPHxPH[k][tj][i]:amin;
bmax = bzC010TaperzPHxPH[k][tj][i]>bmax?bzC010TaperzPHxPH[k][tj][i]:bmax;
bmin = bzC010TaperzPHxPH[k][tj][i]<bmin?bzC010TaperzPHxPH[k][tj][i]:bmin;
      bzC010TaperzPHyPH[k][tj][i] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC010TaperzPHyPH[k][tj][i] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC010TaperzPHyPH[k][tj][i]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC010TaperzPHyPH[k][tj][i]>amax?azC010TaperzPHyPH[k][tj][i]:amax;
amin = azC010TaperzPHyPH[k][tj][i]<amin?azC010TaperzPHyPH[k][tj][i]:amin;
bmax = bzC010TaperzPHyPH[k][tj][i]>bmax?bzC010TaperzPHyPH[k][tj][i]:bmax;
bmin = bzC010TaperzPHyPH[k][tj][i]<bmin?bzC010TaperzPHyPH[k][tj][i]:bmin;
      bzC010TaperyPH[k][tj][i] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC010TaperyPH[k][tj][i] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC010TaperyPH[k][tj][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC010TaperyPH[k][tj][i]>amax?azC010TaperyPH[k][tj][i]:amax;
amin = azC010TaperyPH[k][tj][i]<amin?azC010TaperyPH[k][tj][i]:amin;
bmax = bzC010TaperyPH[k][tj][i]>bmax?bzC010TaperyPH[k][tj][i]:bmax;
bmin = bzC010TaperyPH[k][tj][i]<bmin?bzC010TaperyPH[k][tj][i]:bmin;
      bzC010TaperyPHxPH[k][tj][i] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC010TaperyPHxPH[k][tj][i] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j])*(bzC010TaperyPHxPH[k][tj][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC010TaperyPHxPH[k][tj][i]>amax?azC010TaperyPHxPH[k][tj][i]:amax;
amin = azC010TaperyPHxPH[k][tj][i]<amin?azC010TaperyPHxPH[k][tj][i]:amin;
bmax = bzC010TaperyPHxPH[k][tj][i]>bmax?bzC010TaperyPHxPH[k][tj][i]:bmax;
bmin = bzC010TaperyPHxPH[k][tj][i]<bmin?bzC010TaperyPHxPH[k][tj][i]:bmin;
      bzC010TaperxPH[k][tj][i] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC010TaperxPH[k][tj][i] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC010TaperxPH[k][tj][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC010TaperxPH[k][tj][i]>amax?azC010TaperxPH[k][tj][i]:amax;
amin = azC010TaperxPH[k][tj][i]<amin?azC010TaperxPH[k][tj][i]:amin;
bmax = bzC010TaperxPH[k][tj][i]>bmax?bzC010TaperxPH[k][tj][i]:bmax;
bmin = bzC010TaperxPH[k][tj][i]<bmin?bzC010TaperxPH[k][tj][i]:bmin;
    }
    for(int i=xEnd1;i<NX;++i) {
      int ti = i-xEnd1;
      bxC110Taper[k][tj][ti] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC110Taper[k][tj][ti] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC110Taper[k][tj][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC110Taper[k][tj][ti]>amax?axC110Taper[k][tj][ti]:amax;
amin = axC110Taper[k][tj][ti]<amin?axC110Taper[k][tj][ti]:amin;
bmax = bxC110Taper[k][tj][ti]>bmax?bxC110Taper[k][tj][ti]:bmax;
bmin = bxC110Taper[k][tj][ti]<bmin?bxC110Taper[k][tj][ti]:bmin;
      bxC110TaperzPH[k][tj][ti] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC110TaperzPH[k][tj][ti] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC110TaperzPH[k][tj][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC110TaperzPH[k][tj][ti]>amax?axC110TaperzPH[k][tj][ti]:amax;
amin = axC110TaperzPH[k][tj][ti]<amin?axC110TaperzPH[k][tj][ti]:amin;
bmax = bxC110TaperzPH[k][tj][ti]>bmax?bxC110TaperzPH[k][tj][ti]:bmax;
bmin = bxC110TaperzPH[k][tj][ti]<bmin?bxC110TaperzPH[k][tj][ti]:bmin;
      bxC110TaperzPHxPH[k][tj][ti] = bxTaperPH[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC110TaperzPHxPH[k][tj][ti] = (xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC110TaperzPHxPH[k][tj][ti]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC110TaperzPHxPH[k][tj][ti]>amax?axC110TaperzPHxPH[k][tj][ti]:amax;
amin = axC110TaperzPHxPH[k][tj][ti]<amin?axC110TaperzPHxPH[k][tj][ti]:amin;
bmax = bxC110TaperzPHxPH[k][tj][ti]>bmax?bxC110TaperzPHxPH[k][tj][ti]:bmax;
bmin = bxC110TaperzPHxPH[k][tj][ti]<bmin?bxC110TaperzPHxPH[k][tj][ti]:bmin;
      bxC110TaperzPHyPH[k][tj][ti] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC110TaperzPHyPH[k][tj][ti] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j])*(bxC110TaperzPHyPH[k][tj][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC110TaperzPHyPH[k][tj][ti]>amax?axC110TaperzPHyPH[k][tj][ti]:amax;
amin = axC110TaperzPHyPH[k][tj][ti]<amin?axC110TaperzPHyPH[k][tj][ti]:amin;
bmax = bxC110TaperzPHyPH[k][tj][ti]>bmax?bxC110TaperzPHyPH[k][tj][ti]:bmax;
bmin = bxC110TaperzPHyPH[k][tj][ti]<bmin?bxC110TaperzPHyPH[k][tj][ti]:bmin;
      bxC110TaperyPH[k][tj][ti] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC110TaperyPH[k][tj][ti] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC110TaperyPH[k][tj][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC110TaperyPH[k][tj][ti]>amax?axC110TaperyPH[k][tj][ti]:amax;
amin = axC110TaperyPH[k][tj][ti]<amin?axC110TaperyPH[k][tj][ti]:amin;
bmax = bxC110TaperyPH[k][tj][ti]>bmax?bxC110TaperyPH[k][tj][ti]:bmax;
bmin = bxC110TaperyPH[k][tj][ti]<bmin?bxC110TaperyPH[k][tj][ti]:bmin;
      bxC110TaperyPHxPH[k][tj][ti] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC110TaperyPHxPH[k][tj][ti] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC110TaperyPHxPH[k][tj][ti]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC110TaperyPHxPH[k][tj][ti]>amax?axC110TaperyPHxPH[k][tj][ti]:amax;
amin = axC110TaperyPHxPH[k][tj][ti]<amin?axC110TaperyPHxPH[k][tj][ti]:amin;
bmax = bxC110TaperyPHxPH[k][tj][ti]>bmax?bxC110TaperyPHxPH[k][tj][ti]:bmax;
bmin = bxC110TaperyPHxPH[k][tj][ti]<bmin?bxC110TaperyPHxPH[k][tj][ti]:bmin;
      bxC110TaperxPH[k][tj][ti] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC110TaperxPH[k][tj][ti] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC110TaperxPH[k][tj][ti]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC110TaperxPH[k][tj][ti]>amax?axC110TaperxPH[k][tj][ti]:amax;
amin = axC110TaperxPH[k][tj][ti]<amin?axC110TaperxPH[k][tj][ti]:amin;
bmax = bxC110TaperxPH[k][tj][ti]>bmax?bxC110TaperxPH[k][tj][ti]:bmax;
bmin = bxC110TaperxPH[k][tj][ti]<bmin?bxC110TaperxPH[k][tj][ti]:bmin;
      byC110Taper[k][tj][ti] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC110Taper[k][tj][ti] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC110Taper[k][tj][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC110Taper[k][tj][ti]>amax?ayC110Taper[k][tj][ti]:amax;
amin = ayC110Taper[k][tj][ti]<amin?ayC110Taper[k][tj][ti]:amin;
bmax = byC110Taper[k][tj][ti]>bmax?byC110Taper[k][tj][ti]:bmax;
bmin = byC110Taper[k][tj][ti]<bmin?byC110Taper[k][tj][ti]:bmin;
      byC110TaperzPH[k][tj][ti] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC110TaperzPH[k][tj][ti] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC110TaperzPH[k][tj][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC110TaperzPH[k][tj][ti]>amax?ayC110TaperzPH[k][tj][ti]:amax;
amin = ayC110TaperzPH[k][tj][ti]<amin?ayC110TaperzPH[k][tj][ti]:amin;
bmax = byC110TaperzPH[k][tj][ti]>bmax?byC110TaperzPH[k][tj][ti]:bmax;
bmin = byC110TaperzPH[k][tj][ti]<bmin?byC110TaperzPH[k][tj][ti]:bmin;
      byC110TaperzPHxPH[k][tj][ti] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC110TaperzPHxPH[k][tj][ti] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i])*(byC110TaperzPHxPH[k][tj][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC110TaperzPHxPH[k][tj][ti]>amax?ayC110TaperzPHxPH[k][tj][ti]:amax;
amin = ayC110TaperzPHxPH[k][tj][ti]<amin?ayC110TaperzPHxPH[k][tj][ti]:amin;
bmax = byC110TaperzPHxPH[k][tj][ti]>bmax?byC110TaperzPHxPH[k][tj][ti]:bmax;
bmin = byC110TaperzPHxPH[k][tj][ti]<bmin?byC110TaperzPHxPH[k][tj][ti]:bmin;
      byC110TaperzPHyPH[k][tj][ti] = byTaperPH[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC110TaperzPHyPH[k][tj][ti] = (yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC110TaperzPHyPH[k][tj][ti]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC110TaperzPHyPH[k][tj][ti]>amax?ayC110TaperzPHyPH[k][tj][ti]:amax;
amin = ayC110TaperzPHyPH[k][tj][ti]<amin?ayC110TaperzPHyPH[k][tj][ti]:amin;
bmax = byC110TaperzPHyPH[k][tj][ti]>bmax?byC110TaperzPHyPH[k][tj][ti]:bmax;
bmin = byC110TaperzPHyPH[k][tj][ti]<bmin?byC110TaperzPHyPH[k][tj][ti]:bmin;
      byC110TaperyPH[k][tj][ti] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC110TaperyPH[k][tj][ti] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC110TaperyPH[k][tj][ti]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC110TaperyPH[k][tj][ti]>amax?ayC110TaperyPH[k][tj][ti]:amax;
amin = ayC110TaperyPH[k][tj][ti]<amin?ayC110TaperyPH[k][tj][ti]:amin;
bmax = byC110TaperyPH[k][tj][ti]>bmax?byC110TaperyPH[k][tj][ti]:bmax;
bmin = byC110TaperyPH[k][tj][ti]<bmin?byC110TaperyPH[k][tj][ti]:bmin;
      byC110TaperyPHxPH[k][tj][ti] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC110TaperyPHxPH[k][tj][ti] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC110TaperyPHxPH[k][tj][ti]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC110TaperyPHxPH[k][tj][ti]>amax?ayC110TaperyPHxPH[k][tj][ti]:amax;
amin = ayC110TaperyPHxPH[k][tj][ti]<amin?ayC110TaperyPHxPH[k][tj][ti]:amin;
bmax = byC110TaperyPHxPH[k][tj][ti]>bmax?byC110TaperyPHxPH[k][tj][ti]:bmax;
bmin = byC110TaperyPHxPH[k][tj][ti]<bmin?byC110TaperyPHxPH[k][tj][ti]:bmin;
      byC110TaperxPH[k][tj][ti] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC110TaperxPH[k][tj][ti] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC110TaperxPH[k][tj][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC110TaperxPH[k][tj][ti]>amax?ayC110TaperxPH[k][tj][ti]:amax;
amin = ayC110TaperxPH[k][tj][ti]<amin?ayC110TaperxPH[k][tj][ti]:amin;
bmax = byC110TaperxPH[k][tj][ti]>bmax?byC110TaperxPH[k][tj][ti]:bmax;
bmin = byC110TaperxPH[k][tj][ti]<bmin?byC110TaperxPH[k][tj][ti]:bmin;
      bzC110Taper[k][tj][ti] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC110Taper[k][tj][ti] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC110Taper[k][tj][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC110Taper[k][tj][ti]>amax?azC110Taper[k][tj][ti]:amax;
amin = azC110Taper[k][tj][ti]<amin?azC110Taper[k][tj][ti]:amin;
bmax = bzC110Taper[k][tj][ti]>bmax?bzC110Taper[k][tj][ti]:bmax;
bmin = bzC110Taper[k][tj][ti]<bmin?bzC110Taper[k][tj][ti]:bmin;
      bzC110TaperzPH[k][tj][ti] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC110TaperzPH[k][tj][ti] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC110TaperzPH[k][tj][ti]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC110TaperzPH[k][tj][ti]>amax?azC110TaperzPH[k][tj][ti]:amax;
amin = azC110TaperzPH[k][tj][ti]<amin?azC110TaperzPH[k][tj][ti]:amin;
bmax = bzC110TaperzPH[k][tj][ti]>bmax?bzC110TaperzPH[k][tj][ti]:bmax;
bmin = bzC110TaperzPH[k][tj][ti]<bmin?bzC110TaperzPH[k][tj][ti]:bmin;
      bzC110TaperzPHxPH[k][tj][ti] = bzTaperPH[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC110TaperzPHxPH[k][tj][ti] = (zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC110TaperzPHxPH[k][tj][ti]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC110TaperzPHxPH[k][tj][ti]>amax?azC110TaperzPHxPH[k][tj][ti]:amax;
amin = azC110TaperzPHxPH[k][tj][ti]<amin?azC110TaperzPHxPH[k][tj][ti]:amin;
bmax = bzC110TaperzPHxPH[k][tj][ti]>bmax?bzC110TaperzPHxPH[k][tj][ti]:bmax;
bmin = bzC110TaperzPHxPH[k][tj][ti]<bmin?bzC110TaperzPHxPH[k][tj][ti]:bmin;
      bzC110TaperzPHyPH[k][tj][ti] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC110TaperzPHyPH[k][tj][ti] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC110TaperzPHyPH[k][tj][ti]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC110TaperzPHyPH[k][tj][ti]>amax?azC110TaperzPHyPH[k][tj][ti]:amax;
amin = azC110TaperzPHyPH[k][tj][ti]<amin?azC110TaperzPHyPH[k][tj][ti]:amin;
bmax = bzC110TaperzPHyPH[k][tj][ti]>bmax?bzC110TaperzPHyPH[k][tj][ti]:bmax;
bmin = bzC110TaperzPHyPH[k][tj][ti]<bmin?bzC110TaperzPHyPH[k][tj][ti]:bmin;
      bzC110TaperyPH[k][tj][ti] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC110TaperyPH[k][tj][ti] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC110TaperyPH[k][tj][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC110TaperyPH[k][tj][ti]>amax?azC110TaperyPH[k][tj][ti]:amax;
amin = azC110TaperyPH[k][tj][ti]<amin?azC110TaperyPH[k][tj][ti]:amin;
bmax = bzC110TaperyPH[k][tj][ti]>bmax?bzC110TaperyPH[k][tj][ti]:bmax;
bmin = bzC110TaperyPH[k][tj][ti]<bmin?bzC110TaperyPH[k][tj][ti]:bmin;
      bzC110TaperyPHxPH[k][tj][ti] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC110TaperyPHxPH[k][tj][ti] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j])*(bzC110TaperyPHxPH[k][tj][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC110TaperyPHxPH[k][tj][ti]>amax?azC110TaperyPHxPH[k][tj][ti]:amax;
amin = azC110TaperyPHxPH[k][tj][ti]<amin?azC110TaperyPHxPH[k][tj][ti]:amin;
bmax = bzC110TaperyPHxPH[k][tj][ti]>bmax?bzC110TaperyPHxPH[k][tj][ti]:bmax;
bmin = bzC110TaperyPHxPH[k][tj][ti]<bmin?bzC110TaperyPHxPH[k][tj][ti]:bmin;
      bzC110TaperxPH[k][tj][ti] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC110TaperxPH[k][tj][ti] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC110TaperxPH[k][tj][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC110TaperxPH[k][tj][ti]>amax?azC110TaperxPH[k][tj][ti]:amax;
amin = azC110TaperxPH[k][tj][ti]<amin?azC110TaperxPH[k][tj][ti]:amin;
bmax = bzC110TaperxPH[k][tj][ti]>bmax?bzC110TaperxPH[k][tj][ti]:bmax;
bmin = bzC110TaperxPH[k][tj][ti]<bmin?bzC110TaperxPH[k][tj][ti]:bmin;
    }
  }
}
for(int k=zEnd1;k<NZ;++k) {
  int tk = k-zEnd1;
  for(int j=0;j<_nYmin;++j) {
    for(int i=0;i<_nXmin;++i) {
      bxC001Taper[tk][j][i] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC001Taper[tk][j][i] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC001Taper[tk][j][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC001Taper[tk][j][i]>amax?axC001Taper[tk][j][i]:amax;
amin = axC001Taper[tk][j][i]<amin?axC001Taper[tk][j][i]:amin;
bmax = bxC001Taper[tk][j][i]>bmax?bxC001Taper[tk][j][i]:bmax;
bmin = bxC001Taper[tk][j][i]<bmin?bxC001Taper[tk][j][i]:bmin;
      bxC001TaperzPH[tk][j][i] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC001TaperzPH[tk][j][i] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC001TaperzPH[tk][j][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC001TaperzPH[tk][j][i]>amax?axC001TaperzPH[tk][j][i]:amax;
amin = axC001TaperzPH[tk][j][i]<amin?axC001TaperzPH[tk][j][i]:amin;
bmax = bxC001TaperzPH[tk][j][i]>bmax?bxC001TaperzPH[tk][j][i]:bmax;
bmin = bxC001TaperzPH[tk][j][i]<bmin?bxC001TaperzPH[tk][j][i]:bmin;
      bxC001TaperzPHxPH[tk][j][i] = bxTaperPH[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC001TaperzPHxPH[tk][j][i] = (xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC001TaperzPHxPH[tk][j][i]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC001TaperzPHxPH[tk][j][i]>amax?axC001TaperzPHxPH[tk][j][i]:amax;
amin = axC001TaperzPHxPH[tk][j][i]<amin?axC001TaperzPHxPH[tk][j][i]:amin;
bmax = bxC001TaperzPHxPH[tk][j][i]>bmax?bxC001TaperzPHxPH[tk][j][i]:bmax;
bmin = bxC001TaperzPHxPH[tk][j][i]<bmin?bxC001TaperzPHxPH[tk][j][i]:bmin;
      bxC001TaperzPHyPH[tk][j][i] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC001TaperzPHyPH[tk][j][i] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j])*(bxC001TaperzPHyPH[tk][j][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC001TaperzPHyPH[tk][j][i]>amax?axC001TaperzPHyPH[tk][j][i]:amax;
amin = axC001TaperzPHyPH[tk][j][i]<amin?axC001TaperzPHyPH[tk][j][i]:amin;
bmax = bxC001TaperzPHyPH[tk][j][i]>bmax?bxC001TaperzPHyPH[tk][j][i]:bmax;
bmin = bxC001TaperzPHyPH[tk][j][i]<bmin?bxC001TaperzPHyPH[tk][j][i]:bmin;
      bxC001TaperyPH[tk][j][i] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC001TaperyPH[tk][j][i] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC001TaperyPH[tk][j][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC001TaperyPH[tk][j][i]>amax?axC001TaperyPH[tk][j][i]:amax;
amin = axC001TaperyPH[tk][j][i]<amin?axC001TaperyPH[tk][j][i]:amin;
bmax = bxC001TaperyPH[tk][j][i]>bmax?bxC001TaperyPH[tk][j][i]:bmax;
bmin = bxC001TaperyPH[tk][j][i]<bmin?bxC001TaperyPH[tk][j][i]:bmin;
      bxC001TaperyPHxPH[tk][j][i] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC001TaperyPHxPH[tk][j][i] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC001TaperyPHxPH[tk][j][i]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC001TaperyPHxPH[tk][j][i]>amax?axC001TaperyPHxPH[tk][j][i]:amax;
amin = axC001TaperyPHxPH[tk][j][i]<amin?axC001TaperyPHxPH[tk][j][i]:amin;
bmax = bxC001TaperyPHxPH[tk][j][i]>bmax?bxC001TaperyPHxPH[tk][j][i]:bmax;
bmin = bxC001TaperyPHxPH[tk][j][i]<bmin?bxC001TaperyPHxPH[tk][j][i]:bmin;
      bxC001TaperxPH[tk][j][i] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC001TaperxPH[tk][j][i] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC001TaperxPH[tk][j][i]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC001TaperxPH[tk][j][i]>amax?axC001TaperxPH[tk][j][i]:amax;
amin = axC001TaperxPH[tk][j][i]<amin?axC001TaperxPH[tk][j][i]:amin;
bmax = bxC001TaperxPH[tk][j][i]>bmax?bxC001TaperxPH[tk][j][i]:bmax;
bmin = bxC001TaperxPH[tk][j][i]<bmin?bxC001TaperxPH[tk][j][i]:bmin;
      byC001Taper[tk][j][i] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC001Taper[tk][j][i] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC001Taper[tk][j][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC001Taper[tk][j][i]>amax?ayC001Taper[tk][j][i]:amax;
amin = ayC001Taper[tk][j][i]<amin?ayC001Taper[tk][j][i]:amin;
bmax = byC001Taper[tk][j][i]>bmax?byC001Taper[tk][j][i]:bmax;
bmin = byC001Taper[tk][j][i]<bmin?byC001Taper[tk][j][i]:bmin;
      byC001TaperzPH[tk][j][i] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC001TaperzPH[tk][j][i] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC001TaperzPH[tk][j][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC001TaperzPH[tk][j][i]>amax?ayC001TaperzPH[tk][j][i]:amax;
amin = ayC001TaperzPH[tk][j][i]<amin?ayC001TaperzPH[tk][j][i]:amin;
bmax = byC001TaperzPH[tk][j][i]>bmax?byC001TaperzPH[tk][j][i]:bmax;
bmin = byC001TaperzPH[tk][j][i]<bmin?byC001TaperzPH[tk][j][i]:bmin;
      byC001TaperzPHxPH[tk][j][i] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC001TaperzPHxPH[tk][j][i] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i])*(byC001TaperzPHxPH[tk][j][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC001TaperzPHxPH[tk][j][i]>amax?ayC001TaperzPHxPH[tk][j][i]:amax;
amin = ayC001TaperzPHxPH[tk][j][i]<amin?ayC001TaperzPHxPH[tk][j][i]:amin;
bmax = byC001TaperzPHxPH[tk][j][i]>bmax?byC001TaperzPHxPH[tk][j][i]:bmax;
bmin = byC001TaperzPHxPH[tk][j][i]<bmin?byC001TaperzPHxPH[tk][j][i]:bmin;
      byC001TaperzPHyPH[tk][j][i] = byTaperPH[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC001TaperzPHyPH[tk][j][i] = (yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC001TaperzPHyPH[tk][j][i]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC001TaperzPHyPH[tk][j][i]>amax?ayC001TaperzPHyPH[tk][j][i]:amax;
amin = ayC001TaperzPHyPH[tk][j][i]<amin?ayC001TaperzPHyPH[tk][j][i]:amin;
bmax = byC001TaperzPHyPH[tk][j][i]>bmax?byC001TaperzPHyPH[tk][j][i]:bmax;
bmin = byC001TaperzPHyPH[tk][j][i]<bmin?byC001TaperzPHyPH[tk][j][i]:bmin;
      byC001TaperyPH[tk][j][i] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC001TaperyPH[tk][j][i] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC001TaperyPH[tk][j][i]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC001TaperyPH[tk][j][i]>amax?ayC001TaperyPH[tk][j][i]:amax;
amin = ayC001TaperyPH[tk][j][i]<amin?ayC001TaperyPH[tk][j][i]:amin;
bmax = byC001TaperyPH[tk][j][i]>bmax?byC001TaperyPH[tk][j][i]:bmax;
bmin = byC001TaperyPH[tk][j][i]<bmin?byC001TaperyPH[tk][j][i]:bmin;
      byC001TaperyPHxPH[tk][j][i] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC001TaperyPHxPH[tk][j][i] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC001TaperyPHxPH[tk][j][i]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC001TaperyPHxPH[tk][j][i]>amax?ayC001TaperyPHxPH[tk][j][i]:amax;
amin = ayC001TaperyPHxPH[tk][j][i]<amin?ayC001TaperyPHxPH[tk][j][i]:amin;
bmax = byC001TaperyPHxPH[tk][j][i]>bmax?byC001TaperyPHxPH[tk][j][i]:bmax;
bmin = byC001TaperyPHxPH[tk][j][i]<bmin?byC001TaperyPHxPH[tk][j][i]:bmin;
      byC001TaperxPH[tk][j][i] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC001TaperxPH[tk][j][i] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC001TaperxPH[tk][j][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC001TaperxPH[tk][j][i]>amax?ayC001TaperxPH[tk][j][i]:amax;
amin = ayC001TaperxPH[tk][j][i]<amin?ayC001TaperxPH[tk][j][i]:amin;
bmax = byC001TaperxPH[tk][j][i]>bmax?byC001TaperxPH[tk][j][i]:bmax;
bmin = byC001TaperxPH[tk][j][i]<bmin?byC001TaperxPH[tk][j][i]:bmin;
      bzC001Taper[tk][j][i] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC001Taper[tk][j][i] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC001Taper[tk][j][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC001Taper[tk][j][i]>amax?azC001Taper[tk][j][i]:amax;
amin = azC001Taper[tk][j][i]<amin?azC001Taper[tk][j][i]:amin;
bmax = bzC001Taper[tk][j][i]>bmax?bzC001Taper[tk][j][i]:bmax;
bmin = bzC001Taper[tk][j][i]<bmin?bzC001Taper[tk][j][i]:bmin;
      bzC001TaperzPH[tk][j][i] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC001TaperzPH[tk][j][i] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC001TaperzPH[tk][j][i]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC001TaperzPH[tk][j][i]>amax?azC001TaperzPH[tk][j][i]:amax;
amin = azC001TaperzPH[tk][j][i]<amin?azC001TaperzPH[tk][j][i]:amin;
bmax = bzC001TaperzPH[tk][j][i]>bmax?bzC001TaperzPH[tk][j][i]:bmax;
bmin = bzC001TaperzPH[tk][j][i]<bmin?bzC001TaperzPH[tk][j][i]:bmin;
      bzC001TaperzPHxPH[tk][j][i] = bzTaperPH[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC001TaperzPHxPH[tk][j][i] = (zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC001TaperzPHxPH[tk][j][i]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC001TaperzPHxPH[tk][j][i]>amax?azC001TaperzPHxPH[tk][j][i]:amax;
amin = azC001TaperzPHxPH[tk][j][i]<amin?azC001TaperzPHxPH[tk][j][i]:amin;
bmax = bzC001TaperzPHxPH[tk][j][i]>bmax?bzC001TaperzPHxPH[tk][j][i]:bmax;
bmin = bzC001TaperzPHxPH[tk][j][i]<bmin?bzC001TaperzPHxPH[tk][j][i]:bmin;
      bzC001TaperzPHyPH[tk][j][i] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC001TaperzPHyPH[tk][j][i] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC001TaperzPHyPH[tk][j][i]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC001TaperzPHyPH[tk][j][i]>amax?azC001TaperzPHyPH[tk][j][i]:amax;
amin = azC001TaperzPHyPH[tk][j][i]<amin?azC001TaperzPHyPH[tk][j][i]:amin;
bmax = bzC001TaperzPHyPH[tk][j][i]>bmax?bzC001TaperzPHyPH[tk][j][i]:bmax;
bmin = bzC001TaperzPHyPH[tk][j][i]<bmin?bzC001TaperzPHyPH[tk][j][i]:bmin;
      bzC001TaperyPH[tk][j][i] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC001TaperyPH[tk][j][i] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC001TaperyPH[tk][j][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC001TaperyPH[tk][j][i]>amax?azC001TaperyPH[tk][j][i]:amax;
amin = azC001TaperyPH[tk][j][i]<amin?azC001TaperyPH[tk][j][i]:amin;
bmax = bzC001TaperyPH[tk][j][i]>bmax?bzC001TaperyPH[tk][j][i]:bmax;
bmin = bzC001TaperyPH[tk][j][i]<bmin?bzC001TaperyPH[tk][j][i]:bmin;
      bzC001TaperyPHxPH[tk][j][i] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC001TaperyPHxPH[tk][j][i] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j])*(bzC001TaperyPHxPH[tk][j][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC001TaperyPHxPH[tk][j][i]>amax?azC001TaperyPHxPH[tk][j][i]:amax;
amin = azC001TaperyPHxPH[tk][j][i]<amin?azC001TaperyPHxPH[tk][j][i]:amin;
bmax = bzC001TaperyPHxPH[tk][j][i]>bmax?bzC001TaperyPHxPH[tk][j][i]:bmax;
bmin = bzC001TaperyPHxPH[tk][j][i]<bmin?bzC001TaperyPHxPH[tk][j][i]:bmin;
      bzC001TaperxPH[tk][j][i] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC001TaperxPH[tk][j][i] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC001TaperxPH[tk][j][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC001TaperxPH[tk][j][i]>amax?azC001TaperxPH[tk][j][i]:amax;
amin = azC001TaperxPH[tk][j][i]<amin?azC001TaperxPH[tk][j][i]:amin;
bmax = bzC001TaperxPH[tk][j][i]>bmax?bzC001TaperxPH[tk][j][i]:bmax;
bmin = bzC001TaperxPH[tk][j][i]<bmin?bzC001TaperxPH[tk][j][i]:bmin;
    }
    for(int i=xEnd1;i<NX;++i) {
      int ti = i-xEnd1;
      bxC101Taper[tk][j][ti] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC101Taper[tk][j][ti] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC101Taper[tk][j][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC101Taper[tk][j][ti]>amax?axC101Taper[tk][j][ti]:amax;
amin = axC101Taper[tk][j][ti]<amin?axC101Taper[tk][j][ti]:amin;
bmax = bxC101Taper[tk][j][ti]>bmax?bxC101Taper[tk][j][ti]:bmax;
bmin = bxC101Taper[tk][j][ti]<bmin?bxC101Taper[tk][j][ti]:bmin;
      bxC101TaperzPH[tk][j][ti] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC101TaperzPH[tk][j][ti] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC101TaperzPH[tk][j][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC101TaperzPH[tk][j][ti]>amax?axC101TaperzPH[tk][j][ti]:amax;
amin = axC101TaperzPH[tk][j][ti]<amin?axC101TaperzPH[tk][j][ti]:amin;
bmax = bxC101TaperzPH[tk][j][ti]>bmax?bxC101TaperzPH[tk][j][ti]:bmax;
bmin = bxC101TaperzPH[tk][j][ti]<bmin?bxC101TaperzPH[tk][j][ti]:bmin;
      bxC101TaperzPHxPH[tk][j][ti] = bxTaperPH[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC101TaperzPHxPH[tk][j][ti] = (xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC101TaperzPHxPH[tk][j][ti]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC101TaperzPHxPH[tk][j][ti]>amax?axC101TaperzPHxPH[tk][j][ti]:amax;
amin = axC101TaperzPHxPH[tk][j][ti]<amin?axC101TaperzPHxPH[tk][j][ti]:amin;
bmax = bxC101TaperzPHxPH[tk][j][ti]>bmax?bxC101TaperzPHxPH[tk][j][ti]:bmax;
bmin = bxC101TaperzPHxPH[tk][j][ti]<bmin?bxC101TaperzPHxPH[tk][j][ti]:bmin;
      bxC101TaperzPHyPH[tk][j][ti] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC101TaperzPHyPH[tk][j][ti] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j])*(bxC101TaperzPHyPH[tk][j][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC101TaperzPHyPH[tk][j][ti]>amax?axC101TaperzPHyPH[tk][j][ti]:amax;
amin = axC101TaperzPHyPH[tk][j][ti]<amin?axC101TaperzPHyPH[tk][j][ti]:amin;
bmax = bxC101TaperzPHyPH[tk][j][ti]>bmax?bxC101TaperzPHyPH[tk][j][ti]:bmax;
bmin = bxC101TaperzPHyPH[tk][j][ti]<bmin?bxC101TaperzPHyPH[tk][j][ti]:bmin;
      bxC101TaperyPH[tk][j][ti] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC101TaperyPH[tk][j][ti] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC101TaperyPH[tk][j][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC101TaperyPH[tk][j][ti]>amax?axC101TaperyPH[tk][j][ti]:amax;
amin = axC101TaperyPH[tk][j][ti]<amin?axC101TaperyPH[tk][j][ti]:amin;
bmax = bxC101TaperyPH[tk][j][ti]>bmax?bxC101TaperyPH[tk][j][ti]:bmax;
bmin = bxC101TaperyPH[tk][j][ti]<bmin?bxC101TaperyPH[tk][j][ti]:bmin;
      bxC101TaperyPHxPH[tk][j][ti] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC101TaperyPHxPH[tk][j][ti] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC101TaperyPHxPH[tk][j][ti]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC101TaperyPHxPH[tk][j][ti]>amax?axC101TaperyPHxPH[tk][j][ti]:amax;
amin = axC101TaperyPHxPH[tk][j][ti]<amin?axC101TaperyPHxPH[tk][j][ti]:amin;
bmax = bxC101TaperyPHxPH[tk][j][ti]>bmax?bxC101TaperyPHxPH[tk][j][ti]:bmax;
bmin = bxC101TaperyPHxPH[tk][j][ti]<bmin?bxC101TaperyPHxPH[tk][j][ti]:bmin;
      bxC101TaperxPH[tk][j][ti] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC101TaperxPH[tk][j][ti] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC101TaperxPH[tk][j][ti]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC101TaperxPH[tk][j][ti]>amax?axC101TaperxPH[tk][j][ti]:amax;
amin = axC101TaperxPH[tk][j][ti]<amin?axC101TaperxPH[tk][j][ti]:amin;
bmax = bxC101TaperxPH[tk][j][ti]>bmax?bxC101TaperxPH[tk][j][ti]:bmax;
bmin = bxC101TaperxPH[tk][j][ti]<bmin?bxC101TaperxPH[tk][j][ti]:bmin;
      byC101Taper[tk][j][ti] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC101Taper[tk][j][ti] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC101Taper[tk][j][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC101Taper[tk][j][ti]>amax?ayC101Taper[tk][j][ti]:amax;
amin = ayC101Taper[tk][j][ti]<amin?ayC101Taper[tk][j][ti]:amin;
bmax = byC101Taper[tk][j][ti]>bmax?byC101Taper[tk][j][ti]:bmax;
bmin = byC101Taper[tk][j][ti]<bmin?byC101Taper[tk][j][ti]:bmin;
      byC101TaperzPH[tk][j][ti] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC101TaperzPH[tk][j][ti] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC101TaperzPH[tk][j][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC101TaperzPH[tk][j][ti]>amax?ayC101TaperzPH[tk][j][ti]:amax;
amin = ayC101TaperzPH[tk][j][ti]<amin?ayC101TaperzPH[tk][j][ti]:amin;
bmax = byC101TaperzPH[tk][j][ti]>bmax?byC101TaperzPH[tk][j][ti]:bmax;
bmin = byC101TaperzPH[tk][j][ti]<bmin?byC101TaperzPH[tk][j][ti]:bmin;
      byC101TaperzPHxPH[tk][j][ti] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC101TaperzPHxPH[tk][j][ti] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i])*(byC101TaperzPHxPH[tk][j][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC101TaperzPHxPH[tk][j][ti]>amax?ayC101TaperzPHxPH[tk][j][ti]:amax;
amin = ayC101TaperzPHxPH[tk][j][ti]<amin?ayC101TaperzPHxPH[tk][j][ti]:amin;
bmax = byC101TaperzPHxPH[tk][j][ti]>bmax?byC101TaperzPHxPH[tk][j][ti]:bmax;
bmin = byC101TaperzPHxPH[tk][j][ti]<bmin?byC101TaperzPHxPH[tk][j][ti]:bmin;
      byC101TaperzPHyPH[tk][j][ti] = byTaperPH[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC101TaperzPHyPH[tk][j][ti] = (yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC101TaperzPHyPH[tk][j][ti]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC101TaperzPHyPH[tk][j][ti]>amax?ayC101TaperzPHyPH[tk][j][ti]:amax;
amin = ayC101TaperzPHyPH[tk][j][ti]<amin?ayC101TaperzPHyPH[tk][j][ti]:amin;
bmax = byC101TaperzPHyPH[tk][j][ti]>bmax?byC101TaperzPHyPH[tk][j][ti]:bmax;
bmin = byC101TaperzPHyPH[tk][j][ti]<bmin?byC101TaperzPHyPH[tk][j][ti]:bmin;
      byC101TaperyPH[tk][j][ti] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC101TaperyPH[tk][j][ti] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC101TaperyPH[tk][j][ti]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC101TaperyPH[tk][j][ti]>amax?ayC101TaperyPH[tk][j][ti]:amax;
amin = ayC101TaperyPH[tk][j][ti]<amin?ayC101TaperyPH[tk][j][ti]:amin;
bmax = byC101TaperyPH[tk][j][ti]>bmax?byC101TaperyPH[tk][j][ti]:bmax;
bmin = byC101TaperyPH[tk][j][ti]<bmin?byC101TaperyPH[tk][j][ti]:bmin;
      byC101TaperyPHxPH[tk][j][ti] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC101TaperyPHxPH[tk][j][ti] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC101TaperyPHxPH[tk][j][ti]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC101TaperyPHxPH[tk][j][ti]>amax?ayC101TaperyPHxPH[tk][j][ti]:amax;
amin = ayC101TaperyPHxPH[tk][j][ti]<amin?ayC101TaperyPHxPH[tk][j][ti]:amin;
bmax = byC101TaperyPHxPH[tk][j][ti]>bmax?byC101TaperyPHxPH[tk][j][ti]:bmax;
bmin = byC101TaperyPHxPH[tk][j][ti]<bmin?byC101TaperyPHxPH[tk][j][ti]:bmin;
      byC101TaperxPH[tk][j][ti] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC101TaperxPH[tk][j][ti] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC101TaperxPH[tk][j][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC101TaperxPH[tk][j][ti]>amax?ayC101TaperxPH[tk][j][ti]:amax;
amin = ayC101TaperxPH[tk][j][ti]<amin?ayC101TaperxPH[tk][j][ti]:amin;
bmax = byC101TaperxPH[tk][j][ti]>bmax?byC101TaperxPH[tk][j][ti]:bmax;
bmin = byC101TaperxPH[tk][j][ti]<bmin?byC101TaperxPH[tk][j][ti]:bmin;
      bzC101Taper[tk][j][ti] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC101Taper[tk][j][ti] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC101Taper[tk][j][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC101Taper[tk][j][ti]>amax?azC101Taper[tk][j][ti]:amax;
amin = azC101Taper[tk][j][ti]<amin?azC101Taper[tk][j][ti]:amin;
bmax = bzC101Taper[tk][j][ti]>bmax?bzC101Taper[tk][j][ti]:bmax;
bmin = bzC101Taper[tk][j][ti]<bmin?bzC101Taper[tk][j][ti]:bmin;
      bzC101TaperzPH[tk][j][ti] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC101TaperzPH[tk][j][ti] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC101TaperzPH[tk][j][ti]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC101TaperzPH[tk][j][ti]>amax?azC101TaperzPH[tk][j][ti]:amax;
amin = azC101TaperzPH[tk][j][ti]<amin?azC101TaperzPH[tk][j][ti]:amin;
bmax = bzC101TaperzPH[tk][j][ti]>bmax?bzC101TaperzPH[tk][j][ti]:bmax;
bmin = bzC101TaperzPH[tk][j][ti]<bmin?bzC101TaperzPH[tk][j][ti]:bmin;
      bzC101TaperzPHxPH[tk][j][ti] = bzTaperPH[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC101TaperzPHxPH[tk][j][ti] = (zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC101TaperzPHxPH[tk][j][ti]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC101TaperzPHxPH[tk][j][ti]>amax?azC101TaperzPHxPH[tk][j][ti]:amax;
amin = azC101TaperzPHxPH[tk][j][ti]<amin?azC101TaperzPHxPH[tk][j][ti]:amin;
bmax = bzC101TaperzPHxPH[tk][j][ti]>bmax?bzC101TaperzPHxPH[tk][j][ti]:bmax;
bmin = bzC101TaperzPHxPH[tk][j][ti]<bmin?bzC101TaperzPHxPH[tk][j][ti]:bmin;
      bzC101TaperzPHyPH[tk][j][ti] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC101TaperzPHyPH[tk][j][ti] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC101TaperzPHyPH[tk][j][ti]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC101TaperzPHyPH[tk][j][ti]>amax?azC101TaperzPHyPH[tk][j][ti]:amax;
amin = azC101TaperzPHyPH[tk][j][ti]<amin?azC101TaperzPHyPH[tk][j][ti]:amin;
bmax = bzC101TaperzPHyPH[tk][j][ti]>bmax?bzC101TaperzPHyPH[tk][j][ti]:bmax;
bmin = bzC101TaperzPHyPH[tk][j][ti]<bmin?bzC101TaperzPHyPH[tk][j][ti]:bmin;
      bzC101TaperyPH[tk][j][ti] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC101TaperyPH[tk][j][ti] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC101TaperyPH[tk][j][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC101TaperyPH[tk][j][ti]>amax?azC101TaperyPH[tk][j][ti]:amax;
amin = azC101TaperyPH[tk][j][ti]<amin?azC101TaperyPH[tk][j][ti]:amin;
bmax = bzC101TaperyPH[tk][j][ti]>bmax?bzC101TaperyPH[tk][j][ti]:bmax;
bmin = bzC101TaperyPH[tk][j][ti]<bmin?bzC101TaperyPH[tk][j][ti]:bmin;
      bzC101TaperyPHxPH[tk][j][ti] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC101TaperyPHxPH[tk][j][ti] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j])*(bzC101TaperyPHxPH[tk][j][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC101TaperyPHxPH[tk][j][ti]>amax?azC101TaperyPHxPH[tk][j][ti]:amax;
amin = azC101TaperyPHxPH[tk][j][ti]<amin?azC101TaperyPHxPH[tk][j][ti]:amin;
bmax = bzC101TaperyPHxPH[tk][j][ti]>bmax?bzC101TaperyPHxPH[tk][j][ti]:bmax;
bmin = bzC101TaperyPHxPH[tk][j][ti]<bmin?bzC101TaperyPHxPH[tk][j][ti]:bmin;
      bzC101TaperxPH[tk][j][ti] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC101TaperxPH[tk][j][ti] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC101TaperxPH[tk][j][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC101TaperxPH[tk][j][ti]>amax?azC101TaperxPH[tk][j][ti]:amax;
amin = azC101TaperxPH[tk][j][ti]<amin?azC101TaperxPH[tk][j][ti]:amin;
bmax = bzC101TaperxPH[tk][j][ti]>bmax?bzC101TaperxPH[tk][j][ti]:bmax;
bmin = bzC101TaperxPH[tk][j][ti]<bmin?bzC101TaperxPH[tk][j][ti]:bmin;
    }
  }
  for(int j=yEnd1;j<NY;++j) {
    int tj = j-yEnd1;
    for(int i=0;i<_nXmin;++i) {
      bxC011Taper[tk][tj][i] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC011Taper[tk][tj][i] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC011Taper[tk][tj][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC011Taper[tk][tj][i]>amax?axC011Taper[tk][tj][i]:amax;
amin = axC011Taper[tk][tj][i]<amin?axC011Taper[tk][tj][i]:amin;
bmax = bxC011Taper[tk][tj][i]>bmax?bxC011Taper[tk][tj][i]:bmax;
bmin = bxC011Taper[tk][tj][i]<bmin?bxC011Taper[tk][tj][i]:bmin;
      bxC011TaperzPH[tk][tj][i] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC011TaperzPH[tk][tj][i] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC011TaperzPH[tk][tj][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC011TaperzPH[tk][tj][i]>amax?axC011TaperzPH[tk][tj][i]:amax;
amin = axC011TaperzPH[tk][tj][i]<amin?axC011TaperzPH[tk][tj][i]:amin;
bmax = bxC011TaperzPH[tk][tj][i]>bmax?bxC011TaperzPH[tk][tj][i]:bmax;
bmin = bxC011TaperzPH[tk][tj][i]<bmin?bxC011TaperzPH[tk][tj][i]:bmin;
      bxC011TaperzPHxPH[tk][tj][i] = bxTaperPH[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC011TaperzPHxPH[tk][tj][i] = (xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC011TaperzPHxPH[tk][tj][i]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC011TaperzPHxPH[tk][tj][i]>amax?axC011TaperzPHxPH[tk][tj][i]:amax;
amin = axC011TaperzPHxPH[tk][tj][i]<amin?axC011TaperzPHxPH[tk][tj][i]:amin;
bmax = bxC011TaperzPHxPH[tk][tj][i]>bmax?bxC011TaperzPHxPH[tk][tj][i]:bmax;
bmin = bxC011TaperzPHxPH[tk][tj][i]<bmin?bxC011TaperzPHxPH[tk][tj][i]:bmin;
      bxC011TaperzPHyPH[tk][tj][i] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC011TaperzPHyPH[tk][tj][i] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j])*(bxC011TaperzPHyPH[tk][tj][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC011TaperzPHyPH[tk][tj][i]>amax?axC011TaperzPHyPH[tk][tj][i]:amax;
amin = axC011TaperzPHyPH[tk][tj][i]<amin?axC011TaperzPHyPH[tk][tj][i]:amin;
bmax = bxC011TaperzPHyPH[tk][tj][i]>bmax?bxC011TaperzPHyPH[tk][tj][i]:bmax;
bmin = bxC011TaperzPHyPH[tk][tj][i]<bmin?bxC011TaperzPHyPH[tk][tj][i]:bmin;
      bxC011TaperyPH[tk][tj][i] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC011TaperyPH[tk][tj][i] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC011TaperyPH[tk][tj][i]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC011TaperyPH[tk][tj][i]>amax?axC011TaperyPH[tk][tj][i]:amax;
amin = axC011TaperyPH[tk][tj][i]<amin?axC011TaperyPH[tk][tj][i]:amin;
bmax = bxC011TaperyPH[tk][tj][i]>bmax?bxC011TaperyPH[tk][tj][i]:bmax;
bmin = bxC011TaperyPH[tk][tj][i]<bmin?bxC011TaperyPH[tk][tj][i]:bmin;
      bxC011TaperyPHxPH[tk][tj][i] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC011TaperyPHxPH[tk][tj][i] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC011TaperyPHxPH[tk][tj][i]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC011TaperyPHxPH[tk][tj][i]>amax?axC011TaperyPHxPH[tk][tj][i]:amax;
amin = axC011TaperyPHxPH[tk][tj][i]<amin?axC011TaperyPHxPH[tk][tj][i]:amin;
bmax = bxC011TaperyPHxPH[tk][tj][i]>bmax?bxC011TaperyPHxPH[tk][tj][i]:bmax;
bmin = bxC011TaperyPHxPH[tk][tj][i]<bmin?bxC011TaperyPHxPH[tk][tj][i]:bmin;
      bxC011TaperxPH[tk][tj][i] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC011TaperxPH[tk][tj][i] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC011TaperxPH[tk][tj][i]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC011TaperxPH[tk][tj][i]>amax?axC011TaperxPH[tk][tj][i]:amax;
amin = axC011TaperxPH[tk][tj][i]<amin?axC011TaperxPH[tk][tj][i]:amin;
bmax = bxC011TaperxPH[tk][tj][i]>bmax?bxC011TaperxPH[tk][tj][i]:bmax;
bmin = bxC011TaperxPH[tk][tj][i]<bmin?bxC011TaperxPH[tk][tj][i]:bmin;
      byC011Taper[tk][tj][i] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC011Taper[tk][tj][i] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC011Taper[tk][tj][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC011Taper[tk][tj][i]>amax?ayC011Taper[tk][tj][i]:amax;
amin = ayC011Taper[tk][tj][i]<amin?ayC011Taper[tk][tj][i]:amin;
bmax = byC011Taper[tk][tj][i]>bmax?byC011Taper[tk][tj][i]:bmax;
bmin = byC011Taper[tk][tj][i]<bmin?byC011Taper[tk][tj][i]:bmin;
      byC011TaperzPH[tk][tj][i] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC011TaperzPH[tk][tj][i] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC011TaperzPH[tk][tj][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC011TaperzPH[tk][tj][i]>amax?ayC011TaperzPH[tk][tj][i]:amax;
amin = ayC011TaperzPH[tk][tj][i]<amin?ayC011TaperzPH[tk][tj][i]:amin;
bmax = byC011TaperzPH[tk][tj][i]>bmax?byC011TaperzPH[tk][tj][i]:bmax;
bmin = byC011TaperzPH[tk][tj][i]<bmin?byC011TaperzPH[tk][tj][i]:bmin;
      byC011TaperzPHxPH[tk][tj][i] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC011TaperzPHxPH[tk][tj][i] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i])*(byC011TaperzPHxPH[tk][tj][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC011TaperzPHxPH[tk][tj][i]>amax?ayC011TaperzPHxPH[tk][tj][i]:amax;
amin = ayC011TaperzPHxPH[tk][tj][i]<amin?ayC011TaperzPHxPH[tk][tj][i]:amin;
bmax = byC011TaperzPHxPH[tk][tj][i]>bmax?byC011TaperzPHxPH[tk][tj][i]:bmax;
bmin = byC011TaperzPHxPH[tk][tj][i]<bmin?byC011TaperzPHxPH[tk][tj][i]:bmin;
      byC011TaperzPHyPH[tk][tj][i] = byTaperPH[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC011TaperzPHyPH[tk][tj][i] = (yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC011TaperzPHyPH[tk][tj][i]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC011TaperzPHyPH[tk][tj][i]>amax?ayC011TaperzPHyPH[tk][tj][i]:amax;
amin = ayC011TaperzPHyPH[tk][tj][i]<amin?ayC011TaperzPHyPH[tk][tj][i]:amin;
bmax = byC011TaperzPHyPH[tk][tj][i]>bmax?byC011TaperzPHyPH[tk][tj][i]:bmax;
bmin = byC011TaperzPHyPH[tk][tj][i]<bmin?byC011TaperzPHyPH[tk][tj][i]:bmin;
      byC011TaperyPH[tk][tj][i] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC011TaperyPH[tk][tj][i] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC011TaperyPH[tk][tj][i]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC011TaperyPH[tk][tj][i]>amax?ayC011TaperyPH[tk][tj][i]:amax;
amin = ayC011TaperyPH[tk][tj][i]<amin?ayC011TaperyPH[tk][tj][i]:amin;
bmax = byC011TaperyPH[tk][tj][i]>bmax?byC011TaperyPH[tk][tj][i]:bmax;
bmin = byC011TaperyPH[tk][tj][i]<bmin?byC011TaperyPH[tk][tj][i]:bmin;
      byC011TaperyPHxPH[tk][tj][i] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC011TaperyPHxPH[tk][tj][i] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC011TaperyPHxPH[tk][tj][i]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC011TaperyPHxPH[tk][tj][i]>amax?ayC011TaperyPHxPH[tk][tj][i]:amax;
amin = ayC011TaperyPHxPH[tk][tj][i]<amin?ayC011TaperyPHxPH[tk][tj][i]:amin;
bmax = byC011TaperyPHxPH[tk][tj][i]>bmax?byC011TaperyPHxPH[tk][tj][i]:bmax;
bmin = byC011TaperyPHxPH[tk][tj][i]<bmin?byC011TaperyPHxPH[tk][tj][i]:bmin;
      byC011TaperxPH[tk][tj][i] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC011TaperxPH[tk][tj][i] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC011TaperxPH[tk][tj][i]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC011TaperxPH[tk][tj][i]>amax?ayC011TaperxPH[tk][tj][i]:amax;
amin = ayC011TaperxPH[tk][tj][i]<amin?ayC011TaperxPH[tk][tj][i]:amin;
bmax = byC011TaperxPH[tk][tj][i]>bmax?byC011TaperxPH[tk][tj][i]:bmax;
bmin = byC011TaperxPH[tk][tj][i]<bmin?byC011TaperxPH[tk][tj][i]:bmin;
      bzC011Taper[tk][tj][i] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC011Taper[tk][tj][i] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC011Taper[tk][tj][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC011Taper[tk][tj][i]>amax?azC011Taper[tk][tj][i]:amax;
amin = azC011Taper[tk][tj][i]<amin?azC011Taper[tk][tj][i]:amin;
bmax = bzC011Taper[tk][tj][i]>bmax?bzC011Taper[tk][tj][i]:bmax;
bmin = bzC011Taper[tk][tj][i]<bmin?bzC011Taper[tk][tj][i]:bmin;
      bzC011TaperzPH[tk][tj][i] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC011TaperzPH[tk][tj][i] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC011TaperzPH[tk][tj][i]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC011TaperzPH[tk][tj][i]>amax?azC011TaperzPH[tk][tj][i]:amax;
amin = azC011TaperzPH[tk][tj][i]<amin?azC011TaperzPH[tk][tj][i]:amin;
bmax = bzC011TaperzPH[tk][tj][i]>bmax?bzC011TaperzPH[tk][tj][i]:bmax;
bmin = bzC011TaperzPH[tk][tj][i]<bmin?bzC011TaperzPH[tk][tj][i]:bmin;
      bzC011TaperzPHxPH[tk][tj][i] = bzTaperPH[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC011TaperzPHxPH[tk][tj][i] = (zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC011TaperzPHxPH[tk][tj][i]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC011TaperzPHxPH[tk][tj][i]>amax?azC011TaperzPHxPH[tk][tj][i]:amax;
amin = azC011TaperzPHxPH[tk][tj][i]<amin?azC011TaperzPHxPH[tk][tj][i]:amin;
bmax = bzC011TaperzPHxPH[tk][tj][i]>bmax?bzC011TaperzPHxPH[tk][tj][i]:bmax;
bmin = bzC011TaperzPHxPH[tk][tj][i]<bmin?bzC011TaperzPHxPH[tk][tj][i]:bmin;
      bzC011TaperzPHyPH[tk][tj][i] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC011TaperzPHyPH[tk][tj][i] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC011TaperzPHyPH[tk][tj][i]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC011TaperzPHyPH[tk][tj][i]>amax?azC011TaperzPHyPH[tk][tj][i]:amax;
amin = azC011TaperzPHyPH[tk][tj][i]<amin?azC011TaperzPHyPH[tk][tj][i]:amin;
bmax = bzC011TaperzPHyPH[tk][tj][i]>bmax?bzC011TaperzPHyPH[tk][tj][i]:bmax;
bmin = bzC011TaperzPHyPH[tk][tj][i]<bmin?bzC011TaperzPHyPH[tk][tj][i]:bmin;
      bzC011TaperyPH[tk][tj][i] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC011TaperyPH[tk][tj][i] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC011TaperyPH[tk][tj][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC011TaperyPH[tk][tj][i]>amax?azC011TaperyPH[tk][tj][i]:amax;
amin = azC011TaperyPH[tk][tj][i]<amin?azC011TaperyPH[tk][tj][i]:amin;
bmax = bzC011TaperyPH[tk][tj][i]>bmax?bzC011TaperyPH[tk][tj][i]:bmax;
bmin = bzC011TaperyPH[tk][tj][i]<bmin?bzC011TaperyPH[tk][tj][i]:bmin;
      bzC011TaperyPHxPH[tk][tj][i] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC011TaperyPHxPH[tk][tj][i] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j])*(bzC011TaperyPHxPH[tk][tj][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC011TaperyPHxPH[tk][tj][i]>amax?azC011TaperyPHxPH[tk][tj][i]:amax;
amin = azC011TaperyPHxPH[tk][tj][i]<amin?azC011TaperyPHxPH[tk][tj][i]:amin;
bmax = bzC011TaperyPHxPH[tk][tj][i]>bmax?bzC011TaperyPHxPH[tk][tj][i]:bmax;
bmin = bzC011TaperyPHxPH[tk][tj][i]<bmin?bzC011TaperyPHxPH[tk][tj][i]:bmin;
      bzC011TaperxPH[tk][tj][i] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC011TaperxPH[tk][tj][i] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC011TaperxPH[tk][tj][i]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC011TaperxPH[tk][tj][i]>amax?azC011TaperxPH[tk][tj][i]:amax;
amin = azC011TaperxPH[tk][tj][i]<amin?azC011TaperxPH[tk][tj][i]:amin;
bmax = bzC011TaperxPH[tk][tj][i]>bmax?bzC011TaperxPH[tk][tj][i]:bmax;
bmin = bzC011TaperxPH[tk][tj][i]<bmin?bzC011TaperxPH[tk][tj][i]:bmin;
    }
    for(int i=xEnd1;i<NX;++i) {
      int ti = i-xEnd1;
      bxC111Taper[tk][tj][ti] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC111Taper[tk][tj][ti] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC111Taper[tk][tj][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC111Taper[tk][tj][ti]>amax?axC111Taper[tk][tj][ti]:amax;
amin = axC111Taper[tk][tj][ti]<amin?axC111Taper[tk][tj][ti]:amin;
bmax = bxC111Taper[tk][tj][ti]>bmax?bxC111Taper[tk][tj][ti]:bmax;
bmin = bxC111Taper[tk][tj][ti]<bmin?bxC111Taper[tk][tj][ti]:bmin;
      bxC111TaperzPH[tk][tj][ti] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC111TaperzPH[tk][tj][ti] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC111TaperzPH[tk][tj][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC111TaperzPH[tk][tj][ti]>amax?axC111TaperzPH[tk][tj][ti]:amax;
amin = axC111TaperzPH[tk][tj][ti]<amin?axC111TaperzPH[tk][tj][ti]:amin;
bmax = bxC111TaperzPH[tk][tj][ti]>bmax?bxC111TaperzPH[tk][tj][ti]:bmax;
bmin = bxC111TaperzPH[tk][tj][ti]<bmin?bxC111TaperzPH[tk][tj][ti]:bmin;
      bxC111TaperzPHxPH[tk][tj][ti] = bxTaperPH[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC111TaperzPHxPH[tk][tj][ti] = (xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j])*(bxC111TaperzPHxPH[tk][tj][ti]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaperPH[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC111TaperzPHxPH[tk][tj][ti]>amax?axC111TaperzPHxPH[tk][tj][ti]:amax;
amin = axC111TaperzPHxPH[tk][tj][ti]<amin?axC111TaperzPHxPH[tk][tj][ti]:amin;
bmax = bxC111TaperzPHxPH[tk][tj][ti]>bmax?bxC111TaperzPHxPH[tk][tj][ti]:bmax;
bmin = bxC111TaperzPHxPH[tk][tj][ti]<bmin?bxC111TaperzPHxPH[tk][tj][ti]:bmin;
      bxC111TaperzPHyPH[tk][tj][ti] = bxTaper[i]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC111TaperzPHyPH[tk][tj][ti] = (xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j])*(bxC111TaperzPHyPH[tk][tj][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaperPH[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC111TaperzPHyPH[tk][tj][ti]>amax?axC111TaperzPHyPH[tk][tj][ti]:amax;
amin = axC111TaperzPHyPH[tk][tj][ti]<amin?axC111TaperzPHyPH[tk][tj][ti]:amin;
bmax = bxC111TaperzPHyPH[tk][tj][ti]>bmax?bxC111TaperzPHyPH[tk][tj][ti]:bmax;
bmin = bxC111TaperzPHyPH[tk][tj][ti]<bmin?bxC111TaperzPHyPH[tk][tj][ti]:bmin;
      bxC111TaperyPH[tk][tj][ti] = bxTaper[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC111TaperyPH[tk][tj][ti] = (xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC111TaperyPH[tk][tj][ti]-1.)/(kxTaper[i]*(xTaper[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaper[i]*alphaXTaper[i]));
amax = axC111TaperyPH[tk][tj][ti]>amax?axC111TaperyPH[tk][tj][ti]:amax;
amin = axC111TaperyPH[tk][tj][ti]<amin?axC111TaperyPH[tk][tj][ti]:amin;
bmax = bxC111TaperyPH[tk][tj][ti]>bmax?bxC111TaperyPH[tk][tj][ti]:bmax;
bmin = bxC111TaperyPH[tk][tj][ti]<bmin?bxC111TaperyPH[tk][tj][ti]:bmin;
      bxC111TaperyPHxPH[tk][tj][ti] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaperPH[j]*dt);
      axC111TaperyPHxPH[tk][tj][ti] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j])*(bxC111TaperyPHxPH[tk][tj][ti]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaperPH[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC111TaperyPHxPH[tk][tj][ti]>amax?axC111TaperyPHxPH[tk][tj][ti]:amax;
amin = axC111TaperyPHxPH[tk][tj][ti]<amin?axC111TaperyPHxPH[tk][tj][ti]:amin;
bmax = bxC111TaperyPHxPH[tk][tj][ti]>bmax?bxC111TaperyPHxPH[tk][tj][ti]:bmax;
bmin = bxC111TaperyPHxPH[tk][tj][ti]<bmin?bxC111TaperyPHxPH[tk][tj][ti]:bmin;
      bxC111TaperxPH[tk][tj][ti] = bxTaperPH[i]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*yTaper[j]*dt);
      axC111TaperxPH[tk][tj][ti] = (xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j])*(bxC111TaperxPH[tk][tj][ti]-1.)/(kxTaperPH[i]*(xTaperPH[i]+xtap*zTaper[k]+xtap*yTaper[j]+kxTaperPH[i]*alphaXTaperPH[i]));
amax = axC111TaperxPH[tk][tj][ti]>amax?axC111TaperxPH[tk][tj][ti]:amax;
amin = axC111TaperxPH[tk][tj][ti]<amin?axC111TaperxPH[tk][tj][ti]:amin;
bmax = bxC111TaperxPH[tk][tj][ti]>bmax?bxC111TaperxPH[tk][tj][ti]:bmax;
bmin = bxC111TaperxPH[tk][tj][ti]<bmin?bxC111TaperxPH[tk][tj][ti]:bmin;
      byC111Taper[tk][tj][ti] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC111Taper[tk][tj][ti] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC111Taper[tk][tj][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC111Taper[tk][tj][ti]>amax?ayC111Taper[tk][tj][ti]:amax;
amin = ayC111Taper[tk][tj][ti]<amin?ayC111Taper[tk][tj][ti]:amin;
bmax = byC111Taper[tk][tj][ti]>bmax?byC111Taper[tk][tj][ti]:bmax;
bmin = byC111Taper[tk][tj][ti]<bmin?byC111Taper[tk][tj][ti]:bmin;
      byC111TaperzPH[tk][tj][ti] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC111TaperzPH[tk][tj][ti] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC111TaperzPH[tk][tj][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC111TaperzPH[tk][tj][ti]>amax?ayC111TaperzPH[tk][tj][ti]:amax;
amin = ayC111TaperzPH[tk][tj][ti]<amin?ayC111TaperzPH[tk][tj][ti]:amin;
bmax = byC111TaperzPH[tk][tj][ti]>bmax?byC111TaperzPH[tk][tj][ti]:bmax;
bmin = byC111TaperzPH[tk][tj][ti]<bmin?byC111TaperzPH[tk][tj][ti]:bmin;
      byC111TaperzPHxPH[tk][tj][ti] = byTaper[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC111TaperzPHxPH[tk][tj][ti] = (yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i])*(byC111TaperzPHxPH[tk][tj][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaperPH[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC111TaperzPHxPH[tk][tj][ti]>amax?ayC111TaperzPHxPH[tk][tj][ti]:amax;
amin = ayC111TaperzPHxPH[tk][tj][ti]<amin?ayC111TaperzPHxPH[tk][tj][ti]:amin;
bmax = byC111TaperzPHxPH[tk][tj][ti]>bmax?byC111TaperzPHxPH[tk][tj][ti]:bmax;
bmin = byC111TaperzPHxPH[tk][tj][ti]<bmin?byC111TaperzPHxPH[tk][tj][ti]:bmin;
      byC111TaperzPHyPH[tk][tj][ti] = byTaperPH[j]*exp(-xtap*zTaperPH[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC111TaperzPHyPH[tk][tj][ti] = (yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i])*(byC111TaperzPHyPH[tk][tj][ti]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaperPH[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC111TaperzPHyPH[tk][tj][ti]>amax?ayC111TaperzPHyPH[tk][tj][ti]:amax;
amin = ayC111TaperzPHyPH[tk][tj][ti]<amin?ayC111TaperzPHyPH[tk][tj][ti]:amin;
bmax = byC111TaperzPHyPH[tk][tj][ti]>bmax?byC111TaperzPHyPH[tk][tj][ti]:bmax;
bmin = byC111TaperzPHyPH[tk][tj][ti]<bmin?byC111TaperzPHyPH[tk][tj][ti]:bmin;
      byC111TaperyPH[tk][tj][ti] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaper[i]*dt);
      ayC111TaperyPH[tk][tj][ti] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i])*(byC111TaperyPH[tk][tj][ti]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaper[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC111TaperyPH[tk][tj][ti]>amax?ayC111TaperyPH[tk][tj][ti]:amax;
amin = ayC111TaperyPH[tk][tj][ti]<amin?ayC111TaperyPH[tk][tj][ti]:amin;
bmax = byC111TaperyPH[tk][tj][ti]>bmax?byC111TaperyPH[tk][tj][ti]:bmax;
bmin = byC111TaperyPH[tk][tj][ti]<bmin?byC111TaperyPH[tk][tj][ti]:bmin;
      byC111TaperyPHxPH[tk][tj][ti] = byTaperPH[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC111TaperyPHxPH[tk][tj][ti] = (yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC111TaperyPHxPH[tk][tj][ti]-1.)/(kyTaperPH[j]*(yTaperPH[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaperPH[j]*alphaYTaperPH[j]));
amax = ayC111TaperyPHxPH[tk][tj][ti]>amax?ayC111TaperyPHxPH[tk][tj][ti]:amax;
amin = ayC111TaperyPHxPH[tk][tj][ti]<amin?ayC111TaperyPHxPH[tk][tj][ti]:amin;
bmax = byC111TaperyPHxPH[tk][tj][ti]>bmax?byC111TaperyPHxPH[tk][tj][ti]:bmax;
bmin = byC111TaperyPHxPH[tk][tj][ti]<bmin?byC111TaperyPHxPH[tk][tj][ti]:bmin;
      byC111TaperxPH[tk][tj][ti] = byTaper[j]*exp(-xtap*zTaper[k]*dt)*exp(-xtap*xTaperPH[i]*dt);
      ayC111TaperxPH[tk][tj][ti] = (yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i])*(byC111TaperxPH[tk][tj][ti]-1.)/(kyTaper[j]*(yTaper[j]+xtap*zTaper[k]+xtap*xTaperPH[i]+kyTaper[j]*alphaYTaper[j]));
amax = ayC111TaperxPH[tk][tj][ti]>amax?ayC111TaperxPH[tk][tj][ti]:amax;
amin = ayC111TaperxPH[tk][tj][ti]<amin?ayC111TaperxPH[tk][tj][ti]:amin;
bmax = byC111TaperxPH[tk][tj][ti]>bmax?byC111TaperxPH[tk][tj][ti]:bmax;
bmin = byC111TaperxPH[tk][tj][ti]<bmin?byC111TaperxPH[tk][tj][ti]:bmin;
      bzC111Taper[tk][tj][ti] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC111Taper[tk][tj][ti] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC111Taper[tk][tj][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC111Taper[tk][tj][ti]>amax?azC111Taper[tk][tj][ti]:amax;
amin = azC111Taper[tk][tj][ti]<amin?azC111Taper[tk][tj][ti]:amin;
bmax = bzC111Taper[tk][tj][ti]>bmax?bzC111Taper[tk][tj][ti]:bmax;
bmin = bzC111Taper[tk][tj][ti]<bmin?bzC111Taper[tk][tj][ti]:bmin;
      bzC111TaperzPH[tk][tj][ti] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC111TaperzPH[tk][tj][ti] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j])*(bzC111TaperzPH[tk][tj][ti]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC111TaperzPH[tk][tj][ti]>amax?azC111TaperzPH[tk][tj][ti]:amax;
amin = azC111TaperzPH[tk][tj][ti]<amin?azC111TaperzPH[tk][tj][ti]:amin;
bmax = bzC111TaperzPH[tk][tj][ti]>bmax?bzC111TaperzPH[tk][tj][ti]:bmax;
bmin = bzC111TaperzPH[tk][tj][ti]<bmin?bzC111TaperzPH[tk][tj][ti]:bmin;
      bzC111TaperzPHxPH[tk][tj][ti] = bzTaperPH[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC111TaperzPHxPH[tk][tj][ti] = (zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC111TaperzPHxPH[tk][tj][ti]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC111TaperzPHxPH[tk][tj][ti]>amax?azC111TaperzPHxPH[tk][tj][ti]:amax;
amin = azC111TaperzPHxPH[tk][tj][ti]<amin?azC111TaperzPHxPH[tk][tj][ti]:amin;
bmax = bzC111TaperzPHxPH[tk][tj][ti]>bmax?bzC111TaperzPHxPH[tk][tj][ti]:bmax;
bmin = bzC111TaperzPHxPH[tk][tj][ti]<bmin?bzC111TaperzPHxPH[tk][tj][ti]:bmin;
      bzC111TaperzPHyPH[tk][tj][ti] = bzTaperPH[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC111TaperzPHyPH[tk][tj][ti] = (zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC111TaperzPHyPH[tk][tj][ti]-1.)/(kzTaperPH[k]*(zTaperPH[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaperPH[k]*alphaZTaperPH[k]));
amax = azC111TaperzPHyPH[tk][tj][ti]>amax?azC111TaperzPHyPH[tk][tj][ti]:amax;
amin = azC111TaperzPHyPH[tk][tj][ti]<amin?azC111TaperzPHyPH[tk][tj][ti]:amin;
bmax = bzC111TaperzPHyPH[tk][tj][ti]>bmax?bzC111TaperzPHyPH[tk][tj][ti]:bmax;
bmin = bzC111TaperzPHyPH[tk][tj][ti]<bmin?bzC111TaperzPHyPH[tk][tj][ti]:bmin;
      bzC111TaperyPH[tk][tj][ti] = bzTaper[k]*exp(-xtap*xTaper[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC111TaperyPH[tk][tj][ti] = (zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j])*(bzC111TaperyPH[tk][tj][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaper[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC111TaperyPH[tk][tj][ti]>amax?azC111TaperyPH[tk][tj][ti]:amax;
amin = azC111TaperyPH[tk][tj][ti]<amin?azC111TaperyPH[tk][tj][ti]:amin;
bmax = bzC111TaperyPH[tk][tj][ti]>bmax?bzC111TaperyPH[tk][tj][ti]:bmax;
bmin = bzC111TaperyPH[tk][tj][ti]<bmin?bzC111TaperyPH[tk][tj][ti]:bmin;
      bzC111TaperyPHxPH[tk][tj][ti] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaperPH[j]*dt);
      azC111TaperyPHxPH[tk][tj][ti] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j])*(bzC111TaperyPHxPH[tk][tj][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaperPH[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC111TaperyPHxPH[tk][tj][ti]>amax?azC111TaperyPHxPH[tk][tj][ti]:amax;
amin = azC111TaperyPHxPH[tk][tj][ti]<amin?azC111TaperyPHxPH[tk][tj][ti]:amin;
bmax = bzC111TaperyPHxPH[tk][tj][ti]>bmax?bzC111TaperyPHxPH[tk][tj][ti]:bmax;
bmin = bzC111TaperyPHxPH[tk][tj][ti]<bmin?bzC111TaperyPHxPH[tk][tj][ti]:bmin;
      bzC111TaperxPH[tk][tj][ti] = bzTaper[k]*exp(-xtap*xTaperPH[i]*dt)*exp(-xtap*yTaper[j]*dt);
      azC111TaperxPH[tk][tj][ti] = (zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j])*(bzC111TaperxPH[tk][tj][ti]-1.)/(kzTaper[k]*(zTaper[k]+xtap*xTaperPH[i]+xtap*yTaper[j]+kzTaper[k]*alphaZTaper[k]));
amax = azC111TaperxPH[tk][tj][ti]>amax?azC111TaperxPH[tk][tj][ti]:amax;
amin = azC111TaperxPH[tk][tj][ti]<amin?azC111TaperxPH[tk][tj][ti]:amin;
bmax = bzC111TaperxPH[tk][tj][ti]>bmax?bzC111TaperxPH[tk][tj][ti]:bmax;
bmin = bzC111TaperxPH[tk][tj][ti]<bmin?bzC111TaperxPH[tk][tj][ti]:bmin;
    }
  }
}

}

//Here are the updating equations for Vx.  These are done in exactly
//the order that vx is stored in memory starting with the minZ, minY
//minX corner proceeding to maxZ, maxY, and maxX with x fastest, then
//y.  The variables minXb, fminXb, etc are used to define where this
//domain is situated relative to the MPML zones and the absolute edges
//of the model.  minXb, maxXb, etc. will be true if the MPML zone extends
//into this domain on the specified side.  fminXb, fmaxXb, etc. will only
//be true if this domain is actually on an absolute edge of the model.
//These two sets will the identical except in the cases where the MPML
//zone is wider than the width of the domain on the absolute edge.
//4th order accurate differences are used when they can be, but near the
//flanks, second order differences are used and finally, in some cases,
//second order with the assumption of zero valued variables off grid
//are used right on the flanks.
void updateVxAttenMPMLBounds(modelDefStruct* modelDef,
                        float* __restrict__ vx,
                        float* __restrict__ pressure,
                        float cx[2],float cy[2],float cz[2],
                        float* __restrict__ rho,unsigned char* __restrict__ vxfunc,
                        bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                        bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb) {
  using namespace acousticAtten_MPML;
  DEF_MODEL_SIZE(modelDef);
  
  //  minXb = maxXb = minYb = maxYb = minZb = maxZb = false;
  //  maxZb=false;
  
  //These variables are primarily variables of convenience.
  register int i,j,k;
  int zLim=NZ-2,yLim=NY-2,xLim=NX-2;
  float cx0=cx[0],cx1=cx[1];
  float cy0=cy[0],cy1=cy[1];
  float cz0=cz[0],cz1=cz[1];
  float dcx0=modelDef->dt/modelDef->dx*modelDef->scalarSpeed;
  float dcy0=modelDef->dt/modelDef->dy*modelDef->scalarSpeed;
  float dcz0=modelDef->dt/modelDef->dz*modelDef->scalarSpeed;
  int nxymin = NX*_nYmin, nxymax = NX*_nYmax, nxminy = _nXmin*(_yMaxStart-_nYmin), nxmaxy=_nXmax*(_yMaxStart-_nYmin);
  int mxStart1 = MAX(_nXmin,1);
  int mxEnd1 = MIN(_xMaxStart-1,NX-2);
  
  //printf("vx %d %d %d, %d %d %d %d %d %d\n",NX,NY,NZ,minXb,maxXb,minYb,maxYb,minZb,maxZb);
  //The following variables keep track of the starting and ending indicies
  //for the memory variables and will change depending on which memory
  //variable we are using
  int kk=0, jj=0, ii=0;
  int nx = NX, nxy=NXY;
  //Z min MPML zone
  if(minZb) {
    if(fminZb) {  //if this is the absolute edge of the model
      k=0;
      for(j=2;j<yStart2;j++){
        for(i=1;i<_nXmin;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC000TaperxPH[k][j][i],ayC000TaperxPH[k][j][i],azC000TaperxPH[k][j][i],bxC000TaperxPH[k][j][i],byC000TaperxPH[k][j][i],bzC000TaperxPH[k][j][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axE011Taper[k][j],ayE011Taper[k][j],azE011Taper[k][j],bxE011Taper[k][j],byE011Taper[k][j],bzE011Taper[k][j],i,j,k);
        }
        for(i=mxEnd1;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC100TaperxPH[k][j][ti],ayC100TaperxPH[k][j][ti],azC100TaperxPH[k][j][ti],bxC100TaperxPH[k][j][ti],byC100TaperxPH[k][j][ti],bzC100TaperxPH[k][j][ti],i,j,k);
        }
      }
      for(j=yStart2;j<yEnd;j++){
        for(i=1;i<_nXmin;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axE101TaperxPH[k][i],ayE101TaperxPH[k][i],azE101TaperxPH[k][i],bxE101TaperxPH[k][i],byE101TaperxPH[k][i],bzE101TaperxPH[k][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axZ0Taper[k],ayZ0Taper[k],azTaper[k],bxZ0Taper[k],byZ0Taper[k],bzTaper[k],i,j,k);
        }
        for(i=mxEnd1;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axE201TaperxPH[k][ti],ayE201TaperxPH[k][ti],azE201TaperxPH[k][ti],bxE201TaperxPH[k][ti],byE201TaperxPH[k][ti],bzE201TaperxPH[k][ti],i,j,k);
        }
      }
      for(j=yEnd;j<yLim;j++){
        int tj = j-yEnd1;
        for(i=1;i<_nXmin;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC010TaperxPH[k][tj][i],ayC010TaperxPH[k][tj][i],azC010TaperxPH[k][tj][i],bxC010TaperxPH[k][tj][i],byC010TaperxPH[k][tj][i],bzC010TaperxPH[k][tj][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axE021Taper[k][tj],ayE021Taper[k][tj],azE021Taper[k][tj],bxE021Taper[k][tj],byE021Taper[k][tj],bzE021Taper[k][tj],i,j,k);
        }
        for(i=mxEnd1;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC110TaperxPH[k][tj][ti],ayC110TaperxPH[k][tj][ti],azC110TaperxPH[k][tj][ti],bxC110TaperxPH[k][tj][ti],byC110TaperxPH[k][tj][ti],bzC110TaperxPH[k][tj][ti],i,j,k);
        }
      }
      k=1;
      for(j=2;j<yStart2;j++){
        for(i=1;i<_nXmin;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC000TaperxPH[k][j][i],ayC000TaperxPH[k][j][i],azC000TaperxPH[k][j][i],bxC000TaperxPH[k][j][i],byC000TaperxPH[k][j][i],bzC000TaperxPH[k][j][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axE011Taper[k][j],ayE011Taper[k][j],azE011Taper[k][j],bxE011Taper[k][j],byE011Taper[k][j],bzE011Taper[k][j],i,j,k);
        }
        for(i=mxEnd1;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC100TaperxPH[k][j][ti],ayC100TaperxPH[k][j][ti],azC100TaperxPH[k][j][ti],bxC100TaperxPH[k][j][ti],byC100TaperxPH[k][j][ti],bzC100TaperxPH[k][j][ti],i,j,k);
        }
      }
      for(j=yStart2;j<yEnd;j++){
        for(i=1;i<_nXmin;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axE101TaperxPH[k][i],ayE101TaperxPH[k][i],azE101TaperxPH[k][i],bxE101TaperxPH[k][i],byE101TaperxPH[k][i],bzE101TaperxPH[k][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axZ0Taper[k],ayZ0Taper[k],azTaper[k],bxZ0Taper[k],byZ0Taper[k],bzTaper[k],i,j,k);
        }
        for(i=mxEnd1;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axE201TaperxPH[k][ti],ayE201TaperxPH[k][ti],azE201TaperxPH[k][ti],bxE201TaperxPH[k][ti],byE201TaperxPH[k][ti],bzE201TaperxPH[k][ti],i,j,k);
        }
      }
      for(j=yEnd;j<yLim;j++){
        int tj = j-yEnd1;
        for(i=1;i<_nXmin;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC010TaperxPH[k][tj][i],ayC010TaperxPH[k][tj][i],azC010TaperxPH[k][tj][i],bxC010TaperxPH[k][tj][i],byC010TaperxPH[k][tj][i],bzC010TaperxPH[k][tj][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axE021Taper[k][tj],ayE021Taper[k][tj],azE021Taper[k][tj],bxE021Taper[k][tj],byE021Taper[k][tj],bzE021Taper[k][tj],i,j,k);
        }
        for(i=mxEnd1;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC110TaperxPH[k][tj][ti],ayC110TaperxPH[k][tj][ti],azC110TaperxPH[k][tj][ti],bxC110TaperxPH[k][tj][ti],byC110TaperxPH[k][tj][ti],bzC110TaperxPH[k][tj][ti],i,j,k);
        }
      }
    }
    for(k=2;k<zStart;k++) {
      int pIndex;
      if(fminYb) {
        j=0;
        for(i=1;i<_nXmin;++i) {
          int index=i+j*NX+k*NXY;
          pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC000TaperxPH[k][j][i],ayC000TaperxPH[k][j][i],azC000TaperxPH[k][j][i],bxC000TaperxPH[k][j][i],byC000TaperxPH[k][j][i],bzC000TaperxPH[k][j][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;++i) {
          int index=i+j*NX+k*NXY;
          pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axE011Taper[k][j],ayE011Taper[k][j],azE011Taper[k][j],bxE011Taper[k][j],byE011Taper[k][j],bzE011Taper[k][j],i,j,k);
        }
        for(i=mxEnd1;i<xLim;++i) {
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC100TaperxPH[k][j][ti],ayC100TaperxPH[k][j][ti],azC100TaperxPH[k][j][ti],bxC100TaperxPH[k][j][ti],byC100TaperxPH[k][j][ti],bzC100TaperxPH[k][j][ti],i,j,k);
        }
        j=1;
        for(i=1;i<_nXmin;++i) {
          int index=i+j*NX+k*NXY;
          pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC000TaperxPH[k][j][i],ayC000TaperxPH[k][j][i],azC000TaperxPH[k][j][i],bxC000TaperxPH[k][j][i],byC000TaperxPH[k][j][i],bzC000TaperxPH[k][j][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;++i) {
          int index=i+j*NX+k*NXY;
          pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axE011Taper[k][j],ayE011Taper[k][j],azE011Taper[k][j],bxE011Taper[k][j],byE011Taper[k][j],bzE011Taper[k][j],i,j,k);
        }
        for(i=mxEnd1;i<xLim;++i) {
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC100TaperxPH[k][j][ti],ayC100TaperxPH[k][j][ti],azC100TaperxPH[k][j][ti],bxC100TaperxPH[k][j][ti],byC100TaperxPH[k][j][ti],bzC100TaperxPH[k][j][ti],i,j,k);
        }
      }
      for(j=2;j<yStart2;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC000TaperxPH[k][j][i],ayC000TaperxPH[k][j][i],azC000TaperxPH[k][j][i],bxC000TaperxPH[k][j][i],byC000TaperxPH[k][j][i],bzC000TaperxPH[k][j][i],i,j,k);
        }
        for(i=1;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,index,pxZmin,vx,axC000TaperxPH[k][j][i],ayC000TaperxPH[k][j][i],azC000TaperxPH[k][j][i],bxC000TaperxPH[k][j][i],byC000TaperxPH[k][j][i],bzC000TaperxPH[k][j][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,index,pxZmin,vx,axE011Taper[k][j],ayE011Taper[k][j],azE011Taper[k][j],bxE011Taper[k][j],byE011Taper[k][j],bzE011Taper[k][j],i,j,k);
        }
        for(i=mxEnd1;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,index,pxZmin,vx,axC100TaperxPH[k][j][ti],ayC100TaperxPH[k][j][ti],azC100TaperxPH[k][j][ti],bxC100TaperxPH[k][j][ti],byC100TaperxPH[k][j][ti],bzC100TaperxPH[k][j][ti],i,j,k);
        }
        if (fmaxXb) {
          i=xLim;
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC100TaperxPH[k][j][ti],ayC100TaperxPH[k][j][ti],azC100TaperxPH[k][j][ti],bxC100TaperxPH[k][j][ti],byC100TaperxPH[k][j][ti],bzC100TaperxPH[k][j][ti],i,j,k);
        }
      }
      for(j=yStart2;j<yEnd;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axE101TaperxPH[k][i],ayE101TaperxPH[k][i],azE101TaperxPH[k][i],bxE101TaperxPH[k][i],byE101TaperxPH[k][i],bzE101TaperxPH[k][i],i,j,k);
        }
        for(i=1;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,index,pxZmin,vx,axE101TaperxPH[k][i],ayE101TaperxPH[k][i],azE101TaperxPH[k][i],bxE101TaperxPH[k][i],byE101TaperxPH[k][i],bzE101TaperxPH[k][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,index,pxZmin,vx,axZ0Taper[k],ayZ0Taper[k],azTaper[k],bxZ0Taper[k],byZ0Taper[k],bzTaper[k],i,j,k);
        }
        for(i=mxEnd1;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,index,pxZmin,vx,axE201TaperxPH[k][ti],ayE201TaperxPH[k][ti],azE201TaperxPH[k][ti],bxE201TaperxPH[k][ti],byE201TaperxPH[k][ti],bzE201TaperxPH[k][ti],i,j,k);
        }
        if (fmaxXb) {
          i=xLim;
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axE201TaperxPH[k][ti],ayE201TaperxPH[k][ti],azE201TaperxPH[k][ti],bxE201TaperxPH[k][ti],byE201TaperxPH[k][ti],bzE201TaperxPH[k][ti],i,j,k);
        }
      }
      for(j=yEnd;j<yLim;++j) {
        int tj = j-yEnd1;
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC010TaperxPH[k][tj][i],ayC010TaperxPH[k][tj][i],azC010TaperxPH[k][tj][i],bxC010TaperxPH[k][tj][i],byC010TaperxPH[k][tj][i],bzC010TaperxPH[k][tj][i],i,j,k);
        }
        for(i=1;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,index,pxZmin,vx,axC010TaperxPH[k][tj][i],ayC010TaperxPH[k][tj][i],azC010TaperxPH[k][tj][i],bxC010TaperxPH[k][tj][i],byC010TaperxPH[k][tj][i],bzC010TaperxPH[k][tj][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,index,pxZmin,vx,axE021Taper[k][tj],ayE021Taper[k][tj],azE021Taper[k][tj],bxE021Taper[k][tj],byE021Taper[k][tj],bzE021Taper[k][tj],i,j,k);
        }
        for(i=mxEnd1;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,index,pxZmin,vx,axC110TaperxPH[k][tj][ti],ayC110TaperxPH[k][tj][ti],azC110TaperxPH[k][tj][ti],bxC110TaperxPH[k][tj][ti],byC110TaperxPH[k][tj][ti],bzC110TaperxPH[k][tj][ti],i,j,k);
        }
        if (fmaxXb) {
          i=xLim;
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC110TaperxPH[k][tj][ti],ayC110TaperxPH[k][tj][ti],azC110TaperxPH[k][tj][ti],bxC110TaperxPH[k][tj][ti],byC110TaperxPH[k][tj][ti],bzC110TaperxPH[k][tj][ti],i,j,k);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          int tj = j-yEnd1;
          for(i=1;i<_nXmin;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = index;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC010TaperxPH[k][tj][i],ayC010TaperxPH[k][tj][i],azC010TaperxPH[k][tj][i],bxC010TaperxPH[k][tj][i],byC010TaperxPH[k][tj][i],bzC010TaperxPH[k][tj][i],i,j,k);
          }
          for(i=_nXmin;i<mxEnd1;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = index;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axE021Taper[k][tj],ayE021Taper[k][tj],azE021Taper[k][tj],bxE021Taper[k][tj],byE021Taper[k][tj],bzE021Taper[k][tj],i,j,k);
          }
          for(i=mxEnd1;i<xLim;++i) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = index;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,axC110TaperxPH[k][tj][ti],ayC110TaperxPH[k][tj][ti],azC110TaperxPH[k][tj][ti],bxC110TaperxPH[k][tj][ti],byC110TaperxPH[k][tj][ti],bzC110TaperxPH[k][tj][ti],i,j,k);
          }
        }
      }
    }
  }
  kk = zStart;
  for(k=zStart;k<zEnd;++k) {
    jj=0;
    nxy = NX*_nYmin;
    if(minYb) {
      if(fminYb) {
        j=0;
        for(i=1;i<_nXmin;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmin,vx,axE110TaperxPH[j][i],ayE110TaperxPH[j][i],azE110TaperxPH[j][i],bxE110TaperxPH[j][i],byE110TaperxPH[j][i],bzE110TaperxPH[j][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmin,vx,axY0Taper[j],ayTaper[j],azY0Taper[j],bxY0Taper[j],byTaper[j],bzY0Taper[j],i,j,k);
        }
        for(i=mxEnd1;i<xLim;++i) {
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmin,vx,axE210TaperxPH[j][ti],ayE210TaperxPH[j][ti],azE210TaperxPH[j][ti],bxE210TaperxPH[j][ti],byE210TaperxPH[j][ti],bzE210TaperxPH[j][ti],i,j,k);
        }
        j=1;
        for(i=1;i<_nXmin;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmin,vx,axE110TaperxPH[j][i],ayE110TaperxPH[j][i],azE110TaperxPH[j][i],bxE110TaperxPH[j][i],byE110TaperxPH[j][i],bzE110TaperxPH[j][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmin,vx,axY0Taper[j],ayTaper[j],azY0Taper[j],bxY0Taper[j],byTaper[j],bzY0Taper[j],i,j,k);
        }
        for(i=mxEnd1;i<xLim;++i) {
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmin,vx,axE210TaperxPH[j][ti],ayE210TaperxPH[j][ti],azE210TaperxPH[j][ti],bxE210TaperxPH[j][ti],byE210TaperxPH[j][ti],bzE210TaperxPH[j][ti],i,j,k);
        }
      }
      for(j=2;j<yStart2;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmin,vx,axE110TaperxPH[j][i],ayE110TaperxPH[j][i],azE110TaperxPH[j][i],bxE110TaperxPH[j][i],byE110TaperxPH[j][i],bzE110TaperxPH[j][i],i,j,k);
        }
        for(i=1;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY, pIndex=i+j*NX+(k-zStart)*nxymin;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxYmin,vx,axE110TaperxPH[j][i],ayE110TaperxPH[j][i],azE110TaperxPH[j][i],bxE110TaperxPH[j][i],byE110TaperxPH[j][i],bzE110TaperxPH[j][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;i++) {
          int index = i+j*NX+k*NXY, pIndex=i+j*NX+(k-zStart)*nxymin;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxYmin,vx,axY0Taper[j],ayTaper[j],azY0Taper[j],bxY0Taper[j],byTaper[j],bzY0Taper[j],i,j,k);
        }
        for(i=mxEnd1;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=i+j*NX+(k-zStart)*nxymin;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxYmin,vx,axE210TaperxPH[j][ti],ayE210TaperxPH[j][ti],azE210TaperxPH[j][ti],bxE210TaperxPH[j][ti],byE210TaperxPH[j][ti],bzE210TaperxPH[j][ti],i,j,k);
        }
        if (fmaxXb) {
          i=xLim;
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmin,vx,axE210TaperxPH[j][ti],ayE210TaperxPH[j][ti],azE210TaperxPH[j][ti],bxE210TaperxPH[j][ti],byE210TaperxPH[j][ti],bzE210TaperxPH[j][ti],i,j,k);
        }
      }
    }
    jj = yStart2;
    for(j=yStart2;j<yEnd;++j) {
      if(minXb) {
        nxy = _nXmin*nyInt;
        nx = _nXmin;
        ii=0;
        if(fminXb) {
          i=0;
          int index = i+j*NX+k*NXY, pIndex=i+(j-yStart2)*_nXmin+(k-zStart)*nxminy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxXmin,vx,axTaperPH[i],ayX0TaperPH[i],azX0TaperPH[i],bxTaperPH[i],byX0TaperPH[i],bzX0TaperPH[i],i,j,k);
        }
        for(i=1;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY, pIndex=i+(j-yStart2)*_nXmin+(k-zStart)*nxminy;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxXmin,vx,axTaperPH[i],ayX0TaperPH[i],azX0TaperPH[i],bxTaperPH[i],byX0TaperPH[i],bzX0TaperPH[i],i,j,k);
        }
      }
      for(i=mxStart1;i<mxEnd1;++i) {
        int index = i+j*NX+k*NXY;
        updateAcousticVxNoPML(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,vx);
      }
      if (maxXb) {
        ii=_xMaxStart-1;
        nxy = _nXmax*nyInt;
        nx = _nXmax;
        for(i=xEnd1;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=i-_xMaxStart+1+(j-yStart2)*_nXmax+(k-zStart)*nxmaxy;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxXmax,vx,axTaperPH[i],ayX1TaperPH[ti],azX1TaperPH[ti],bxTaperPH[i],byX1TaperPH[ti],bzX1TaperPH[ti],i,j,k);
        }
        if(fmaxXb) {
          i=xLim;
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=i-_xMaxStart+1+(j-yStart2)*_nXmax+(k-zStart)*nxmaxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxXmax,vx,axTaperPH[i],ayX1TaperPH[ti],azX1TaperPH[ti],bxTaperPH[i],byX1TaperPH[ti],bzX1TaperPH[ti],i,j,k);
        }
      }
    }
    ii=0;
    jj = _yMaxStart;
    nxy = NX*_nYmax;
    nx = NX;
    if(maxYb) {
      for(j=yEnd;j<yLim;++j) {
        int tj = j-yEnd1;
        if (fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmax,vx,axE120TaperxPH[tj][i],ayE120TaperxPH[tj][i],azE120TaperxPH[tj][i],bxE120TaperxPH[tj][i],byE120TaperxPH[tj][i],bzE120TaperxPH[tj][i],i,j,k);
        }
        for(i=1;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY, pIndex=i+(j-_yMaxStart)*NX+(k-zStart)*nxymax;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxYmax,vx,axE120TaperxPH[tj][i],ayE120TaperxPH[tj][i],azE120TaperxPH[tj][i],bxE120TaperxPH[tj][i],byE120TaperxPH[tj][i],bzE120TaperxPH[tj][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;i++) {
          int index = i+j*NX+k*NXY, pIndex=i+(j-_yMaxStart)*NX+(k-zStart)*nxymax;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxYmax,vx,axY1Taper[tj],ayTaper[j],azY1Taper[tj],bxY1Taper[tj],byTaper[j],bzY1Taper[tj],i,j,k);
        }
        for(i=mxEnd1;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=i+(j-_yMaxStart)*NX+(k-zStart)*nxymax;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxYmax,vx,axE220TaperxPH[tj][ti],ayE220TaperxPH[tj][ti],azE220TaperxPH[tj][ti],bxE220TaperxPH[tj][ti],byE220TaperxPH[tj][ti],bzE220TaperxPH[tj][ti],i,j,k);
        }
        if (fmaxXb) {
          i=xLim;
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmax,vx,axE220TaperxPH[tj][ti],ayE220TaperxPH[tj][ti],azE220TaperxPH[tj][ti],bxE220TaperxPH[tj][ti],byE220TaperxPH[tj][ti],bzE220TaperxPH[tj][ti],i,j,k);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          int tj = j-yEnd1;
          for(i=1;i<_nXmin;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmax,vx,axE120TaperxPH[tj][i],ayE120TaperxPH[tj][i],azE120TaperxPH[tj][i],bxE120TaperxPH[tj][i],byE120TaperxPH[tj][i],bzE120TaperxPH[tj][i],i,j,k);
          }
          for(i=_nXmin;i<mxEnd1;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmax,vx,axY1Taper[tj],ayTaper[j],azY1Taper[tj],bxY1Taper[tj],byTaper[j],bzY1Taper[tj],i,j,k);
          }
          for(i=mxEnd1;i<xLim;++i) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmax,vx,axE220TaperxPH[tj][ti],ayE220TaperxPH[tj][ti],azE220TaperxPH[tj][ti],bxE220TaperxPH[tj][ti],byE220TaperxPH[tj][ti],bzE220TaperxPH[tj][ti],i,j,k);
          }
        }
      }
    }
  }
  jj=0;
  kk = _zMaxStart;
  nxy = NXY;
  if(maxZb) {
    for(k=zEnd;k<zLim;++k) {
      int tk = k-zEnd1;
      if(fminYb) {
        j=0;
        for(i=1;i<_nXmin;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axC001TaperxPH[tk][j][i],ayC001TaperxPH[tk][j][i],azC001TaperxPH[tk][j][i],bxC001TaperxPH[tk][j][i],byC001TaperxPH[tk][j][i],bzC001TaperxPH[tk][j][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axE012Taper[tk][j],ayE012Taper[tk][j],azE012Taper[tk][j],bxE012Taper[tk][j],byE012Taper[tk][j],bzE012Taper[tk][j],i,j,k);
        }
        for(i=mxEnd1;i<xLim;++i) {
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axC101TaperxPH[tk][j][ti],ayC101TaperxPH[tk][j][ti],azC101TaperxPH[tk][j][ti],bxC101TaperxPH[tk][j][ti],byC101TaperxPH[tk][j][ti],bzC101TaperxPH[tk][j][ti],i,j,k);
        }
        j=1;
        for(i=1;i<_nXmin;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axC001TaperxPH[tk][j][i],ayC001TaperxPH[tk][j][i],azC001TaperxPH[tk][j][i],bxC001TaperxPH[tk][j][i],byC001TaperxPH[tk][j][i],bzC001TaperxPH[tk][j][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axE012Taper[tk][j],ayE012Taper[tk][j],azE012Taper[tk][j],bxE012Taper[tk][j],byE012Taper[tk][j],bzE012Taper[tk][j],i,j,k);
        }
        for(i=mxEnd1;i<xLim;++i) {
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axC101TaperxPH[tk][j][ti],ayC101TaperxPH[tk][j][ti],azC101TaperxPH[tk][j][ti],bxC101TaperxPH[tk][j][ti],byC101TaperxPH[tk][j][ti],bzC101TaperxPH[tk][j][ti],i,j,k);
        }
      }
      for(j=2;j<yStart2;++j) {
        if (fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axC001TaperxPH[tk][j][i],ayC001TaperxPH[tk][j][i],azC001TaperxPH[tk][j][i],bxC001TaperxPH[tk][j][i],byC001TaperxPH[tk][j][i],bzC001TaperxPH[tk][j][i],i,j,k);
        }
        for(i=1;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxZmax,vx,axC001TaperxPH[tk][j][i],ayC001TaperxPH[tk][j][i],azC001TaperxPH[tk][j][i],bxC001TaperxPH[tk][j][i],byC001TaperxPH[tk][j][i],bzC001TaperxPH[tk][j][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxZmax,vx,axE012Taper[tk][j],ayE012Taper[tk][j],azE012Taper[tk][j],bxE012Taper[tk][j],byE012Taper[tk][j],bzE012Taper[tk][j],i,j,k);
        }
        for(i=mxEnd1;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxZmax,vx,axC101TaperxPH[tk][j][ti],ayC101TaperxPH[tk][j][ti],azC101TaperxPH[tk][j][ti],bxC101TaperxPH[tk][j][ti],byC101TaperxPH[tk][j][ti],bzC101TaperxPH[tk][j][ti],i,j,k);
        }
        if (fmaxXb) {
          i=xLim;
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axC101TaperxPH[tk][j][ti],ayC101TaperxPH[tk][j][ti],azC101TaperxPH[tk][j][ti],bxC101TaperxPH[tk][j][ti],byC101TaperxPH[tk][j][ti],bzC101TaperxPH[tk][j][ti],i,j,k);
        }
      }
      for(j=yStart2;j<yEnd;++j) {
        if (fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axE102TaperxPH[tk][i],ayE102TaperxPH[tk][i],azE102TaperxPH[tk][i],bxE102TaperxPH[tk][i],byE102TaperxPH[tk][i],bzE102TaperxPH[tk][i],i,j,k);
        }
        for(i=1;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxZmax,vx,axE102TaperxPH[tk][i],ayE102TaperxPH[tk][i],azE102TaperxPH[tk][i],bxE102TaperxPH[tk][i],byE102TaperxPH[tk][i],bzE102TaperxPH[tk][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxZmax,vx,axZ1Taper[tk],ayZ1Taper[tk],azTaper[k],bxZ1Taper[tk],byZ1Taper[tk],bzTaper[k],i,j,k);
        }
        for(i=mxEnd1;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxZmax,vx,axE202TaperxPH[tk][ti],ayE202TaperxPH[tk][ti],azE202TaperxPH[tk][ti],bxE202TaperxPH[tk][ti],byE202TaperxPH[tk][ti],bzE202TaperxPH[tk][ti],i,j,k);
        }
        if (fmaxXb) {
          i=xLim;
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axE202TaperxPH[tk][ti],ayE202TaperxPH[tk][ti],azE202TaperxPH[tk][ti],bxE202TaperxPH[tk][ti],byE202TaperxPH[tk][ti],bzE202TaperxPH[tk][ti],i,j,k);
        }
      }
      for(j=yEnd;j<yLim;++j) {
        int tj = j-yEnd1;
        if (fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axC011TaperxPH[tk][tj][i],ayC011TaperxPH[tk][tj][i],azC011TaperxPH[tk][tj][i],bxC011TaperxPH[tk][tj][i],byC011TaperxPH[tk][tj][i],bzC011TaperxPH[tk][tj][i],i,j,k);
        }
        for(i=1;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxZmax,vx,axC011TaperxPH[tk][tj][i],ayC011TaperxPH[tk][tj][i],azC011TaperxPH[tk][tj][i],bxC011TaperxPH[tk][tj][i],byC011TaperxPH[tk][tj][i],bzC011TaperxPH[tk][tj][i],i,j,k);
        }
        for(i=_nXmin;i<mxEnd1;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxZmax,vx,axE022Taper[tk][tj],ayE022Taper[tk][tj],azE022Taper[tk][tj],bxE022Taper[tk][tj],byE022Taper[tk][tj],bzE022Taper[tk][tj],i,j,k);
        }
        for(i=mxEnd1;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,pIndex,index,pxZmax,vx,axC111TaperxPH[tk][tj][ti],ayC111TaperxPH[tk][tj][ti],azC111TaperxPH[tk][tj][ti],bxC111TaperxPH[tk][tj][ti],byC111TaperxPH[tk][tj][ti],bzC111TaperxPH[tk][tj][ti],i,j,k);
        }
        if (fmaxXb) {
          i=xLim;
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axC111TaperxPH[tk][tj][ti],ayC111TaperxPH[tk][tj][ti],azC111TaperxPH[tk][tj][ti],bxC111TaperxPH[tk][tj][ti],byC111TaperxPH[tk][tj][ti],bzC111TaperxPH[tk][tj][ti],i,j,k);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          int tj = j-yEnd1;
          for(i=1;i<_nXmin;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axC011TaperxPH[tk][tj][i],ayC011TaperxPH[tk][tj][i],azC011TaperxPH[tk][tj][i],bxC011TaperxPH[tk][tj][i],byC011TaperxPH[tk][tj][i],bzC011TaperxPH[tk][tj][i],i,j,k);
          }
          for(i=_nXmin;i<mxEnd1;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axE022Taper[tk][tj],ayE022Taper[tk][tj],azE022Taper[tk][tj],bxE022Taper[tk][tj],byE022Taper[tk][tj],bzE022Taper[tk][tj],i,j,k);
          }
          for(i=mxEnd1;i<xLim;++i) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axC111TaperxPH[tk][tj][ti],ayC111TaperxPH[tk][tj][ti],azC111TaperxPH[tk][tj][ti],bxC111TaperxPH[tk][tj][ti],byC111TaperxPH[tk][tj][ti],bzC111TaperxPH[tk][tj][ti],i,j,k);
          }
        }
      }
    }
    if(fmaxZb) {
      for(k=zLim;k<NZ;k++) {
        int tk = k-zEnd1;
        for(j=2;j<yStart2;j++){
          for(i=1;i<_nXmin;i++){
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axC001TaperxPH[tk][j][i],ayC001TaperxPH[tk][j][i],azC001TaperxPH[tk][j][i],bxC001TaperxPH[tk][j][i],byC001TaperxPH[tk][j][i],bzC001TaperxPH[tk][j][i],i,j,k);
          }
          for(i=_nXmin;i<mxEnd1;i++){
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axE012Taper[tk][j],ayE012Taper[tk][j],azE012Taper[tk][j],bxE012Taper[tk][j],byE012Taper[tk][j],bzE012Taper[tk][j],i,j,k);
          }
          for(i=mxEnd1;i<xLim;i++){
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axC101TaperxPH[tk][j][ti],ayC101TaperxPH[tk][j][ti],azC101TaperxPH[tk][j][ti],bxC101TaperxPH[tk][j][ti],byC101TaperxPH[tk][j][ti],bzC101TaperxPH[tk][j][ti],i,j,k);
          }
        }
        for(j=yStart2;j<yEnd;j++){
          for(i=1;i<_nXmin;i++){
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axE102TaperxPH[tk][i],ayE102TaperxPH[tk][i],azE102TaperxPH[tk][i],bxE102TaperxPH[tk][i],byE102TaperxPH[tk][i],bzE102TaperxPH[tk][i],i,j,k);
          }
          for(i=_nXmin;i<mxEnd1;i++){
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axZ1Taper[tk],ayZ1Taper[tk],azTaper[k],bxZ1Taper[tk],byZ1Taper[tk],bzTaper[k],i,j,k);
          }
          for(i=mxEnd1;i<xLim;i++){
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axE202TaperxPH[tk][ti],ayE202TaperxPH[tk][ti],azE202TaperxPH[tk][ti],bxE202TaperxPH[tk][ti],byE202TaperxPH[tk][ti],bzE202TaperxPH[tk][ti],i,j,k);
          }
        }
        for(j=yEnd;j<yLim;j++){
          int tj = j-yEnd1;
          for(i=1;i<_nXmin;i++){
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axC011TaperxPH[tk][tj][i],ayC011TaperxPH[tk][tj][i],azC011TaperxPH[tk][tj][i],bxC011TaperxPH[tk][tj][i],byC011TaperxPH[tk][tj][i],bzC011TaperxPH[tk][tj][i],i,j,k);
          }
          for(i=_nXmin;i<mxEnd1;i++){
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axE022Taper[tk][tj],ayE022Taper[tk][tj],azE022Taper[tk][tj],bxE022Taper[tk][tj],byE022Taper[tk][tj],bzE022Taper[tk][tj],i,j,k);
          }
          for(i=mxEnd1;i<xLim;i++){
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,axC111TaperxPH[tk][tj][ti],ayC111TaperxPH[tk][tj][ti],azC111TaperxPH[tk][tj][ti],bxC111TaperxPH[tk][tj][ti],byC111TaperxPH[tk][tj][ti],bzC111TaperxPH[tk][tj][ti],i,j,k);
          }
        }
      }
    }
  }
}

//vy updating. see vx updating for description
void updateVyAttenMPMLBounds(modelDefStruct* modelDef,
                        float* __restrict__ vy,
                        float* __restrict__ pressure,
                        float cx[2],float cy[2],float cz[2],
                        float* __restrict__ rho,
                        unsigned char* __restrict__ vyfunc,
                        bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                        bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb) {
  using namespace acousticAtten_MPML;
  DEF_MODEL_SIZE(modelDef);
  
  //  minXb = maxXb = minYb = maxYb = minZb = maxZb = false;
  //  maxZb=false;
  
  register int i,j,k;
  int zLim=NZ-2,yLim=NY-2,xLim=NX-2;
  float cx0=cx[0],cx1=cx[1];
  float cy0=cy[0],cy1=cy[1];
  float cz0=cz[0],cz1=cz[1];
  float dcx0=modelDef->dt/modelDef->dx*modelDef->scalarSpeed;
  float dcy0=modelDef->dt/modelDef->dy*modelDef->scalarSpeed;
  float dcz0=modelDef->dt/modelDef->dz*modelDef->scalarSpeed;
  int nxymin = NX*_nYmin, nxymax = NX*_nYmax, nxminy = _nXmin*(_yMaxStart-_nYmin), nxmaxy=_nXmax*(_yMaxStart-_nYmin);
  int xStart2 = MAX(_nXmin,2);
  int xEnd = MIN(_xMaxStart,NX-2);
  
  //printf("vy %d %d %d, %d %d %d %d %d %d\n",NX,NY,NZ,minXb,maxXb,minYb,maxYb,minZb,maxZb);
  int ii=0, jj=0, kk=0;
  int nxy = NXY, nx = NX;
  if(minZb) {
    if(fminZb) {
      k=0;
      for(j=1;j<yStart1;j++){
        for(i=2;i<xStart2;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC000TaperyPH[k][j][i],ayC000TaperyPH[k][j][i],azC000TaperyPH[k][j][i],bxC000TaperyPH[k][j][i],byC000TaperyPH[k][j][i],bzC000TaperyPH[k][j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axE011TaperyPH[k][j],ayE011TaperyPH[k][j],azE011TaperyPH[k][j],bxE011TaperyPH[k][j],byE011TaperyPH[k][j],bzE011TaperyPH[k][j],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC100TaperyPH[k][j][ti],ayC100TaperyPH[k][j][ti],azC100TaperyPH[k][j][ti],bxC100TaperyPH[k][j][ti],byC100TaperyPH[k][j][ti],bzC100TaperyPH[k][j][ti],i,j,k);
        }
      }
      for(j=yStart1;j<yEnd1;j++){
        for(i=2;i<xStart2;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axE101Taper[k][i],ayE101Taper[k][i],azE101Taper[k][i],bxE101Taper[k][i],byE101Taper[k][i],bzE101Taper[k][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axZ0Taper[k],ayZ0Taper[k],azTaper[k],bxZ0Taper[k],byZ0Taper[k],bzTaper[k],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axE201Taper[k][ti],ayE201Taper[k][ti],azE201Taper[k][ti],bxE201Taper[k][ti],byE201Taper[k][ti],bzE201Taper[k][ti],i,j,k);
        }
      }
      for(j=yEnd1;j<yLim;j++){
        int tj = j-yEnd1;
        for(i=2;i<xStart2;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC010TaperyPH[k][tj][i],ayC010TaperyPH[k][tj][i],azC010TaperyPH[k][tj][i],bxC010TaperyPH[k][tj][i],byC010TaperyPH[k][tj][i],bzC010TaperyPH[k][tj][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axE021TaperyPH[k][tj],ayE021TaperyPH[k][tj],azE021TaperyPH[k][tj],bxE021TaperyPH[k][tj],byE021TaperyPH[k][tj],bzE021TaperyPH[k][tj],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC110TaperyPH[k][tj][ti],ayC110TaperyPH[k][tj][ti],azC110TaperyPH[k][tj][ti],bxC110TaperyPH[k][tj][ti],byC110TaperyPH[k][tj][ti],bzC110TaperyPH[k][tj][ti],i,j,k);
        }
      }
      k=1;
      for(j=1;j<yStart1;j++){
        for(i=2;i<xStart2;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC000TaperyPH[k][j][i],ayC000TaperyPH[k][j][i],azC000TaperyPH[k][j][i],bxC000TaperyPH[k][j][i],byC000TaperyPH[k][j][i],bzC000TaperyPH[k][j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axE011TaperyPH[k][j],ayE011TaperyPH[k][j],azE011TaperyPH[k][j],bxE011TaperyPH[k][j],byE011TaperyPH[k][j],bzE011TaperyPH[k][j],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC100TaperyPH[k][j][ti],ayC100TaperyPH[k][j][ti],azC100TaperyPH[k][j][ti],bxC100TaperyPH[k][j][ti],byC100TaperyPH[k][j][ti],bzC100TaperyPH[k][j][ti],i,j,k);
        }
      }
      for(j=yStart1;j<yEnd1;j++){
        for(i=2;i<xStart2;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axE101Taper[k][i],ayE101Taper[k][i],azE101Taper[k][i],bxE101Taper[k][i],byE101Taper[k][i],bzE101Taper[k][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axZ0Taper[k],ayZ0Taper[k],azTaper[k],bxZ0Taper[k],byZ0Taper[k],bzTaper[k],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axE201Taper[k][ti],ayE201Taper[k][ti],azE201Taper[k][ti],bxE201Taper[k][ti],byE201Taper[k][ti],bzE201Taper[k][ti],i,j,k);
        }
      }
      for(j=yEnd1;j<yLim;j++){
        int tj = j-yEnd1;
        for(i=2;i<xStart2;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC010TaperyPH[k][tj][i],ayC010TaperyPH[k][tj][i],azC010TaperyPH[k][tj][i],bxC010TaperyPH[k][tj][i],byC010TaperyPH[k][tj][i],bzC010TaperyPH[k][tj][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axE021TaperyPH[k][tj],ayE021TaperyPH[k][tj],azE021TaperyPH[k][tj],bxE021TaperyPH[k][tj],byE021TaperyPH[k][tj],bzE021TaperyPH[k][tj],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC110TaperyPH[k][tj][ti],ayC110TaperyPH[k][tj][ti],azC110TaperyPH[k][tj][ti],bxC110TaperyPH[k][tj][ti],byC110TaperyPH[k][tj][ti],bzC110TaperyPH[k][tj][ti],i,j,k);
        }
      }
    }
    for(k=2;k<zStart;++k) {
      if(fminYb) {
        j=0;
        for(i=2;i<xStart2;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC000TaperyPH[k][j][i],ayC000TaperyPH[k][j][i],azC000TaperyPH[k][j][i],bxC000TaperyPH[k][j][i],byC000TaperyPH[k][j][i],bzC000TaperyPH[k][j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axE011TaperyPH[k][j],ayE011TaperyPH[k][j],azE011TaperyPH[k][j],bxE011TaperyPH[k][j],byE011TaperyPH[k][j],bzE011TaperyPH[k][j],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC100TaperyPH[k][j][ti],ayC100TaperyPH[k][j][ti],azC100TaperyPH[k][j][ti],bxC100TaperyPH[k][j][ti],byC100TaperyPH[k][j][ti],bzC100TaperyPH[k][j][ti],i,j,k);
        }
      }
      for(j=1;j<yStart1;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC000TaperyPH[k][j][i],ayC000TaperyPH[k][j][i],azC000TaperyPH[k][j][i],bxC000TaperyPH[k][j][i],byC000TaperyPH[k][j][i],bzC000TaperyPH[k][j][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC000TaperyPH[k][j][i],ayC000TaperyPH[k][j][i],azC000TaperyPH[k][j][i],bxC000TaperyPH[k][j][i],byC000TaperyPH[k][j][i],bzC000TaperyPH[k][j][i],i,j,k);
        }
        for(i=2;i<xStart2;i++) {
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmin,vy,axC000TaperyPH[k][j][i],ayC000TaperyPH[k][j][i],azC000TaperyPH[k][j][i],bxC000TaperyPH[k][j][i],byC000TaperyPH[k][j][i],bzC000TaperyPH[k][j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmin,vy,axE011TaperyPH[k][j],ayE011TaperyPH[k][j],azE011TaperyPH[k][j],bxE011TaperyPH[k][j],byE011TaperyPH[k][j],bzE011TaperyPH[k][j],i,j,k);
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmin,vy,axC100TaperyPH[k][j][ti],ayC100TaperyPH[k][j][ti],azC100TaperyPH[k][j][ti],bxC100TaperyPH[k][j][ti],byC100TaperyPH[k][j][ti],bzC100TaperyPH[k][j][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = index;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC100TaperyPH[k][j][ti],ayC100TaperyPH[k][j][ti],azC100TaperyPH[k][j][ti],bxC100TaperyPH[k][j][ti],byC100TaperyPH[k][j][ti],bzC100TaperyPH[k][j][ti],i,j,k);
          }
        }
      }
      for(j=yStart1;j<yEnd1;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axE101Taper[k][i],ayE101Taper[k][i],azE101Taper[k][i],bxE101Taper[k][i],byE101Taper[k][i],bzE101Taper[k][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axE101Taper[k][i],ayE101Taper[k][i],azE101Taper[k][i],bxE101Taper[k][i],byE101Taper[k][i],bzE101Taper[k][i],i,j,k);
        }
        for(i=2;i<xStart2;i++) {
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmin,vy,axE101Taper[k][i],ayE101Taper[k][i],azE101Taper[k][i],bxE101Taper[k][i],byE101Taper[k][i],bzE101Taper[k][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmin,vy,axZ0Taper[k],ayZ0Taper[k],azTaper[k],bxZ0Taper[k],byZ0Taper[k],bzTaper[k],i,j,k);
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmin,vy,axE201Taper[k][ti],ayE201Taper[k][ti],azE201Taper[k][ti],bxE201Taper[k][ti],byE201Taper[k][ti],bzE201Taper[k][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = index;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axE201Taper[k][ti],ayE201Taper[k][ti],azE201Taper[k][ti],bxE201Taper[k][ti],byE201Taper[k][ti],bzE201Taper[k][ti],i,j,k);
          }
        }
      }
      for(j=yEnd1;j<yLim;++j) {
        int tj = j-yEnd1;
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC010TaperyPH[k][tj][i],ayC010TaperyPH[k][tj][i],azC010TaperyPH[k][tj][i],bxC010TaperyPH[k][tj][i],byC010TaperyPH[k][tj][i],bzC010TaperyPH[k][tj][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC010TaperyPH[k][tj][i],ayC010TaperyPH[k][tj][i],azC010TaperyPH[k][tj][i],bxC010TaperyPH[k][tj][i],byC010TaperyPH[k][tj][i],bzC010TaperyPH[k][tj][i],i,j,k);
        }
        for(i=2;i<xStart2;i++) {
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmin,vy,axC010TaperyPH[k][tj][i],ayC010TaperyPH[k][tj][i],azC010TaperyPH[k][tj][i],bxC010TaperyPH[k][tj][i],byC010TaperyPH[k][tj][i],bzC010TaperyPH[k][tj][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmin,vy,axE021TaperyPH[k][tj],ayE021TaperyPH[k][tj],azE021TaperyPH[k][tj],bxE021TaperyPH[k][tj],byE021TaperyPH[k][tj],bzE021TaperyPH[k][tj],i,j,k);
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmin,vy,axC110TaperyPH[k][tj][ti],ayC110TaperyPH[k][tj][ti],azC110TaperyPH[k][tj][ti],bxC110TaperyPH[k][tj][ti],byC110TaperyPH[k][tj][ti],bzC110TaperyPH[k][tj][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = index;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC110TaperyPH[k][tj][ti],ayC110TaperyPH[k][tj][ti],azC110TaperyPH[k][tj][ti],bxC110TaperyPH[k][tj][ti],byC110TaperyPH[k][tj][ti],bzC110TaperyPH[k][tj][ti],i,j,k);
          }
        }
      }
      if(fmaxYb) {
        j=yLim;
        int tj = j-yEnd1;
        for(i=2;i<xStart2;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC010TaperyPH[k][tj][i],ayC010TaperyPH[k][tj][i],azC010TaperyPH[k][tj][i],bxC010TaperyPH[k][tj][i],byC010TaperyPH[k][tj][i],bzC010TaperyPH[k][tj][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axE021TaperyPH[k][tj],ayE021TaperyPH[k][tj],azE021TaperyPH[k][tj],bxE021TaperyPH[k][tj],byE021TaperyPH[k][tj],bzE021TaperyPH[k][tj],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,axC110TaperyPH[k][tj][ti],ayC110TaperyPH[k][tj][ti],azC110TaperyPH[k][tj][ti],bxC110TaperyPH[k][tj][ti],byC110TaperyPH[k][tj][ti],bzC110TaperyPH[k][tj][ti],i,j,k);
        }
      }
    }
  }
  kk = zStart;
  for(k=zStart;k<zEnd;++k) {
    jj=0;
    nxy = NX*_nYmin;
    if(minYb) {
      if(fminYb) {
        j=0;
        for(i=2;i<xStart2;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmin,vy,axE110TaperyPH[j][i],ayE110TaperyPH[j][i],azE110TaperyPH[j][i],bxE110TaperyPH[j][i],byE110TaperyPH[j][i],bzE110TaperyPH[j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmin,vy,axY0TaperPH[j],ayTaperPH[j],azY0TaperPH[j],bxY0TaperPH[j],byTaperPH[j],bzY0TaperPH[j],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmin,vy,axE210TaperyPH[j][ti],ayE210TaperyPH[j][ti],azE210TaperyPH[j][ti],bxE210TaperyPH[j][ti],byE210TaperyPH[j][ti],bzE210TaperyPH[j][ti],i,j,k);
        }
      }
      for(j=1;j<yStart1;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmin,vy,axE110TaperyPH[j][i],ayE110TaperyPH[j][i],azE110TaperyPH[j][i],bxE110TaperyPH[j][i],byE110TaperyPH[j][i],bzE110TaperyPH[j][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmin,vy,axE110TaperyPH[j][i],ayE110TaperyPH[j][i],azE110TaperyPH[j][i],bxE110TaperyPH[j][i],byE110TaperyPH[j][i],bzE110TaperyPH[j][i],i,j,k);
        }
        for(i=2;i<xStart2;++i) {
          int index = i+j*NX+k*NXY, pIndex=i+j*NX+(k-zStart)*nxymin;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyYmin,vy,axE110TaperyPH[j][i],ayE110TaperyPH[j][i],azE110TaperyPH[j][i],bxE110TaperyPH[j][i],byE110TaperyPH[j][i],bzE110TaperyPH[j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;++i) {
          int index = i+j*NX+k*NXY, pIndex=i+j*NX+(k-zStart)*nxymin;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyYmin,vy,axY0TaperPH[j],ayTaperPH[j],azY0TaperPH[j],bxY0TaperPH[j],byTaperPH[j],bzY0TaperPH[j],i,j,k);
        }
        for(i=xEnd;i<xLim;++i) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=i+j*NX+(k-zStart)*nxymin;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyYmin,vy,axE210TaperyPH[j][ti],ayE210TaperyPH[j][ti],azE210TaperyPH[j][ti],bxE210TaperyPH[j][ti],byE210TaperyPH[j][ti],bzE210TaperyPH[j][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmin,vy,axE210TaperyPH[j][ti],ayE210TaperyPH[j][ti],azE210TaperyPH[j][ti],bxE210TaperyPH[j][ti],byE210TaperyPH[j][ti],bzE210TaperyPH[j][ti],i,j,k);
          }
        }
      }
    }
    jj = yStart1;
    for(j=yStart1;j<yEnd1;++j) {
      if(minXb) {
        ii=0;
        nxy = _nXmin*nyInt;
        nx = _nXmin;
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyXmin,vy,axTaper[i],ayX0Taper[i],azX0Taper[i],bxTaper[i],byX0Taper[i],bzX0Taper[i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyXmin,vy,axTaper[i],ayX0Taper[i],azX0Taper[i],bxTaper[i],byX0Taper[i],bzX0Taper[i],i,j,k);
        }
        for(i=2;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY, pIndex=i+(j-yStart1)*_nXmin+(k-zStart)*nxminy;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyXmin,vy,axTaper[i],ayX0Taper[i],azX0Taper[i],bxTaper[i],byX0Taper[i],bzX0Taper[i],i,j,k);
        }
      }
      for(i=xStart2;i<xEnd;++i) {
        int index = i+j*NX+k*NXY;
        updateAcousticVyNoPML(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,index,vy);
      }
      if(maxXb) {
        ii=_xMaxStart;
        nxy = _nXmax*nyInt;
        nx = _nXmax;
        for(i=xEnd;i<xLim;++i) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=i-_xMaxStart+(j-yStart1)*_nXmax+(k-zStart)*nxmaxy;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyXmax,vy,axTaper[i],ayX1Taper[ti],azX1Taper[ti],bxTaper[i],byX1Taper[ti],bzX1Taper[ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyXmax,vy,axTaper[i],ayX1Taper[ti],azX1Taper[ti],bxTaper[i],byX1Taper[ti],bzX1Taper[ti],i,j,k);
          }
        }
      }
    }
    ii=0;
    jj = _yMaxStart-1;
    nxy = NX*_nYmax;
    nx = NX;
    if(maxYb) {
      for(j=yEnd1;j<yLim;++j) {
        int tj = j-yEnd1;
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmax,vy,axE120TaperyPH[tj][i],ayE120TaperyPH[tj][i],azE120TaperyPH[tj][i],bxE120TaperyPH[tj][i],byE120TaperyPH[tj][i],bzE120TaperyPH[tj][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmax,vy,axE120TaperyPH[tj][i],ayE120TaperyPH[tj][i],azE120TaperyPH[tj][i],bxE120TaperyPH[tj][i],byE120TaperyPH[tj][i],bzE120TaperyPH[tj][i],i,j,k);
        }
        for(i=2;i<xStart2;i++) {
          int index = i+j*NX+k*NXY, pIndex=i+(j-_yMaxStart+1)*NX+(k-zStart)*nxymax;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyYmax,vy,axE120TaperyPH[tj][i],ayE120TaperyPH[tj][i],azE120TaperyPH[tj][i],bxE120TaperyPH[tj][i],byE120TaperyPH[tj][i],bzE120TaperyPH[tj][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+j*NX+k*NXY, pIndex=i+(j-_yMaxStart+1)*NX+(k-zStart)*nxymax;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyYmax,vy,axY1TaperPH[tj],ayTaperPH[j],azY1TaperPH[tj],bxY1TaperPH[tj],byTaperPH[j],bzY1TaperPH[tj],i,j,k);
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=i+(j-_yMaxStart+1)*NX+(k-zStart)*nxymax;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyYmax,vy,axE220TaperyPH[tj][ti],ayE220TaperyPH[tj][ti],azE220TaperyPH[tj][ti],bxE220TaperyPH[tj][ti],byE220TaperyPH[tj][ti],bzE220TaperyPH[tj][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmax,vy,axE220TaperyPH[tj][ti],ayE220TaperyPH[tj][ti],azE220TaperyPH[tj][ti],bxE220TaperyPH[tj][ti],byE220TaperyPH[tj][ti],bzE220TaperyPH[tj][ti],i,j,k);
          }
        }
      }
      if(fmaxYb) {
        j=yLim;
        int tj = j-yEnd1;
        for(i=2;i<xStart2;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmax,vy,axE120TaperyPH[tj][i],ayE120TaperyPH[tj][i],azE120TaperyPH[tj][i],bxE120TaperyPH[tj][i],byE120TaperyPH[tj][i],bzE120TaperyPH[tj][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmax,vy,axY1TaperPH[tj],ayTaperPH[j],azY1TaperPH[tj],bxY1TaperPH[tj],byTaperPH[j],bzY1TaperPH[tj],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmax,vy,axE220TaperyPH[tj][ti],ayE220TaperyPH[tj][ti],azE220TaperyPH[tj][ti],bxE220TaperyPH[tj][ti],byE220TaperyPH[tj][ti],bzE220TaperyPH[tj][ti],i,j,k);
        }
      }
    }
  }
  jj=0;
  nxy=NXY;
  kk = _zMaxStart;
  if(maxZb) {
    for(k=zEnd;k<zLim;++k) {
      int tk = k-zEnd1;
      if(fminYb) {
        j=0;
        for(i=2;i<xStart2;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axC001TaperyPH[tk][j][i],ayC001TaperyPH[tk][j][i],azC001TaperyPH[tk][j][i],bxC001TaperyPH[tk][j][i],byC001TaperyPH[tk][j][i],bzC001TaperyPH[tk][j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axE012TaperyPH[tk][j],ayE012TaperyPH[tk][j],azE012TaperyPH[tk][j],bxE012TaperyPH[tk][j],byE012TaperyPH[tk][j],bzE012TaperyPH[tk][j],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axC101TaperyPH[tk][j][ti],ayC101TaperyPH[tk][j][ti],azC101TaperyPH[tk][j][ti],bxC101TaperyPH[tk][j][ti],byC101TaperyPH[tk][j][ti],bzC101TaperyPH[tk][j][ti],i,j,k);
        }
      }
      for(j=1;j<yStart1;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axC001TaperyPH[tk][j][i],ayC001TaperyPH[tk][j][i],azC001TaperyPH[tk][j][i],bxC001TaperyPH[tk][j][i],byC001TaperyPH[tk][j][i],bzC001TaperyPH[tk][j][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axC001TaperyPH[tk][j][i],ayC001TaperyPH[tk][j][i],azC001TaperyPH[tk][j][i],bxC001TaperyPH[tk][j][i],byC001TaperyPH[tk][j][i],bzC001TaperyPH[tk][j][i],i,j,k);
        }
        for(i=2;i<xStart2;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmax,vy,axC001TaperyPH[tk][j][i],ayC001TaperyPH[tk][j][i],azC001TaperyPH[tk][j][i],bxC001TaperyPH[tk][j][i],byC001TaperyPH[tk][j][i],bzC001TaperyPH[tk][j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmax,vy,axE012TaperyPH[tk][j],ayE012TaperyPH[tk][j],azE012TaperyPH[tk][j],bxE012TaperyPH[tk][j],byE012TaperyPH[tk][j],bzE012TaperyPH[tk][j],i,j,k);
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmax,vy,axC101TaperyPH[tk][j][ti],ayC101TaperyPH[tk][j][ti],azC101TaperyPH[tk][j][ti],bxC101TaperyPH[tk][j][ti],byC101TaperyPH[tk][j][ti],bzC101TaperyPH[tk][j][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axC101TaperyPH[tk][j][ti],ayC101TaperyPH[tk][j][ti],azC101TaperyPH[tk][j][ti],bxC101TaperyPH[tk][j][ti],byC101TaperyPH[tk][j][ti],bzC101TaperyPH[tk][j][ti],i,j,k);
          }
        }
      }
      for(j=yStart1;j<yEnd1;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axE102Taper[tk][i],ayE102Taper[tk][i],azE102Taper[tk][i],bxE102Taper[tk][i],byE102Taper[tk][i],bzE102Taper[tk][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axE102Taper[tk][i],ayE102Taper[tk][i],azE102Taper[tk][i],bxE102Taper[tk][i],byE102Taper[tk][i],bzE102Taper[tk][i],i,j,k);
        }
        for(i=2;i<xStart2;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmax,vy,axE102Taper[tk][i],ayE102Taper[tk][i],azE102Taper[tk][i],bxE102Taper[tk][i],byE102Taper[tk][i],bzE102Taper[tk][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmax,vy,axZ1Taper[tk],ayZ1Taper[tk],azTaper[k],bxZ1Taper[tk],byZ1Taper[tk],bzTaper[k],i,j,k);
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmax,vy,axE202Taper[tk][ti],ayE202Taper[tk][ti],azE202Taper[tk][ti],bxE202Taper[tk][ti],byE202Taper[tk][ti],bzE202Taper[tk][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axE202Taper[tk][ti],ayE202Taper[tk][ti],azE202Taper[tk][ti],bxE202Taper[tk][ti],byE202Taper[tk][ti],bzE202Taper[tk][ti],i,j,k);
          }
        }
      }
      for(j=yEnd1;j<yLim;++j) {
        int tj = j-yEnd1;
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axC011TaperyPH[tk][tj][i],ayC011TaperyPH[tk][tj][i],azC011TaperyPH[tk][tj][i],bxC011TaperyPH[tk][tj][i],byC011TaperyPH[tk][tj][i],bzC011TaperyPH[tk][tj][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axC011TaperyPH[tk][tj][i],ayC011TaperyPH[tk][tj][i],azC011TaperyPH[tk][tj][i],bxC011TaperyPH[tk][tj][i],byC011TaperyPH[tk][tj][i],bzC011TaperyPH[tk][tj][i],i,j,k);
        }
        for(i=2;i<xStart2;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmax,vy,axC011TaperyPH[tk][tj][i],ayC011TaperyPH[tk][tj][i],azC011TaperyPH[tk][tj][i],bxC011TaperyPH[tk][tj][i],byC011TaperyPH[tk][tj][i],bzC011TaperyPH[tk][tj][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmax,vy,axE022TaperyPH[tk][tj],ayE022TaperyPH[tk][tj],azE022TaperyPH[tk][tj],bxE022TaperyPH[tk][tj],byE022TaperyPH[tk][tj],bzE022TaperyPH[tk][tj],i,j,k);
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmax,vy,axC111TaperyPH[tk][tj][ti],ayC111TaperyPH[tk][tj][ti],azC111TaperyPH[tk][tj][ti],bxC111TaperyPH[tk][tj][ti],byC111TaperyPH[tk][tj][ti],bzC111TaperyPH[tk][tj][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axC111TaperyPH[tk][tj][ti],ayC111TaperyPH[tk][tj][ti],azC111TaperyPH[tk][tj][ti],bxC111TaperyPH[tk][tj][ti],byC111TaperyPH[tk][tj][ti],bzC111TaperyPH[tk][tj][ti],i,j,k);
          }
        }
      }
      if(fmaxYb) {
        j=yLim;
        int tj = j-yEnd1;
        for(i=2;i<xStart2;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axC011TaperyPH[tk][tj][i],ayC011TaperyPH[tk][tj][i],azC011TaperyPH[tk][tj][i],bxC011TaperyPH[tk][tj][i],byC011TaperyPH[tk][tj][i],bzC011TaperyPH[tk][tj][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axE022TaperyPH[tk][tj],ayE022TaperyPH[tk][tj],azE022TaperyPH[tk][tj],bxE022TaperyPH[tk][tj],byE022TaperyPH[tk][tj],bzE022TaperyPH[tk][tj],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axC111TaperyPH[tk][tj][ti],ayC111TaperyPH[tk][tj][ti],azC111TaperyPH[tk][tj][ti],bxC111TaperyPH[tk][tj][ti],byC111TaperyPH[tk][tj][ti],bzC111TaperyPH[tk][tj][ti],i,j,k);
        }
      }
    }
    if(fmaxZb) {
      for(k=zLim;k<NZ;k++) {
        int tk = k-zEnd1;
        for(j=1;j<yStart1;j++){
          for(i=2;i<xStart2;i++){
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axC001TaperyPH[tk][j][i],ayC001TaperyPH[tk][j][i],azC001TaperyPH[tk][j][i],bxC001TaperyPH[tk][j][i],byC001TaperyPH[tk][j][i],bzC001TaperyPH[tk][j][i],i,j,k);
          }
          for(i=xStart2;i<xEnd;i++){
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axE012TaperyPH[tk][j],ayE012TaperyPH[tk][j],azE012TaperyPH[tk][j],bxE012TaperyPH[tk][j],byE012TaperyPH[tk][j],bzE012TaperyPH[tk][j],i,j,k);
          }
          for(i=xEnd;i<xLim;i++){
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axC101TaperyPH[tk][j][ti],ayC101TaperyPH[tk][j][ti],azC101TaperyPH[tk][j][ti],bxC101TaperyPH[tk][j][ti],byC101TaperyPH[tk][j][ti],bzC101TaperyPH[tk][j][ti],i,j,k);
          }
        }
        for(j=yStart1;j<yEnd1;j++){
          for(i=2;i<xStart2;i++){
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axE102Taper[tk][i],ayE102Taper[tk][i],azE102Taper[tk][i],bxE102Taper[tk][i],byE102Taper[tk][i],bzE102Taper[tk][i],i,j,k);
          }
          for(i=xStart2;i<xEnd;i++){
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axZ1Taper[tk],ayZ1Taper[tk],azTaper[k],bxZ1Taper[tk],byZ1Taper[tk],bzTaper[k],i,j,k);
          }
          for(i=xEnd;i<xLim;i++){
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axE202Taper[tk][ti],ayE202Taper[tk][ti],azE202Taper[tk][ti],bxE202Taper[tk][ti],byE202Taper[tk][ti],bzE202Taper[tk][ti],i,j,k);
          }
        }
        for(j=yEnd1;j<yLim;j++){
          int tj = j-yEnd1;
          for(i=2;i<xStart2;i++){
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axC011TaperyPH[tk][tj][i],ayC011TaperyPH[tk][tj][i],azC011TaperyPH[tk][tj][i],bxC011TaperyPH[tk][tj][i],byC011TaperyPH[tk][tj][i],bzC011TaperyPH[tk][tj][i],i,j,k);
          }
          for(i=xStart2;i<xEnd;i++){
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axE022TaperyPH[tk][tj],ayE022TaperyPH[tk][tj],azE022TaperyPH[tk][tj],bxE022TaperyPH[tk][tj],byE022TaperyPH[tk][tj],bzE022TaperyPH[tk][tj],i,j,k);
          }
          for(i=xEnd;i<xLim;i++){
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,axC111TaperyPH[tk][tj][ti],ayC111TaperyPH[tk][tj][ti],azC111TaperyPH[tk][tj][ti],bxC111TaperyPH[tk][tj][ti],byC111TaperyPH[tk][tj][ti],bzC111TaperyPH[tk][tj][ti],i,j,k);
          }
        }
      }
    }
  }
}

void updateVzAttenMPMLBounds(modelDefStruct* modelDef,
                        float* __restrict__ vz,
                        float* __restrict__ pressure,
                        float cx[2],float cy[2],float cz[2],
                        float* __restrict__ rho,
                        unsigned char* __restrict__ vzfunc,
                        bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                        bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb)  {
  using namespace acousticAtten_MPML;
  DEF_MODEL_SIZE(modelDef);
  
  register int i,j,k;
  int zLim=NZ-2,yLim=NY-2,xLim=NX-2;
  float cx0=cx[0],cx1=cx[1];
  float cy0=cy[0],cy1=cy[1];
  float cz0=cz[0],cz1=cz[1];
  float dcx0=modelDef->dt/modelDef->dx*modelDef->scalarSpeed;
  float dcy0=modelDef->dt/modelDef->dy*modelDef->scalarSpeed;
  float dcz0=modelDef->dt/modelDef->dz*modelDef->scalarSpeed;
  int nxymin = NX*_nYmin, nxymax = NX*_nYmax, nxminy = _nXmin*(_yMaxStart-_nYmin), nxmaxy=_nXmax*(_yMaxStart-_nYmin);
  int xStart2 = MAX(_nXmin,2);
  int xEnd = MIN(_xMaxStart,NX-2);
  
  //printf("vz %d %d %d, %d %d %d %d %d %d\n",NX,NY,NZ,minXb,maxXb,minYb,maxYb,minZb,maxZb);
  int ii=0, jj=0, kk=0;
  int nxy = NXY, nx = NX;
  if(minZb) {
    if(fminZb) {
      k=0;
      for(j=2;j<yStart2;j++){
        for(i=2;i<xStart2;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC000TaperzPH[k][j][i],ayC000TaperzPH[k][j][i],azC000TaperzPH[k][j][i],bxC000TaperzPH[k][j][i],byC000TaperzPH[k][j][i],bzC000TaperzPH[k][j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axE011TaperzPH[k][j],ayE011TaperzPH[k][j],azE011TaperzPH[k][j],bxE011TaperzPH[k][j],byE011TaperzPH[k][j],bzE011TaperzPH[k][j],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC100TaperzPH[k][j][ti],ayC100TaperzPH[k][j][ti],azC100TaperzPH[k][j][ti],bxC100TaperzPH[k][j][ti],byC100TaperzPH[k][j][ti],bzC100TaperzPH[k][j][ti],i,j,k);
        }
      }
      for(j=yStart2;j<yEnd;j++){
        for(i=2;i<xStart2;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axE101TaperzPH[k][i],ayE101TaperzPH[k][i],azE101TaperzPH[k][i],bxE101TaperzPH[k][i],byE101TaperzPH[k][i],bzE101TaperzPH[k][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axZ0TaperPH[k],ayZ0TaperPH[k],azTaperPH[k],bxZ0TaperPH[k],byZ0TaperPH[k],bzTaperPH[k],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axE201TaperzPH[k][ti],ayE201TaperzPH[k][ti],azE201TaperzPH[k][ti],bxE201TaperzPH[k][ti],byE201TaperzPH[k][ti],bzE201TaperzPH[k][ti],i,j,k);
        }
      }
      for(j=yEnd;j<yLim;j++){
        int tj = j-yEnd1;
        for(i=2;i<xStart2;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC010TaperzPH[k][tj][i],ayC010TaperzPH[k][tj][i],azC010TaperzPH[k][tj][i],bxC010TaperzPH[k][tj][i],byC010TaperzPH[k][tj][i],bzC010TaperzPH[k][tj][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axE021TaperzPH[k][tj],ayE021TaperzPH[k][tj],azE021TaperzPH[k][tj],bxE021TaperzPH[k][tj],byE021TaperzPH[k][tj],bzE021TaperzPH[k][tj],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC110TaperzPH[k][tj][ti],ayC110TaperzPH[k][tj][ti],azC110TaperzPH[k][tj][ti],bxC110TaperzPH[k][tj][ti],byC110TaperzPH[k][tj][ti],bzC110TaperzPH[k][tj][ti],i,j,k);
        }
      }
    }
    for(k=1;k<zStart_z;++k) {
      if(fminYb) {
        j=0;
        for(i=2;i<xStart2;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC000TaperzPH[k][j][i],ayC000TaperzPH[k][j][i],azC000TaperzPH[k][j][i],bxC000TaperzPH[k][j][i],byC000TaperzPH[k][j][i],bzC000TaperzPH[k][j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axE011TaperzPH[k][j],ayE011TaperzPH[k][j],azE011TaperzPH[k][j],bxE011TaperzPH[k][j],byE011TaperzPH[k][j],bzE011TaperzPH[k][j],i,j,k);
        }
        for(i=xEnd;i<xLim;++i) {
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC100TaperzPH[k][j][ti],ayC100TaperzPH[k][j][ti],azC100TaperzPH[k][j][ti],bxC100TaperzPH[k][j][ti],byC100TaperzPH[k][j][ti],bzC100TaperzPH[k][j][ti],i,j,k);
        }
        j=1;
        for(i=2;i<xStart2;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC000TaperzPH[k][j][i],ayC000TaperzPH[k][j][i],azC000TaperzPH[k][j][i],bxC000TaperzPH[k][j][i],byC000TaperzPH[k][j][i],bzC000TaperzPH[k][j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axE011TaperzPH[k][j],ayE011TaperzPH[k][j],azE011TaperzPH[k][j],bxE011TaperzPH[k][j],byE011TaperzPH[k][j],bzE011TaperzPH[k][j],i,j,k);
        }
        for(i=xEnd;i<xLim;++i) {
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC100TaperzPH[k][j][ti],ayC100TaperzPH[k][j][ti],azC100TaperzPH[k][j][ti],bxC100TaperzPH[k][j][ti],byC100TaperzPH[k][j][ti],bzC100TaperzPH[k][j][ti],i,j,k);
        }
      }
      for(j=2;j<yStart2;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC000TaperzPH[k][j][i],ayC000TaperzPH[k][j][i],azC000TaperzPH[k][j][i],bxC000TaperzPH[k][j][i],byC000TaperzPH[k][j][i],bzC000TaperzPH[k][j][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC000TaperzPH[k][j][i],ayC000TaperzPH[k][j][i],azC000TaperzPH[k][j][i],bxC000TaperzPH[k][j][i],byC000TaperzPH[k][j][i],bzC000TaperzPH[k][j][i],i,j,k);
        }
        for(i=2;i<xStart2;i++) {
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmin,vz,axC000TaperzPH[k][j][i],ayC000TaperzPH[k][j][i],azC000TaperzPH[k][j][i],bxC000TaperzPH[k][j][i],byC000TaperzPH[k][j][i],bzC000TaperzPH[k][j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmin,vz,axE011TaperzPH[k][j],ayE011TaperzPH[k][j],azE011TaperzPH[k][j],bxE011TaperzPH[k][j],byE011TaperzPH[k][j],bzE011TaperzPH[k][j],i,j,k);
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmin,vz,axC100TaperzPH[k][j][ti],ayC100TaperzPH[k][j][ti],azC100TaperzPH[k][j][ti],bxC100TaperzPH[k][j][ti],byC100TaperzPH[k][j][ti],bzC100TaperzPH[k][j][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC100TaperzPH[k][j][ti],ayC100TaperzPH[k][j][ti],azC100TaperzPH[k][j][ti],bxC100TaperzPH[k][j][ti],byC100TaperzPH[k][j][ti],bzC100TaperzPH[k][j][ti],i,j,k);
          }
        }
      }
      for(j=yStart2;j<yEnd;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axE101TaperzPH[k][i],ayE101TaperzPH[k][i],azE101TaperzPH[k][i],bxE101TaperzPH[k][i],byE101TaperzPH[k][i],bzE101TaperzPH[k][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axE101TaperzPH[k][i],ayE101TaperzPH[k][i],azE101TaperzPH[k][i],bxE101TaperzPH[k][i],byE101TaperzPH[k][i],bzE101TaperzPH[k][i],i,j,k);
        }
        for(i=2;i<xStart2;i++) {
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmin,vz,axE101TaperzPH[k][i],ayE101TaperzPH[k][i],azE101TaperzPH[k][i],bxE101TaperzPH[k][i],byE101TaperzPH[k][i],bzE101TaperzPH[k][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmin,vz,axZ0TaperPH[k],ayZ0TaperPH[k],azTaperPH[k],bxZ0TaperPH[k],byZ0TaperPH[k],bzTaperPH[k],i,j,k);
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmin,vz,axE201TaperzPH[k][ti],ayE201TaperzPH[k][ti],azE201TaperzPH[k][ti],bxE201TaperzPH[k][ti],byE201TaperzPH[k][ti],bzE201TaperzPH[k][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axE201TaperzPH[k][ti],ayE201TaperzPH[k][ti],azE201TaperzPH[k][ti],bxE201TaperzPH[k][ti],byE201TaperzPH[k][ti],bzE201TaperzPH[k][ti],i,j,k);
          }
        }
      }
      for(j=yEnd;j<yLim;++j) {
        int tj = j-yEnd1;
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC010TaperzPH[k][tj][i],ayC010TaperzPH[k][tj][i],azC010TaperzPH[k][tj][i],bxC010TaperzPH[k][tj][i],byC010TaperzPH[k][tj][i],bzC010TaperzPH[k][tj][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC010TaperzPH[k][tj][i],ayC010TaperzPH[k][tj][i],azC010TaperzPH[k][tj][i],bxC010TaperzPH[k][tj][i],byC010TaperzPH[k][tj][i],bzC010TaperzPH[k][tj][i],i,j,k);
        }
        for(i=2;i<xStart2;i++) {
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmin,vz,axC010TaperzPH[k][tj][i],ayC010TaperzPH[k][tj][i],azC010TaperzPH[k][tj][i],bxC010TaperzPH[k][tj][i],byC010TaperzPH[k][tj][i],bzC010TaperzPH[k][tj][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmin,vz,axE021TaperzPH[k][tj],ayE021TaperzPH[k][tj],azE021TaperzPH[k][tj],bxE021TaperzPH[k][tj],byE021TaperzPH[k][tj],bzE021TaperzPH[k][tj],i,j,k);
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmin,vz,axC110TaperzPH[k][tj][ti],ayC110TaperzPH[k][tj][ti],azC110TaperzPH[k][tj][ti],bxC110TaperzPH[k][tj][ti],byC110TaperzPH[k][tj][ti],bzC110TaperzPH[k][tj][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC110TaperzPH[k][tj][ti],ayC110TaperzPH[k][tj][ti],azC110TaperzPH[k][tj][ti],bxC110TaperzPH[k][tj][ti],byC110TaperzPH[k][tj][ti],bzC110TaperzPH[k][tj][ti],i,j,k);
          }
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          int tj = j-yEnd1;
          for(i=2;i<xStart2;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC010TaperzPH[k][tj][i],ayC010TaperzPH[k][tj][i],azC010TaperzPH[k][tj][i],bxC010TaperzPH[k][tj][i],byC010TaperzPH[k][tj][i],bzC010TaperzPH[k][tj][i],i,j,k);
          }
          for(i=xStart2;i<xEnd;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axE021TaperzPH[k][tj],ayE021TaperzPH[k][tj],azE021TaperzPH[k][tj],bxE021TaperzPH[k][tj],byE021TaperzPH[k][tj],bzE021TaperzPH[k][tj],i,j,k);
          }
          for(i=xEnd;i<xLim;++i) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,axC110TaperzPH[k][tj][ti],ayC110TaperzPH[k][tj][ti],azC110TaperzPH[k][tj][ti],bxC110TaperzPH[k][tj][ti],byC110TaperzPH[k][tj][ti],bzC110TaperzPH[k][tj][ti],i,j,k);
          }
        }
      }
    }
  }
  kk = zStart_z;
  for(k=zStart_z;k<zEnd1;++k) {
    jj=0;
    nxy = NX*_nYmin;
    if(minYb) {
      if(fminYb) {
        j=0;
        for(i=2;i<xStart2;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmin,vz,axE110Taper[j][i],ayE110Taper[j][i],azE110Taper[j][i],bxE110Taper[j][i],byE110Taper[j][i],bzE110Taper[j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmin,vz,axY0Taper[j],ayTaper[j],azY0Taper[j],bxY0Taper[j],byTaper[j],bzY0Taper[j],i,j,k);
        }
        for(i=xEnd;i<xLim;++i) {
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmin,vz,axE210Taper[j][ti],ayE210Taper[j][ti],azE210Taper[j][ti],bxE210Taper[j][ti],byE210Taper[j][ti],bzE210Taper[j][ti],i,j,k);
        }
        j=1;
        for(i=2;i<xStart2;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmin,vz,axE110Taper[j][i],ayE110Taper[j][i],azE110Taper[j][i],bxE110Taper[j][i],byE110Taper[j][i],bzE110Taper[j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmin,vz,axY0Taper[j],ayTaper[j],azY0Taper[j],bxY0Taper[j],byTaper[j],bzY0Taper[j],i,j,k);
        }
        for(i=xEnd;i<xLim;++i) {
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmin,vz,axE210Taper[j][ti],ayE210Taper[j][ti],azE210Taper[j][ti],bxE210Taper[j][ti],byE210Taper[j][ti],bzE210Taper[j][ti],i,j,k);
        }
      }
      for(j=2;j<yStart2;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmin,vz,axE110Taper[j][i],ayE110Taper[j][i],azE110Taper[j][i],bxE110Taper[j][i],byE110Taper[j][i],bzE110Taper[j][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmin,vz,axE110Taper[j][i],ayE110Taper[j][i],azE110Taper[j][i],bxE110Taper[j][i],byE110Taper[j][i],bzE110Taper[j][i],i,j,k);
        }
        for(i=2;i<xStart2;i++) {
          int index = i+j*NX+k*NXY, pIndex=i+j*NX+(k-zStart_z)*nxymin;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzYmin,vz,axE110Taper[j][i],ayE110Taper[j][i],azE110Taper[j][i],bxE110Taper[j][i],byE110Taper[j][i],bzE110Taper[j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+j*NX+k*NXY, pIndex=i+j*NX+(k-zStart_z)*nxymin;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzYmin,vz,axY0Taper[j],ayTaper[j],azY0Taper[j],bxY0Taper[j],byTaper[j],bzY0Taper[j],i,j,k);
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=i+j*NX+(k-zStart_z)*nxymin;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzYmin,vz,axE210Taper[j][ti],ayE210Taper[j][ti],azE210Taper[j][ti],bxE210Taper[j][ti],byE210Taper[j][ti],bzE210Taper[j][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmin,vz,axE210Taper[j][ti],ayE210Taper[j][ti],azE210Taper[j][ti],bxE210Taper[j][ti],byE210Taper[j][ti],bzE210Taper[j][ti],i,j,k);
          }
        }
      }
    }
    jj = yStart2;
    for(j=yStart2;j<yEnd;++j) {
      if(minXb) {
        ii=0;
        nxy = _nXmin*nyInt;
        nx = _nXmin;
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzXmin,vz,axTaper[i],ayX0Taper[i],azX0Taper[i],bxTaper[i],byX0Taper[i],bzX0Taper[i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzXmin,vz,axTaper[i],ayX0Taper[i],azX0Taper[i],bxTaper[i],byX0Taper[i],bzX0Taper[i],i,j,k);
        }
        for(i=2;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY, pIndex=i+(j-jj)*nx+(k-kk)*nxminy;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzXmin,vz,axTaper[i],ayX0Taper[i],azX0Taper[i],bxTaper[i],byX0Taper[i],bzX0Taper[i],i,j,k);
        }
      }
      for(i=xStart2;i<xEnd;++i) {
        int index = i+j*NX+k*NXY;
        updateAcousticVzNoPML(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,index,vz);
      }
      if(maxXb) {
        ii=xEnd;
        nxy = _nXmax*nyInt;
        nx = _nXmax;
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=i-_xMaxStart+(j-yStart2)*_nXmax+(k-zStart_z)*nxmaxy;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzXmax,vz,axTaper[i],ayX1Taper[ti],azX1Taper[ti],bxTaper[i],byX1Taper[ti],bzX1Taper[ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzXmax,vz,axTaper[i],ayX1Taper[ti],azX1Taper[ti],bxTaper[i],byX1Taper[ti],bzX1Taper[ti],i,j,k);
          }
        }
      }
    }
    ii=0;
    nx=NX;
    nxy = NX*_nYmax;
    jj = _yMaxStart;
    if(maxYb) {
      for(j=yEnd;j<yLim;++j) {
        int tj = j-yEnd1;
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmax,vz,axE120Taper[tj][i],ayE120Taper[tj][i],azE120Taper[tj][i],bxE120Taper[tj][i],byE120Taper[tj][i],bzE120Taper[tj][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmax,vz,axE120Taper[tj][i],ayE120Taper[tj][i],azE120Taper[tj][i],bxE120Taper[tj][i],byE120Taper[tj][i],bzE120Taper[tj][i],i,j,k);
        }
        for(i=2;i<xStart2;++i) {
          int index = i+j*NX+k*NXY, pIndex=i+(j-_yMaxStart)*NX+(k-zStart_z)*nxymax;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzYmax,vz,axE120Taper[tj][i],ayE120Taper[tj][i],azE120Taper[tj][i],bxE120Taper[tj][i],byE120Taper[tj][i],bzE120Taper[tj][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;++i) {
          int index = i+j*NX+k*NXY, pIndex=i+(j-_yMaxStart)*NX+(k-zStart_z)*nxymax;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzYmax,vz,axY1Taper[tj],ayTaper[j],azY1Taper[tj],bxY1Taper[tj],byTaper[j],bzY1Taper[tj],i,j,k);
        }
        for(i=xEnd;i<xLim;++i) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=i+(j-_yMaxStart)*NX+(k-zStart_z)*nxymax;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzYmax,vz,axE220Taper[tj][ti],ayE220Taper[tj][ti],azE220Taper[tj][ti],bxE220Taper[tj][ti],byE220Taper[tj][ti],bzE220Taper[tj][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmax,vz,axE220Taper[tj][ti],ayE220Taper[tj][ti],azE220Taper[tj][ti],bxE220Taper[tj][ti],byE220Taper[tj][ti],bzE220Taper[tj][ti],i,j,k);
          }
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          int tj = j-yEnd1;
          for(i=2;i<xStart2;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmax,vz,axE120Taper[tj][i],ayE120Taper[tj][i],azE120Taper[tj][i],bxE120Taper[tj][i],byE120Taper[tj][i],bzE120Taper[tj][i],i,j,k);
          }
          for(i=xStart2;i<xEnd;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmax,vz,axY1Taper[tj],ayTaper[j],azY1Taper[tj],bxY1Taper[tj],byTaper[j],bzY1Taper[tj],i,j,k);
          }
          for(i=xEnd;i<xLim;++i) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmax,vz,axE220Taper[tj][ti],ayE220Taper[tj][ti],azE220Taper[tj][ti],bxE220Taper[tj][ti],byE220Taper[tj][ti],bzE220Taper[tj][ti],i,j,k);
          }
        }
      }
    }
  }
  jj=0;
  kk = _zMaxStart-1;
  nxy = NXY;
  if(maxZb) {
    for(k=zEnd1;k<zLim;++k) {
      int tk = k-zEnd1;
      if(fminYb) {
        j=0;
        for(i=2;i<xStart2;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC001TaperzPH[tk][j][i],ayC001TaperzPH[tk][j][i],azC001TaperzPH[tk][j][i],bxC001TaperzPH[tk][j][i],byC001TaperzPH[tk][j][i],bzC001TaperzPH[tk][j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axE012TaperzPH[tk][j],ayE012TaperzPH[tk][j],azE012TaperzPH[tk][j],bxE012TaperzPH[tk][j],byE012TaperzPH[tk][j],bzE012TaperzPH[tk][j],i,j,k);
        }
        for(i=xEnd;i<xLim;++i) {
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC101TaperzPH[tk][j][ti],ayC101TaperzPH[tk][j][ti],azC101TaperzPH[tk][j][ti],bxC101TaperzPH[tk][j][ti],byC101TaperzPH[tk][j][ti],bzC101TaperzPH[tk][j][ti],i,j,k);
        }
        j=1;
        for(i=2;i<xStart2;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC001TaperzPH[tk][j][i],ayC001TaperzPH[tk][j][i],azC001TaperzPH[tk][j][i],bxC001TaperzPH[tk][j][i],byC001TaperzPH[tk][j][i],bzC001TaperzPH[tk][j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axE012TaperzPH[tk][j],ayE012TaperzPH[tk][j],azE012TaperzPH[tk][j],bxE012TaperzPH[tk][j],byE012TaperzPH[tk][j],bzE012TaperzPH[tk][j],i,j,k);
        }
        for(i=xEnd;i<xLim;++i) {
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC101TaperzPH[tk][j][ti],ayC101TaperzPH[tk][j][ti],azC101TaperzPH[tk][j][ti],bxC101TaperzPH[tk][j][ti],byC101TaperzPH[tk][j][ti],bzC101TaperzPH[tk][j][ti],i,j,k);
        }
      }
      for(j=2;j<yStart2;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC001TaperzPH[tk][j][i],ayC001TaperzPH[tk][j][i],azC001TaperzPH[tk][j][i],bxC001TaperzPH[tk][j][i],byC001TaperzPH[tk][j][i],bzC001TaperzPH[tk][j][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC001TaperzPH[tk][j][i],ayC001TaperzPH[tk][j][i],azC001TaperzPH[tk][j][i],bxC001TaperzPH[tk][j][i],byC001TaperzPH[tk][j][i],bzC001TaperzPH[tk][j][i],i,j,k);
        }
        for(i=2;i<xStart2;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmax,vz,axC001TaperzPH[tk][j][i],ayC001TaperzPH[tk][j][i],azC001TaperzPH[tk][j][i],bxC001TaperzPH[tk][j][i],byC001TaperzPH[tk][j][i],bzC001TaperzPH[tk][j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmax,vz,axE012TaperzPH[tk][j],ayE012TaperzPH[tk][j],azE012TaperzPH[tk][j],bxE012TaperzPH[tk][j],byE012TaperzPH[tk][j],bzE012TaperzPH[tk][j],i,j,k);
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmax,vz,axC101TaperzPH[tk][j][ti],ayC101TaperzPH[tk][j][ti],azC101TaperzPH[tk][j][ti],bxC101TaperzPH[tk][j][ti],byC101TaperzPH[tk][j][ti],bzC101TaperzPH[tk][j][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC101TaperzPH[tk][j][ti],ayC101TaperzPH[tk][j][ti],azC101TaperzPH[tk][j][ti],bxC101TaperzPH[tk][j][ti],byC101TaperzPH[tk][j][ti],bzC101TaperzPH[tk][j][ti],i,j,k);
          }
        }
      }
      for(j=yStart2;j<yEnd;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axE102TaperzPH[tk][i],ayE102TaperzPH[tk][i],azE102TaperzPH[tk][i],bxE102TaperzPH[tk][i],byE102TaperzPH[tk][i],bzE102TaperzPH[tk][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axE102TaperzPH[tk][i],ayE102TaperzPH[tk][i],azE102TaperzPH[tk][i],bxE102TaperzPH[tk][i],byE102TaperzPH[tk][i],bzE102TaperzPH[tk][i],i,j,k);
        }
        for(i=2;i<xStart2;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmax,vz,axE102TaperzPH[tk][i],ayE102TaperzPH[tk][i],azE102TaperzPH[tk][i],bxE102TaperzPH[tk][i],byE102TaperzPH[tk][i],bzE102TaperzPH[tk][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmax,vz,axZ1TaperPH[tk],ayZ1TaperPH[tk],azTaperPH[k],bxZ1TaperPH[tk],byZ1TaperPH[tk],bzTaperPH[k],i,j,k);
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmax,vz,axE202TaperzPH[tk][ti],ayE202TaperzPH[tk][ti],azE202TaperzPH[tk][ti],bxE202TaperzPH[tk][ti],byE202TaperzPH[tk][ti],bzE202TaperzPH[tk][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axE202TaperzPH[tk][ti],ayE202TaperzPH[tk][ti],azE202TaperzPH[tk][ti],bxE202TaperzPH[tk][ti],byE202TaperzPH[tk][ti],bzE202TaperzPH[tk][ti],i,j,k);
          }
        }
      }
      for(j=yEnd;j<yLim;++j) {
        int tj = j-yEnd1;
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC011TaperzPH[tk][tj][i],ayC011TaperzPH[tk][tj][i],azC011TaperzPH[tk][tj][i],bxC011TaperzPH[tk][tj][i],byC011TaperzPH[tk][tj][i],bzC011TaperzPH[tk][tj][i],i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC011TaperzPH[tk][tj][i],ayC011TaperzPH[tk][tj][i],azC011TaperzPH[tk][tj][i],bxC011TaperzPH[tk][tj][i],byC011TaperzPH[tk][tj][i],bzC011TaperzPH[tk][tj][i],i,j,k);
        }
        for(i=2;i<xStart2;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmax,vz,axC011TaperzPH[tk][tj][i],ayC011TaperzPH[tk][tj][i],azC011TaperzPH[tk][tj][i],bxC011TaperzPH[tk][tj][i],byC011TaperzPH[tk][tj][i],bzC011TaperzPH[tk][tj][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmax,vz,axE022TaperzPH[tk][tj],ayE022TaperzPH[tk][tj],azE022TaperzPH[tk][tj],bxE022TaperzPH[tk][tj],byE022TaperzPH[tk][tj],bzE022TaperzPH[tk][tj],i,j,k);
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmax,vz,axC111TaperzPH[tk][tj][ti],ayC111TaperzPH[tk][tj][ti],azC111TaperzPH[tk][tj][ti],bxC111TaperzPH[tk][tj][ti],byC111TaperzPH[tk][tj][ti],bzC111TaperzPH[tk][tj][ti],i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC111TaperzPH[tk][tj][ti],ayC111TaperzPH[tk][tj][ti],azC111TaperzPH[tk][tj][ti],bxC111TaperzPH[tk][tj][ti],byC111TaperzPH[tk][tj][ti],bzC111TaperzPH[tk][tj][ti],i,j,k);
          }
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          int tj = j-yEnd1;
          for(i=2;i<xStart2;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC011TaperzPH[tk][tj][i],ayC011TaperzPH[tk][tj][i],azC011TaperzPH[tk][tj][i],bxC011TaperzPH[tk][tj][i],byC011TaperzPH[tk][tj][i],bzC011TaperzPH[tk][tj][i],i,j,k);
          }
          for(i=xStart2;i<xEnd;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axE022TaperzPH[tk][tj],ayE022TaperzPH[tk][tj],azE022TaperzPH[tk][tj],bxE022TaperzPH[tk][tj],byE022TaperzPH[tk][tj],bzE022TaperzPH[tk][tj],i,j,k);
          }
          for(i=xEnd;i<xLim;++i) {
            int ti = i-xEnd1;
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC111TaperzPH[tk][tj][ti],ayC111TaperzPH[tk][tj][ti],azC111TaperzPH[tk][tj][ti],bxC111TaperzPH[tk][tj][ti],byC111TaperzPH[tk][tj][ti],bzC111TaperzPH[tk][tj][ti],i,j,k);
          }
        }
      }
    }
    if(fmaxZb) {
      k=zLim;
      int tk = k-zEnd1;
      
      for(j=2;j<yStart2;j++){
        for(i=2;i<xStart2;i++){
          int index=i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC001TaperzPH[tk][j][i],ayC001TaperzPH[tk][j][i],azC001TaperzPH[tk][j][i],bxC001TaperzPH[tk][j][i],byC001TaperzPH[tk][j][i],bzC001TaperzPH[tk][j][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          int index=i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axE012TaperzPH[tk][j],ayE012TaperzPH[tk][j],azE012TaperzPH[tk][j],bxE012TaperzPH[tk][j],byE012TaperzPH[tk][j],bzE012TaperzPH[tk][j],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC101TaperzPH[tk][j][ti],ayC101TaperzPH[tk][j][ti],azC101TaperzPH[tk][j][ti],bxC101TaperzPH[tk][j][ti],byC101TaperzPH[tk][j][ti],bzC101TaperzPH[tk][j][ti],i,j,k);
        }
      }
      for(j=yStart2;j<yEnd;j++){
        for(i=2;i<xStart2;i++){
          int index=i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axE102TaperzPH[tk][i],ayE102TaperzPH[tk][i],azE102TaperzPH[tk][i],bxE102TaperzPH[tk][i],byE102TaperzPH[tk][i],bzE102TaperzPH[tk][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          int index=i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axZ1TaperPH[tk],ayZ1TaperPH[tk],azTaperPH[k],bxZ1TaperPH[tk],byZ1TaperPH[tk],bzTaperPH[k],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axE202TaperzPH[tk][ti],ayE202TaperzPH[tk][ti],azE202TaperzPH[tk][ti],bxE202TaperzPH[tk][ti],byE202TaperzPH[tk][ti],bzE202TaperzPH[tk][ti],i,j,k);
        }
      }
      for(j=yEnd;j<yLim;j++){
        int tj = j-yEnd1;
        for(i=2;i<xStart2;i++){
          int index=i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC011TaperzPH[tk][tj][i],ayC011TaperzPH[tk][tj][i],azC011TaperzPH[tk][tj][i],bxC011TaperzPH[tk][tj][i],byC011TaperzPH[tk][tj][i],bzC011TaperzPH[tk][tj][i],i,j,k);
        }
        for(i=xStart2;i<xEnd;i++){
          int index=i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axE022TaperzPH[tk][tj],ayE022TaperzPH[tk][tj],azE022TaperzPH[tk][tj],bxE022TaperzPH[tk][tj],byE022TaperzPH[tk][tj],bzE022TaperzPH[tk][tj],i,j,k);
        }
        for(i=xEnd;i<xLim;i++){
          int ti = i-xEnd1;
          int index=i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,axC111TaperzPH[tk][tj][ti],ayC111TaperzPH[tk][tj][ti],azC111TaperzPH[tk][tj][ti],bxC111TaperzPH[tk][tj][ti],byC111TaperzPH[tk][tj][ti],bzC111TaperzPH[tk][tj][ti],i,j,k);
        }
      }
    }
  }
}

//This is a special updating formula for vz directly above the free
//surface.  It extrapolates from vz below the free surface using
//4th order extrapolation and using dvz/dz=0 on the free surface
void updateVzPressFreeAttenMPMLBounds(modelDefStruct* modelDef,
                        float* __restrict__ vz,
                        float* __restrict__ pressure,
                        float cx[2],float cy[2],float cz[2],
                        float* __restrict__ rho,
                        unsigned char* __restrict__ vzfunc,
                        bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb)  {
  using namespace acousticAtten_MPML;
  DEF_MODEL_SIZE(modelDef);
  
  register int i,j,k;
  const float dvzC0 = 7./8., dvzC1 = 1./8., dvzC2 = -1./24., dvzC3 = -24./23.;
  int xStart2 = MAX(_nXmin,2);
  int xEnd = MIN(_xMaxStart,NX-2);
  
  //printf("vz %d %d %d, %d %d %d %d %d %d\n",NX,NY,NZ,minXb,maxXb,minYb,maxYb,minZb,maxZb);
  k=1;
  for(j=0;j<yStart2;++j) {
    for(i=0;i<xStart2;++i) {
      int index = i+j*NX+k*NXY;
      //vz[index] = vz[index+NXY];  //second order extrapolator
      vz[index] = dvzC3*(-dvzC0*vz[index+NXY]-dvzC1*vz[index+2*NXY]-dvzC2*vz[index+3*NXY]);  //4th order extrapolator
    }
    for(i=xStart2;i<xEnd;++i) {
      int index = i+j*NX+k*NXY;
      //vz[index] = vz[index+NXY];  //second order extrapolator
      vz[index] = dvzC3*(-dvzC0*vz[index+NXY]-dvzC1*vz[index+2*NXY]-dvzC2*vz[index+3*NXY]);  //4th order extrapolator
    }
    for(i=xEnd;i<NX;++i) {
      int index = i+j*NX+k*NXY;
      //vz[index] = vz[index+NXY];  //second order extrapolator
      vz[index] = dvzC3*(-dvzC0*vz[index+NXY]-dvzC1*vz[index+2*NXY]-dvzC2*vz[index+3*NXY]);  //4th order extrapolator
    }
  }
  for(j=yStart2;j<yEnd;++j) {
    for(i=0;i<xStart2;++i) {
      int index = i+j*NX+k*NXY;
      //vz[index] = vz[index+NXY];  //second order extrapolator
      vz[index] = dvzC3*(-dvzC0*vz[index+NXY]-dvzC1*vz[index+2*NXY]-dvzC2*vz[index+3*NXY]);  //4th order extrapolator
    }
    for(i=xStart2;i<xEnd;++i) {
      int index = i+j*NX+k*NXY;
      //vz[index] = vz[index+NXY];  //second order extrapolator
      vz[index] = dvzC3*(-dvzC0*vz[index+NXY]-dvzC1*vz[index+2*NXY]-dvzC2*vz[index+3*NXY]);  //4th order extrapolator
    }
    for(i=xEnd;i<NX;++i) {
      int index = i+j*NX+k*NXY;
      //vz[index] = vz[index+NXY];  //second order extrapolator
      vz[index] = dvzC3*(-dvzC0*vz[index+NXY]-dvzC1*vz[index+2*NXY]-dvzC2*vz[index+3*NXY]);  //4th order extrapolator
    }
  }
  for(j=yEnd;j<NY;++j) {
    for(i=0;i<xStart2;++i) {
      int index = i+j*NX+k*NXY;
      //vz[index] = vz[index+NXY];  //second order extrapolator
      vz[index] = dvzC3*(-dvzC0*vz[index+NXY]-dvzC1*vz[index+2*NXY]-dvzC2*vz[index+3*NXY]);  //4th order extrapolator
    }
    for(i=xStart2;i<xEnd;++i) {
      int index = i+j*NX+k*NXY;
      //vz[index] = vz[index+NXY];  //second order extrapolator
      vz[index] = dvzC3*(-dvzC0*vz[index+NXY]-dvzC1*vz[index+2*NXY]-dvzC2*vz[index+3*NXY]);  //4th order extrapolator
    }
    for(i=xEnd;i<NX;++i) {
      int index = i+j*NX+k*NXY;
      //vz[index] = vz[index+NXY];  //second order extrapolator
      vz[index] = dvzC3*(-dvzC0*vz[index+NXY]-dvzC1*vz[index+2*NXY]-dvzC2*vz[index+3*NXY]);  //4th order extrapolator
    }
  }
}

//Update the pressure.  See vx above for a description.
void updateAcousticAttenPressureMPML(modelDefStruct* modelDef,
                                   float*__restrict__ pressure,
                                float*__restrict__ vx,float*__restrict__ vy,float*__restrict__ vz,  int*__restrict__ Qindex, int nMechs,
                                float**__restrict__ decayRates, float**__restrict__ ampP,
                                float**__restrict__ rp,
                                   float cx[2],float cy[2],float cz[2],
                                   float*__restrict__ bulk,
                                   unsigned char*__restrict__ ssfunc,
                                   bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                                   bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb) {
  using namespace acousticAtten_MPML;
  DEF_MODEL_SIZE(modelDef);	   
  
//  minXb = maxXb = minYb = maxYb = minZb = maxZb = false;
//  maxZb=false;
  
  register int i,j,k;
  int zLim=NZ-2,yLim=NY-2,xLim=NX-2;
  float cx0=cx[0],cx1=cx[1];
  float cy0=cy[0],cy1=cy[1];
  float cz0=cz[0],cz1=cz[1];
  float dcz0=modelDef->dt/modelDef->dz*modelDef->scalarSpeed;
  float dcx0=dcz0, dcy0=dcz0;
  int xStart2 = MAX(_nXmin,2);
  int xEnd = MIN(_xMaxStart,NX-2);
  int kk=0, jj=0, iii=0;
  int nx = NX, nxy=NXY;
  
  if(minZb) {
    if(fminZb) {
      k=0;
      for(j=2;j<yStart2;++j) {
        int jkind = j*NX;
        //atten adds the following lines
        int ii=0;
        float* __restrict__ currRp=rp[ii];
        float* __restrict__ cRates=decayRates[ii];
        float* __restrict__ cAmpP=ampP[ii];
        for(i=2;i<xLim;++i) {
          int qind = Qindex[i+jkind];
          tcr[i] = cRates[qind];
          tap[i] = cAmpP[qind];
        }
        //atten adds the above lines
        
        //the inline function below is the only difference compared to the no pml case
        for(i=2;i<xStart2;i++) {
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*vz[index];
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axC000Taper[k][j][i],ayC000Taper[k][j][i],azC000Taper[k][j][i],
                              bxC000Taper[k][j][i],byC000Taper[k][j][i],bzC000Taper[k][j][i],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          //atten adds the lines to the next atten comment
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          //atten added the above lines plus the totalRp next line
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*vz[index];
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axE011Taper[k][j],ayE011Taper[k][j],azE011Taper[k][j],
                              bxE011Taper[k][j],byE011Taper[k][j],bzE011Taper[k][j],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*vz[index];
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axC100Taper[k][j][ti],ayC100Taper[k][j][ti],azC100Taper[k][j][ti],
                              bxC100Taper[k][j][ti],byC100Taper[k][j][ti],bzC100Taper[k][j][ti],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(ii=1;ii<nMechs;++ii) {
          currRp=rp[ii];
          cRates=decayRates[ii];
          cAmpP=ampP[ii];
          for(i=2;i<xLim;++i) {
            int qind = Qindex[i+jkind];
            tcr[i] = cRates[qind];
            tap[i] = cAmpP[qind];
          }
          for(i=2;i<xLim;i++) {
            int index = i+jkind;
            updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
          }
        }
      }
      for(j=yStart2;j<yEnd;++j) {
        int jkind = j*NX;
        int ii=0;
        float* __restrict__ currRp=rp[ii];
        float* __restrict__ cRates=decayRates[ii];
        float* __restrict__ cAmpP=ampP[ii];
        for(i=2;i<xLim;++i) {
          int qind = Qindex[i+jkind];
          tcr[i] = cRates[qind];
          tap[i] = cAmpP[qind];
        }
        //the inline function below is the only difference compared to the no pml case
        for(i=2;i<xStart2;i++) {
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*vz[index];
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axE101Taper[k][i],ayE101Taper[k][i],azE101Taper[k][i],
                              bxE101Taper[k][i],byE101Taper[k][i],bzE101Taper[k][i],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*vz[index];
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axZ0Taper[k],ayZ0Taper[k],azTaper[k],
                              bxZ0Taper[k],byZ0Taper[k],bzTaper[k],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*vz[index];
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axE201Taper[k][ti],ayE201Taper[k][ti],azE201Taper[k][ti],
                              bxE201Taper[k][ti],byE201Taper[k][ti],bzE201Taper[k][ti],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(ii=1;ii<nMechs;++ii) {
          currRp=rp[ii];
          cRates=decayRates[ii];
          cAmpP=ampP[ii];
          for(i=2;i<xLim;++i) {
            int qind = Qindex[i+jkind];
            tcr[i] = cRates[qind];
            tap[i] = cAmpP[qind];
          }
          for(i=2;i<xLim;i++) {
            int index = i+jkind;
            updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
          }
        }
      }
      for(j=yEnd;j<yLim;++j) {
        int tj = j-yEnd1;
        int jkind = j*NX;
        int ii=0;
        float* __restrict__ currRp=rp[ii];
        float* __restrict__ cRates=decayRates[ii];
        float* __restrict__ cAmpP=ampP[ii];
        for(i=2;i<xLim;++i) {
          int qind = Qindex[i+jkind];
          tcr[i] = cRates[qind];
          tap[i] = cAmpP[qind];
        }
        //the inline function below is the only difference compared to the no pml case
        for(i=2;i<xStart2;i++) {
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*vz[index];
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axC010Taper[k][tj][i],ayC010Taper[k][tj][i],azC010Taper[k][tj][i],
                              bxC010Taper[k][tj][i],byC010Taper[k][tj][i],bzC010Taper[k][tj][i],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*vz[index];
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axE021Taper[k][tj],ayE021Taper[k][tj],azE021Taper[k][tj],
                              bxE021Taper[k][tj],byE021Taper[k][tj],bzE021Taper[k][tj],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*vz[index];
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axC110Taper[k][tj][ti],ayC110Taper[k][tj][ti],azC110Taper[k][tj][ti],
                              bxC110Taper[k][tj][ti],byC110Taper[k][tj][ti],bzC110Taper[k][tj][ti],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(ii=1;ii<nMechs;++ii) {
          currRp=rp[ii];
          cRates=decayRates[ii];
          cAmpP=ampP[ii];
          for(i=2;i<xLim;++i) {
            int qind = Qindex[i+jkind];
            tcr[i] = cRates[qind];
            tap[i] = cAmpP[qind];
          }
          for(i=2;i<xLim;i++) {
            int index = i+jkind;
            updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
          }
        }
      }
      k=1;
      for(j=2;j<yStart2;++j) {
        int jkind = j*NX+NXY;
        //below function call assumes ssfunc is has been set to 2nd order for k=1
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axC000Taper[k][j],ayC000Taper[k][j],azC000Taper[k][j],bxC000Taper[k][j],byC000Taper[k][j],bzC000Taper[k][j],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axE011Taper[k][j],ayE011Taper[k][j],azE011Taper[k][j],bxE011Taper[k][j],byE011Taper[k][j],bzE011Taper[k][j],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axC100Taper[k][j],ayC100Taper[k][j],azC100Taper[k][j],bxC100Taper[k][j],byC100Taper[k][j],bzC100Taper[k][j],kyTaper[j],kzTaper[k]);
      }
      for(j=yStart2;j<yEnd;++j) {
        int jkind = j*NX+NXY;
        //below function call assumes ssfunc is has been set to 2nd order for k=1
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axE101Taper[k],ayE101Taper[k],azE101Taper[k],bxE101Taper[k],byE101Taper[k],bzE101Taper[k],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axZ0Taper[k],ayZ0Taper[k],azTaper[k],bxZ0Taper[k],byZ0Taper[k],bzTaper[k],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axE201Taper[k],ayE201Taper[k],azE201Taper[k],bxE201Taper[k],byE201Taper[k],bzE201Taper[k],kyTaper[j],kzTaper[k]);
      }
      for(j=yEnd;j<yLim;++j) {
        int tj = j-yEnd1;
        int jkind = j*NX+NXY;
        //below function call assumes ssfunc is has been set to 2nd order for k=1
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axC010Taper[k][tj],ayC010Taper[k][tj],azC010Taper[k][tj],bxC010Taper[k][tj],byC010Taper[k][tj],bzC010Taper[k][tj],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axE021Taper[k][tj],ayE021Taper[k][tj],azE021Taper[k][tj],bxE021Taper[k][tj],byE021Taper[k][tj],bzE021Taper[k][tj],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axC110Taper[k][tj],ayC110Taper[k][tj],azC110Taper[k][tj],bxC110Taper[k][tj],byC110Taper[k][tj],bzC110Taper[k][tj],kyTaper[j],kzTaper[k]);
      }
    }
    for(k=2;k<zStart;++k) {
            int kind = k*NXY;
      if(fminYb) {
        j=0;
        int jkind = kind;
        int ii=0;
        float* __restrict__ currRp=rp[ii];
        float* __restrict__ cRates=decayRates[ii];
        float* __restrict__ cAmpP=ampP[ii];
        for(i=2;i<xLim;++i) {
          int qind = Qindex[i+jkind];
          tcr[i] = cRates[qind];
          tap[i] = cAmpP[qind];
        }
        //the inline function below is the only difference compared to the no pml case
        for(i=2;i<xStart2;i++) {
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axC000Taper[k][j][i],ayC000Taper[k][j][i],azC000Taper[k][j][i],
                              bxC000Taper[k][j][i],byC000Taper[k][j][i],bzC000Taper[k][j][i],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axE011Taper[k][j],ayE011Taper[k][j],azE011Taper[k][j],
                              bxE011Taper[k][j],byE011Taper[k][j],bzE011Taper[k][j],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axC100Taper[k][j][ti],ayC100Taper[k][j][ti],azC100Taper[k][j][ti],
                              bxC100Taper[k][j][ti],byC100Taper[k][j][ti],bzC100Taper[k][j][ti],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(ii=1;ii<nMechs;++ii) {
          currRp=rp[ii];
          cRates=decayRates[ii];
          cAmpP=ampP[ii];
          for(i=2;i<xLim;++i) {
            int qind = Qindex[i+jkind];
            tcr[i] = cRates[qind];
            tap[i] = cAmpP[qind];
          }
          for(i=2;i<xLim;i++) {
            int index = i+jkind;
            updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
          }
        }
        j=1;
        jkind = NX+kind;
        //assumes ssfunc is O2 for j=1
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axC000Taper[k][j],ayC000Taper[k][j],azC000Taper[k][j],bxC000Taper[k][j],byC000Taper[k][j],bzC000Taper[k][j],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axE011Taper[k][j],ayE011Taper[k][j],azE011Taper[k][j],bxE011Taper[k][j],byE011Taper[k][j],bzE011Taper[k][j],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axC100Taper[k][j],ayC100Taper[k][j],azC100Taper[k][j],bxC100Taper[k][j],byC100Taper[k][j],bzC100Taper[k][j],kyTaper[j],kzTaper[k]);
      }
      for(j=2;j<yStart2;++j) {
        int jkind = j*NX+kind;
        if(fminXb) {
          i=0;
          int ii=0;
          float* __restrict__ currRp=rp[ii];
          float* __restrict__ cRates=decayRates[ii];
          float* __restrict__ cAmpP=ampP[ii];
          int qind = Qindex[i+jkind];
          tcr[i] = cRates[qind];
          tap[i] = cAmpP[qind];
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axC000Taper[k][j][i],ayC000Taper[k][j][i],azC000Taper[k][j][i],
                              bxC000Taper[k][j][i],byC000Taper[k][j][i],bzC000Taper[k][j][i],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
          for(ii=1;ii<nMechs;++ii) {
            currRp=rp[ii];
            cRates=decayRates[ii];
            cAmpP=ampP[ii];
            int qind = Qindex[i+jkind];
            tcr[i] = cRates[qind];
            tap[i] = cAmpP[qind];
            int index = i+jkind;
            updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
          }
          //assumes ssfunc is O2 for x=1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,1,2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axC000Taper[k][j],ayC000Taper[k][j],azC000Taper[k][j],bxC000Taper[k][j],byC000Taper[k][j],bzC000Taper[k][j],kyTaper[j],kzTaper[k]);
        }
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,jkind,pressure,axC000Taper[k][j],ayC000Taper[k][j],azC000Taper[k][j],bxC000Taper[k][j],byC000Taper[k][j],bzC000Taper[k][j],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,jkind,pressure,axE011Taper[k][j],ayE011Taper[k][j],azE011Taper[k][j],bxE011Taper[k][j],byE011Taper[k][j],bzE011Taper[k][j],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,jkind,pressure,axC100Taper[k][j],ayC100Taper[k][j],azC100Taper[k][j],bxC100Taper[k][j],byC100Taper[k][j],bzC100Taper[k][j],kyTaper[j],kzTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xLim,NX,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axC100Taper[k][j],ayC100Taper[k][j],azC100Taper[k][j],bxC100Taper[k][j],byC100Taper[k][j],bzC100Taper[k][j],kyTaper[j],kzTaper[k]);
        }
      }
      for(j=yStart2;j<yEnd;++j) {
        int jkind = j*NX+kind;
        if(fminXb) {
          i=0;
          int ii=0;
          float* __restrict__ currRp=rp[ii];
          float* __restrict__ cRates=decayRates[ii];
          float* __restrict__ cAmpP=ampP[ii];
          int qind = Qindex[i+jkind];
          tcr[i] = cRates[qind];
          tap[i] = cAmpP[qind];
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axE101Taper[k][i],ayE101Taper[k][i],azE101Taper[k][i],
                              bxE101Taper[k][i],byE101Taper[k][i],bzE101Taper[k][i],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
          for(ii=1;ii<nMechs;++ii) {
            currRp=rp[ii];
            cRates=decayRates[ii];
            cAmpP=ampP[ii];
            int qind = Qindex[i+jkind];
            tcr[i] = cRates[qind];
            tap[i] = cAmpP[qind];
            int index = i+jkind;
            updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
          }
          //assumes ssfunc is O2 for x=1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,1,2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axE101Taper[k],ayE101Taper[k],azE101Taper[k],bxE101Taper[k],byE101Taper[k],bzE101Taper[k],kyTaper[j],kzTaper[k]);
        }
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,jkind,pressure,axE101Taper[k],ayE101Taper[k],azE101Taper[k],bxE101Taper[k],byE101Taper[k],bzE101Taper[k],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,jkind,pressure,axZ0Taper[k],ayZ0Taper[k],azTaper[k],bxZ0Taper[k],byZ0Taper[k],bzTaper[k],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,jkind,pressure,axE201Taper[k],ayE201Taper[k],azE201Taper[k],bxE201Taper[k],byE201Taper[k],bzE201Taper[k],kyTaper[j],kzTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xLim,NX,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axE201Taper[k],ayE201Taper[k],azE201Taper[k],bxE201Taper[k],byE201Taper[k],bzE201Taper[k],kyTaper[j],kzTaper[k]);
        }
      }
      for(j=yEnd;j<yLim;++j) {
        int tj = j-yEnd1;
        int jkind = j*NX+kind;
        if(fminXb) {
          i=0;
          int ii=0;
          float* __restrict__ currRp=rp[ii];
          float* __restrict__ cRates=decayRates[ii];
          float* __restrict__ cAmpP=ampP[ii];
          int qind = Qindex[i+jkind];
          tcr[i] = cRates[qind];
          tap[i] = cAmpP[qind];
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axC010Taper[k][tj][i],ayC010Taper[k][tj][i],azC010Taper[k][tj][i],
                              bxC010Taper[k][tj][i],byC010Taper[k][tj][i],bzC010Taper[k][tj][i],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
          for(ii=1;ii<nMechs;++ii) {
            currRp=rp[ii];
            cRates=decayRates[ii];
            cAmpP=ampP[ii];
            int qind = Qindex[i+jkind];
            tcr[i] = cRates[qind];
            tap[i] = cAmpP[qind];
            int index = i+jkind;
            updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
          }
          //assumes ssfunc is O2 for x=1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,1,2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axC010Taper[k][tj],ayC010Taper[k][tj],azC010Taper[k][tj],bxC010Taper[k][tj],byC010Taper[k][tj],bzC010Taper[k][tj],kyTaper[j],kzTaper[k]);
        }
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,jkind,pressure,axC010Taper[k][tj],ayC010Taper[k][tj],azC010Taper[k][tj],bxC010Taper[k][tj],byC010Taper[k][tj],bzC010Taper[k][tj],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,jkind,pressure,axE021Taper[k][tj],ayE021Taper[k][tj],azE021Taper[k][tj],bxE021Taper[k][tj],byE021Taper[k][tj],bzE021Taper[k][tj],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,jkind,pressure,axC110Taper[k][tj],ayC110Taper[k][tj],azC110Taper[k][tj],bxC110Taper[k][tj],byC110Taper[k][tj],bzC110Taper[k][tj],kyTaper[j],kzTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xLim,NX,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axC110Taper[k][tj],ayC110Taper[k][tj],azC110Taper[k][tj],bxC110Taper[k][tj],byC110Taper[k][tj],bzC110Taper[k][tj],kyTaper[j],kzTaper[k]);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          int tj = j-yEnd1;
          int jkind = j*NX+kind;
          //assumes ssfunc is O2 for y=yLim and NY-1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axC010Taper[k][tj],ayC010Taper[k][tj],azC010Taper[k][tj],bxC010Taper[k][tj],byC010Taper[k][tj],bzC010Taper[k][tj],kyTaper[j],kzTaper[k]);
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axE021Taper[k][tj],ayE021Taper[k][tj],azE021Taper[k][tj],bxE021Taper[k][tj],byE021Taper[k][tj],bzE021Taper[k][tj],kyTaper[j],kzTaper[k]);
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,axC110Taper[k][tj],ayC110Taper[k][tj],azC110Taper[k][tj],bxC110Taper[k][tj],byC110Taper[k][tj],bzC110Taper[k][tj],kyTaper[j],kzTaper[k]);
        }
      }
    }
  }    
  kk = zStart;
  for(k=zStart;k<zEnd;++k) {
        int kind = k*NXY;
    jj=0;
    nxy = NX*_nYmin;
    if(minYb) {
      if(fminYb) {
        j=0;
        int jkind = kind;
        int ii=0;
        float* __restrict__ currRp=rp[ii];
        float* __restrict__ cRates=decayRates[ii];
        float* __restrict__ cAmpP=ampP[ii];
        for(i=2;i<xLim;++i) {
          int qind = Qindex[i+jkind];
          tcr[i] = cRates[qind];
          tap[i] = cAmpP[qind];
        }
        //the inline function below is the only difference compared to the no pml case
        for(i=2;i<xStart2;i++) {
          int index = i+jkind;
          int pIndex = i+j*NX+(k-kk)*nxy;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxYmin,vyyYmin,vzzYmin,pIndex,axE110Taper[j][i],ayE110Taper[j][i],azE110Taper[j][i],
                              bxE110Taper[j][i],byE110Taper[j][i],bzE110Taper[j][i],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxYmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyYmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzYmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+jkind;
          int pIndex = i+j*NX+(k-kk)*nxy;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxYmin,vyyYmin,vzzYmin,pIndex,axY0Taper[j],ayTaper[j],azY0Taper[j],
                              bxY0Taper[j],byTaper[j],bzY0Taper[j],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxYmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyYmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzYmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+jkind;
          int pIndex = i+j*NX+(k-kk)*nxy;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxYmin,vyyYmin,vzzYmin,pIndex,axE210Taper[j][ti],ayE210Taper[j][ti],azE210Taper[j][ti],
                              bxE210Taper[j][ti],byE210Taper[j][ti],bzE210Taper[j][ti],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxYmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyYmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzYmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(ii=1;ii<nMechs;++ii) {
          currRp=rp[ii];
          cRates=decayRates[ii];
          cAmpP=ampP[ii];
          for(i=2;i<xLim;++i) {
            int qind = Qindex[i+jkind];
            tcr[i] = cRates[qind];
            tap[i] = cAmpP[qind];
          }
          for(i=2;i<xLim;i++) {
            int index = i+jkind;
            updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
          }
        }
        j=1;
        jkind = NX+kind;
        //assumes ssfunc is O2 for j=1
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxYmin,vyyYmin,vzzYmin,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,j*NX+(k-kk)*nxy,pressure,axE110Taper[j],ayE110Taper[j],azE110Taper[j],bxE110Taper[j],byE110Taper[j],bzE110Taper[j],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxYmin,vyyYmin,vzzYmin,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,j*NX+(k-kk)*nxy,pressure,axY0Taper[j],ayTaper[j],azY0Taper[j],bxY0Taper[j],byTaper[j],bzY0Taper[j],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxYmin,vyyYmin,vzzYmin,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,j*NX+(k-kk)*nxy,pressure,axE210Taper[j],ayE210Taper[j],azE210Taper[j],bxE210Taper[j],byE210Taper[j],bzE210Taper[j],kyTaper[j],kzTaper[k]);
      }
      for(j=2;j<yStart2;++j) {
        int jkind = j*NX+kind;
        if(fminXb) {
          i=0;
          int ii=0;
          float* __restrict__ currRp=rp[ii];
          float* __restrict__ cRates=decayRates[ii];
          float* __restrict__ cAmpP=ampP[ii];
          int qind = Qindex[i+jkind];
          tcr[i] = cRates[qind];
          tap[i] = cAmpP[qind];
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          int pIndex = i+j*NX+(k-kk)*nxy;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxYmin,vyyYmin,vzzYmin,pIndex,axE110Taper[j][i],ayE110Taper[j][i],azE110Taper[j][i],
                              bxE110Taper[j][i],byE110Taper[j][i],bzE110Taper[j][i],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxYmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyYmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzYmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
          //assumes ssfunc is O2 for x=1
          for(ii=1;ii<nMechs;++ii) {
            currRp=rp[ii];
            cRates=decayRates[ii];
            cAmpP=ampP[ii];
            int qind = Qindex[i+jkind];
            tcr[i] = cRates[qind];
            tap[i] = cAmpP[qind];
            int index = i+jkind;
            updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
          }
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxYmin,vyyYmin,vzzYmin,rp,decayRates,ampP,Qindex,nMechs,1,2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,j*NX+(k-kk)*nxy,pressure,axE110Taper[j],ayE110Taper[j],azE110Taper[j],bxE110Taper[j],byE110Taper[j],bzE110Taper[j],kyTaper[j],kzTaper[k]);
        }
        updateAcousticAttenPressurePML(vx,vy,vz,vxxYmin,vyyYmin,vzzYmin,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,j*NX+(k-kk)*nxy,pressure,axE110Taper[j],ayE110Taper[j],azE110Taper[j],bxE110Taper[j],byE110Taper[j],bzE110Taper[j],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxYmin,vyyYmin,vzzYmin,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,j*NX+(k-kk)*nxy,pressure,axY0Taper[j],ayTaper[j],azY0Taper[j],bxY0Taper[j],byTaper[j],bzY0Taper[j],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxYmin,vyyYmin,vzzYmin,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,j*NX+(k-kk)*nxy,pressure,axE210Taper[j],ayE210Taper[j],azE210Taper[j],bxE210Taper[j],byE210Taper[j],bzE210Taper[j],kyTaper[j],kzTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxYmin,vyyYmin,vzzYmin,rp,decayRates,ampP,Qindex,nMechs,xLim,NX,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,j*NX+(k-kk)*nxy,pressure,axE210Taper[j],ayE210Taper[j],azE210Taper[j],bxE210Taper[j],byE210Taper[j],bzE210Taper[j],kyTaper[j],kzTaper[k]);
        }
      }
    }
    jj = yStart2;
    for(j=yStart2;j<yEnd;++j) {
      int jkind = j*NX+kind;
      if(minXb) {
        nxy = _nXmin*nyInt;
        nx = _nXmin;
        if(fminXb) {
          i=0;
          int ii=0;
          float* __restrict__ currRp=rp[ii];
          float* __restrict__ cRates=decayRates[ii];
          float* __restrict__ cAmpP=ampP[ii];
          int qind = Qindex[i+jkind];
          tcr[i] = cRates[qind];
          tap[i] = cAmpP[qind];
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxXmin,vyyXmin,vzzXmin,pIndex,axTaper[i],ayX0Taper[i],azX0Taper[i],
                              bxTaper[i],byX0Taper[i],bzX0Taper[i],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxXmin[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyXmin[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzXmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
          for(ii=1;ii<nMechs;++ii) {
            currRp=rp[ii];
            cRates=decayRates[ii];
            cAmpP=ampP[ii];
            int qind = Qindex[i+jkind];
            tcr[i] = cRates[qind];
            tap[i] = cAmpP[qind];
            int index = i+jkind;
            updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
          }
          //assumes ssfunc is O2 for x=1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxXmin,vyyXmin,vzzXmin,rp,decayRates,ampP,Qindex,nMechs,1,2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,(j-jj)*nx+(k-kk)*nxy,pressure,axTaper,ayX0Taper,azX0Taper,bxTaper,byX0Taper,bzX0Taper,kyTaper[j],kzTaper[k]);
        }
        updateAcousticAttenPressurePML(vx,vy,vz,vxxXmin,vyyXmin,vzzXmin,rp,decayRates,ampP,Qindex,nMechs,2,_nXmin,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,(j-jj)*nx+(k-kk)*nxy,pressure,axTaper,ayX0Taper,azX0Taper,bxTaper,byX0Taper,bzX0Taper,kyTaper[j],kzTaper[k]);
      }
      updateAcousticAttenPressureNoPML(vx,vy,vz,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pressure);
      if(maxXb) {
        iii=_xMaxStart;
        nxy = _nXmax*nyInt;
        nx = _nXmax;
        //assumes ssfunc is O2 for x=xLim and NX-1
        updateAcousticAttenPressurePML(vx,vy,vz,vxxXmax,vyyXmax,vzzXmax,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,(j-jj)*nx+(k-kk)*nxy-iii,pressure,axTaper+xEnd1,ayX1Taper,azX1Taper,bxTaper+xEnd1,byX1Taper,bzX1Taper,kyTaper[j],kzTaper[k]);
        if(fmaxXb) {
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxXmax,vyyXmax,vzzXmax,rp,decayRates,ampP,Qindex,nMechs,xLim,NX,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,(j-jj)*nx+(k-kk)*nxy-iii,pressure,axTaper+xEnd1,ayX1Taper,azX1Taper,bxTaper+xEnd1,byX1Taper,bzX1Taper,kyTaper[j],kzTaper[k]);
        }
      }
    }
    iii=0;
    jj=_yMaxStart;
    nxy = NX*_nYmax;
    nx = NX;
    if(maxYb) {
      for(j=yEnd;j<yLim;++j) {
        int tj = j-yEnd1;
        int jkind = j*NX+kind;
        if(fminXb) {
          i=0;
          int ii=0;
          float* __restrict__ currRp=rp[ii];
          float* __restrict__ cRates=decayRates[ii];
          float* __restrict__ cAmpP=ampP[ii];
          int qind = Qindex[i+jkind];
          tcr[i] = cRates[qind];
          tap[i] = cAmpP[qind];
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxYmax,vyyYmax,vzzYmax,pIndex,axE120Taper[tj][i],ayE120Taper[tj][i],azE120Taper[tj][i],
                              bxE120Taper[tj][i],byE120Taper[tj][i],bzE120Taper[tj][i],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxYmax[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyYmax[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzYmax[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
          for(ii=1;ii<nMechs;++ii) {
            currRp=rp[ii];
            cRates=decayRates[ii];
            cAmpP=ampP[ii];
            int qind = Qindex[i+jkind];
            tcr[i] = cRates[qind];
            tap[i] = cAmpP[qind];
            int index = i+jkind;
            updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
          }
          //assumes ssfunc is O2 for x=1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxYmax,vyyYmax,vzzYmax,rp,decayRates,ampP,Qindex,nMechs,1,2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,(j-jj)*NX+(k-kk)*nxy,pressure,axE120Taper[tj],ayE120Taper[tj],azE120Taper[tj],bxE120Taper[tj],byE120Taper[tj],bzE120Taper[tj],kyTaper[j],kzTaper[k]);
        }
        updateAcousticAttenPressurePML(vx,vy,vz,vxxYmax,vyyYmax,vzzYmax,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,(j-jj)*NX+(k-kk)*nxy,pressure,axE120Taper[tj],ayE120Taper[tj],azE120Taper[tj],bxE120Taper[tj],byE120Taper[tj],bzE120Taper[tj],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxYmax,vyyYmax,vzzYmax,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,(j-jj)*NX+(k-kk)*nxy,pressure,axY1Taper[tj],ayTaper[j],azY1Taper[tj],bxY1Taper[tj],byTaper[j],bzY1Taper[tj],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxYmax,vyyYmax,vzzYmax,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,(j-jj)*NX+(k-kk)*nxy,pressure,axE220Taper[tj],ayE220Taper[tj],azE220Taper[tj],bxE220Taper[tj],byE220Taper[tj],bzE220Taper[tj],kyTaper[j],kzTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxYmax,vyyYmax,vzzYmax,rp,decayRates,ampP,Qindex,nMechs,xLim,NX,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,(j-jj)*NX+(k-kk)*nxy,pressure,axE220Taper[tj],ayE220Taper[tj],azE220Taper[tj],bxE220Taper[tj],byE220Taper[tj],bzE220Taper[tj],kyTaper[j],kzTaper[k]);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          int tj = j-yEnd1;
          int jkind = j*NX+kind;
          //assumes ssfunc is O2 for y=yLim and NY-1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxYmax,vyyYmax,vzzYmax,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,(j-jj)*NX+(k-kk)*nxy,pressure,axE120Taper[tj],ayE120Taper[tj],azE120Taper[tj],bxE120Taper[tj],byE120Taper[tj],bzE120Taper[tj],kyTaper[j],kzTaper[k]);
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxYmax,vyyYmax,vzzYmax,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,(j-jj)*NX+(k-kk)*nxy,pressure,axY1Taper[tj],ayTaper[j],azY1Taper[tj],bxY1Taper[tj],byTaper[j],bzY1Taper[tj],kyTaper[j],kzTaper[k]);
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxYmax,vyyYmax,vzzYmax,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,(j-jj)*NX+(k-kk)*nxy,pressure,axE220Taper[tj],ayE220Taper[tj],azE220Taper[tj],bxE220Taper[tj],byE220Taper[tj],bzE220Taper[tj],kyTaper[j],kzTaper[k]);
        }
      }
    }
  }
  if(maxZb) {
    jj=0;
    kk = _zMaxStart;
    nxy = NXY;
    for(k=zEnd;k<zLim;++k) {
      int tk = k-zEnd1;
      int kind = k*NXY;
      int pkind = (k-kk)*NXY;
      if(fminYb) {
        j=0;
        int jkind = kind;
        int ii=0;
        float* __restrict__ currRp=rp[ii];
        float* __restrict__ cRates=decayRates[ii];
        float* __restrict__ cAmpP=ampP[ii];
        for(i=2;i<xLim;++i) {
          int qind = Qindex[i+jkind];
          tcr[i] = cRates[qind];
          tap[i] = cAmpP[qind];
        }
        int pjkind = pkind;
        //the inline function below is the only difference compared to the no pml case
        for(i=2;i<xStart2;i++) {
          int index = i+jkind;
          int pIndex = i+pjkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmax,vyyZmax,vzzZmax,pIndex,axC001Taper[tk][j][i],ayC001Taper[tk][j][i],azC001Taper[tk][j][i],
                              bxC001Taper[tk][j][i],byC001Taper[tk][j][i],bzC001Taper[tk][j][i],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmax[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmax[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmax[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(i=xStart2;i<xEnd;i++) {
          int index = i+jkind;
          int pIndex = i+pjkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmax,vyyZmax,vzzZmax,pIndex,axE012Taper[tk][j],ayE012Taper[tk][j],azE012Taper[tk][j],
                              bxE012Taper[tk][j],byE012Taper[tk][j],bzE012Taper[tk][j],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmax[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmax[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmax[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(i=xEnd;i<xLim;i++) {
          int ti = i-xEnd1;
          int index = i+jkind;
          int pIndex = i+pjkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmax,vyyZmax,vzzZmax,pIndex,axC101Taper[tk][j][ti],ayC101Taper[tk][j][ti],azC101Taper[tk][j][ti],
                              bxC101Taper[tk][j][ti],byC101Taper[tk][j][ti],bzC101Taper[tk][j][ti],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmax[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmax[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmax[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
        }
        for(ii=1;ii<nMechs;++ii) {
          currRp=rp[ii];
          cRates=decayRates[ii];
          cAmpP=ampP[ii];
          for(i=2;i<xLim;++i) {
            int qind = Qindex[i+jkind];
            tcr[i] = cRates[qind];
            tap[i] = cAmpP[qind];
          }
          for(i=2;i<xLim;i++) {
            int index = i+jkind;
            updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
          }
        }
        j=1;
        jkind = NX+kind;
        pjkind = NX+pkind;
        //assumes ssfunc is O2 for j=1
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axC001Taper[tk][j],ayC001Taper[tk][j],azC001Taper[tk][j],bxC001Taper[tk][j],byC001Taper[tk][j],bzC001Taper[tk][j],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axE012Taper[tk][j],ayE012Taper[tk][j],azE012Taper[tk][j],bxE012Taper[tk][j],byE012Taper[tk][j],bzE012Taper[tk][j],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axC101Taper[tk][j],ayC101Taper[tk][j],azC101Taper[tk][j],bxC101Taper[tk][j],byC101Taper[tk][j],bzC101Taper[tk][j],kyTaper[j],kzTaper[k]);
      }
      for(j=2;j<yStart2;++j) {
        int jkind = j*NX+kind;
        int pjkind = j*NX+pkind;
        if(fminXb) {
          i=0;
          int ii=0;
          float* __restrict__ currRp=rp[ii];
          float* __restrict__ cRates=decayRates[ii];
          float* __restrict__ cAmpP=ampP[ii];
          int qind = Qindex[i+jkind];
          tcr[i] = cRates[qind];
          tap[i] = cAmpP[qind];
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          int pIndex = i+pjkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmax,vyyZmax,vzzZmax,pIndex,axC001Taper[tk][j][i],ayC001Taper[tk][j][i],azC001Taper[tk][j][i],
                              bxC001Taper[tk][j][i],byC001Taper[tk][j][i],bzC001Taper[tk][j][i],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmax[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmax[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmax[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
          for(ii=1;ii<nMechs;++ii) {
            currRp=rp[ii];
            cRates=decayRates[ii];
            cAmpP=ampP[ii];
            int qind = Qindex[i+jkind];
            tcr[i] = cRates[qind];
            tap[i] = cAmpP[qind];
            int index = i+jkind;
            updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
          }
          //assumes ssfunc is O2 for x=1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,1,2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axC001Taper[tk][j],ayC001Taper[tk][j],azC001Taper[tk][j],bxC001Taper[tk][j],byC001Taper[tk][j],bzC001Taper[tk][j],kyTaper[j],kzTaper[k]);
        }
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pjkind,pressure,axC001Taper[tk][j],ayC001Taper[tk][j],azC001Taper[tk][j],bxC001Taper[tk][j],byC001Taper[tk][j],bzC001Taper[tk][j],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pjkind,pressure,axE012Taper[tk][j],ayE012Taper[tk][j],azE012Taper[tk][j],bxE012Taper[tk][j],byE012Taper[tk][j],bzE012Taper[tk][j],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pjkind,pressure,axC101Taper[tk][j],ayC101Taper[tk][j],azC101Taper[tk][j],bxC101Taper[tk][j],byC101Taper[tk][j],bzC101Taper[tk][j],kyTaper[j],kzTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xLim,NX,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axC101Taper[tk][j],ayC101Taper[tk][j],azC101Taper[tk][j],bxC101Taper[tk][j],byC101Taper[tk][j],bzC101Taper[tk][j],kyTaper[j],kzTaper[k]);
        }
      }
      for(j=yStart2;j<yEnd;++j) {
        int jkind = j*NX+kind;
        int pjkind = j*NX+pkind;
        if(fminXb) {
          i=0;
          int ii=0;
          float* __restrict__ currRp=rp[ii];
          float* __restrict__ cRates=decayRates[ii];
          float* __restrict__ cAmpP=ampP[ii];
          int qind = Qindex[i+jkind];
          tcr[i] = cRates[qind];
          tap[i] = cAmpP[qind];
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          int pIndex = i+pjkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmax,vyyZmax,vzzZmax,pIndex,axE102Taper[tk][i],ayE102Taper[tk][i],azE102Taper[tk][i],
                              bxE102Taper[tk][i],byE102Taper[tk][i],bzE102Taper[tk][i],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmax[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmax[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmax[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
          for(ii=1;ii<nMechs;++ii) {
            currRp=rp[ii];
            cRates=decayRates[ii];
            cAmpP=ampP[ii];
            int qind = Qindex[i+jkind];
            tcr[i] = cRates[qind];
            tap[i] = cAmpP[qind];
            int index = i+jkind;
            updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
          }
          //assumes ssfunc is O2 for x=1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,1,2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axE102Taper[tk],ayE102Taper[tk],azE102Taper[tk],bxE102Taper[tk],byE102Taper[tk],bzE102Taper[tk],kyTaper[j],kzTaper[k]);
        }
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pjkind,pressure,axE102Taper[tk],ayE102Taper[tk],azE102Taper[tk],bxE102Taper[tk],byE102Taper[tk],bzE102Taper[tk],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pjkind,pressure,axZ1Taper[tk],ayZ1Taper[tk],azTaper[k],bxZ1Taper[tk],byZ1Taper[tk],bzTaper[k],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pjkind,pressure,axE202Taper[tk],ayE202Taper[tk],azE202Taper[tk],bxE202Taper[tk],byE202Taper[tk],bzE202Taper[tk],kyTaper[j],kzTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xLim,NX,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axE202Taper[tk],ayE202Taper[tk],azE202Taper[tk],bxE202Taper[tk],byE202Taper[tk],bzE202Taper[tk],kyTaper[j],kzTaper[k]);
        }
      }
      for(j=yEnd;j<yLim;++j) {
        int tj = j-yEnd1;
        int jkind = j*NX+kind;
        int pjkind = j*NX+pkind;
        if(fminXb) {
          i=0;
          int ii=0;
          float* __restrict__ currRp=rp[ii];
          float* __restrict__ cRates=decayRates[ii];
          float* __restrict__ cAmpP=ampP[ii];
          int qind = Qindex[i+jkind];
          tcr[i] = cRates[qind];
          tap[i] = cAmpP[qind];
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          int pIndex = i+pjkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmax,vyyZmax,vzzZmax,pIndex,axC011Taper[tk][tj][i],ayC011Taper[tk][tj][i],azC011Taper[tk][tj][i],
                              bxC011Taper[tk][tj][i],byC011Taper[tk][tj][i],bzC011Taper[tk][tj][i],

                              dvx,dvy,dvz);
          tdvx[i] = dvx = dvx*kxTaper[i]+vxxZmax[pIndex];
          tdvy[i] = dvy = dvy*kyTaper[j]+vyyZmax[pIndex];
          tdvz[i] = dvz = dvz*kzTaper[k]+vzzZmax[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          float omega=tcr[i];
#if USE_ALT_MEM_VAR
          float ccAmpP=tap[i]*bb;
#else
          float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
          
          float tmoo2po = (2.0f-omega)/(2.0f+omega);
          float to2po = 2.0f/(2.0f+omega);
          float dVa = to2po*ccAmpP*dvTot;
          
          float prevRp = currRp[index];
          currRp[index]=tmoo2po*prevRp-dVa;
          
          //update the totals for the pressure update
          float totalRp = omega*(currRp[index]+prevRp)*0.5f;
          
#if !USE_ALT_MEM_VAR
          totalRp*=bb;
#endif
          pressure[index]-=dvTot*bb+totalRp;
          for(ii=1;ii<nMechs;++ii) {
            currRp=rp[ii];
            cRates=decayRates[ii];
            cAmpP=ampP[ii];
            int qind = Qindex[i+jkind];
            tcr[i] = cRates[qind];
            tap[i] = cAmpP[qind];
            int index = i+jkind;
            updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
          }
          //assumes ssfunc is O2 for x=1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,1,2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axC011Taper[tk][tj],ayC011Taper[tk][tj],azC011Taper[tk][tj],bxC011Taper[tk][tj],byC011Taper[tk][tj],bzC011Taper[tk][tj],kyTaper[j],kzTaper[k]);
        }
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pjkind,pressure,axC011Taper[tk][tj],ayC011Taper[tk][tj],azC011Taper[tk][tj],bxC011Taper[tk][tj],byC011Taper[tk][tj],bzC011Taper[tk][tj],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pjkind,pressure,axE022Taper[tk][tj],ayE022Taper[tk][tj],azE022Taper[tk][tj],bxE022Taper[tk][tj],byE022Taper[tk][tj],bzE022Taper[tk][tj],kyTaper[j],kzTaper[k]);
        updateAcousticAttenPressurePML(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pjkind,pressure,axC111Taper[tk][tj],ayC111Taper[tk][tj],azC111Taper[tk][tj],bxC111Taper[tk][tj],byC111Taper[tk][tj],bzC111Taper[tk][tj],kyTaper[j],kzTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xLim,NX,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axC111Taper[tk][tj],ayC111Taper[tk][tj],azC111Taper[tk][tj],bxC111Taper[tk][tj],byC111Taper[tk][tj],bzC111Taper[tk][tj],kyTaper[j],kzTaper[k]);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          int tj = j-yEnd1;
          int jkind = j*NX+kind;
          int pjkind = j*NX+pkind;
          //assumes ssfunc is O2 for y=yLim and NY-1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axC011Taper[tk][tj],ayC011Taper[tk][tj],azC011Taper[tk][tj],bxC011Taper[tk][tj],byC011Taper[tk][tj],bzC011Taper[tk][tj],kyTaper[j],kzTaper[k]);
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axE022Taper[tk][tj],ayE022Taper[tk][tj],azE022Taper[tk][tj],bxE022Taper[tk][tj],byE022Taper[tk][tj],bzE022Taper[tk][tj],kyTaper[j],kzTaper[k]);
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axC111Taper[tk][tj],ayC111Taper[tk][tj],azC111Taper[tk][tj],bxC111Taper[tk][tj],byC111Taper[tk][tj],bzC111Taper[tk][tj],kyTaper[j],kzTaper[k]);
        }
      }
    }
    if(fmaxZb) {
      for(k=zLim;k<NZ;++k) {
        int tk = k-zEnd1;
        int kind = k*NXY;
        int pkind = (k-kk)*NXY;
        for(j=2;j<yStart2;++j) {
          int jkind = j*NX+kind;
          int pjkind = j*NX+pkind;
          //assumes ssfunc is O2 for z=zLim and NZ-1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axC001Taper[tk][j],ayC001Taper[tk][j],azC001Taper[tk][j],bxC001Taper[tk][j],byC001Taper[tk][j],bzC001Taper[tk][j],kyTaper[j],kzTaper[k]);
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axE012Taper[tk][j],ayE012Taper[tk][j],azE012Taper[tk][j],bxE012Taper[tk][j],byE012Taper[tk][j],bzE012Taper[tk][j],kyTaper[j],kzTaper[k]);
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axC101Taper[tk][j],ayC101Taper[tk][j],azC101Taper[tk][j],bxC101Taper[tk][j],byC101Taper[tk][j],bzC101Taper[tk][j],kyTaper[j],kzTaper[k]);
        }
        for(j=yStart2;j<yEnd;++j) {
          int jkind = j*NX+kind;
          int pjkind = j*NX+pkind;
          //assumes ssfunc is O2 for z=zLim and NZ-1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axE102Taper[tk],ayE102Taper[tk],azE102Taper[tk],bxE102Taper[tk],byE102Taper[tk],bzE102Taper[tk],kyTaper[j],kzTaper[k]);
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axZ1Taper[tk],ayZ1Taper[tk],azTaper[k],bxZ1Taper[tk],byZ1Taper[tk],bzTaper[k],kyTaper[j],kzTaper[k]);
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axE202Taper[tk],ayE202Taper[tk],azE202Taper[tk],bxE202Taper[tk],byE202Taper[tk],bzE202Taper[tk],kyTaper[j],kzTaper[k]);
        }
        for(j=yEnd;j<yLim;++j) {
          int tj = j-yEnd1;
          int jkind = j*NX+kind;
          int pjkind = j*NX+pkind;
          //assumes ssfunc is O2 for z=zLim and NZ-1
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,2,xStart2,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axC011Taper[tk][tj],ayC011Taper[tk][tj],azC011Taper[tk][tj],bxC011Taper[tk][tj],byC011Taper[tk][tj],bzC011Taper[tk][tj],kyTaper[j],kzTaper[k]);
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xStart2,xEnd,0,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axE022Taper[tk][tj],ayE022Taper[tk][tj],azE022Taper[tk][tj],bxE022Taper[tk][tj],byE022Taper[tk][tj],bzE022Taper[tk][tj],kyTaper[j],kzTaper[k]);
          updateAcousticAttenPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,rp,decayRates,ampP,Qindex,nMechs,xEnd,xLim,xEnd1,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,axC111Taper[tk][tj],ayC111Taper[tk][tj],azC111Taper[tk][tj],bxC111Taper[tk][tj],byC111Taper[tk][tj],bzC111Taper[tk][tj],kyTaper[j],kzTaper[k]);
        }
      }
    }
  }
}

//Extrapolate the pressure one node above the free surface using 4th
//order extrapolation using pressure=0 on the free surface
void updateAcousticAttenPressurePressFreeMPML(modelDefStruct* modelDef,
                                float*__restrict__ pressure,
                                float*__restrict__ vx,float*__restrict__ vy,float*__restrict__ vz,
                                float cx[2],float cy[2],float cz[2],
                                float*__restrict__ bulk,
                                unsigned char*__restrict__ ssfunc,
                                bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb) {
  using namespace acousticAtten_MPML;
  DEF_MODEL_SIZE(modelDef);
  
  //  minXb = maxXb = minYb = maxYb = minZb = maxZb = false;
  //  maxZb=false;
  
  register int i,j,k;
  int xStart2 = MAX(_nXmin,2);
  int xEnd = MIN(_xMaxStart,NX-2);
  
  k=1;
  for(j=0;j<yStart2;++j) {
    for(i=0;i<xStart2;++i) {
      int index = i+j*NX+k*NXY;
      pressure[index] = 4.*pressure[index+3*NXY]-(6.*pressure[index+2*NXY]+pressure[index+4*NXY]);
    }
    for(i=xStart2;i<xEnd;++i) {
      int index = i+j*NX+k*NXY;
      pressure[index] = 4.*pressure[index+3*NXY]-(6.*pressure[index+2*NXY]+pressure[index+4*NXY]);
    }
    for(i=xEnd;i<NX;++i) {
      int index = i+j*NX+k*NXY;
      pressure[index] = 4.*pressure[index+3*NXY]-(6.*pressure[index+2*NXY]+pressure[index+4*NXY]);
    }
  }
  for(j=yStart2;j<yEnd;++j) {
    for(i=0;i<xStart2;++i) {
      int index = i+j*NX+k*NXY;
      pressure[index] = 4.*pressure[index+3*NXY]-(6.*pressure[index+2*NXY]+pressure[index+4*NXY]);
    }
    for(i=xStart2;i<xEnd;++i) {
      int index = i+j*NX+k*NXY;
      pressure[index] = 4.*pressure[index+3*NXY]-(6.*pressure[index+2*NXY]+pressure[index+4*NXY]);
    }
    for(i=xEnd;i<NX;++i) {
      int index = i+j*NX+k*NXY;
      pressure[index] = 4.*pressure[index+3*NXY]-(6.*pressure[index+2*NXY]+pressure[index+4*NXY]);
    }
  }
  for(j=yEnd;j<NY;++j) {
    for(i=0;i<xStart2;++i) {
      int index = i+j*NX+k*NXY;
      pressure[index] = 4.*pressure[index+3*NXY]-(6.*pressure[index+2*NXY]+pressure[index+4*NXY]);
    }
    for(i=xStart2;i<xEnd;++i) {
      int index = i+j*NX+k*NXY;
      pressure[index] = 4.*pressure[index+3*NXY]-(6.*pressure[index+2*NXY]+pressure[index+4*NXY]);
    }
    for(i=xEnd;i<NX;++i) {
      int index = i+j*NX+k*NXY;
      pressure[index] = 4.*pressure[index+3*NXY]-(6.*pressure[index+2*NXY]+pressure[index+4*NXY]);
    }
  }
}
