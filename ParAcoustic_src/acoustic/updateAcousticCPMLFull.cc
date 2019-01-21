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
 *  updateAcousticCPMLFull.cpp
 *  
 *
 *  This file contains the setup and updating functions for all the
 *  dependent variables for the case where CPMLs serve as the absorbing
 *  BCs.  It also includes the specialized functions it uses if a pressure
 *  free surface is used.  For the updating functions, updating proceeds
 *  in exactly the order that the dependent variables are stored in
 *  memory.  It gives a big speed boost compared to updating interior
 *  points then positive x flanks, then negative x flanks, etc.
 *
 *  Defines the following functions:
 *  setupAcousticCPMLBoundsFull
 *  updateVxCPMLBounds
 *  updateVyCPMLBounds
 *  updateVzCPMLBounds
 *  updateVzPressFreeCPMLBounds
 *  updateAcousticPressureCPML
 *  updateAcousticPressurePressFreeCPML
 *
 */

#include "updateAcousticCPMLFull.h"
#include "acousticCPMLNameSpace.h"
#include "sgfd.hh"
#include "constants.h"

#define ELASTIC_FS_2_0 3
#define SURFACE_PRESS_FREE 2

//  This function sets up the extra variables and arrays needed for the
//  CPML BCs and fills in the values of the CPML sigmas, kappas and
//  alphas as a function of distance into the CPML zone
void setupAcousticCPMLBoundsFull(int nXmin, float xMinVal, float xMinAVal, float xMinKVal,int nXmax, float xMaxVal, float xMaxAval, float xMaxKVal, int nYmin, float yMinVal, float yMinAVal, float yMinKVal,
                               int nYmax, float yMaxVal, float yMaxAVal, float yMaxKVal, int nZmin, float zMinVal, float zMinAVal, float zMinKVal,int nZmax, float zMaxVal, float zMaxAVal, float zMaxKVal,int kstart, int kstart_z, sgfdModel* model, int sbcmode) {
  using namespace acoustic_CPML;
  DEF_MODEL_SIZE(model->modelDef());
  DEF_MODEL_LIMITS(model->modelDef());
  DEF_PARALLEL(model->parallelDef());
  
  //define the thickness of the CPML zones for this domain
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
  
  //temporary holding variable
  del = new float[NX];
  
  //Allocate space for all the 1-D profile functions needed for the CPMLs
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
  //These parameters define the sizes for the CPML memory variables
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
  
  //Initialize the memory variables to 0
  for(int i=0;i<NXY*_nZmin;i++) {vxxZmin[i]=vyyZmin[i]=vzzZmin[i]=pxZmin[i]=pyZmin[i]=pzZmin[i]=0.;}
  for(int i=0;i<NXY*_nZmax;i++) {vxxZmax[i]=vyyZmax[i]=vzzZmax[i]=pxZmax[i]=pyZmax[i]=pzZmax[i]=0.;}
  for(int i=0;i<NX*_nYmin*nzInt;i++) {vxxYmin[i]=vyyYmin[i]=vzzYmin[i]=pxYmin[i]=pyYmin[i]=pzYmin[i]=0.;}
  for(int i=0;i<NX*_nYmax*nzInt;i++) {vxxYmax[i]=vyyYmax[i]=vzzYmax[i]=pxYmax[i]=pyYmax[i]=pzYmax[i]=0.;}
  for(int i=0;i<_nXmin*nyInt*nzInt;i++) {vxxXmin[i]=vyyXmin[i]=vzzXmin[i]=pxXmin[i]=pyXmin[i]=pzXmin[i]=0.;}
  for(int i=0;i<_nXmax*nyInt*nzInt;i++) {vxxXmax[i]=vyyXmax[i]=vzzXmax[i]=pxXmax[i]=pyXmax[i]=pzXmax[i]=0.;}
  
  //Fill in the flank profile variables
  //fromStart and fromEnd are maximum at the min side or max side, respectively.  Appropriate values are filled in when outside the CPML layers; however, these values are not used.
  //sigma and kappa tapers vary quadratically.  sigma starts at zero at the interior edge of the CPML zone and increases to the flank sigma computed above at the flank of the model.  kappa starts at one at the interior and will decrease to the flank value at the edge.
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
}

//Here are the updating equations for Vx.  These are done in exactly
//the order that vx is stored in memory starting with the minZ, minY
//minX corner proceeding to maxZ, maxY, and maxX with x fastest, then
//y.  The variables minXb, fminXb, etc are used to define where this
//domain is situated relative to the CPML zones and the absolute edges
//of the model.  minXb, maxXb, etc. will be true if the CPML zone extends
//into this domain on the specified side.  fminXb, fmaxXb, etc. will only
//be true if this domain is actually on an absolute edge of the model.
//These two sets will the identical except in the cases where the CPML
//zone is wider than the width of the domain on the absolute edge.
//4th order accurate differences are used when they can be, but near the
//flanks, second order differences are used and finally, in some cases,
//second order with the assumption of zero valued variables off grid
//are used right on the flanks.
void updateVxCPMLBounds(modelDefStruct* modelDef,
                       float* __restrict__ vx,
                       float* __restrict__ pressure,
                       float cx[2],float cy[2],float cz[2],
                       float* __restrict__ rho,unsigned char* __restrict__ vxfunc,
                       bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                       bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb) {
  using namespace acoustic_CPML;
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
  //Z min CPML zone
  if(minZb) {
    if(fminZb) {  //if this is the absolute edge of the model
      k=0;
      for(j=2;j<yLim;j++){
        for(i=1;i<xLim;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,i,j,k);
        }
      }
      k=1;
      for(j=2;j<yLim;j++){
        for(i=1;i<xLim;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,i,j,k);
        }
      }
    }
    for(k=2;k<zStart;k++) {
      int pIndex;
      if(fminYb) {
        j=0;
        for(i=1;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,i,j,k);
        }
        j=1;
        for(i=1;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,i,j,k);
        }
      }
      for(j=2;j<yLim;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,i,j,k);
        }
        //for(i=1;i<xLim;i++) {
        //  int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,j*NX+k*NXY,j*NX+k*NXY,pxZmin,vx,1,xLim,j,k);
        //}
        if (fmaxXb) {
          i=xLim;
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,i,j,k);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          for(i=1;i<xLim;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = index;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmin,vx,i,j,k);
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
        for(i=1;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmin,vx,i,j,k);
        }
        j=1;
        for(i=1;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmin,vx,i,j,k);
        }
      }
      for(j=2;j<yStart2;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmin,vx,i,j,k);
        }
        //for(i=1;i<xLim;i++) {
        //  int index = i+j*NX+k*NXY, pIndex=i+j*NX+(k-zStart)*nxymin;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,j*NX+(k-zStart)*nxymin,j*NX+k*NXY,pxYmin,vx,1,xLim,j,k);
        //}
        if (fmaxXb) {
          i=xLim;
          int index=i+j*NX+k*NXY;
          int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmin,vx,i,j,k);
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
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxXmin,vx,i,j,k);
        }
        //for(i=1;i<_nXmin;i++) {
        //  int index = i+j*NX+k*NXY, pIndex=i+(j-yStart2)*_nXmin+(k-zStart)*nxminy;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,(j-yStart2)*_nXmin+(k-zStart)*nxminy,j*NX+k*NXY,pxXmin,vx,1,_nXmin,j,k);
        //}
      }
      for(i=mxStart1;i<mxEnd1;++i) {
        int index = i+j*NX+k*NXY;
        updateAcousticVxNoPML(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,vx);
      }
      if (maxXb) {
        ii=_xMaxStart-1;
        nxy = _nXmax*nyInt;
        nx = _nXmax;
        //for(i=xEnd1;i<xLim;i++) {
        //  int index = i+j*NX+k*NXY, pIndex=i-_xMaxStart+1+(j-yStart2)*_nXmax+(k-zStart)*nxmaxy;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,-_xMaxStart+1+(j-yStart2)*_nXmax+(k-zStart)*nxmaxy,j*NX+k*NXY,pxXmax,vx,xEnd1,xLim,j,k);
        //}
        if(fmaxXb) {
          i=xLim;
          int index = i+j*NX+k*NXY, pIndex=i-_xMaxStart+1+(j-yStart2)*_nXmax+(k-zStart)*nxmaxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxXmax,vx,i,j,k);
        }
      }
    }
    ii=0;
    jj = _yMaxStart;
    nxy = NX*_nYmax;
    nx = NX;
    if(maxYb) {
      for(j=yEnd;j<yLim;++j) {
        if (fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmax,vx,i,j,k);
        }
        //for(i=1;i<xLim;i++) {
        //  int index = i+j*NX+k*NXY, pIndex=i+(j-_yMaxStart)*NX+(k-zStart)*nxymax;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,(j-_yMaxStart)*NX+(k-zStart)*nxymax,j*NX+k*NXY,pxYmax,vx,1,xLim,j,k);
        //}
        if (fmaxXb) {
          i=xLim;
          int index=i+j*NX+k*NXY;
          int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmax,vx,i,j,k);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          for(i=1;i<xLim;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxYmax,vx,i,j,k);
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
      if(fminYb) {
        j=0;
        for(i=1;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,i,j,k);
        }
        j=1;
        for(i=1;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,i,j,k);
        }
      }
      for(j=2;j<yLim;++j) {
        if (fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,i,j,k);
        }
        //for(i=1;i<xLim;i++) {
        //  int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,j*NX+(k-_zMaxStart)*NXY,j*NX+k*NXY,pxZmax,vx,1,xLim,j,k);
        //}
        if (fmaxXb) {
          i=xLim;
          int index=i+j*NX+k*NXY;
          int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,i,j,k);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          for(i=1;i<xLim;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,i,j,k);
          }
        }
      }
    }
    if(fmaxZb) {
      for(k=zLim;k<NZ;k++) {
        for(j=2;j<yLim;j++){
          for(i=1;i<xLim;i++){
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pxZmax,vx,i,j,k);
          }
        }
      }
    }
  }
}

//vy updating. see vx updating for description
void updateVyCPMLBounds(modelDefStruct* modelDef,
                       float* __restrict__ vy,
                       float* __restrict__ pressure,
                       float cx[2],float cy[2],float cz[2],
                       float* __restrict__ rho,
                       unsigned char* __restrict__ vyfunc,
                       bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                       bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb) {
  using namespace acoustic_CPML;
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
      for(j=1;j<yLim;j++){
        for(i=2;i<xLim;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,i,j,k);
        }
      }
      k=1;
      for(j=1;j<yLim;j++){
        for(i=2;i<xLim;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,i,j,k);
        }
      }
    }
    for(k=2;k<zStart;++k) {
      if(fminYb) {
        j=0;
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,i,j,k);
        }
      }
      for(j=1;j<yLim;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = index;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,i,j,k);
        }
        for(i=2;i<xLim;i++) {
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmin,vy,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = index;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,i,j,k);
          }
        }
      }  
      if(fmaxYb) {
        j=yLim;
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmin,vy,i,j,k);
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
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmin,vy,i,j,k);
        }
      }
      for(j=1;j<yStart1;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmin,vy,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmin,vy,i,j,k);
        }
        for(i=2;i<xLim;++i) {
          int index = i+j*NX+k*NXY, pIndex=i+j*NX+(k-zStart)*nxymin;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyYmin,vy,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmin,vy,i,j,k);
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
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyXmin,vy,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyXmin,vy,i,j,k);
        }
        for(i=2;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY, pIndex=i+(j-yStart1)*_nXmin+(k-zStart)*nxminy;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyXmin,vy,i,j,k);
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
          int index = i+j*NX+k*NXY, pIndex=i-_xMaxStart+(j-yStart1)*_nXmax+(k-zStart)*nxmaxy;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyXmax,vy,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyXmax,vy,i,j,k);
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
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmax,vy,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmax,vy,i,j,k);
        }
        for(i=2;i<xLim;i++) {
          int index = i+j*NX+k*NXY, pIndex=i+(j-_yMaxStart+1)*NX+(k-zStart)*nxymax;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyYmax,vy,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmax,vy,i,j,k);
          }
        }
      }
      if(fmaxYb) {
        j=yLim;
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyYmax,vy,i,j,k);
        }
      }
    }
  }
  jj=0;
  nxy=NXY;
  kk = _zMaxStart;
  if(maxZb) {
    for(k=zEnd;k<zLim;++k) {
      if(fminYb) {
        j=0;
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,i,j,k);
        }
      }
      for(j=1;j<yLim;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,i,j,k);
        }
        for(i=2;i<xLim;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,pIndex,index,pyZmax,vy,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,i,j,k);
          }
        }
      }      
      if(fmaxYb) {
        j=yLim;
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,i,j,k);
        }
      }
    }
    if(fmaxZb) {
      for(k=zLim;k<NZ;k++) {
        for(j=1;j<yLim;j++){
          for(i=2;i<xLim;i++){
            int index=i+j*NX+k*NXY, pIndex=index-NXY*_zMaxStart;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pyZmax,vy,i,j,k);
          }
        }
      }
    }
  }
}

void updateVzCPMLBounds(modelDefStruct* modelDef,
                       float* __restrict__ vz,
                       float* __restrict__ pressure,
                       float cx[2],float cy[2],float cz[2],
                       float* __restrict__ rho,
                       unsigned char* __restrict__ vzfunc,
                       bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                       bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb)  {
  using namespace acoustic_CPML;
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
      for(j=2;j<yLim;j++){
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY, pIndex=index;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,i,j,k);
        }
      }
    }
    for(k=1;k<zStart_z;++k) {
      if(fminYb) {
        j=0;
        for(i=2;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,i,j,k);
        }
        j=1;
        for(i=2;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,i,j,k);
        }
      }
      for(j=2;j<yLim;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,i,j,k);
        }
        for(i=2;i<xLim;i++) {
          int index = i+j*NX+k*NXY, pIndex = index;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmin,vz,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,i,j,k);
          }
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          for(i=2;i<xLim;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmin,vz,i,j,k);
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
        for(i=2;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmin,vz,i,j,k);
        }
        j=1;
        for(i=2;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmin,vz,i,j,k);
        }
      }
      for(j=2;j<yStart2;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmin,vz,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmin,vz,i,j,k);
        }
        for(i=2;i<xLim;i++) {
          int index = i+j*NX+k*NXY, pIndex=i+j*NX+(k-zStart_z)*nxymin;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzYmin,vz,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmin,vz,i,j,k);
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
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzXmin,vz,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzXmin,vz,i,j,k);
        }
        for(i=2;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY, pIndex=i+(j-jj)*nx+(k-kk)*nxminy;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzXmin,vz,i,j,k);
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
          int index = i+j*NX+k*NXY, pIndex=i-_xMaxStart+(j-yStart2)*_nXmax+(k-zStart_z)*nxmaxy;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzXmax,vz,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzXmax,vz,i,j,k);
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
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmax,vz,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmax,vz,i,j,k);
        }
        for(i=2;i<xLim;++i) {
          int index = i+j*NX+k*NXY, pIndex=i+(j-_yMaxStart)*NX+(k-zStart_z)*nxymax;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzYmax,vz,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmax,vz,i,j,k);
          }
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          for(i=2;i<xLim;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzYmax,vz,i,j,k);
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
      if(fminYb) {
        j=0;
        for(i=2;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,i,j,k);
        }
        j=1;
        for(i=2;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          int pIndex = i+j*NX+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,i,j,k);
        }
      }      
      for(j=2;j<yLim;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,i,j,k);
        }
        for(i=2;i<xLim;i++) {
          int index = i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,pIndex,index,pzZmax,vz,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int index=i+j*NX+k*NXY;
            int pIndex = i-ii+(j-jj)*nx+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,i,j,k);
          }
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          for(i=2;i<xLim;++i) {
            int index=i+j*NX+k*NXY;
            int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,i,j,k);
          }
        }
      }      
    }
    if(fmaxZb) {
      k=zLim;
      for(j=2;j<yLim;j++){
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY, pIndex=index-NXY*(_zMaxStart-1);
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,pIndex,index,pzZmax,vz,i,j,k);
        }
      }
    }
  }
}

//This is a special updating formula for vz directly above the free
//surface.  It extrapolates from vz below the free surface using
//4th order extrapolation and using dvz/dz=0 on the free surface
void updateVzPressFreeCPMLBounds(modelDefStruct* modelDef,
                        float* __restrict__ vz,
                        float* __restrict__ pressure,
                        float cx[2],float cy[2],float cz[2],
                        float* __restrict__ rho,
                        unsigned char* __restrict__ vzfunc,
                        bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb)  {
  using namespace acoustic_CPML;
  DEF_MODEL_SIZE(modelDef);
  
  register int i,j,k;
  const float dvzC0 = 7./8., dvzC1 = 1./8., dvzC2 = -1./24., dvzC3 = -24./23.;
  
  //printf("vz %d %d %d, %d %d %d %d %d %d\n",NX,NY,NZ,minXb,maxXb,minYb,maxYb,minZb,maxZb);
  k=1;
  for(j=0;j<NY;++j) {
    for(i=0;i<NX;++i) {
      int index = i+j*NX+k*NXY;
      //vz[index] = vz[index+NXY];  //second order extrapolator
      vz[index] = dvzC3*(-dvzC0*vz[index+NXY]-dvzC1*vz[index+2*NXY]-dvzC2*vz[index+3*NXY]);  //4th order extrapolator
    }
  }
}

//Update the pressure.  See vx above for a description.
void updateAcousticPressureCPML(modelDefStruct* modelDef,
                                   float*__restrict__ pressure,
                                   float*__restrict__ vx,float*__restrict__ vy,float*__restrict__ vz,
                                   float cx[2],float cy[2],float cz[2],
                                   float*__restrict__ bulk,
                                   unsigned char*__restrict__ ssfunc,
                                   bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                                   bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb) {
  using namespace acoustic_CPML;
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
      for(j=2;j<yLim;++j) {
        int jkind = j*NX;
        
        //the inline function below is the only difference compared to the no pml case
        for(i=2;i<xLim;i++) {
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*vz[index];
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axTaper[i],ayTaper[j],azTaper[k],
                              bxTaper[i],byTaper[j],bzTaper[k],
                              dvx,dvy,dvz);
          dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
        }
      }
      k=1;
      for(j=2;j<yLim;++j) {
        int jkind = j*NX+NXY;
        //below function call assumes ssfunc is has been set to 2nd order for k=1
        updateAcousticPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
      }
    }
    for(k=2;k<zStart;++k) {
      int kind = k*NXY;
      if(fminYb) {
        j=0;
        int jkind = kind;
        //the inline function below is the only difference compared to the no pml case
        for(i=2;i<xLim;i++) {
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axTaper[i],ayTaper[j],azTaper[k],
                              bxTaper[i],byTaper[j],bzTaper[k],
                              dvx,dvy,dvz);
          dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
        }
        j=1;
        jkind = NX+kind;
        //assumes ssfunc is O2 for j=1
        updateAcousticPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
      }
      for(j=2;j<yLim;++j) {
        int jkind = j*NX+kind;
        if(fminXb) {
          i=0;
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          int pIndex = index;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmin,vyyZmin,vzzZmin,pIndex,axTaper[i],ayTaper[j],azTaper[k],
                              bxTaper[i],byTaper[j],bzTaper[k],
                              dvx,dvy,dvz);
          dvx = dvx*kxTaper[i]+vxxZmin[pIndex];
          dvy = dvy*kyTaper[j]+vyyZmin[pIndex];
          dvz = dvz*kzTaper[k]+vzzZmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
          //assumes ssfunc is O2 for x=1
          updateAcousticPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,1,2,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        }
        updateAcousticPressurePML(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,2,xLim,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,jkind,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,xLim,NX,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          int jkind = j*NX+kind;
          //assumes ssfunc is O2 for y=yLim and NY-1
          updateAcousticPressurePMLO2(vx,vy,vz,vxxZmin,vyyZmin,vzzZmin,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,jkind,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
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
        //the inline function below is the only difference compared to the no pml case
        for(i=2;i<xLim;i++) {
          int index = i+jkind;
          int pIndex = i+j*NX+(k-kk)*nxy;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxYmin,vyyYmin,vzzYmin,pIndex,axTaper[i],ayTaper[j],azTaper[k],
                              bxTaper[i],byTaper[j],bzTaper[k],
                              dvx,dvy,dvz);
          dvx = dvx*kxTaper[i]+vxxYmin[pIndex];
          dvy = dvy*kyTaper[j]+vyyYmin[pIndex];
          dvz = dvz*kzTaper[k]+vzzYmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
        }
        j=1;
        jkind = NX+kind;
        //assumes ssfunc is O2 for j=1
        updateAcousticPressurePMLO2(vx,vy,vz,vxxYmin,vyyYmin,vzzYmin,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,j*NX+(k-kk)*nxy,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
      }
      for(j=2;j<yStart2;++j) {
        int jkind = j*NX+kind;
        if(fminXb) {
          i=0;
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          int pIndex = i+j*NX+(k-kk)*nxy;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxYmin,vyyYmin,vzzYmin,pIndex,axTaper[i],ayTaper[j],azTaper[k],
                              bxTaper[i],byTaper[j],bzTaper[k],
                              dvx,dvy,dvz);
          dvx = dvx*kxTaper[i]+vxxYmin[pIndex];
          dvy = dvy*kyTaper[j]+vyyYmin[pIndex];
          dvz = dvz*kzTaper[k]+vzzYmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
          //assumes ssfunc is O2 for x=1
          updateAcousticPressurePMLO2(vx,vy,vz,vxxYmin,vyyYmin,vzzYmin,1,2,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,j*NX+(k-kk)*nxy,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        }
        updateAcousticPressurePML(vx,vy,vz,vxxYmin,vyyYmin,vzzYmin,2,xLim,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,j*NX+(k-kk)*nxy,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticPressurePMLO2(vx,vy,vz,vxxYmin,vyyYmin,vzzYmin,xLim,NX,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,j*NX+(k-kk)*nxy,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
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
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          int pIndex = i+(j-jj)*nx+(k-kk)*nxy;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxXmin,vyyXmin,vzzXmin,pIndex,axTaper[i],ayTaper[j],azTaper[k],
                              bxTaper[i],byTaper[j],bzTaper[k],
                              dvx,dvy,dvz);
          dvx = dvx*kxTaper[i]+vxxXmin[pIndex];
          dvy = dvy*kyTaper[j]+vyyXmin[pIndex];
          dvz = dvz*kzTaper[k]+vzzXmin[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
          //assumes ssfunc is O2 for x=1
          updateAcousticPressurePMLO2(vx,vy,vz,vxxXmin,vyyXmin,vzzXmin,1,2,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,(j-jj)*nx+(k-kk)*nxy,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        }
        updateAcousticPressurePML(vx,vy,vz,vxxXmin,vyyXmin,vzzXmin,2,_nXmin,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,(j-jj)*nx+(k-kk)*nxy,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
      }
      updateAcousticPressureNoPML(vx,vy,vz,xStart2,xEnd,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pressure);
      if(maxXb) {
        iii=_xMaxStart;
        nxy = _nXmax*nyInt;
        nx = _nXmax;
        //assumes ssfunc is O2 for x=xLim and NX-1
        updateAcousticPressurePML(vx,vy,vz,vxxXmax,vyyXmax,vzzXmax,xEnd,xLim,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,(j-jj)*nx+(k-kk)*nxy-iii,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        if(fmaxXb) {
          updateAcousticPressurePMLO2(vx,vy,vz,vxxXmax,vyyXmax,vzzXmax,xLim,NX,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,(j-jj)*nx+(k-kk)*nxy-iii,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        }
      }
    }
    iii=0;
    jj=_yMaxStart;
    nxy = NX*_nYmax;
    nx = NX;
    if(maxYb) {
      for(j=yEnd;j<yLim;++j) {
        int jkind = j*NX+kind;
        if(fminXb) {
          i=0;
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          int pIndex = i+(j-jj)*NX+(k-kk)*nxy;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxYmax,vyyYmax,vzzYmax,pIndex,axTaper[i],ayTaper[j],azTaper[k],
                              bxTaper[i],byTaper[j],bzTaper[k],
                              dvx,dvy,dvz);
          dvx = dvx*kxTaper[i]+vxxYmax[pIndex];
          dvy = dvy*kyTaper[j]+vyyYmax[pIndex];
          dvz = dvz*kzTaper[k]+vzzYmax[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
          //assumes ssfunc is O2 for x=1
          updateAcousticPressurePMLO2(vx,vy,vz,vxxYmax,vyyYmax,vzzYmax,1,2,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,(j-jj)*NX+(k-kk)*nxy,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        }
        updateAcousticPressurePML(vx,vy,vz,vxxYmax,vyyYmax,vzzYmax,2,xLim,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,(j-jj)*NX+(k-kk)*nxy,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticPressurePMLO2(vx,vy,vz,vxxYmax,vyyYmax,vzzYmax,xLim,NX,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,(j-jj)*NX+(k-kk)*nxy,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          int jkind = j*NX+kind;
          //assumes ssfunc is O2 for y=yLim and NY-1
          updateAcousticPressurePMLO2(vx,vy,vz,vxxYmax,vyyYmax,vzzYmax,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,(j-jj)*NX+(k-kk)*nxy,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        }
      }
    }
  }
  if(maxZb) {
    jj=0;
    kk = _zMaxStart;
    nxy = NXY;
    for(k=zEnd;k<zLim;++k) {
      int kind = k*NXY;
      int pkind = (k-kk)*NXY;
      if(fminYb) {
        j=0;
        int jkind = kind;
        int pjkind = pkind;
        //the inline function below is the only difference compared to the no pml case
        for(i=2;i<xLim;i++) {
          int index = i+jkind;
          int pIndex = i+pjkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmax,vyyZmax,vzzZmax,pIndex,axTaper[i],ayTaper[j],azTaper[k],
                              bxTaper[i],byTaper[j],bzTaper[k],
                              dvx,dvy,dvz);
          dvx = dvx*kxTaper[i]+vxxZmax[pIndex];
          dvy = dvy*kyTaper[j]+vyyZmax[pIndex];
          dvz = dvz*kzTaper[k]+vzzZmax[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
        }
        j=1;
        jkind = NX+kind;
        pjkind = NX+pkind;
        //assumes ssfunc is O2 for j=1
        updateAcousticPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
      }
      for(j=2;j<yLim;++j) {
        int jkind = j*NX+kind;
        int pjkind = j*NX+pkind;
        if(fminXb) {
          i=0;
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          int pIndex = i+pjkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          //PML updating in group of lines below
          updatePMLMemVarComp(vxxZmax,vyyZmax,vzzZmax,pIndex,axTaper[i],ayTaper[j],azTaper[k],
                              bxTaper[i],byTaper[j],bzTaper[k],
                              dvx,dvy,dvz);
          dvx = dvx*kxTaper[i]+vxxZmax[pIndex];
          dvy = dvy*kyTaper[j]+vyyZmax[pIndex];
          dvz = dvz*kzTaper[k]+vzzZmax[pIndex];
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
          //assumes ssfunc is O2 for x=1
          updateAcousticPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,1,2,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        }
        updateAcousticPressurePML(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,2,xLim,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pjkind,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,xLim,NX,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          int jkind = j*NX+kind;
          int pjkind = j*NX+pkind;
          //assumes ssfunc is O2 for y=yLim and NY-1
          updateAcousticPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        }
      }
    }
    if(fmaxZb) {
      for(k=zLim;k<NZ;++k) {
        int kind = k*NXY;
        int pkind = (k-kk)*NXY;
        for(j=2;j<yLim;++j) {
          int jkind = j*NX+kind;
          int pjkind = j*NX+pkind;
          //assumes ssfunc is O2 for z=zLim and NZ-1
          updateAcousticPressurePMLO2(vx,vy,vz,vxxZmax,vyyZmax,vzzZmax,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pjkind,pressure,ayTaper[j],azTaper[k],byTaper[j],bzTaper[k],kyTaper[j],kzTaper[k]);
        }
      }
    }
  }
}

//Extrapolate the pressure one node above the free surface using 4th
//order extrapolation using pressure=0 on the free surface
void updateAcousticPressurePressFreeCPML(modelDefStruct* modelDef,
                                float*__restrict__ pressure,
                                float*__restrict__ vx,float*__restrict__ vy,float*__restrict__ vz,
                                float cx[2],float cy[2],float cz[2],
                                float*__restrict__ bulk,
                                unsigned char*__restrict__ ssfunc,
                                bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb) {
  using namespace acoustic_CPML;
  DEF_MODEL_SIZE(modelDef);
  
  //  minXb = maxXb = minYb = maxYb = minZb = maxZb = false;
  //  maxZb=false;
  
  register int i,j,k;
  
  k=1;
  for(j=0;j<NY;++j) {
    for(i=0;i<NX;++i) {
      int index = i+j*NX+k*NXY;
      pressure[index] = 4.*pressure[index+3*NXY]-(6.*pressure[index+2*NXY]+pressure[index+4*NXY]);
    }
  }
}
