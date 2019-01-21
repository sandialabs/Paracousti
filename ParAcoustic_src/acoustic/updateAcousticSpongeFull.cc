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
 *  updateAnelasticSpongeFull.cpp
 *  
 *
 *  This file contains the setup and updating functions for all the
 *  dependent variables for the case where sponges serve as the absorbing
 *  BCs.  It also includes the specialized functions it uses if a pressure
 *  free surface is used.  For the updating functions, updating proceeds
 *  in exactly the order that the dependent variables are stored in
 *  memory.  It gives a big speed boost compared to updating interior
 *  points then positive x flanks, then negative x flanks, etc.
 *
 *  Defines the following functions:
 *  setupAcousticSpongeBoundsFull
 *  updateVxSpongeBounds
 *  updateVySpongeBounds
 *  updateVzSpongeBounds
 *  updateVzPressFreeSpongeBounds
 *  updateAcousticPressureSponge
 *  updateAcousticPressurePressFreeSponge
 *
 */

#include "updateAcousticSpongeFull.h"
#include "acousticSpongeNameSpace.h"
#include "sgfd.hh"
#include "constants.h"

#define ELASTIC_FS_2_0 3
#define SURFACE_PRESS_FREE 2

const float pio2 = 0.5*acos(-1.);

//  This function sets up the extra variables and arrays needed for the
//  sponge BCs and fills in the values of the sponge taper zones relative
//  to the flanks of the model.
void setupAcousticSpongeBoundsFull(int nXmin, float xMinVal, int nXmax, float xMaxVal, int nYmin, float yMinVal,
                               int nYmax, float yMaxVal, int nZmin, float zMinVal, int nZmax, float zMaxVal, int kstart, int kstart_z, sgfdModel* model, int sbcmode) {
  using namespace acoustic_Sponge;
  DEF_MODEL_SIZE(model->modelDef());
  DEF_PARALLEL(model->parallelDef());
  
  //define the thickness of the sponge zone in this domain
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
  
  //declare taper variables
  xTaper = new float[NX];
  yTaper = new float[NY];
  zTaper = new float[NZ];
  xTaperPH = new float[NX];
  yTaperPH= new float[NY];
  zTaperPH = new float[NZ];
  
  //Fill in the flank profile variables
  //fromStart and fromEnd are maximum at the min side or max side, respectively.  Appropriate values are filled in when outside the sponge layers; however, these values are not used.
  //cosine tapers are used with a value of 1 in equivalent to no sponge BC and tapers to the minVal at the flank.
  //X-profile
  for(int i=0;i<NX;i++) {
    float fromStart = nXmin-i-procLim[0];
    float fromEnd = nXmax-globalNX+procLim[1]-(NX-1-i);
    if(fromStart>0.f) {
      xTaper[i] = xMinVal+(1.-xMinVal)*cos(fromStart/nXmin*pio2);
      xTaperPH[i] = xMinVal+(1.-xMinVal)*cos((fromStart-0.5)/nXmin*pio2);
    } else if(fromEnd>-0.5f) {
      xTaper[i] = xMaxVal+(1.-xMaxVal)*cos(fromEnd/nXmax*pio2);
      xTaperPH[i] = xMaxVal+(1.-xMaxVal)*cos((fromEnd+0.5)/nXmax*pio2);
    } else {
      xTaper[i] = xTaperPH[i] = 1.;
    }
  }
  //Y-profile
  for(int i=0;i<NY;i++) {
    float fromStart = nYmin-i-procLim[2];
    float fromEnd = nYmax-globalNY+procLim[3]-(NY-1-i);
    if(fromStart>0.f) {
      yTaper[i] = yMinVal+(1.-yMinVal)*cos(fromStart/nYmin*pio2);
      yTaperPH[i] = yMinVal+(1.-yMinVal)*cos((fromStart-0.5)/nYmin*pio2);
    } else if(fromEnd>-0.5f){
      yTaper[i] = yMaxVal+(1.-yMaxVal)*cos(fromEnd/nYmax*pio2);
      yTaperPH[i] = yMaxVal+(1.-yMaxVal)*cos((fromEnd+0.5)/nYmax*pio2);
    } else {
      yTaper[i] = yTaperPH[i] = 1.;
    }
  }
  //Z-profile
  for(int i=0;i<NZ;i++) {
    float fromStart = nZmin-i-procLim[4];
    float fromEnd = nZmax-globalNZ+procLim[5]-(NZ-1-i);
    if(fromStart>0.f) {
      zTaper[i] = zMinVal+(1.-zMinVal)*cos(fromStart/nZmin*pio2);
      zTaperPH[i] = zMinVal+(1.-zMinVal)*cos((fromStart-0.5)/nZmin*pio2);
    } else if(fromEnd>-0.5f) {
      zTaper[i] = zMaxVal+(1.-zMaxVal)*cos(fromEnd/nZmax*pio2);
      zTaperPH[i] = zMaxVal+(1.-zMaxVal)*cos((fromEnd+0.5)/nZmax*pio2);
    } else {
      zTaper[i] = zTaperPH[i] = 1.;
    }
  }
}

//Here are the updating equations for Vx.  These are done in exactly
//the order that vx is stored in memory starting with the minZ, minY
//minX corner proceeding to maxZ, maxY, and maxX with x fastest, then
//y.  The variables minXb, fminXb, etc are used to define where this
//domain is situated relative to the sponge zones and the absolute edges
//of the model.  minXb, maxXb, etc. will be true if the sponge zone extends
//into this domain on the specified side.  fminXb, fmaxXb, etc. will only
//be true if this domain is actually on an absolute edge of the model.
//These two sets will the identical except in the cases where the sponge
//zone is wider than the width of the domain on the absolute edge.
//4th order accurate differences are used when they can be, but near the
//flanks, second order differences are used and finally, in some cases,
//second order with the assumption of zero valued variables off grid
//are used right on the flanks.
void updateVxSpongeBounds(modelDefStruct* modelDef,
                       float* __restrict__ vx,
                       float* __restrict__ pressure,
                       float cx[2],float cy[2],float cz[2],
                       float* __restrict__ rho,unsigned char* __restrict__ vxfunc,
                       bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                       bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb) {
  using namespace acoustic_Sponge;
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
  int mxStart1 = MAX(_nXmin,1);
  int mxEnd1 = MIN(_xMaxStart-1,NX-2);
  
  //printf("vx %d %d %d, %d %d %d %d %d %d\n",NX,NY,NZ,minXb,maxXb,minYb,maxYb,minZb,maxZb);
  if(minZb) {
    if(fminZb) {
      k=0;
      for(j=2;j<yLim;j++){
        for(i=1;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
      }
      k=1;
      for(j=2;j<yLim;j++){
        for(i=1;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
      }
    }
    for(k=2;k<zStart;k++) {
      if(fminYb) {
        j=0;
        for(i=1;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
        j=1;
        for(i=1;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
      }
      for(j=2;j<yLim;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
        for(i=1;i<xLim;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,vx,i,j,k);
        }
        if (fmaxXb) {
          i=xLim;
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          for(i=1;i<xLim;++i) {
            int index=i+j*NX+k*NXY;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
          }
        }
      }
    }
  }
  for(k=zStart;k<zEnd;++k) {
    if(minYb) {
      if(fminYb) {
        j=0;
        for(i=1;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
        j=1;
        for(i=1;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
      }
      for(j=2;j<yStart2;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
        for(i=1;i<xLim;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,vx,i,j,k);
        }
        if (fmaxXb) {
          i=xLim;
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
      }
    }
    for(j=yStart2;j<yEnd;++j) {
      if(minXb) {
        if(fminXb) {
          i=0;
          int index = i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
        for(i=1;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,vx,i,j,k);
        }
      }
      for(i=mxStart1;i<mxEnd1;++i) {
        int index = i+j*NX+k*NXY;
        updateAcousticVxNoSponge(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,vx);
      }
      if (maxXb) {
        for(i=xEnd1;i<xLim;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,vx,i,j,k);
        }
        if(fmaxXb) {
          i=xLim;
          int index = i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
      }
    }
    if(maxYb) {
      for(j=yEnd;j<yLim;++j) {
        if (fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
        for(i=1;i<xLim;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,vx,i,j,k);
        }
        if (fmaxXb) {
          i=xLim;
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          for(i=1;i<xLim;++i) {
            int index=i+j*NX+k*NXY;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
          }
        }
      }
    }
  }
  if(maxZb) {
    for(k=zEnd;k<zLim;++k) {
      if(fminYb) {
        j=0;
        for(i=1;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
        j=1;
        for(i=1;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
      }
      for(j=2;j<yLim;++j) {
        if (fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
        for(i=1;i<xLim;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVx(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vxfunc,index,vx,i,j,k);
        }
        if (fmaxXb) {
          i=xLim;
          int index=i+j*NX+k*NXY;
          updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          for(i=1;i<xLim;++i) {
            int index=i+j*NX+k*NXY;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
          }
        }
      }
    }
    if(fmaxZb) {
      for(k=zLim;k<NZ;k++) {
        for(j=2;j<yLim;j++){
          for(i=1;i<xLim;i++){
            int index=i+j*NX+k*NXY;
            updateAcousticVxO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vx,i,j,k);
          }
        }
      }
    }
  }
}

//vy updating. see vx updating for description
void updateVySpongeBounds(modelDefStruct* modelDef,
                       float* __restrict__ vy,
                       float* __restrict__ pressure,
                       float cx[2],float cy[2],float cz[2],
                       float* __restrict__ rho,
                       unsigned char* __restrict__ vyfunc,
                       bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                       bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb) {
  using namespace acoustic_Sponge;
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
  int xStart2 = MAX(_nXmin,2);
  int xEnd = MIN(_xMaxStart,NX-2);
  
  //printf("vy %d %d %d, %d %d %d %d %d %d\n",NX,NY,NZ,minXb,maxXb,minYb,maxYb,minZb,maxZb);
  if(minZb) {
    if(fminZb) {
      k=0;
      for(j=1;j<yLim;j++){
        for(i=2;i<xLim;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
        }
      }
      k=1;
      for(j=1;j<yLim;j++){
        for(i=2;i<xLim;i++){
          //define indices to save later multiplications
          int index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
        }
      }
    }
    for(k=2;k<zStart;++k) {
      if(fminYb) {
        j=0;
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
        }
      }
      for(j=1;j<yLim;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
        }
        for(i=2;i<xLim;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,index,vy,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int index=i+j*NX+k*NXY;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
          }
        }
      }  
      if(fmaxYb) {
        j=yLim;
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
        }
      }
    }
  }
  for(k=zStart;k<zEnd;++k) {
    if(minYb) {
      if(fminYb) {
        j=0;
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
        }
      }
      for(j=1;j<yStart1;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
        }
        for(i=2;i<xLim;++i) {
          int index = i+j*NX+k*NXY;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,index,vy,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int index=i+j*NX+k*NXY;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
          }
        }
      }
    }
    for(j=yStart1;j<yEnd1;++j) {
      if(minXb) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
        }
        for(i=2;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,index,vy,i,j,k);
        }
      }
      for(i=xStart2;i<xEnd;++i) {
        int index = i+j*NX+k*NXY;
        updateAcousticVyNoSponge(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,index,vy);
      }
      if(maxXb) {
        for(i=xEnd;i<xLim;++i) {
          int index = i+j*NX+k*NXY;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,index,vy,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int index=i+j*NX+k*NXY;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
          }
        }
      }
    }
    if(maxYb) {
      for(j=yEnd1;j<yLim;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
        }
        for(i=2;i<xLim;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,index,vy,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int index=i+j*NX+k*NXY;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
          }
        }
      }
      if(fmaxYb) {
        j=yLim;
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
        }
      }
    }
  }
  if(maxZb) {
    for(k=zEnd;k<zLim;++k) {
      if(fminYb) {
        j=0;
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
        }
      }
      for(j=1;j<yLim;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
        }
        for(i=2;i<xLim;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVy(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vyfunc,index,vy,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;++i) {
            int index=i+j*NX+k*NXY;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
          }
        }
      }      
      if(fmaxYb) {
        j=yLim;
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
        }
      }
    }
    if(fmaxZb) {
      for(k=zLim;k<NZ;k++) {
        for(j=1;j<yLim;j++){
          for(i=2;i<xLim;i++){
            int index=i+j*NX+k*NXY;
            updateAcousticVyO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vy,i,j,k);
          }
        }
      }
    }
  }
}

void updateVzSpongeBounds(modelDefStruct* modelDef,
                       float* __restrict__ vz,
                       float* __restrict__ pressure,
                       float cx[2],float cy[2],float cz[2],
                       float* __restrict__ rho,
                       unsigned char* __restrict__ vzfunc,
                       bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                       bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb)  {
  using namespace acoustic_Sponge;
  DEF_MODEL_SIZE(modelDef);
  
  register int i,j,k;
  int zLim=NZ-2,yLim=NY-2,xLim=NX-2;
  float cx0=cx[0],cx1=cx[1];
  float cy0=cy[0],cy1=cy[1];
  float cz0=cz[0],cz1=cz[1];
  float dcx0=modelDef->dt/modelDef->dx*modelDef->scalarSpeed;
  float dcy0=modelDef->dt/modelDef->dy*modelDef->scalarSpeed;
  float dcz0=modelDef->dt/modelDef->dz*modelDef->scalarSpeed;
  int xStart2 = MAX(_nXmin,2);
  int xEnd = MIN(_xMaxStart,NX-2);
  
  //printf("vz %d %d %d, %d %d %d %d %d %d\n",NX,NY,NZ,minXb,maxXb,minYb,maxYb,minZb,maxZb);
  if(minZb) {
    if(fminZb) {
      k=0;
      for(j=2;j<yLim;j++){
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
        }
      }
    }
    for(k=1;k<zStart_z;++k) {
      if(fminYb) {
        j=0;
        for(i=2;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
        }
        j=1;
        for(i=2;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
        }
      }
      for(j=2;j<yLim;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
        }
        for(i=2;i<xLim;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,index,vz,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int index=i+j*NX+k*NXY;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
          }
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          for(i=2;i<xLim;++i) {
            int index=i+j*NX+k*NXY;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
          }
        }
      }
    }
  }
  for(k=zStart_z;k<zEnd1;++k) {
    if(minYb) {
      if(fminYb) {
        j=0;
        for(i=2;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
        }
        j=1;
        for(i=2;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
        }
      }
      for(j=2;j<yStart2;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
        }
        for(i=2;i<xLim;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,index,vz,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int index=i+j*NX+k*NXY;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
          }
        }
      }
    }
    for(j=yStart2;j<yEnd;++j) {
      if(minXb) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
        }
        for(i=2;i<_nXmin;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,index,vz,i,j,k);
        }
      }
      for(i=xStart2;i<xEnd;++i) {
        int index = i+j*NX+k*NXY;
        updateAcousticVzNoSponge(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,index,vz);
      }        
      if(maxXb) {
        for(i=xEnd;i<xLim;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,index,vz,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int index=i+j*NX+k*NXY;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
          }
        }
      }
    }
    if(maxYb) {
      for(j=yEnd;j<yLim;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
        }
        for(i=2;i<xLim;++i) {
          int index = i+j*NX+k*NXY;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,index,vz,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int index=i+j*NX+k*NXY;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
          }
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          for(i=2;i<xLim;++i) {
            int index=i+j*NX+k*NXY;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
          }
        }
      }
    }
  }
  if(maxZb) {
    for(k=zEnd1;k<zLim;++k) {
      if(fminYb) {
        j=0;
        for(i=2;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
        }
        j=1;
        for(i=2;i<xLim;++i) {
          int index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
        }
      }      
      for(j=2;j<yLim;++j) {
        if(fminXb) {
          i=0;
          int index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
          
          i=1;
          index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
        }
        for(i=2;i<xLim;i++) {
          int index = i+j*NX+k*NXY;
          updateAcousticVz(pressure,rho,NX,NXY,cx0,cx1,dcx0,cy0,cy1,dcy0,cz0,cz1,dcz0,vzfunc,index,vz,i,j,k);
        }
        if(fmaxXb) {
          for(i=xLim;i<NX;i++) {
            int index=i+j*NX+k*NXY;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
          }
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          for(i=2;i<xLim;++i) {
            int index=i+j*NX+k*NXY;
            updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
          }
        }
      }      
    }
    if(fmaxZb) {
      k=zLim;
      for(j=2;j<yLim;j++){
        for(i=2;i<xLim;i++){
          int index=i+j*NX+k*NXY;
          updateAcousticVzO2(pressure,rho,NX,NXY,dcx0,dcy0,dcz0,index,vz,i,j,k);
        }
      }
    }
  }
}

//This is a special updating formula for vz directly above the free
//surface.  It extrapolates from vz below the free surface using
//4th order extrapolation and using dvz/dz=0 on the free surface
void updateVzPressFreeSpongeBounds(modelDefStruct* modelDef,
                        float* __restrict__ vz,
                        float* __restrict__ pressure,
                        float cx[2],float cy[2],float cz[2],
                        float* __restrict__ rho,
                        unsigned char* __restrict__ vzfunc,
                        bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb)  {
  using namespace acoustic_Sponge;
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
void updateAcousticPressureSponge(modelDefStruct* modelDef,
                                   float*__restrict__ pressure,
                                   float*__restrict__ vx,float*__restrict__ vy,float*__restrict__ vz,
                                   float cx[2],float cy[2],float cz[2],
                                   float*__restrict__ bulk,
                                   unsigned char*__restrict__ ssfunc,
                                   bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                                   bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb) {
  using namespace acoustic_Sponge;
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
  
  if(minZb) {
    if(fminZb) {
      k=0;
      for(j=2;j<yLim;++j) {
        int jkind = j*NX;
        
        //the inline function below is the only difference compared to the no pml case
        for(i=2;i<xLim;i++) {
          int index = i+jkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*vz[index];
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
          pressure[index]*=xTaper[i]*yTaper[j]*zTaper[k];
        }
      }
      k=1;
      for(j=2;j<yLim;++j) {
        int jkind = j*NX+NXY;
        //below function call assumes ssfunc is has been set to 2nd order for k=1
        updateAcousticPressureSpongeO2(vx,vy,vz,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
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
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);

          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
          pressure[index]*=xTaper[i]*yTaper[j]*zTaper[k];
        }
        j=1;
        jkind = NX+kind;
        //assumes ssfunc is O2 for j=1
        updateAcousticPressureSpongeO2(vx,vy,vz,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
      }
      for(j=2;j<yLim;++j) {
        int jkind = j*NX+kind;
        if(fminXb) {
          i=0;
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
          pressure[index]*=xTaper[i]*yTaper[j]*zTaper[k];
          //assumes ssfunc is O2 for x=1
          updateAcousticPressureSpongeO2(vx,vy,vz,1,2,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
        }
        updateAcousticPressureSponge(vx,vy,vz,2,xLim,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pressure,yTaper[j],zTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticPressureSpongeO2(vx,vy,vz,xLim,NX,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          int jkind = j*NX+kind;
          //assumes ssfunc is O2 for y=yLim and NY-1
          updateAcousticPressureSpongeO2(vx,vy,vz,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
        }
      }
    }
  }    
  for(k=zStart;k<zEnd;++k) {
    int kind = k*NXY;
    if(minYb) {
      if(fminYb) {
        j=0;
        int jkind = kind;
        //the inline function below is the only difference compared to the no pml case
        for(i=2;i<xLim;i++) {
          int index = i+jkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
          pressure[index]*=xTaper[i]*yTaper[j]*zTaper[k];
        }
        j=1;
        jkind = NX+kind;
        //assumes ssfunc is O2 for j=1
        updateAcousticPressureSpongeO2(vx,vy,vz,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
      }
      for(j=2;j<yStart2;++j) {
        int jkind = j*NX+kind;
        if(fminXb) {
          i=0;
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
          pressure[index]*=xTaper[i]*yTaper[j]*zTaper[k];
          //assumes ssfunc is O2 for x=1
          updateAcousticPressureSpongeO2(vx,vy,vz,1,2,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
        }
        updateAcousticPressureSponge(vx,vy,vz,2,xLim,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pressure,yTaper[j],zTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticPressureSpongeO2(vx,vy,vz,xLim,NX,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
        }
      }
    }
    for(j=yStart2;j<yEnd;++j) {
      int jkind = j*NX+kind;
      if(minXb) {
        if(fminXb) {
          i=0;
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
          pressure[index]*=xTaper[i]*yTaper[j]*zTaper[k];
          //assumes ssfunc is O2 for x=1
          updateAcousticPressureSpongeO2(vx,vy,vz,1,2,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
        }
        updateAcousticPressureSponge(vx,vy,vz,2,_nXmin,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pressure,yTaper[j],zTaper[k]);
      }
      updateAcousticPressureNoSponge(vx,vy,vz,xStart2,xEnd,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pressure);
      if(maxXb) {
        //assumes ssfunc is O2 for x=xLim and NX-1
        updateAcousticPressureSponge(vx,vy,vz,xEnd,xLim,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pressure,yTaper[j],zTaper[k]);
        if(fmaxXb) {
          updateAcousticPressureSpongeO2(vx,vy,vz,xLim,NX,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
        }
      }
    }
    if(maxYb) {
      for(j=yEnd;j<yLim;++j) {
        int jkind = j*NX+kind;
        if(fminXb) {
          i=0;
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);

          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
          pressure[index]*=xTaper[i]*yTaper[j]*zTaper[k];
          //assumes ssfunc is O2 for x=1
          updateAcousticPressureSpongeO2(vx,vy,vz,1,2,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
        }
        updateAcousticPressureSponge(vx,vy,vz,2,xLim,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pressure,yTaper[j],zTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticPressureSpongeO2(vx,vy,vz,xLim,NX,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          int jkind = j*NX+kind;
          //assumes ssfunc is O2 for y=yLim and NY-1
          updateAcousticPressureSpongeO2(vx,vy,vz,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
        }
      }
    }
  }
  if(maxZb) {
    for(k=zEnd;k<zLim;++k) {
      int kind = k*NXY;
      if(fminYb) {
        j=0;
        int jkind = kind;
        //the inline function below is the only difference compared to the no pml case
        for(i=2;i<xLim;i++) {
          int index = i+jkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]-vx[index-1]);
          dvy = dcy0*(vy[index]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
          pressure[index]*=xTaper[i]*yTaper[j]*zTaper[k];
        }
        j=1;
        jkind = NX+kind;
        //assumes ssfunc is O2 for j=1
        updateAcousticPressureSpongeO2(vx,vy,vz,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
      }
      for(j=2;j<yLim;++j) {
        int jkind = j*NX+kind;
        if(fminXb) {
          i=0;
          //the inline function below is the only difference compared to the no pml case
          int index = i+jkind;
          float dvx, dvy, dvz;
          float bb = bulk[index];
          dvx = dcx0*(vx[index]);
          dvy = dcy0*(vy[index]-vy[index-NX]);
          dvz = dcz0*(vz[index]-vz[index-NXY]);
          
          float dvTot = dvx+dvy+dvz;
          pressure[index]-=dvTot*bb;
          pressure[index]*=xTaper[i]*yTaper[j]*zTaper[k];
          //assumes ssfunc is O2 for x=1
          updateAcousticPressureSpongeO2(vx,vy,vz,1,2,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
        }
        updateAcousticPressureSponge(vx,vy,vz,2,xLim,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,jkind,pressure,yTaper[j],zTaper[k]);
        if(fmaxXb) {
          //assumes ssfunc is O2 for x=xLim and NX-1
          updateAcousticPressureSpongeO2(vx,vy,vz,xLim,NX,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
        }
      }
      if(fmaxYb) {
        for(j=yLim;j<NY;++j) {
          int jkind = j*NX+kind;
          //assumes ssfunc is O2 for y=yLim and NY-1
          updateAcousticPressureSpongeO2(vx,vy,vz,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
        }
      }
    }
    if(fmaxZb) {
      for(k=zLim;k<NZ;++k) {
        int kind = k*NXY;
        for(j=2;j<yLim;++j) {
          int jkind = j*NX+kind;
          //assumes ssfunc is O2 for z=zLim and NZ-1
          updateAcousticPressureSpongeO2(vx,vy,vz,2,xLim,bulk,NX,NXY,dcx0,dcy0,dcz0,jkind,pressure,yTaper[j],zTaper[k]);
        }
      }
    }
  }
}

//Extrapolate the pressure one node above the free surface using 4th
//order extrapolation using pressure=0 on the free surface
void updateAcousticPressurePressFreeSponge(modelDefStruct* modelDef,
                                float*__restrict__ pressure,
                                float*__restrict__ vx,float*__restrict__ vy,float*__restrict__ vz,
                                float cx[2],float cy[2],float cz[2],
                                float*__restrict__ bulk,
                                unsigned char*__restrict__ ssfunc,
                                bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb) {
  using namespace acoustic_Sponge;
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
