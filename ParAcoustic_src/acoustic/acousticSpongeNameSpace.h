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
 *  acousticSpongeNameSpace.h
 *  
 *
 *  This file contains the updating kernels that actually form the core
 *  of the code for the sponge boundary.  Both O2 and O4 updating for the
 *  sponge zones and O4 for the interior zone.  These are called from
 *  updateAcousticSpongeFull.cc
 *
 *  The first section includes the varibles that are actually defined
 *  in acousticSpongeNameSpace.cc.
 *
 *  Defines several inline functions used in updateAcousticSpongeFull.cc
 *  or internally:
 *    sqr
 *    updateSpongeOne
 *    updateAcousticVx
 *    updateAcousticVxO2
 *    updateAcousticVy
 *    updateAcousticVyO2
 *    updateAcousticVz
 *    updateAcousticVzO2
 *    updateAcousticVxNoSponge
 *    updateAcousticVyNoSponge
 *    updateAcousticVzNoSponge
 *    updateAcousticPressureNoSpongeCore
 *    updateAcousticPressureNoSpongeCoreO2
 *    updateAcousticPressureNoSponge
 *    updateAcousticPressureNoSpongeO2
 *    updateAcousticPressureSpongeCore
 *    updateAcousticPressureSpongeCoreO2
 *    updateAcousticPressureSponge
 *    updateAcousticPressureSpongeO2
 */
#ifndef _ACOUSTIC_SPONGE_NS
#define _ACOUSTIC_SPONGE_NS

#include<math.h>

#include<stdio.h>

namespace acoustic_Sponge {
  extern float *__restrict__ xTaper, *__restrict__ yTaper, *__restrict__ zTaper;
  extern float *__restrict__ xTaperPH, *__restrict__ yTaperPH, *__restrict__ zTaperPH;
    
  extern int _nXmin, _nXmax, _nYmin, _nYmax, _nZmin, _nZmax;
  extern int _xMaxStart, _yMaxStart, _zMaxStart;
  extern int xStart1, xStart2, yStart1, yStart2, zStart, zStart_z;
  extern int xEnd, yEnd, zEnd, xEnd1, yEnd1, zEnd1;
  extern int _surfaceBCMode;
  extern int _kstart;
  
  //calculates the square of a number, saving some typing
  inline float sqr(float x) {return x*x;}
  
  //Acutally update a velocity
  //component (vx) using sponge parameters xtaper, ytaper, and ztaper
  //at index in the vx array (which corresponds to pindex in the memory
  //variable array and the density (rho) and spatial derivatives of
  //the appropriate stresses (derX)
  inline void updateSpongeOne(float*__restrict__ vx, int index,
                           float xtaper,float ytaper,
                           float ztaper,float derX, float rho) {
    vx[index] -= derX/rho;  //note: kxtaper is one over the original
    vx[index] *= xtaper*ytaper*ztaper;
  }
  
  //Update the vx for O4 sponge.  The c variables refer to O4 FD coefficients
  //while the d variables refer to O2 coefficients.  High contrast
  //interfaces trigger the change from O4 to O2.  This is determined
  //in acoustic_control.cc and a flag array vxfunc holds this info.  px
  //is the sponge memory variable and ii, jj and kk refer to the 3-D indices
  //within vx
  inline void updateAcousticVx(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float cx0, float cx1, float dcx0,
                               float cy0, float cy1, float dcy0, float cz0, float cz1, float dcz0, unsigned char*__restrict__ vxfunc, int index,
                               float*__restrict__ vx, int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+1]);
    float dfx;
    unsigned char dd = vxfunc[index];
    unsigned char dd1 = 1-dd;
    float cx0d = cx0*dd1+dcx0*dd;
    float cx1d = cx1*dd1;
    dfx = cx0d*(pressure[index+1]-pressure[index   ]) +
    cx1d*(pressure[index+2]-pressure[index-1]);
    updateSpongeOne(vx,index,xTaperPH[ii],yTaper[jj],
                 zTaper[kk],dfx,rhoave);
  }    
  
  //Update the vx for O2 sponge.  This is only used near the very flanks of
  //the model
  inline void updateAcousticVxO2(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float dcx0,
                                 float dcy0,float dcz0, int index,
                                 float*__restrict__ vx, int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+1]);
    float dfx;
    dfx = dcx0*(pressure[index+1]-pressure[index   ]);
    updateSpongeOne(vx,index,xTaperPH[ii],yTaper[jj],
                 zTaper[kk],dfx,rhoave);
  }
  
  //Update vy for O4 sponge.
  inline void updateAcousticVy(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float cx0, float cx1, float dcx0,
                               float cy0, float cy1, float dcy0, float cz0, float cz1, float dcz0, unsigned char*__restrict__ vyfunc,int index,
                               float*__restrict__ vy, int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+NX]);
    float dfy;
    unsigned char dd = vyfunc[index];
    unsigned char dd1 = 1-dd;
    float cy0d = cy0*dd1+dcy0*dd;
    float cy1d = cy1*dd1;
    dfy = cy0d*(pressure[index+NX]-pressure[index   ]) +
    cy1d*(pressure[index+2*NX]-pressure[index-NX]);
    updateSpongeOne(vy,index,xTaper[ii],yTaperPH[jj],
                 zTaper[kk],dfy,rhoave);
  }
  
  //update vy for O2 sponge.  Only used on the very flanks of the model
  inline void updateAcousticVyO2(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float dcx0,
                                 float dcy0,float dcz0, int index,
                                 float*__restrict__ vy, int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+NX]);
    float dfy;
    dfy = dcy0*(pressure[index+NX]-pressure[index   ]);
    updateSpongeOne(vy,index,xTaper[ii],yTaperPH[jj],
                 zTaper[kk],dfy,rhoave);
  }
  
  //update vz for O4 sponge.
  inline void updateAcousticVz(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY,
                               float cx0, float cx1, float dcx0, 
                               float cy0, float cy1, float dcy0, float cz0, float cz1, float dcz0, unsigned char*__restrict__ vzfunc,int index,
                               float*__restrict__ vz, int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+NXY]);
    float dfz;
    unsigned char dd = vzfunc[index];
    unsigned char dd1 = 1-dd;
    float cz0d = cz0*dd1+dcz0*dd;
    float cz1d = cz1*dd1;
    dfz=cz0d*(pressure[index+NXY]-pressure[index   ]) +
    cz1d*(pressure[index+2*NXY]-pressure[index-NXY]);
    updateSpongeOne(vz,index,xTaper[ii],yTaper[jj],
                 zTaperPH[kk],dfz,rhoave);
  }
  
  //update vz for O2 sponge.  Only used on the very flanks of the model.
  inline void updateAcousticVzO2(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY,
                                 float dcx0, 
                                 float dcy0,float dcz0,int index,
                                 float*__restrict__ vz, int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+NXY]);
    float dfz;
    dfz=dcz0*(pressure[index+NXY]-pressure[index   ]);
    updateSpongeOne(vz,index,xTaper[ii],yTaper[jj],
                 zTaperPH[kk],dfz,rhoave);
  }
  
  //update vx in the interior (no sponge)
  inline void updateAcousticVxNoSponge(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float cx0, float cx1, float dcx0,
                                    float cy0, float cy1, float dcy0, float cz0, float cz1, float dcz0, unsigned char*__restrict__ vxfunc, int index,
                                    float*__restrict__ vx) {
    float rhoave=(rho[index]+rho[index+1]);
    unsigned char dd = vxfunc[index];
    unsigned char dd1 = 1-dd;
    float cx0d = cx0*dd1+dcx0*dd;
    float cx1d = cx1*dd1;
    vx[index]-=(cx0d*(pressure[index+1]-pressure[index   ]) +
                cx1d*(pressure[index+2]-pressure[index-1]))/rhoave;
  }    
  
  //update vy in the interior (no sponge)
  inline void updateAcousticVyNoSponge(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float cx0, float cx1, float dcx0,
                                    float cy0, float cy1, float dcy0, float cz0, float cz1, float dcz0, unsigned char*__restrict__ vyfunc, int index,
                                    float*__restrict__ vy) {
    float rhoave=(rho[index]+rho[index+NX]);
    unsigned char dd = vyfunc[index];
    unsigned char dd1 = 1-dd;
    float cy0d = cy0*dd1+dcy0*dd;
    float cy1d = cy1*dd1;
    vy[index]-=(cy0d*(pressure[index+NX]-pressure[index   ]) +
                cy1d*(pressure[index+2*NX]-pressure[index-NX]))/rhoave;
  }    
  
  //update vz in the interior (no sponge)
  inline void updateAcousticVzNoSponge(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY,
                                    float cx0, float cx1, float dcx0, 
                                    float cy0, float cy1, float dcy0, float cz0, float cz1, float dcz0, unsigned char*__restrict__ vzfunc,int index,
                                    float*__restrict__ vz) {
    float rhoave=(rho[index]+rho[index+NXY]);
    unsigned char dd = vzfunc[index];
    unsigned char dd1 = 1-dd;
    float cz0d = cz0*dd1+dcz0*dd;
    float cz1d = cz1*dd1;
    vz[index]-=(cz0d*(pressure[index+NXY]-pressure[index   ]) +
                cz1d*(pressure[index+2*NXY]-pressure[index-NXY]))/rhoave;
  }    
  
  //Update the pressure in the interior (no sponge).  ssfunc fills the same
  //role as vxfunc, etc above.  This is just the kernel.
  inline void updateAcousticPressureNoSpongeCore(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,
                                              float *__restrict__  bulk,int NX,int NXY,float cx0,float cx1, 
                                              float cy0,float cy1,float cz0,float cz1,float dcx0,float dcy0,float dcz0,unsigned char *__restrict__  ssfunc, int index, int i,
                                              float *__restrict__  pressure) {
    float dvx, dvy, dvz;
    float bb = bulk[index];
    unsigned char dd = ssfunc[index];
    unsigned char dd1= 1-dd;
    float cx0d = cx0*dd1+dcx0*dd;
    float cx1d = cx1*dd1;
    float cy0d = cy0*dd1+dcy0*dd;
    float cy1d = cy1*dd1;
    float cz0d = cz0*dd1+dcz0*dd;
    float cz1d = cz1*dd1;
    dvx = cx0d*(vx[index]-vx[index-1])+
    cx1d*(vx[index+1]-vx[index-2]);
    dvy = cy0d*(vy[index]-vy[index-NX])+
    cy1d*(vy[index+NX]-vy[index-2*NX]);
    dvz = cz0d*(vz[index]-vz[index-NXY])+
    cz1d*(vz[index+NXY]-vz[index-2*NXY]);
    
    float dvTot = dvx+dvy+dvz;
    pressure[index]-=dvTot*bb;
  }
  
  //Update the pressure in the interior (no sponge) for O2.  This is not actually ever called since O4 is used everywhere interior.  This is just the kernel.
  inline void updateAcousticPressureNoSpongeCoreO2(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,
                                                float *__restrict__  bulk,int NX,int NXY,float dcx0,float dcy0,float dcz0,
                                                int index, int i,
                                                float *__restrict__  pressure) {
    float dvx, dvy, dvz;
    float bb = bulk[index];
    dvx = dcx0*(vx[index]-vx[index-1]);
    dvy = dcy0*(vy[index]-vy[index-NX]);
    dvz = dcz0*(vz[index]-vz[index-NXY]);
    
    float dvTot = dvx+dvy+dvz;
    pressure[index]-=dvTot*bb;
  }
  
  //This is what is actually called from updateAcousticSpongeFull.cc.  This performs the interior (no sponge) pressure update
  inline void updateAcousticPressureNoSponge(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,
                                          int i0, int ie,
                                          float *__restrict__  bulk,int NX,int NXY,float cx0,float cx1, 
                                          float cy0,float cy1,float cz0,float cz1,float dcx0,float dcy0,float dcz0,unsigned char *__restrict__  ssfunc, int jkind,
                                          float *__restrict__  pressure) {
    for(int i=i0;i<ie;i++) {
      int index = i+jkind;
      updateAcousticPressureNoSpongeCore(vx,vy,vz,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,index,i,pressure);
    }
    
  }
  
  //This is not ever called, but is analogous to updateAcousticPressureNoSponge but for O2 updating
  inline void updateAcousticPressureNoSpongeO2(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,
                                            int i0, int ie,
                                            float *__restrict__  bulk,int NX,int NXY,float dcx0,float dcy0,float dcz0,int jkind,
                                            float *__restrict__  pressure) {
    for(int i=i0;i<ie;i++) {
      int index = i+jkind;
      updateAcousticPressureNoSpongeCoreO2(vx,vy,vz,bulk,NX,NXY,dcx0,dcy0,dcz0,index,i,pressure);
    }
    
  }
  
  //Update pressure with O4 updating in the sponge zone. This is just the kernel.
  inline void updateAcousticPressureSpongeCore(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,
                                            float *__restrict__  bulk,int NX,int NXY,float cx0,float cx1,
                                            float cy0,float cy1,float cz0,float cz1,float dcx0,float dcy0,float dcz0,unsigned char *__restrict__  ssfunc, int index, int i,
                                            float *__restrict__  pressure,
                                            float xtaper, float ytaper, float ztaper) {
    float dvx, dvy, dvz;
    float bb = bulk[index];
    unsigned char dd = ssfunc[index];
    unsigned char dd1= 1-dd;
    float cx0d = cx0*dd1+dcx0*dd;
    float cx1d = cx1*dd1;
    float cy0d = cy0*dd1+dcy0*dd;
    float cy1d = cy1*dd1;
    float cz0d = cz0*dd1+dcz0*dd;
    float cz1d = cz1*dd1;
    dvx = cx0d*(vx[index]-vx[index-1])+
    cx1d*(vx[index+1]-vx[index-2]);
    dvy = cy0d*(vy[index]-vy[index-NX])+
    cy1d*(vy[index+NX]-vy[index-2*NX]);
    dvz = cz0d*(vz[index]-vz[index-NXY])+
    cz1d*(vz[index+NXY]-vz[index-2*NXY]);

    float dvTot = dvx+dvy+dvz;
    pressure[index]-=dvTot*bb;
    pressure[index]*=xtaper*ytaper*ztaper;
  }
  
  //Update the pressure with O2 accuracy in the sponge.  This is only used on the very flank of the model.  This is just the kernel.
  inline void updateAcousticPressureSpongeCoreO2(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,
                                              float *__restrict__  bulk,int NX,int NXY,
                                              float dcx0,float dcy0,float dcz0,int index, int i,
                                              float *__restrict__  pressure,
                                              float xtaper, float ytaper, float ztaper) {
    float dvx, dvy, dvz;
    float bb = bulk[index];
    dvx = dcx0*(vx[index]-vx[index-1]);
    dvy = dcy0*(vy[index]-vy[index-NX]);
    dvz = dcz0*(vz[index]-vz[index-NXY]);
    
    float dvTot = dvx+dvy+dvz;
    pressure[index]-=dvTot*bb;
    pressure[index]*=xtaper*ytaper*ztaper;
  }
  
  //This is what is actually called from updateAcousticSpongeFull.cc for updating pressure in the sponge zone with O4 accuracy
  inline void updateAcousticPressureSponge(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,
                                        int i0, int ie,
                                        float *__restrict__  bulk,int NX,int NXY,float cx0,float cx1, 
                                        float cy0,float cy1,float cz0,float cz1,float dcx0,float dcy0,float dcz0,unsigned char *__restrict__  ssfunc, int jkind,
                                        float *__restrict__  pressure,
                                        float ytaper, float ztaper) {
    //the inline function below is the only difference compared to the no pml case
    for(int i=i0;i<ie;i++) {
      int index = i+jkind;
      updateAcousticPressureSpongeCore(vx,vy,vz,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,index,i,pressure,xTaper[i],ytaper,ztaper);
    }
  }
  
  //This is what is actually called from updateAcousticSpongeFull.cc for updating pressure in the sponge zone with O2 accuracy.  This is only used on the very flanks of the model.
  inline void updateAcousticPressureSpongeO2(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,
                                          int i0, int ie,
                                          float *__restrict__  bulk,int NX,int NXY, 
                                          float dcx0,float dcy0,float dcz0,int jkind,
                                          float *__restrict__  pressure,
                                          float ytaper, float ztaper) {
    //the inline function below is the only difference compared to the no pml case
    for(int i=i0;i<ie;i++) {
      int index = i+jkind;
      updateAcousticPressureSpongeCoreO2(vx,vy,vz,bulk,NX,NXY,dcx0,dcy0,dcz0,index,i,pressure,xTaper[i],ytaper,ztaper);
    }
  }
}

#endif
