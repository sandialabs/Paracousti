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
 *  This file contains the updating kernels that actually form the core
 *  of the code for the CPML boundary.  Both O2 and O4 updating for the 
 *  CPML zones and O4 for the interior zone.  These are called from 
 *  updateAcousticCPMLFull.cc
 *
 *  The first section includes the varibles that are actually defined 
 *  in acousticCPMLNameSpace.cc.
 *
 *  Defines several inline functions used in updateAcousticCPMLFull.cc
 *  or internally:
 *    sqr
 *    reflToD
 *    updatePMLMemVar
 *    updatePMLOne
 *    updatePMLMemVarComp
 *    updateAcousticVx
 *    updateAcousticVxO2
 *    updateAcousticVy
 *    updateAcousticVyO2
 *    updateAcousticVz
 *    updateAcousticVzO2
 *    updateAcousticVxNoPML
 *    updateAcousticVyNoPML
 *    updateAcousticVzNoPML
 *    updateAcousticPressureNoPMLCore
 *    updateAcousticPressureNoPMLCoreO2
 *    updateAcousticPressureNoPML
 *    updateAcousticPressureNoPMLO2
 *    updateAcousticPressurePMLCore
 *    updateAcousticPressurePMLCoreO2
 *    updateAcousticPressurePML
 *    updateAcousticPressurePMLO2
 */
#ifndef _ACOUSTIC_CPML_NS
#define _ACOUSTIC_CPML_NS

#include<math.h>

#include<stdio.h>

namespace acoustic_CPML {
  extern float *__restrict__ kxTaper, *__restrict__ kyTaper, *__restrict__ kzTaper;
  extern float *__restrict__ alphaXTaper, *__restrict__ alphaYTaper, *__restrict__ alphaZTaper;
  extern float *__restrict__ xTaper, *__restrict__ yTaper, *__restrict__ zTaper;
  extern float *__restrict__ axTaper, *__restrict__ ayTaper, *__restrict__ azTaper;
  extern float *__restrict__ bxTaper, *__restrict__ byTaper, *__restrict__ bzTaper;
  extern float *__restrict__ kxTaperPH, *__restrict__ kyTaperPH, *__restrict__ kzTaperPH;
  extern float *__restrict__ alphaXTaperPH, *__restrict__ alphaYTaperPH, *__restrict__ alphaZTaperPH;
  extern float *__restrict__ xTaperPH, *__restrict__ yTaperPH, *__restrict__ zTaperPH;
  extern float *__restrict__ axTaperPH, *__restrict__ ayTaperPH, *__restrict__ azTaperPH;
  extern float *__restrict__ bxTaperPH, *__restrict__ byTaperPH, *__restrict__ bzTaperPH;
  extern float *__restrict__ xTaperRel, *__restrict__ yTaperRel;
  extern float *__restrict__ xTaperRelPH, *__restrict__ yTaperRelPH;
  extern float *__restrict__ del;
  
  extern float *__restrict__ vxxXmin, *__restrict__ vxxXmax, *__restrict__ vxxYmin, *__restrict__ vxxYmax, *__restrict__ vxxZmin, *__restrict__ vxxZmax;
  extern float *__restrict__ vyyXmin, *__restrict__ vyyXmax, *__restrict__ vyyYmin, *__restrict__ vyyYmax, *__restrict__ vyyZmin, *__restrict__ vyyZmax;
  extern float *__restrict__ vzzXmin, *__restrict__ vzzXmax, *__restrict__ vzzYmin, *__restrict__ vzzYmax, *__restrict__ vzzZmin, *__restrict__ vzzZmax;
  
  extern float *__restrict__ pxXmin, *__restrict__ pxXmax, *__restrict__ pxYmin, *__restrict__ pxYmax, *__restrict__ pxZmin, *__restrict__ pxZmax;
  extern float *__restrict__ pyXmin, *__restrict__ pyXmax, *__restrict__ pyYmin, *__restrict__ pyYmax, *__restrict__ pyZmin, *__restrict__ pyZmax;
  extern float *__restrict__ pzXmin, *__restrict__ pzXmax, *__restrict__ pzYmin, *__restrict__ pzYmax, *__restrict__ pzZmin, *__restrict__ pzZmax;
  
  extern int nyInt, nzInt;
  
  extern int xMinSize, xMaxSize, yMinSize, yMaxSize, zMinSize, zMaxSize;
  extern int _nXmin, _nXmax, _nYmin, _nYmax, _nZmin, _nZmax;
  extern int _xMaxStart, _yMaxStart, _zMaxStart;
  extern int xStart1, xStart2, yStart1, yStart2, zStart, zStart_z;
  extern int xEnd, yEnd, zEnd, xEnd1, yEnd1, zEnd1;
  extern int _surfaceBCMode;
  extern int _kstart;
  
  //calculates the square of a number, saving some typing
  inline float sqr(float x) {return x*x;}
  //this function converts the input desired reflection coefficient (r)
  //to the sigma at the flank of the model using the theoretical
  //reflection coefficient for CPMLs with the width (w, in meters) and
  //the maximum vp in the model anywhere (vMax)
  inline float reflToD(float r, float w, float vMax) {
    return -1.5f*logf(r)*vMax/w;
  }
  
  //This returns the new value of CPML memory variable at this time step
  //using the previous memory varibale value (prev), the CPML parameters
  //at this point (ataper and btaper) and the spatial derivative of the
  //corresponding dependent variable (derivative)
  inline float updatePMLMemVar(float prev,float ataper,float btaper,float derivative) {
    return btaper*prev+ataper*derivative;
  }
  
  //Acutally update the memory variable (xxx) and velocities for velocity
  //components (vx) using CPML parameters axtaper, bxtaper, and kxtaper
  //at index in the vx array (which corresponds to pindex in the memory
  //variable array and the density (rho) and spatial derivatives of
  //the appropriate stresses (derX)
  inline void updatePMLOne(float*__restrict__ vx, int index,float*__restrict__ xxx, int pindex,float axtaper,float bxtaper, 
                           float kxtaper,float derX, float rho) {
    xxx[pindex]=updatePMLMemVar(xxx[pindex],axtaper,bxtaper,derX);
    vx[index] -= (derX*kxtaper+xxx[pindex])/rho;  //note: kxtaper is one over the original
  }
  
  //Actually update the memory variables (xxx, xyy and xzz) at pindex
  //using CPML parameters axtaper, aytaper, aztaper, etc. and spatial
  //derivatives of velocties (derX, derY, derZ).  Used for pressure
  //updates.
  inline void updatePMLMemVarComp(float*__restrict__ xxx, float*__restrict__ xyy,float*__restrict__ xzz,int pindex,float axtaper,float aytaper, float aztaper, float bxtaper, float bytaper,
                                  float bztaper, float derX, float derY, float derZ) {
    xxx[pindex]=updatePMLMemVar(xxx[pindex],axtaper,bxtaper,derX);
    xyy[pindex]=updatePMLMemVar(xyy[pindex],aytaper,bytaper,derY);
    xzz[pindex]=updatePMLMemVar(xzz[pindex],aztaper,bztaper,derZ);
  }    
  
  //Update the vx for O4 CPML.  The c variables refer to O4 FD coefficients
  //while the d variables refer to O2 coefficients.  High contrast
  //interfaces trigger the change from O4 to O2.  This is determined
  //in acoustic_control.cc and a flag array vxfunc holds this info.  px
  //is the CPML memory variable and ii, jj and kk refer to the 3-D indices
  //within vx
  inline void updateAcousticVx(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float cx0, float cx1, float dcx0, 
                               float cy0, float cy1, float dcy0, float cz0, float cz1, float dcz0, unsigned char*__restrict__ vxfunc, int pIndex0, int index0,
                               float*__restrict__ px, float*__restrict__ vx, int i0, int ie, int jj, int kk) {
    for(int i=i0;i<ie;++i) {
      int index = index0+i;
      unsigned char dd = vxfunc[index];
      unsigned char dd1 = 1-dd;
      float cx0d = cx0*dd1+dcx0*dd;
      float cx1d = cx1*dd1;
      del[i] = cx0d*(pressure[index+1]-pressure[index   ]) +
        cx1d*(pressure[index+2]-pressure[index-1]);
    }
    for(int i=i0;i<ie;++i) {
      int pIndex = pIndex0+i;
      px[pIndex] = updatePMLMemVar(px[pIndex],axTaperPH[i],bxTaperPH[i],del[i]);
    }
    for(int i=i0;i<ie;++i) {
      int index = index0+i;
      int pIndex = pIndex0+i;
      float rhoave=(rho[index]+rho[index+1]);
      vx[index] -= (del[i]*kxTaperPH[i]+px[pIndex])/rhoave;
    }
  }
  
  //Update the vx for O2 CPML.  This is only used near the very flanks of
  //the model
  inline void updateAcousticVxO2(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float dcx0, 
                                 float dcy0,float dcz0, int pIndex, int index,
                                 float*__restrict__ px, float*__restrict__ vx, int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+1]);
    float dfx;
    dfx = dcx0*(pressure[index+1]-pressure[index   ]);
    updatePMLOne(vx,index,px,pIndex,axTaperPH[ii],bxTaperPH[ii],
                 kxTaperPH[ii],dfx,rhoave);
  }    
  
  //Update vy for O4 CPML.
  inline void updateAcousticVy(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float cx0, float cx1, float dcx0, 
                               float cy0, float cy1, float dcy0, float cz0, float cz1, float dcz0, unsigned char*__restrict__ vyfunc,int pIndex, int index,
                               float*__restrict__ py, float*__restrict__ vy, int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+NX]);
    float dfy;
    unsigned char dd = vyfunc[index];
    unsigned char dd1 = 1-dd;
    float cy0d = cy0*dd1+dcy0*dd;
    float cy1d = cy1*dd1;
    dfy = cy0d*(pressure[index+NX]-pressure[index   ]) +
    cy1d*(pressure[index+2*NX]-pressure[index-NX]);
    updatePMLOne(vy,index,py,pIndex,ayTaperPH[jj],byTaperPH[jj],
                 kyTaperPH[jj],dfy,rhoave);
  }    
  
  //update vy for O2 CPML.  Only used on the very flanks of the model
  inline void updateAcousticVyO2(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float dcx0, 
                                 float dcy0,float dcz0, int pIndex, int index,
                                 float*__restrict__ py, float*__restrict__ vy, int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+NX]);
    float dfy;
    dfy = dcy0*(pressure[index+NX]-pressure[index   ]);
    updatePMLOne(vy,index,py,pIndex,ayTaperPH[jj],byTaperPH[jj],
                 kyTaperPH[jj],dfy,rhoave);
  }    
  
  //update vz for O4 CPML.
  inline void updateAcousticVz(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, 
                               float cx0, float cx1, float dcx0, 
                               float cy0, float cy1, float dcy0, float cz0, float cz1, float dcz0, unsigned char*__restrict__ vzfunc,int pIndex, int index,
                               float*__restrict__ pz,float*__restrict__ vz, int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+NXY]);
    float dfz;
    unsigned char dd = vzfunc[index];
    unsigned char dd1 = 1-dd;
    float cz0d = cz0*dd1+dcz0*dd;
    float cz1d = cz1*dd1;
    dfz=cz0d*(pressure[index+NXY]-pressure[index   ]) +
    cz1d*(pressure[index+2*NXY]-pressure[index-NXY]);
    updatePMLOne(vz,index,pz,pIndex,azTaperPH[kk],bzTaperPH[kk],
                 kzTaperPH[kk],dfz,rhoave);
  }    
  
  //update vz for O2 CPML.  Only used on the very flanks of the model.
  inline void updateAcousticVzO2(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, 
                                 float dcx0, 
                                 float dcy0,float dcz0,int pIndex, int index,
                                 float*__restrict__ pz,float*__restrict__ vz, int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+NXY]);
    float dfz;
    dfz=dcz0*(pressure[index+NXY]-pressure[index   ]);
    updatePMLOne(vz,index,pz,pIndex,azTaperPH[kk],bzTaperPH[kk],
                 kzTaperPH[kk],dfz,rhoave);
  }   
  
  //update vx in the interior (no CPML)
  inline void updateAcousticVxNoPML(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float cx0, float cx1, float dcx0, 
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
  
  //update vy in the interior (no CPML)
  inline void updateAcousticVyNoPML(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float cx0, float cx1, float dcx0, 
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
  
  //update vz in the interior (no CPML)
  inline void updateAcousticVzNoPML(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, 
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
  
  //Update the pressure in the interior (no CPML).  ssfunc fills the same
  //role as vxfunc, etc above.  This is just the kernel.
  inline void updateAcousticPressureNoPMLCore(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,                                             
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
  
  //Update the pressure in the interior (no CPML) for O2.  This is not actually ever called since O4 is used everywhere interior.  This is just the kernel.
  inline void updateAcousticPressureNoPMLCoreO2(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,                                             
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
  
  //This is what is actually called from updateAcousticCPMLFull.cc.  This performs the interior (no CPML) pressure update
  inline void updateAcousticPressureNoPML(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz, 
                                          int i0, int ie,
                                          float *__restrict__  bulk,int NX,int NXY,float cx0,float cx1, 
                                          float cy0,float cy1,float cz0,float cz1,float dcx0,float dcy0,float dcz0,unsigned char *__restrict__  ssfunc, int jkind,
                                          float *__restrict__  pressure) {
    for(int i=i0;i<ie;i++) {
      int index = i+jkind;
      updateAcousticPressureNoPMLCore(vx,vy,vz,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,index,i,pressure);
    }
    
  }
  
  //This is not ever called, but is analogous to updateAcousticPressureNoPML but for O2 updating
  inline void updateAcousticPressureNoPMLO2(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz, 
                                            int i0, int ie,
                                            float *__restrict__  bulk,int NX,int NXY,float dcx0,float dcy0,float dcz0,int jkind,
                                            float *__restrict__  pressure) {
    for(int i=i0;i<ie;i++) {
      int index = i+jkind;
      updateAcousticPressureNoPMLCoreO2(vx,vy,vz,bulk,NX,NXY,dcx0,dcy0,dcz0,index,i,pressure);
    }
    
  }
  
  //Update pressure with O4 updating in the CPML zone. CPML memory
  //variables are vxx, vyy and vzz, otherwise the variables have similar
  //meaning to above calls.  This is just the kernel.
  inline void updateAcousticPressurePMLCore(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,
                                            float *__restrict__ vxx,float *__restrict__ vyy,float *__restrict__ vzz,
                                            float *__restrict__  bulk,int NX,int NXY,float cx0,float cx1, 
                                            float cy0,float cy1,float cz0,float cz1,float dcx0,float dcy0,float dcz0,unsigned char *__restrict__  ssfunc, int index, int pIndex, int i,
                                            float *__restrict__  pressure,
                                            float axtaper, float aytaper, float aztaper,
                                            float bxtaper, float bytaper, float bztaper,
                                            float kxtaper, float kytaper, float kztaper) {
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

    //PML updating in group of lines below
    updatePMLMemVarComp(vxx,vyy,vzz,pIndex,axtaper,aytaper,aztaper,
                        bxtaper,bytaper,bztaper,
                        dvx,dvy,dvz);
    dvx = dvx*kxtaper+vxx[pIndex];
    dvy = dvy*kytaper+vyy[pIndex];
    dvz = dvz*kztaper+vzz[pIndex];
    
    float dvTot = dvx+dvy+dvz;
    pressure[index]-=dvTot*bb;
  }
  
  //Update the pressure with O2 accuracy in the CPML.  This is only used on the very flank of the model.  This is just the kernel.
  inline void updateAcousticPressurePMLCoreO2(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,
                                              float *__restrict__ vxx,float *__restrict__ vyy,float *__restrict__ vzz,
                                              float *__restrict__  bulk,int NX,int NXY, 
                                              float dcx0,float dcy0,float dcz0,int index, int pIndex, int i,
                                              float *__restrict__  pressure,
                                              float axtaper, float aytaper, float aztaper,
                                              float bxtaper, float bytaper, float bztaper,
                                              float kxtaper, float kytaper, float kztaper) {
    float dvx, dvy, dvz;
    float bb = bulk[index];
    dvx = dcx0*(vx[index]-vx[index-1]);
    dvy = dcy0*(vy[index]-vy[index-NX]);
    dvz = dcz0*(vz[index]-vz[index-NXY]);
    
    //PML updating in group of lines below
    updatePMLMemVarComp(vxx,vyy,vzz,pIndex,axtaper,aytaper,aztaper,
                        bxtaper,bytaper,bztaper,
                        dvx,dvy,dvz);
    dvx = dvx*kxtaper+vxx[pIndex];
    dvy = dvy*kytaper+vyy[pIndex];
    dvz = dvz*kztaper+vzz[pIndex];
    
    float dvTot = dvx+dvy+dvz;
    pressure[index]-=dvTot*bb;
  }
  
  //This is what is actually called from updateAcousticCPMLFull.cc for updating pressure in the CPML zone with O4 accuracy
  inline void updateAcousticPressurePML(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz, 
                                        float *__restrict__ vxx, float *__restrict__ vyy, float* __restrict__ vzz,
                                        int i0, int ie,
                                        float *__restrict__  bulk,int NX,int NXY,float cx0,float cx1, 
                                        float cy0,float cy1,float cz0,float cz1,float dcx0,float dcy0,float dcz0,unsigned char *__restrict__  ssfunc, int jkind, int jkpIndex,
                                        float *__restrict__  pressure,
                                        float aytaper, float aztaper, float bytaper, float bztaper, float kytaper, float kztaper) {
    //the inline function below is the only difference compared to the no pml case
    for(int i=i0;i<ie;i++) {
      int index = i+jkind;
      int pIndex = i+jkpIndex;
      updateAcousticPressurePMLCore(vx,vy,vz,vxx,vyy,vzz,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,index,pIndex,i,pressure,axTaper[i],aytaper,aztaper,bxTaper[i],bytaper,bztaper,kxTaper[i],kytaper,kztaper);
    }
  }
  
  //This is what is actually called from updateAcousticCPMLFull.cc for updating pressure in the CPML zone with O2 accuracy.  This is only used on the very flanks of the model.
  inline void updateAcousticPressurePMLO2(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz, 
                                          float *__restrict__ vxx, float *__restrict__ vyy, float* __restrict__ vzz,
                                          int i0, int ie,
                                          float *__restrict__  bulk,int NX,int NXY, 
                                          float dcx0,float dcy0,float dcz0,int jkind, int jkpIndex,
                                          float *__restrict__  pressure,
                                          float aytaper, float aztaper, float bytaper, float bztaper, float kytaper, float kztaper) {
    //the inline function below is the only difference compared to the no pml case
    for(int i=i0;i<ie;i++) {
      int index = i+jkind;
      int pIndex = i+jkpIndex;
      updateAcousticPressurePMLCoreO2(vx,vy,vz,vxx,vyy,vzz,bulk,NX,NXY,dcx0,dcy0,dcz0,index,pIndex,i,pressure,axTaper[i],aytaper,aztaper,bxTaper[i],bytaper,bztaper,kxTaper[i],kytaper,kztaper);
    }
  }
}

#endif
