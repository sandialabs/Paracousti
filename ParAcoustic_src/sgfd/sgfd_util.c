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
 *  sgfd_util.c
 *
 *
 *  Define interpolation and extrapolation routines used primarily
 *  by source and receiver classes.
 *
 *  Defines the following functions:
 *  trilinCoeff
 *  trilinInterp
 *  trilinInterpOffset
 *  trilinExtrap
 *  triCubicCoeff
 *  triCubicInterp
 *  triCubicExtrap
 *
 */

#include <math.h>
#include <stdio.h>

#include "nstdutil.h"
#include "sgfd.h"

/*Subroutine TRILIN_COEFF computes the eight coefficients of a  */
/* trilinear interpolator appropriate for a specified receiver  */
/* location. */
void trilinCoeff(float xr,float yr,float zr,
                 float xmin,float dx,float ymin,float dy,float zmin,float dz,
                 interpStruct* interpCoeff){
  //Compute spatial indeces of cell that surrounds receiver.
  int i=(int)floor((xr-xmin)/dx);
  int j=(int)floor((yr-ymin)/dy);
  int k=(int)floor((zr-zmin)/dz);
  
  //Compute coordinates of corner of cell.
  float xi=xmin+i*dx;
  float yj=ymin+j*dy;
  float zk=zmin+k*dz;
  float p,q,r;
  
  //Store spatial indeces of the eight gridpoints surrounding receiver.
  interpCoeff->iptr[0]=i;
  interpCoeff->jptr[0]=j;
  interpCoeff->kptr[0]=k;
  
  interpCoeff->iptr[1]=i+1;
  interpCoeff->jptr[1]=j;
  interpCoeff->kptr[1]=k;
  
  interpCoeff->iptr[2]=i;
  interpCoeff->jptr[2]=j+1;
  interpCoeff->kptr[2]=k;
  
  interpCoeff->iptr[3]=i+1;
  interpCoeff->jptr[3]=j+1;
  interpCoeff->kptr[3]=k;
  
  interpCoeff->iptr[4]=i;
  interpCoeff->jptr[4]=j;
  interpCoeff->kptr[4]=k+1;
  
  interpCoeff->iptr[5]=i+1;
  interpCoeff->jptr[5]=j;
  interpCoeff->kptr[5]=k+1;
  
  interpCoeff->iptr[6]=i;
  interpCoeff->jptr[6]=j+1;
  interpCoeff->kptr[6]=k+1;
  
  interpCoeff->iptr[7]=i+1;
  interpCoeff->jptr[7]=j+1;
  interpCoeff->kptr[7]=k+1;
  
  //Compute and store trilinear interpolator coefficients.
  p=(xr-xi)/dx;
  q=(yr-yj)/dy;
  r=(zr-zk)/dz;
  
  interpCoeff->coeff[0]=(1.0-p)*(1.0-q)*(1.0-r);
  interpCoeff->coeff[1]= p     *(1.0-q)*(1.0-r);
  interpCoeff->coeff[2]=(1.0-p)* q     *(1.0-r);
  interpCoeff->coeff[3]= p     * q     *(1.0-r);
  interpCoeff->coeff[4]=(1.0-p)*(1.0-q)* r;
  interpCoeff->coeff[5]= p     *(1.0-q)* r;
  interpCoeff->coeff[6]=(1.0-p)* q     * r;
  interpCoeff->coeff[7]= p     * q     * r;
}

/*Subroutine TRILIN_INTERP computes an approximate value of an input */
/* 3D array using trilinear interpolation. */
float trilinInterp(modelDefStruct* modelDef,
                   interpStruct* interpCoeff,
                   float* vx){
  DEF_MODEL_SIZE(modelDef);
  
  //Compute trilinear approximation to the array value.
  float value=0.0,totalCoeff=0.0;
  int n;
  for(n=0;n<8;n++){
    int i=interpCoeff->iptr[n];
    int j=interpCoeff->jptr[n];
    int k=interpCoeff->kptr[n];
    if(ISMID(0,i,NX-1) &&
       ISMID(0,j,NY-1) &&
       ISMID(0,k,NZ-1)){
      value+=interpCoeff->coeff[n]*VX(i,j,k);
      totalCoeff+=interpCoeff->coeff[n];
    }
  }
  return totalCoeff>0.0?value/totalCoeff:0.0;
}
/*Subroutine TRILIN_INTERP_Offset computes an approximate value of an input */
/* 3D array using trilinear interpolation. */
float trilinInterpOffset(modelDefStruct* modelDef,
                         interpStruct* interpCoeff,
                         int xoff,int yoff,int zoff,
                         float* vx){
  DEF_MODEL_SIZE(modelDef);
  
  //Compute trilinear approximation to the array value.
  float value=0.0,totalCoeff=0.0;
  int n;
  for(n=0;n<8;n++){
    int i=interpCoeff->iptr[n]+xoff;
    int j=interpCoeff->jptr[n]+yoff;
    int k=interpCoeff->kptr[n]+zoff;
    if(ISMID(0,i,NX-1) &&
       ISMID(0,j,NY-1) &&
       ISMID(0,k,NZ-1)){
      value+=interpCoeff->coeff[n]*VX(i,j,k);
      totalCoeff+=interpCoeff->coeff[n];
    }
  }
  return totalCoeff>0.0?value/totalCoeff:0.0;
}

/*Subroutine TRILIN_EXTRAP updates a 3D array by extrapolating an */
/* input scalar value to the eight surrounding gridpoints. */
float trilinExtrap(modelDefStruct* modelDef,
                   float* vx,
                   interpStruct* interpCoeff,
                   float value){
  DEF_MODEL_SIZE(modelDef);
  
  //Extrapolate scalar value to the eight surrounding gridpoints.
  int n;
  float totalAdded=0.0;
  for(n=0;n<8;n++){
    int i=interpCoeff->iptr[n];
    int j=interpCoeff->jptr[n];
    int k=interpCoeff->kptr[n];
    if(ISMID(0,i,NX-1) &&
       ISMID(0,j,NY-1) &&
       ISMID(0,k,NZ-1)){
      VX(i,j,k)+=value*interpCoeff->coeff[n];
      totalAdded+=fabs(value*interpCoeff->coeff[n]);
    }
  }
  return totalAdded;
}

//
//Following tricubic routines are not used for anything (yet)
// Try them on the disimalar grid problem
//

//Subroutine TRICUB_COEFF computes the sixty-four coefficients of
// a tricubic polynomial interpolator appropriate for a specified
// location.
void triCubicCoeff(float xr,float yr,float zr,
                   float xmin,float dx,float ymin,float dy,float zmin,float dz,
                   cubicInterpStruct *cubicInterp){
  //temp variables to save memory dereferencing
  float* coeff=cubicInterp->coeff;
  
  //Compute spatial indeces of cell that surrounds receiver.
  int i=(int)floor((xr-xmin)/dx);
  int j=(int)floor((yr-ymin)/dy);
  int k=(int)floor((zr-zmin)/dz);
  
  //Compute coordinates of corner of cell.
  float xi=xmin+i*dx;
  float yj=ymin+j*dy;
  float zk=zmin+k*dz;
  float p,q,r;
  float a[4],b[4],c[4];
  int kk,jj,ii;
  
  //Compute and store tricubic polynomial interpolator coefficients
  cubicInterp->shortCircut=FALSE;
  
  p=(xr-xi)/dx;
  q=(yr-yj)/dy;
  r=(zr-zk)/dz;
  
  a[0]=-(p-2.0)*(p-1.0)*(p-0.0)/6.0;
  a[1]= (p-2.0)*(p-1.0)*(p+1.0)/2.0;
  a[2]=-(p-2.0)*(p-0.0)*(p+1.0)/2.0;
  a[3]= (p-0.0)*(p-1.0)*(p+1.0)/6.0;
  
  b[0]=-(q-2.0)*(q-1.0)*(q-0.0)/6.0;
  b[1]= (q-2.0)*(q-1.0)*(q+1.0)/2.0;
  b[2]=-(q-2.0)*(q-0.0)*(q+1.0)/2.0;
  b[3]= (q-0.0)*(q-1.0)*(q+1.0)/6.0;
  
  c[0]=-(r-2.0)*(r-1.0)*(r-0.0)/6.0;
  c[1]= (r-2.0)*(r-1.0)*(r+1.0)/2.0;
  c[2]=-(r-2.0)*(r-0.0)*(r+1.0)/2.0;
  c[3]= (r-0.0)*(r-1.0)*(r+1.0)/6.0;
  
  //fill in
  cubicInterp->iptr=i-1;
  cubicInterp->jptr=j-1;
  cubicInterp->kptr=k-1;
  
  for(kk=0;kk<4;kk++){
    for(jj=0;jj<4;jj++){
      for(ii=0;ii<4;ii++){
        int index=ii+4*jj+16*kk;
        coeff[index]=a[ii]*b[jj]*c[kk];
        
        if(coeff[index]>=0.9999999){//Close enough to 1.0.
          cubicInterp->shortCircut=TRUE;
          cubicInterp->iptr+=ii;
          cubicInterp->jptr+=jj;
          cubicInterp->kptr+=kk;
          return;
        }
      }
    }
  }
}

float triCubicInterp(modelDefStruct* modelDef,
                     cubicInterpStruct* interpCoeff,
                     float* vx){
  DEF_MODEL_SIZE(modelDef);
  if(interpCoeff->shortCircut){
    return vx[interpCoeff->iptr+interpCoeff->jptr*NX+interpCoeff->kptr*NXY];
  }else{
    float* coeff=interpCoeff->coeff;
    int i,j,k,ii,jj,kk;
    float val=0.0;
    
    for(k=interpCoeff->kptr,kk=0;kk<4;k++,kk++){
      for(j=interpCoeff->jptr,jj=0;jj<4;j++,jj++){
        for(i=interpCoeff->iptr,ii=0;ii<4;i++,ii++){
          int cI=ii+4*jj+16*kk,mI=i+j*NX+k*NXY;
          val+=vx[mI]*coeff[cI];
        }
      }
    }
    return val;
  }
}
float triCubicExtrap(modelDefStruct* modelDef,
                     cubicInterpStruct* interpCoeff,
                     float* vx,
                     float value){
  DEF_MODEL_SIZE(modelDef);
  
  int i=interpCoeff->iptr;
  int j=interpCoeff->jptr;
  int k=interpCoeff->kptr;
  if(interpCoeff->shortCircut){
    if(ISMID(0,i,NX-1) &&
       ISMID(0,j,NY-1) &&
       ISMID(0,k,NZ-1)){
      return
      vx[i+j*NX+k*NXY]=value;
    }else{
      return value;
    }
  }else{
    float* coeff=interpCoeff->coeff;
    int ii,jj,kk;
    int cI=0;
    for(kk=0;kk<4;kk++){
      for(jj=0;jj<4;jj++){
        for(ii=0;ii<4;ii++){
          if(ISMID(0,i+ii,NX-1) &&
             ISMID(0,j+jj,NY-1) &&
             ISMID(0,k+kk,NZ-1)){
            int mI=(i+ii)+(j+jj)*NX+(k+kk)*NXY;
            vx[mI]+= value*coeff[cI++];
          }
        }
      }
    }
    return value;
  }
}

