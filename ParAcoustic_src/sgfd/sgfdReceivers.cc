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
 *  sgfdReceivers.cc
 *
 *
 *  Defines class functions used for receivers.  Each object in these classes represent a
 *  single receiver and are divided into the basic receiver types.  These classes handle
 *  interpolating data onto receiver location and time and any conversions such as
 *  integration.
 *
 *  Defines the following class functions:
 *  receiver::receiver
 *  receiver::doCheckpoint
 *  receiver::readCheckpoint
 *  receiver::setActive
 *  receiver::packInit
 *  receiver::packData
 *  receiver::allocateData
 *  velocityReceiver::velocityReceiver
 *  velocityReceiver::velocityReceiver
 *  velocityReceiver::packInit
 *  velocityReceiver::fill
 *  velocityReceiver::scaleData
 *  velocityReceiver::allocateData
 *  velocityReceiver::interpolate
 *  velocityReceiver::integrate
 *  velocityReceiver::differentiate
 *  pressureReceiver::pressureReceiver
 *  pressureReceiver::pressureReceiver
 *  pressureReceiver::packInit
 *  pressureReceiver::fill
 *  pressureReceiver::scaleData
 *  pressureReceiver::allocateData
 *
 */
//These utility functions are for receivers in the SGFD problem

#include "sgfdReceivers.hh"
#include <time.h>

#include "netcdf.h"
#include "nstdutil.hh"
#include "io_procs.h"
#include "message_passing.h"

//basic constructor with location read from message
receiver::receiver(modelDefStruct* modelDef){
  unpackMessage("ffff",&_x,&_y,&_z,&_amp);
  _bx=_by=_bz=0.0; //may not be used but must be present in all forms of receiver
  _data=NULL;
}

///Method for checkpoint writing of receiver data.
void receiver::doCheckpoint(int iteration,FILE* cpFile){
  //Save the current data.
  fwrite(_data,sizeof(float),iteration,cpFile);
}

///Method for checkpoint reading of receiver data.
void receiver::readCheckpoint(int iteration,FILE* cpFile){
  //Read the saved data.
  fread(_data,sizeof(float),iteration,cpFile);
}
int receiver::setActive(modelDefStruct* modelDef){
  if(!modelDef)
    return (_data=NULL)!=NULL; //This evaluates to FALSE without a conversion

  DEF_MODEL_SIZE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  if(!ISMID_CO(minX+dx,_x,minX+dx*(NX-3)) ||
     !ISMID_CO(minY+dy,_y,minY+dy*(NY-3)) ||
     !ISMID_CO(minZ+dz,_z,minZ+dz*(NZ-3))){
    return (_data=NULL)!=NULL; //use this so the active() function from parent class still works
  }else{
    return TRUE;
  }
}

//Pack a message with the parameters for this receiver. Assumes the buffer is ready to 
// be filled.
void receiver::packInit(){
  //These features are passed for all receiver types
  packMessage("i",type());
  packMessage("ffff",_x,_y,_z,_amp);
}

//this does whatever work is required
// to put the output into real units. This includes interpolation for the velocity types
// where the data is on a half integer time raster.
int receiver::packData(modelDefStruct* modelDef,float* wrk,
         int startI,int stopI,int decimation){
  if(decimation==1){
    //This is the easy case, pack the complete receiver.

    //First do the required scaling to physical units. This must be done on the 
    // complete data so interpolation can be done correctly.
    scaleData(modelDef,startI,stopI,_data,wrk);

    packMessage("F",wrk,stopI-startI+1);
  }else{
    //This is more complicated. Determine the number of samples that will actually
    // be sent, then shift these samples to the front of the wrk array.
    int firstI=(int)((float)(startI)/decimation);
    int nSamps=(int)((float)(stopI)/decimation)-firstI+1;
    startI=firstI*decimation; //Reset the starting index to make sure we get the
    //                           lowest needed sample.
    stopI=startI+(nSamps-1)*decimation+1; //Reset the stop index
    //                           to make sure that the last 
    //                           good sample is the same as the
    //                           next start sample (critical if
    //                           the trace is being integrated).

    //First do the required scaling to physical units. This must be done on the 
    // complete data so interpolation can be done correctly.
    scaleData(modelDef,startI,stopI,_data,wrk);

    //Repack the data.
    for(int i=0;i<nSamps;i++){ 
      wrk[i]=wrk[i*decimation];
    }
    //Now pack the required number of samples only.
    packMessage("F",wrk,nSamps);
  }
  return stopI-startI;
}

//Allocate and intialize the data, also fill in any required interpStructs
void receiver::allocateData(modelDefStruct* modelDef){
  DEF_MODEL_SIZE(modelDef);
  assert((_data=(float*)malloc((NT+2)*sizeof(float)))!=NULL,
   "receiver::receiver(%f,%f,%f)--unable to allocate %i floats for data",
   _x,_y,_z,NT);
  for(int j=0;j<NT+2;_data[j++]=0.0);
}

//Add the include for receiverGroups here.
#include "sgfdReceiverGroups.hh"

//Constructor with parameters read from message
velocityReceiver::velocityReceiver(modelDefStruct* modelDef,int allocate):
  receiver(modelDef){
  //read the additional parameters
  unpackMessage("ifffii",&_integrate,&_bx,&_by,&_bz,&_useCubicInterp,&_surfaceReceiver);
    
  if(_surfaceReceiver) _useCubicInterp = TRUE;
    
  if(setActive(modelDef)){
    if(allocate){
      allocateData(modelDef); //this allocates, initializes and fills in any required interp
    }else{
      _data=(float*)0x1; //Flag so this receiver still registers as active.
    }
  }
}

//Constructor with parameters as arguments
velocityReceiver::velocityReceiver(modelDefStruct* modelDef,int allocate,
     float x,float y,float z,
     float bx,float by,float bz,float amp,int useCubic,int surfaceReceiver):
  receiver(modelDef,x,y,z,amp){
  //fill in the additional parameters
  _integrate=FALSE;
  _useCubicInterp = useCubic;
  _surfaceReceiver = surfaceReceiver;
  if(_surfaceReceiver) _useCubicInterp = TRUE;

  _bx=bx;
  _by=by;
  _bz=bz;

  if(setActive(modelDef)){
    if(allocate){
allocateData(modelDef); //this allocates, initializes and fills in any required interp
    }else{
_data=(float*)0x1;
    }
  }
}

//Pack a message with the parameters for this receiver. Assumes the buffer is ready to
// be filled.
void velocityReceiver::packInit(){
  receiver::packInit();
  //And pack the integrate flag and the direction cosines
  packMessage("ifffii",_integrate,_bx,_by,_bz,_useCubicInterp,_surfaceReceiver);
}

//calculate the value of this receiver for the current iteration
float velocityReceiver::fill(modelDefStruct* modelDef,int iteration,
       float* vx,float* vy,float* vz,
       float* xx,float* yy,float* zz,
       float* xy,float* xz,float* yz){
  if(!active()) return -10;

  float ampx, ampy, ampz;
  //Interpolate components of velocity vector onto receiver location.
  if(_useCubicInterp) {
    ampx=triCubicInterp(modelDef,_xCubicInterp,vx);
    ampy=triCubicInterp(modelDef,_yCubicInterp,vy);
    ampz=triCubicInterp(modelDef,_zCubicInterp,vz);
  } else {
    ampx=trilinInterp(modelDef,_xInterp,vx);
    ampy=trilinInterp(modelDef,_yInterp,vy);
    ampz=trilinInterp(modelDef,_zInterp,vz);
  }

  return _data[iteration]=_amp*(_bx*ampx+_by*ampy+_bz*ampz);
}

int velocityReceiver::scaleData(modelDefStruct* modelDef,int startI,int stopI,
    float* data,float* wrk){
  if(!_integrate){
    interpolate(modelDef,data,wrk,startI,stopI);
  }else if(_integrate > 0){
    integrate(modelDef,data,wrk,startI,stopI);
  }else{
    differentiate(modelDef,data,wrk,startI,stopI);
  }
  return stopI-startI;
}
//Fill in any required interpStructs
void velocityReceiver::allocateData(modelDefStruct* modelDef){
  receiver::allocateData(modelDef);

  DEF_MODEL_LIMITS(modelDef);
  
  if(_useCubicInterp) {
    _xCubicInterp = new cubicInterpStruct;
    _yCubicInterp = new cubicInterpStruct;
    _zCubicInterp = new cubicInterpStruct;
    if(!_surfaceReceiver) {
      triCubicCoeff(_x,_y,_z,
                    minX+0.5*dx,dx,minY,dy,minZ,dz,
                    _xCubicInterp);
      triCubicCoeff(_x,_y,_z,
                  minX,dx,minY+0.5*dy,dy,minZ,dz,
                  _yCubicInterp);
      triCubicCoeff(_x,_y,_z,
                  minX,dx,minY,dy,minZ+0.5*dz,dz,
                  _zCubicInterp);
    }
  } else {
    _xInterp = new interpStruct;
    _yInterp = new interpStruct;
    _zInterp = new interpStruct;
    //fill the interpolators for vx, vy, and vz components
    trilinCoeff(_x,_y,_z,
                minX+0.5*dx,dx,minY,dy,minZ,dz,
                _xInterp);
    trilinCoeff(_x,_y,_z,
                minX,dx,minY+0.5*dy,dy,minZ,dz,
                _yInterp);
    trilinCoeff(_x,_y,_z,
                minX,dx,minY,dy,minZ+0.5*dz,dz,
                _zInterp);
  }

}

//Velocities are recorded on the 1/2 integer time raster; move to the integer
// raster by interpolation
float* velocityReceiver::interpolate(modelDefStruct* modelDef,float* data,float* out,
       int startI,int stopI){
  DEF_MODEL_SCALE(modelDef);

  //First value is obtained by linear approximation to the first
  // two points.
  int outIndex=0;
  
  if(!startI)
    out[outIndex++]=scalarVel*
      0.5*data[0];
//	(3.0*data[0]+6.0*data[1]-data[2])/8.0;

  int NT=modelDef->NT;
  int start=MAX(1,startI),stop=MIN(NT-3,stopI);
  for(int l=start;l<=stop;l++){
    //Interior values are determined by linear interpolation.
    out[outIndex++]=scalarVel*0.5*(data[l]+data[l-1]);
//      out[outIndex++]=scalarVel*(9.0*(data[l  ]+data[l+1]) -
//				 (data[l-1]+data[l+2]))/16.0;
  }

  //Last values are obtained by linear approximation to the final points.
  for(int i=stop+1;i<stopI;i++,outIndex++){
    out[outIndex]=0.5*scalarVel*(data[i]+data[i-1]);
//      out[outIndex]=2.0*out[outIndex-1]-out[outIndex-2];
  }

  return out;
}
  
//Particle displacement values are calculated by integrating
// the particle velocity samples over time.  A local parabolic
// polynomial is used to approximate the velocity functon.
// Calculated displacement samples are defined on the integer 
// time raster t=tmin+(l-1)*dt for l=1,2,...,nt.
float* velocityReceiver::integrate(modelDefStruct* modelDef,float* data,float* out,
     int startI,int stopI){
  DEF_MODEL_LIMITS(modelDef);
  DEF_MODEL_SCALE(modelDef);

  //First value (t=tmin) is zeroed.
  int outIndex=0;
  if(!startI){
    out[outIndex++]=0.; //scalarVel*dt*(119.0*data[0]+107.0*data[1]-
          //43.0 *data[2]+9.0  *data[3])/384.0; old way
  }else{
    out[outIndex++]=_lastOut;
  }

  //Interior values are obtained with a centered operator.  
  int start=MAX(1,startI+1),stop=MIN(modelDef->NT,stopI);
  for(int i=start;i<stop;++i,++outIndex){
    float dd=scalarVel*dt*(data[i-1]+22.0*data[i]+data[i+1])/24.0;
    out[outIndex]=out[outIndex-1]+dd;
  }

  //Last value is obtained with a non-centered operator
  out[outIndex]=out[outIndex-1];
  _lastOut=out[outIndex-1];

  return out;
}

float* velocityReceiver::differentiate(modelDefStruct* modelDef,float* data,float* out,
         int startI,int stopI){
  DEF_MODEL_LIMITS(modelDef);
  DEF_MODEL_SCALE(modelDef);

  //First value is obtained with a non-centered differentiator.
  int outIndex=0;
  if(!startI)
    out[outIndex++]=scalarVel/dt*(-23.0*data[0]+21.0*data[1]+
          3.0*data[2]-data[3])/24.0;

  //Interior values are obtained with a centered differentiator.
  int start=MAX(1,startI),stop=MIN(modelDef->NT-1,stopI);
  for(int i=start;i<stop-1;i++)
    out[outIndex++]=scalarVel/dt*(27.0*(data[i+1]-data[i  ])-
          (data[i+2]-data[i-1]))/24.0;

  //Last value is obtained with a non-centered differentiator.
  if(stopI==modelDef->NT)
    out[outIndex++]=scalarVel/dt*
(-23.0*data[stopI  ]+21.0*data[stopI-1]+
 3.0*  data[stopI-2]-     data[stopI-3])/24.0;
       
  return out;
}

//Constructor with parameters read from message
pressureReceiver::pressureReceiver(modelDefStruct* modelDef,int allocate)
  :receiver(modelDef){
  unpackMessage("i",&_useCubicInterp);

  if(setActive(modelDef)){
    if(allocate){
allocateData(modelDef); //this allocates, initializes and fills in any required interp
    }else{
_data=(float*)0x1;
    }
  }
}

//Constructor with parameters passed as arguments
pressureReceiver::pressureReceiver(modelDefStruct* modelDef,int allocate,
     float x,float y,float z,float amp,
     int useCubicInterp):
  receiver(modelDef,x,y,z,amp){
  _useCubicInterp=useCubicInterp;

  if(setActive(modelDef)){
    if(allocate){
allocateData(modelDef); //this allocates, initializes and fills in any required interp
    }else{
_data=(float*)0x1;
    }
  }
}

void pressureReceiver::packInit(){
  receiver::packInit();
  packMessage("i",_useCubicInterp);
}

//calculate the value of this receiver for the current iteration
float pressureReceiver::fill(modelDefStruct* modelDef,int iteration,
       float* vx,float* vy,float* vz,
       float* xx,float* yy,float* zz,
       float* xy,float* xz,float* yz){
  if(!active() ||
     iteration<0 || iteration>modelDef->NT) 
    return -10;
 
  if(!yy){
    //This has been called from the acoustic problem, xx is already pressure.
    if(_useCubicInterp){
      return _data[iteration]=_amp*
      triCubicInterp(modelDef,&_cubicInterp,xx);
    }
      return _data[iteration]=_amp*
                trilinInterp(modelDef,&_interp,xx);
  }

  //Interpolate diagonal components of stress tensor onto receiver location.
  float ampx, ampy, ampz;
  if(_useCubicInterp) {
    ampx = triCubicInterp(modelDef,&_cubicInterp,xx);
    ampy = triCubicInterp(modelDef,&_cubicInterp,yy);
    ampz = triCubicInterp(modelDef,&_cubicInterp,zz);
  } else {
    ampx=trilinInterp(modelDef,&_interp,xx);
    ampy=trilinInterp(modelDef,&_interp,yy);
    ampz=trilinInterp(modelDef,&_interp,zz);
  }

  float val=-1.0/3.0*_amp*(ampx+ampy+ampz);
  return _data[iteration]=val;
}

int pressureReceiver::scaleData(modelDefStruct* modelDef,int startI,int stopI,
    float* data,float* wrk){
  DEF_MODEL_SCALE(modelDef);
  for(int i=startI,j=0;i<stopI;j++,i++)
    wrk[j]=data[i]*scalarStress;
  return stopI-startI;
}

void pressureReceiver::allocateData(modelDefStruct* modelDef){
  receiver::allocateData(modelDef);

  DEF_MODEL_LIMITS(modelDef);
  if(_useCubicInterp){
    triCubicCoeff(_x,_y,_z,
                  minX,dx,minY,dy,minZ,dz,
                  &_cubicInterp);
  } else {
    trilinCoeff(_x,_y,_z,
    minX,dx,minY,dy,minZ,dz,
    &_interp);
  }
}
