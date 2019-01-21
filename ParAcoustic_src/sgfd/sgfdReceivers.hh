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
 *  sgfdReceivers.hh
 *
 *
 *  Declares classes used to represent individual receivers.  These handle their own
 *  interpolation, integration, etc.
 *
 *  Declares the following classes:
 *  receiver
 *  velocityReceiver
 *  pressureReceiver
 *
 *  Defines the following macros:
 *  MAX_R_PER_SEND
 *  VELOCITY_RECEIVER
 *  PRESSURE_RECEIVER
 *
 */
//These utility functions are for receivers in the SGFD problem

#ifndef _sgfdReceivers_hh_
#define _sgfdReceivers_hh_

#include <stdio.h>
#include "array.hh"
#include "sgfd.h"

#define MAX_R_PER_SEND 500

#define VELOCITY_RECEIVER 1
#define PRESSURE_RECEIVER 2

extern const int AllProcesses;

class receiver{
public:
  float _x,_y,_z; //position variables
  float _bx,_by,_bz; //direction variables
  float _amp;   //amplitude for this receiver

  float* _data; //completed information

  //basic constructor with location read from message
  receiver(modelDefStruct* modelDef);

  //basic constructor with location passed as argument
  receiver(modelDefStruct* modelDef,float x,float y,float z,float amp){
    _x=x;
    _y=y;
    _z=z;

    _amp=amp;

    _bx=_by=_bz=0.0; //may not be used but must be present in all forms of receiver
    // since they are written into the cdf file
    _data=NULL;
  }

  virtual ~receiver(){
    if(_data && _data!=(float*)0x1) free(_data);
  }

  ///Method for checkpoint writing of receiver data.
  virtual void doCheckpoint(int iteration,FILE* cpFile);

  ///Method for checkpoint reading of receiver data.
  virtual void readCheckpoint(int iteration,FILE* cpFile);
  virtual int setActive(modelDefStruct* modelDef);

  int active(){return _data != NULL;}
  float& data(int index){return _data[index];}
  virtual int integrate(){return FALSE;}

  virtual int type()=0;
  virtual const char* typeName()=0;

  //Pack a message with the parameters for this receiver. Assumes the buffer is ready to 
  // be filled.
  virtual void packInit();


  //Calculate the value of this receiver for the current iteration, there are spaces here
  // for elastic type codes with 3 velocities and 6 stresses and for acoustic type codes
  // with 3 velocities and pressure.
  virtual float fill(modelDefStruct* modelDef,int iteration,
		     float* vx,float* vy,float* vz,
		     float* xx,float* yy,float* zz,
		     float* xy,float* xz,float* yz)=0;

  //this does whatever work is required
  // to put the output into real units. This includes interpolation for the velocity types
  // where the data is on a half integer time raster.
  virtual int packData(modelDefStruct* modelDef,float* wrk,
                       int startI,int stopI,int decimation);
  virtual int scaleData(modelDefStruct* modelDef,int startI,int stopI,
			float* data,float* wrk)=0;

  //Allocate and intialize the data, also fill in any required interpStructs
  virtual void allocateData(modelDefStruct* modelDef);
  
};
typedef receiver* receiverPtr;
typedef arrayI<receiverPtr> receiverArray;
typedef receiverArray* receiverArrayPtr;

//Now define the various types of receivers that can actually be put into a run.
///velocityReceiver--record a single component of velocity.
class velocityReceiver:public receiver{
public:
  int _integrate; //may want displacement or acceleration (1 or -1 respectivily)
  float _lastOut; //Need to keep for -Rt option.
  int _useCubicInterp, _surfaceReceiver;
  interpStruct *_xInterp,*_yInterp,*_zInterp; //to interpolate from the surronding grid nodes to the
                                           // reciever location
  cubicInterpStruct *_xCubicInterp,*_yCubicInterp,*_zCubicInterp;
  //Constructor with parameters read from message
  velocityReceiver(modelDefStruct* modelDef,int allocate);

  //Constructor with parameters as arguments
  velocityReceiver(modelDefStruct* modelDef,int allocate,
		   float x,float y,float z,
		   float bx,float by,float bz,float amp,int useCubic=TRUE,int surfaceReceiver=FALSE);
  
  virtual ~velocityReceiver() {
    if(_useCubicInterp) {
      delete _xCubicInterp;
      delete _yCubicInterp;
      delete _zCubicInterp;
    } else {
      delete _xInterp;
      delete _yInterp;
      delete _zInterp;
    }
  }

  //Pack a message with the parameters for this receiver. Assumes the buffer is ready to 
  // be filled.
  virtual void packInit();

  virtual int type(){return VELOCITY_RECEIVER;}
  virtual const char* typeName(){return "Velocity";}
  
  int setDisplacement(){return _integrate=1;}
  int setAcceleration(){return _integrate=-1;}
  virtual int integrate(){return _integrate;}

  //calculate the value of this receiver for the current iteration
  virtual float fill(modelDefStruct* modelDef,int iteration,
		     float* vx,float* vy,float* vz,
		     float* xx,float* yy,float* zz,
                     float* xy,float* xz,float* yz);

  virtual int scaleData(modelDefStruct* modelDef,int startI,int stopI,
                        float* data,float* wrk);
  //Fill in any required interpStructs
  virtual void allocateData(modelDefStruct* modelDef);
private:
  //Velocities are recorded on the 1/2 integer time raster; move to the integer
  // raster by interpolation
  float* interpolate(modelDefStruct* modelDef,float* data,float* out,
                     int startI,int stopI);
    
  //Particle displacement values are calculated by integrating
  // the particle velocity samples over time.  A local parabolic
  // polynomial is used to approximate the velocity functon.
  // Calculated displacement samples are defined on the integer 
  // time raster t=tmin+(l-1)*dt for l=1,2,...,nt.
  float* integrate(modelDefStruct* modelDef,float* data,float* out,
                   int startI,int stopI);

  float* differentiate(modelDefStruct* modelDef,float* data,float* out,
                       int startI,int stopI);
};

class pressureReceiver:public receiver{
public:
  interpStruct _interp; //to interpolate from the surronding grid nodes to the reciever location
  
  //Testing the cubic instead of the linear interpolators.
  int _useCubicInterp;
  cubicInterpStruct _cubicInterp;

  //Constructor with parameters read from message
  pressureReceiver(modelDefStruct* modelDef,int allocate);

  //Constructor with parameters passed as arguments
  pressureReceiver(modelDefStruct* modelDef,int allocate,
		   float x,float y,float z,float amp,
		   int useCubicInterp=TRUE);

  virtual void packInit();
  virtual int type(){return PRESSURE_RECEIVER;}
  virtual const char* typeName(){return "Pressure";}

  //calculate the value of this receiver for the current iteration
  virtual float fill(modelDefStruct* modelDef,int iteration,
		     float* vx,float* vy,float* vz,
		     float* xx,float* yy,float* zz,
                     float* xy,float* xz,float* yz);

  virtual int scaleData(modelDefStruct* modelDef,int startI,int stopI,
                        float* data,float* wrk);

  virtual void allocateData(modelDefStruct* modelDef);
protected:
};

#endif //#define _sgfdReceivers_hh_
