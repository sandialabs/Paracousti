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
 *  sgfdReceiverGroups.hh
 *
 *
 *  Declares classes used in groups of receivers, mainly
 *  of class receiverGrid.  These handle common types of receivers
 *  that vary in location, but otherwise are the same.
 *
 *  Declares the following classes:
 *  receiverGrid
 *  pressureReceiverGrid
 *  vxReceiverGrid
 *  vyReceiverGrid
 *  vzReceiverGrid
 *
 */
//These utility functions are for receiver groups in the SGFD problem.
//IMPORTANT NOTE: this is not a standalone header, it must be included within
// sgfdReceivers.hh after the receiver parent class is defined and before the
// receiverGroup container is defined.

#ifndef _sgfdReceiverGroups_hh_
#define _sgfdReceiverGroups_hh_

#include "sgfdReceivers.hh"

#define PRESSURE_RECEIVER_GRID 11
#define VX_RECEIVER_GRID 12
#define VY_RECEIVER_GRID 13
#define VZ_RECEIVER_GRID 14

/// Here is another type of receiver, this is optimised for a dense grid. Start
///  with a virtual class for all receiverGrids.
class receiverGrid:public receiver{
public:
  int _gridIndex,_gridTimeSubsample,_nt,_rawOutput;
  char _rawOutputDir[512];
  interpStruct _interp; //to interpolate from the surronding grid nodes to the reciever location
  int _n,_nx,_ny,_nz,_dx,_dy,_dz;
  int _ix,_iy,_iz;

  float *_xcorr;

  ///Here is the slave constructor, get information from a message that has
  /// already been received.
  receiverGrid(modelDefStruct* modelDef,int gridIndex);

  ///Here is the master version, get the info from arguments that are passed in.
  receiverGrid(modelDefStruct* modelDef,int allocate,float amp,
	       float x,float rdx,int nx,
	       float y,float rdy,int ny,
	       float z,float rdz,int nz,
	       int timeSubsample=1,
	       int rawOutput=FALSE,char* rawOutputDir=NULL);

  ///Method to determine if output is raw.
  int isRaw(){return _rawOutput;}

  ///Pack a message with the parameters for this receiver. Assumes the buffer is 
  /// ready to be filled.
  virtual void packInit();

  ///Read the data in preperation for performing RTM cross-correlation.
  virtual int readData(modelDefStruct* modelDef);

  ///Pack the data, since this one does multiple sends do our own buffer init.
  virtual int packData(modelDefStruct* modelDef,float* wrk,
                       int startI,int stopI,int decimation);

  virtual void allocateData(modelDefStruct* modelDef){
    _n=_nx*_ny*_nz*_nt;
    if(!_n){
      _data=NULL;
    }else{
      assert((_data=(float*)malloc(_n*sizeof(float)))!=NULL,
	     "receiverGrid--unable to allocate %i*%i*%i*%i floats for data",
	     _nx,_ny,_nz,_nt);
      for(int i=0;i<_n;_data[i++]=0.0);
    }
  }
  ///Define the fill method as an iterator, then the real subclasses just need
  /// to define a method to actually fill the data.
  virtual float fill(modelDefStruct* modelDef,int iteration,
		     float* vx,float* vy,float* vz,
		     float* xx,float* yy,float* zz,
                     float* xy,float* xz,float* yz);
  
protected:
  virtual float fillvalue(modelDefStruct* modelDef,int i,int j,int k,
			  float* vx,float* vy,float* vz,
			  float* xx,float* yy,float* zz,
			  float* xy,float* xz,float* yz)=0;
};
typedef receiverGrid* receiverGridPtr;
typedef arrayI<receiverGridPtr> receiverGridArray;
typedef receiverGridArray* receiverGridArrayPtr;

///Here is a real receiver grid type that records pressure.
class pressureReceiverGrid:public receiverGrid{
public:
  interpStruct _interp; //To interpolate from the surronding grid nodes to the 
  //                       reciever location.

  ///Constructor with parameters read from message
  pressureReceiverGrid(modelDefStruct* modelDef,int gridIndex,int allocate):
    receiverGrid(modelDef,gridIndex){	
    if(allocate){
      allocateData(modelDef); //this allocates, initializes and fills in any required interp
      if(_rawOutput==2)
	readData(modelDef);
    }
  }

  //Constructor with parameters passed as arguments
  pressureReceiverGrid(modelDefStruct* modelDef,int allocate,float amp,
		       float x,float rdx,int nx,
		       float y,float rdy,int ny,
		       float z,float rdz,int nz,
		       int timeSubsample=1,
		       int rawOutput=FALSE,char* rawOutputDir=NULL):
    receiverGrid(modelDef,allocate,amp,x,rdx,nx,y,rdy,ny,z,rdz,nz,
		 timeSubsample,rawOutput,rawOutputDir){
    if(setActive(modelDef)){
      if(allocate){
	allocateData(modelDef); //this allocates, initializes and fills in any required interp
      }else{
	_data=(float*)0x1;
      }
    }
  }

  virtual int type(){return PRESSURE_RECEIVER_GRID;}
  virtual const char* typeName(){return "PressureGrid";}

  virtual int scaleData(modelDefStruct* modelDef,int startI,int stopI,
			float* data,float* wrk){
    DEF_MODEL_SCALE(modelDef);
    for(int i=0;i<_nt*_nx*_ny*_nz;_data[i++]*=scalarStress);
    return _nt*_nx*_ny*_nz;
  }

  virtual void allocateData(modelDefStruct* modelDef){
    DEF_MODEL_LIMITS(modelDef);
    trilinCoeff(_x,_y,_z,
		minX,dx,minY,dy,minZ,dz,
		&_interp);

    return receiverGrid::allocateData(modelDef);
  }
protected:
  virtual float fillvalue(modelDefStruct* modelDef,int i,int j,int k,
			  float* vx,float* vy,float* vz,
			  float* xx,float* yy,float* zz,
                          float* xy,float* xz,float* yz);
};

///Here is a real receiver grid type that records vx.
class vxReceiverGrid:public receiverGrid{
public:
  interpStruct _interp; //To interpolate from the surronding grid nodes to the 
  //                       reciever location.

  ///Constructor with parameters read from message
  vxReceiverGrid(modelDefStruct* modelDef,int gridIndex,int allocate):
    receiverGrid(modelDef,gridIndex){	
    if(allocate){
      allocateData(modelDef); //this allocates, initializes and fills in any required interp
      if(_rawOutput==2)
	readData(modelDef);
    }
  }

  //Constructor with parameters passed as arguments
  vxReceiverGrid(modelDefStruct* modelDef,int allocate,float amp,
		       float x,float rdx,int nx,
		       float y,float rdy,int ny,
		       float z,float rdz,int nz,
		       int timeSubsample=1,
		       int rawOutput=FALSE,char* rawOutputDir=NULL):
    receiverGrid(modelDef,allocate,amp,x,rdx,nx,y,rdy,ny,z,rdz,nz,
		 timeSubsample,rawOutput,rawOutputDir){
    if(setActive(modelDef)){
      if(allocate){
	allocateData(modelDef); //this allocates, initializes and fills in any required interp
      }else{
	_data=(float*)0x1;
      }
    }
  }

  virtual int type(){return VX_RECEIVER_GRID;}
  virtual const char* typeName(){return "VxGrid";}

  virtual int scaleData(modelDefStruct* modelDef,int startI,int stopI,
			float* data,float* wrk){
    DEF_MODEL_SCALE(modelDef);
    for(int i=0;i<_nt*_nx*_ny*_nz;_data[i++]*=scalarVel);
    return _nt*_nx*_ny*_nz;
  }

  virtual void allocateData(modelDefStruct* modelDef){
    DEF_MODEL_LIMITS(modelDef);
    trilinCoeff(_x,_y,_z,
		minX+0.5*dx,dx,minY,dy,minZ,dz,
		&_interp);

    return receiverGrid::allocateData(modelDef);
  }
protected:
  virtual float fillvalue(modelDefStruct* modelDef,int i,int j,int k,
			  float* vx,float* vy,float* vz,
			  float* xx,float* yy,float* zz,
			  float* xy,float* xz,float* yz){
    return _amp*trilinInterpOffset(modelDef,&_interp,i*_dx,j*_dy,k*_dz,vx);
  }
};

///Here is a real receiver grid type that records vy.
class vyReceiverGrid:public receiverGrid{
public:
  interpStruct _interp; //To interpolate from the surronding grid nodes to the 
  //                       reciever location.

  ///Constructor with parameters read from message
  vyReceiverGrid(modelDefStruct* modelDef,int gridIndex,int allocate):
    receiverGrid(modelDef,gridIndex){	
    if(allocate){
      allocateData(modelDef); //this allocates, initializes and fills in any required interp
      if(_rawOutput==2)
	readData(modelDef);
    }
  }

  //Constructor with parameters passed as arguments
  vyReceiverGrid(modelDefStruct* modelDef,int allocate,float amp,
		       float x,float rdx,int nx,
		       float y,float rdy,int ny,
		       float z,float rdz,int nz,
		       int timeSubsample=1,
		       int rawOutput=FALSE,char* rawOutputDir=NULL):
    receiverGrid(modelDef,allocate,amp,x,rdx,nx,y,rdy,ny,z,rdz,nz,
		 timeSubsample,rawOutput,rawOutputDir){
    if(setActive(modelDef)){
      if(allocate){
	allocateData(modelDef); //this allocates, initializes and fills in any required interp
      }else{
	_data=(float*)0x1;
      }
    }
  }

  virtual int type(){return VY_RECEIVER_GRID;}
  virtual const char* typeName(){return "VyGrid";}

  virtual int scaleData(modelDefStruct* modelDef,int startI,int stopI,
			float* data,float* wrk){
    DEF_MODEL_SCALE(modelDef);
    for(int i=0;i<_nt*_nx*_ny*_nz;_data[i++]*=scalarVel);
    return _nt*_nx*_ny*_nz;
  }

  virtual void allocateData(modelDefStruct* modelDef){
    DEF_MODEL_LIMITS(modelDef);
    trilinCoeff(_x,_y,_z,
		minX,dx,minY+0.5*dy,dy,minZ,dz,
		&_interp);

    return receiverGrid::allocateData(modelDef);
  }
protected:
  virtual float fillvalue(modelDefStruct* modelDef,int i,int j,int k,
			  float* vx,float* vy,float* vz,
			  float* xx,float* yy,float* zz,
			  float* xy,float* xz,float* yz){
    return _amp*trilinInterpOffset(modelDef,&_interp,i*_dx,j*_dy,k*_dz,vy);
  }
};

///Here is a real receiver grid type that records vz.
class vzReceiverGrid:public receiverGrid{
public:
  interpStruct _interp; //To interpolate from the surronding grid nodes to the 
  //                       reciever location.

  ///Constructor with parameters read from message
  vzReceiverGrid(modelDefStruct* modelDef,int gridIndex,int allocate):
    receiverGrid(modelDef,gridIndex){	
    if(allocate){
      allocateData(modelDef); //this allocates, initializes and fills in any required interp
      if(_rawOutput==2)
	readData(modelDef);
    }
  }

  //Constructor with parameters passed as arguments
  vzReceiverGrid(modelDefStruct* modelDef,int allocate,float amp,
		       float x,float rdx,int nx,
		       float y,float rdy,int ny,
		       float z,float rdz,int nz,
		       int timeSubsample=1,
		       int rawOutput=FALSE,char* rawOutputDir=NULL):
    receiverGrid(modelDef,allocate,amp,x,rdx,nx,y,rdy,ny,z,rdz,nz,
		 timeSubsample,rawOutput,rawOutputDir){
    if(setActive(modelDef)){
      if(allocate){
	allocateData(modelDef); //this allocates, initializes and fills in any required interp
      }else{
	_data=(float*)0x1;
      }
    }
  }

  virtual int type(){return VZ_RECEIVER_GRID;}
  virtual const char* typeName(){return "VzGrid";}

  virtual int scaleData(modelDefStruct* modelDef,int startI,int stopI,
			float* data,float* wrk){
    DEF_MODEL_SCALE(modelDef);
    for(int i=0;i<_nt*_nx*_ny*_nz;_data[i++]*=scalarVel);
    return _nt*_nx*_ny*_nz;
  }

  virtual void allocateData(modelDefStruct* modelDef){
    DEF_MODEL_LIMITS(modelDef);
    trilinCoeff(_x,_y,_z,
		minX,dx,minY,dy,minZ+0.5*dz,dz,
		&_interp);

    return receiverGrid::allocateData(modelDef);
  }
protected:
  virtual float fillvalue(modelDefStruct* modelDef,int i,int j,int k,
			  float* vx,float* vy,float* vz,
			  float* xx,float* yy,float* zz,
			  float* xy,float* xz,float* yz){
    return _amp*trilinInterpOffset(modelDef,&_interp,i*_dx,j*_dy,k*_dz,vz);
  }
};

#endif //#ifndef _sgfdReceiverGroups_hh_
