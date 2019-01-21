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
 *  sgfdBoundary.cc
 *
 *
 *  Declares several base classes for controlling
 *  boundary conditions.  Mostly these classes hold information concerning
 *  the boundaries and pass the tasks to other routines.  Defines the
 *  static function arrays for calling the correct routines for
 *  extrapolating to the very edge of the model.  Also, some of these
 *  classes control message passing among domains and define
 *  special MPI types to efficiently handle the messages.
 *
 *  Declares the following classes:
 *  boundaryCondition
 *  neighborProcess
 *  boundary
 *  sgfdBoundaries
 *
 */

#ifndef _sgfdBoundary_hh_
#define _sgfdBoundary_hh_

class sgfdDependent;
class sgfdModel;

#include "nstdutil.hh"
#include "sgfd.h"

#include "message_passing.h"

#define IS_FACE(i,j,k)((ABS(i)+ABS(j)+ABS(k))==1)
#define IS_EDGE(i,j,k)((ABS(i)+ABS(j)+ABS(k))==2)
#define IS_CORNER(i,j,k)((ABS(i)+ABS(j)+ABS(k))==3)

//Allow for multiple types of surface boundary condtions. So far I
// am interested in: my normal free surface, a pressure free surface for the 
// acoustic problem, and a velocity free surface for the acoustic problem.
#define SURFACE_NO_CONDITION 0
#define SURFACE_STRESS_FREE  1
#define SURFACE_PRESS_FREE   2
#define SURFACE_VEL_FREE     3

///The top level boundaryCondition class
/// It knows how to apply itself to velocity, stress, and pressure updates.
/// This is the base class it does nothing. This is used for edges and corners
/// which should not do anything
class boundaryCondition{
public:
  int _sideI,_sideJ,_sideK; //these are my index
  ///Inialize, needs modelDef, parallelDef and what side am I.
  boundaryCondition(int sideI,int sideJ,int sideK){
    _sideI=sideI; _sideJ=sideJ; _sideK=sideK;
  }

  //do the deletion in two stages
  virtual ~boundaryCondition(){}

  //boolean functions
  ///Is this a real boundary or is a subdomain to trace messages with.
  virtual int isRealBoundary(){return FALSE;}

  ///What is the ID of my neighbor.
  virtual int isNeighbor(){return FALSE;}

  ///Is there a special surface condition
  virtual int surfaceBC(){return SURFACE_NO_CONDITION;}

  //it is easier to break up into neighbor passing and boundary condition implementation
  // subroutines. In this base class none of the subroutines do anything at all
  //Neighbor passing; there are two of each to allow the possiblity of bulk passing

  ///Start passing velocity info to my neighbor.
  virtual int startPassVel(sgfdDependent* dep,
			   int dir,int proc){return FALSE;}
  virtual int startPassVelXY(sgfdDependent* dep,
                           int dir,int proc){return FALSE;}
  virtual int startPassVelZ(sgfdDependent* dep,
                           int dir,int proc){return FALSE;}
  ///Finish passing velocity info to my neighbor.
  virtual int finishPassVel(sgfdDependent* dep,
			    int dir,int proc){return FALSE;}
  virtual int finishPassVelZ(sgfdDependent* dep,
                            int dir,int proc){return FALSE;}
  virtual int finishPassVelXY(sgfdDependent* dep,
                            int dir,int proc){return FALSE;}

  ///Start passing stress info to my neighbor.
  virtual int startPassStress(sgfdDependent* dep,
			      int dir,int proc){return FALSE;}
  ///Finish passing stress info to my neighbor.
  virtual int finishPassStress(sgfdDependent* dep,
			       int dir,int proc){return FALSE;}

#if USE_MPI_SEND & USE_IMMEDIATE_MPI_SEND
  virtual void waitForImmediateSend(){return;}
#endif 
};
typedef boundaryCondition* boundaryConditionPtr;
    
/// The neighbor class; these are ghost cells that must be shared.

///Finally figured out the indexing scheme to pass the minimum number of velocity components,
/// the trick is to skip one.
///  x-5  x-5  x-4  x-4  x-3  x-3  x-2  x-2  x-1  x-1
///   s    v    s    v    s    v    s    v    s    v
///   C    C    C    C    C    C    P    P    P
///             P    P    P    C    C    C    C    C    C    C
///             s    v    s    v    s    v    s    v    s    v
///             0    0    1    1    2    2    3    3    4    4
class neighborProcess: public boundaryCondition{
protected:
  int _tid; //What process am I doing the swap with.

  //These are the indices of the nodes to be sent to the neighbor
  int _sendXStart,_sendXStop,_sendXStartM1,_sendXStopM1;
  int _sendYStart,_sendYStop,_sendYStartM1,_sendYStopM1;
  int _sendZStart,_sendZStop,_sendZStartM1,_sendZStopM1;
  //and the indices to be received from the neighbor
  int _receiveXStart,_receiveXStop,_receiveXStartM1,_receiveXStopM1;
  int _receiveYStart,_receiveYStop,_receiveYStartM1,_receiveYStopM1;
  int _receiveZStart,_receiveZStop,_receiveZStartM1,_receiveZStopM1;

  //These are special data types only used if the messages are being passed with MPI derived
  // types (more complicated but faster).
#if USE_MPI_SEND
  MPI_Datatype _vReceiveBlockTypes[3],_sReceiveBlockTypes[6];
  MPI_Datatype _vSendBlockTypes[3],_sSendBlockTypes[6];

  MPI_Datatype _vReceiveGroupType,_vSendGroupType,
    _vReceiveGroupTypeZ,_vSendGroupTypeZ,
    _sReceiveGroupType,_sSendGroupType;

  float *_vSendAll[3],*_sSendAll[6];
  float *_vReceiveAll[3],*_sReceiveAll[6];
  bool _splitXY_Z;

  int _vSendOffset[3],_vReceiveOffset[3];
  int _sSendOffset[6],_sReceiveOffset[6];
#if USE_IMMEDIATE_MPI_SEND
  MPI_Request _sendRequest,_getRequest;
#endif //#if USE_IMMEDIATE_MPI_SEND
#endif //#if USE_MPI_SEND

public:
  ///Dumb constructor for derived classes.
  neighborProcess(sgfdModel* model,
		  int sideI,int sideJ,int sideK);

  ///Redefine is neighbor to TRUE.
  virtual int isNeighbor(){return TRUE;}

  ///Reimplement to actually pass.
  virtual int startPassVel(sgfdDependent* dependent,
                           int proc,int dir);

  //pass only vx and vy
  virtual int startPassVelXY(sgfdDependent* dependent,
                             int proc,int dir);
  
  virtual int startPassVelZ(sgfdDependent* dependent,
                            int proc,int dir);

  ///Reimplement to actually pass.
  virtual int finishPassVel(sgfdDependent* dependent,
                            int proc,int dir);

  virtual int finishPassVelZ(sgfdDependent* dependent,
                             int proc,int dir);

  virtual int finishPassVelXY(sgfdDependent* dependent,
                              int proc,int dir);

  ///Reimplement to actually pass.
  virtual int startPassStress(sgfdDependent* dependent,
                              int proc,int dir);

  ///Reimplement to actually pass.
  virtual int finishPassStress(sgfdDependent* dependent,
                               int proc,int dir);

#if USE_MPI_SEND & USE_IMMEDIATE_MPI_SEND
  virtual void waitForImmediateSend(){
    MPI_Status status;
    MPI_Wait(&_sendRequest,&status);
//     MPI_Wait(&_getRequest,&status);
  }
#endif 
protected:
  //Internal functions
  virtual void receiveNeighborVel2(sgfdDependent* dependent);
  virtual void sendNeighborVel2(sgfdDependent* dependent);

  virtual void receiveNeighborVel2Z(sgfdDependent* dependent);
  virtual void sendNeighborVel2Z(sgfdDependent* dependent);

  virtual void receiveNeighborStress2(sgfdDependent* dependent);
  virtual void sendNeighborStress2(sgfdDependent* dependent);

  //Internal functions
  virtual void receiveNeighborVel(modelDefStruct* modelDef,
				  float* vx,float* vy,float* vz,
                                  int messageTagInc=0);

  virtual void receiveNeighborVelZ(modelDefStruct* modelDef,
                                  float* vz,
                                  int messageTagInc=0){
#if USE_MPI_SEND
    getBlocks(_tid,&_vReceiveGroupTypeZ,_vReceiveAll+2);
#endif //#if USE_MPI_SEND
  }

  virtual void sendNeighborVelZ(modelDefStruct* modelDef,
                               float* vz,
                               int messageTagInc=0){
#if USE_MPI_SEND
#if USE_IMMEDIATE_MPI_SEND
    immediateSendBlocks(&_sendRequest,_tid,
                        &_vSendGroupTypeZ,_vSendAll+2);
#else
    sendBlocks(_tid,&_vSendGroupTypeZ,_vSendAll+2);
#endif //#if USE_IMMEDIATE_MPI_SEND
#endif //USE_MPI_SEND
  }

  virtual void sendNeighborVel(modelDefStruct* modelDef,
			       float* vx,float* vy,float* vz,
                               int messageTagInc=0);

  virtual void sendNeighborStress(modelDefStruct* modelDef,
				  float* xx,float* yy,float* zz,
				  float* xy,float* xz,float* yz,
                                  int messageTagInc=0);

  virtual void receiveNeighborStress(modelDefStruct* modelDef,
				     float* xx,float* yy,float* zz,
				     float* xy,float* xz,float* yz,
                                     int messageTagInc=0);

  //And for the acoustic problem with a 2,2 passing stencil.
  void setReceiveIndicesMD(modelDefStruct* modelDef,int widthF);
  void setSendIndicesMD(modelDefStruct* modelDef,int widthF);

#if USE_MPI_SEND
  virtual void buildMPITypes(sgfdDependent* dependent);
#endif //#if USE_MPI_SEND
};

///This boundaryCondition is a real model edge boundary.
/// Since this is now being used by both
/// the elastic and the acoustic codes this is another abstract class for the
/// real classes elasticBoundary and acousticBoundary.
class boundary: public boundaryCondition{
protected:
  int _surfaceBCMode;

public:
  ///This initializer needs some extra info.
  boundary(int surfaceBCMode,
	   int sideI,int sideJ,int sideK):
    boundaryCondition(sideI,sideJ,sideK){
    _surfaceBCMode=surfaceBCMode;
  }

  //boolean test functions
  virtual int surfaceBC(){return _surfaceBCMode;};
  virtual int isRealBoundary(){return TRUE;}

protected:
};

/// Boundary conditions will be defined with a 3x3x3 grid.
/// efine a new class for the grid. This is the class that will
/// actually be used by the main processes
class sgfdBoundaries{
protected:
  //there are 26 possible neighbors for the rectangle
  // 6 faces, 12 edges, 8 corners; define a full box where
  // the center is a dummy for the live process
  boundaryCondition* _boundarys[27];
public:
  ///Initializer; use an array with 27 values; non-zero values are
  /// neighbors. Zero values are the edges of the grid
  sgfdBoundaries(){
    for(int i=0;i<27;_boundarys[i++]=NULL);
  }

  ///Access the elements
  boundaryConditionPtr& side(int i,int j,int k){return _boundarys[(i+1)+3*(j+1)+9*(k+1)];}
  ///Access an array of i start indicies for all dependent variables.
  virtual int* fullIstart(){
    assert(FALSE,"sgfdBoundaries::fullIstart function not defined for this class");
    return NULL;
  }
  ///Access an array of i stop indicies for all dependent variables.
  virtual int* fullIstop(){
    assert(FALSE,"sgfdBoundaries::fullIstop function not defined for this class");
    return NULL;
  }
  ///Access an array of j start indicies for all dependent variables.
  virtual int* fullJstart(){
    assert(FALSE,"sgfdBoundaries::fullJstart function not defined for this class");
    return NULL;
  }
  ///Access an array of j stop indicies for all dependent variables.
  virtual int* fullJstop(){
    assert(FALSE,"sgfdBoundaries::fullJstop function not defined for this class");
    return NULL;
  }
  ///Access an array of k start indicies for all dependent variables.
  virtual int* fullKstart(){
    assert(FALSE,"sgfdBoundaries::fullKstart function not defined for this class");
    return 0;
  }
  ///Access an array of k stop indicies for all dependent variables.
  virtual int* fullKstop(){
    assert(FALSE,"sgfdBoundaries::fullKstop function not defined for this class");
    return 0;
  }

  virtual int kstart(modelDefStruct* modelDef){
    assert(FALSE,"sgfdBoundaries::kstart function not defined for this class");
    return 0;
  }
  virtual int kstart_wz(modelDefStruct* modelDef){
    assert(FALSE,"sgfdBoundaries::kstart_wz function not defined for this class");
    return 0;
  }
  virtual int kstop(modelDefStruct* modelDef){
    assert(FALSE,"sgfdBoundaries::kstop function not defined for this class");
    return 0;
  }
  virtual int kstop_wz(modelDefStruct* modelDef){
    assert(FALSE,"sgfdBoundaries::kstop_wz function not defined for this class");
    return 0;
  }

  //Access info about the elements.
  ///Is this run using PML boundaries (width is arbitrary).
  virtual int pmlW(){return 0;}

  ///Is this side a real boundary?
  int isRealBoundary(int i,int j,int k){
    return _boundarys[(i+1)+3*(j+1)+9*(k+1)]->isRealBoundary();
  }

  ///Does this side use a special boundary condition.
  int surfaceBC(int i,int j,int k){
    return _boundarys[(i+1)+3*(j+1)+9*(k+1)]->surfaceBC();
  }

  ///Apply updateVel to all members.
  virtual int updateVel(int iteration,sgfdDependent* dep);

  ///Apply updateStress to all members.
  virtual int updateStress(int iteration,sgfdDependent* dep);

#if USE_MPI_SEND & USE_IMMEDIATE_MPI_SEND
  void checkPassComplete(){
    for(int i=-1;i<=1;i+=2)
      side(i,0,0)->waitForImmediateSend();
    for(int j=-1;j<=1;j+=2)
      side(0,j,0)->waitForImmediateSend();
    for(int k=-1;k<=1;k+=2)
      side(0,0,k)->waitForImmediateSend();

  }
#endif
};

#endif //_sgfdBoundary_hh_
