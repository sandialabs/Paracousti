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
 *  acousticBoundary.cc
 *
 *
 *  Define acoustic boundary related classes.  These are mainly
 *  classes that hold information about the boundaries instead of actually
 *  doing the boundary updates.  The neighbor process classes do the 
 *  message passing among domains.
 *
 *  Defines the following class functions:
 *    acousticNeighborProcess::acousticNeighborProcess
 *    acousticNeighborProcess::receiveNeighborVel
 *    acousticNeighborProcess::sendNeighborVel
 *    acousticNeighborProcess::sendNeighborStress
 *    acousticNeighborProcess::receiveNeighborStress
 *    acousticNeighborProcess::buildMPITypes
 *    acousticSgfdBoundaries::acousticSgfdBoundaries
 *
 */

#include "acousticBoundary.hh"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "message_passing.h"
#include "sgfd.hh"

//Setup message passing classes for specifics for acoustic problem
acousticNeighborProcess::acousticNeighborProcess(sgfdModel* model,sgfdDependent* dependent,
                        int sideI,int sideJ,int sideK):
neighborProcess(model,sideI,sideJ,sideK){
  
  //Set the indices to send and receive
  //Use the new functions which provide more flexability and minimize the velocity layers
  // that are passed (this was the cause of the asymetry in the message passing I observed
  // when profiling previous versions of pe).
  setReceiveIndicesMD(model->modelDef(),2);
  setSendIndicesMD(model->modelDef(),2);
  
#if USE_MPI_SEND
  buildMPITypes(dependent);
#endif
}
//use specific acoustic MPI types for the message passing
void acousticNeighborProcess::receiveNeighborVel(modelDefStruct* modelDef,
                        float* vx,float* vy,float* vz,
                        int messageTagInc){
#if USE_MPI_SEND
  getBlocks(_tid,
            &_vReceiveGroupType,_vReceiveAll);
#else
  neighborProcess::receiveNeighborVel(modelDef,vx,vy,vz);
#endif //#if USE_MPI_SEND
}
//use specific acoustic MPI types for the message passing
void acousticNeighborProcess::sendNeighborVel(modelDefStruct* modelDef,
                     float* vx,float* vy,float* vz,
                     int messageTagInc){
#if USE_MPI_SEND
#if USE_IMMEDIATE_MPI_SEND
  immediateSendBlocks(&_sendRequest,_tid,
                      &_vSendGroupType,_vSendAll);
#else //#if USE_IMMEDIATE_MPI_SEND
  sendBlocks(_tid,
             &_vSendGroupType,_vSendAll);
#endif //#if USE_IMMEDIATE_MPI_SEND
#else //#if USE_MPI_SEND
  neighborProcess::sendNeighborVel(modelDef,vx,vy,vz);
#endif //#if USE_MPI_SEND
}
//use specific acoustic MPI types for the message passing
void acousticNeighborProcess::sendNeighborStress(modelDefStruct* modelDef,
                        float* xx,float* yy,float* zz,
                        float* xy,float* xz,float* yz,
                        int messageTagInc){
#if USE_MPI_SEND
#if USE_IMMEDIATE_MPI_SEND
  immediateSendBlocks(&_sendRequest,
                      _tid,&_sSendGroupType,_sSendAll);
#else //#if USE_IMMEDIATE_MPI_SEND
  sendBlocks(_tid,&_sSendGroupType,_sSendAll);
#endif //#if USE_IMMEDIATE_MPI_SEND
#else //#if USE_MPI_SEND
  neighborProcess::sendNeighborStress(modelDef,xx,yy,zz,xy,xz,yz);
#endif //#if USE_MPI_SEND
}
//use specific acoustic MPI types for the message passing
void acousticNeighborProcess::receiveNeighborStress(modelDefStruct* modelDef,
                           float* xx,float* yy,float* zz,
                           float* xy,float* xz,float* yz,
                           int messageTagInc){
#if USE_MPI_SEND
  getBlocks(_tid,&_sReceiveGroupType,_sReceiveAll);
#else
  neighborProcess::receiveNeighborStress(modelDef,xx,yy,zz,xy,xz,yz);
#endif
}

//here is where all the specific MPI types are built for message passing
//in the acoustic problem
#if USE_MPI_SEND
void acousticNeighborProcess::buildMPITypes(sgfdDependent* dependent){
  DEF_MODEL_SIZE(dependent->modelDef());
  
  //Build locals.
  float *vx=dependent->vx(),*vy=dependent->vy(),*vz=dependent->vz();
  float *p=dependent->P();
  
  //Receive velocity blocks.
  setBlockType(&_vReceiveBlockTypes[0],&_vReceiveOffset[0],
               NX,NY,NZ,
               _receiveXStartM1,_receiveXStopM1,
               _receiveYStart,_receiveYStop,
               _receiveZStart,_receiveZStop);
  _vReceiveAll[0]=&vx[_vReceiveOffset[0]];
  setBlockType(&_vReceiveBlockTypes[1],&_vReceiveOffset[1],
               NX,NY,NZ,
               _receiveXStart,_receiveXStop,
               _receiveYStartM1,_receiveYStopM1,
               _receiveZStart,_receiveZStop);
  _vReceiveAll[1]=&vy[_vReceiveOffset[1]];
  setBlockType(&_vReceiveBlockTypes[2],&_vReceiveOffset[2],
               NX,NY,NZ,
               _receiveXStart,_receiveXStop,
               _receiveYStart,_receiveYStop,
               _receiveZStartM1,_receiveZStopM1);
  _vReceiveAll[2]=&vz[_vReceiveOffset[2]];
  
  //Send velocity blocks.
  setBlockType(&_vSendBlockTypes[0],&_vSendOffset[0],
               NX,NY,NZ,
               _sendXStartM1,_sendXStopM1,
               _sendYStart,_sendYStop,
               _sendZStart,_sendZStop);
  _vSendAll[0]=&vx[_vSendOffset[0]];
  setBlockType(&_vSendBlockTypes[1],&_vSendOffset[1],
               NX,NY,NZ,
               _sendXStart,_sendXStop,
               _sendYStartM1,_sendYStopM1,
               _sendZStart,_sendZStop);
  _vSendAll[1]=&vy[_vSendOffset[1]];
  setBlockType(&_vSendBlockTypes[2],&_vSendOffset[2],
               NX,NY,NZ,
               _sendXStart,_sendXStop,
               _sendYStart,_sendYStop,
               _sendZStartM1,_sendZStopM1);
  _vSendAll[2]=&vz[_vSendOffset[2]];
  
  //Build groups of the blocks for send and receiving velocity.
  if(_sideI){
    setBlocksType(&_vReceiveGroupType,1,
                  _vReceiveAll,&_vReceiveBlockTypes[0]);
    setBlocksType(&_vSendGroupType,1,
                  _vSendAll,&_vSendBlockTypes[0]);
  }else if(_sideJ){
    _vReceiveAll[0]=_vReceiveAll[1];
    setBlocksType(&_vReceiveGroupType,1,
                  _vReceiveAll,&_vReceiveBlockTypes[1]);
    _vSendAll[0]=_vSendAll[1];
    setBlocksType(&_vSendGroupType,1,
                  _vSendAll,&_vSendBlockTypes[1]);
  }else{
    _vReceiveAll[0]=_vReceiveAll[2];
    setBlocksType(&_vReceiveGroupType,1,
                  _vReceiveAll,&_vReceiveBlockTypes[2]);
    _vSendAll[0]=_vSendAll[2];
    setBlocksType(&_vSendGroupType,1,
                  _vSendAll,&_vSendBlockTypes[2]);
  }
  
  //Receive stress blocks.
  setBlockType(&_sReceiveBlockTypes[0],&_sReceiveOffset[0],
               NX,NY,NZ,
               _receiveXStart,_receiveXStop,
               _receiveYStart,_receiveYStop,
               _receiveZStart,_receiveZStop);
  _sReceiveAll[0]=&p[_sReceiveOffset[0]];
  
  //Send stress blocks.
  setBlockType(&_sSendBlockTypes[0],&_sSendOffset[0],
               NX,NY,NZ,
               _sendXStart,_sendXStop,
               _sendYStart,_sendYStop,
               _sendZStart,_sendZStop);
  _sSendAll[0]=&p[_sSendOffset[0]];
  
  //Build groups of the blocks for send and receiving stress.
  setBlocksType(&_sReceiveGroupType,1,
                _sReceiveAll,_sReceiveBlockTypes);
  setBlocksType(&_sSendGroupType,1,
                _sSendAll,_sSendBlockTypes);
}
#endif //#if USE_MPI_SEND

//This class defines what types of boundaries there are on each side
//of this domain.  We can pass information to the other domains (neighbor
//processes), have pressure-free surfaces, use CPML or sponge boundaries.
//It does little itself.  Mainly it defines these edge types and holds
//that information
acousticSgfdBoundaries::acousticSgfdBoundaries(sgfdModel* model,sgfdDependent* dependent,
                       int freeSurfaceMode,int usePML):sgfdBoundaries(){
  _freeSurfaceMode=freeSurfaceMode;
  DEF_PARALLEL(model->parallelDef());
  
  _usePML = usePML;
  //The following variables are used in CPML and Sponge boundaries to
  //determine if we are in a domain that is in the CPML or sponge zone
  //and if we are right on the edge of the model (i.e., no domain adjacent
  //on this side (the f* variables))
  minXb = maxXb = minYb = maxYb = minZb = maxZb = false;
  fminXb = fmaxXb = fminYb = fmaxYb = fminZb = fmaxZb = false;
  
  int numBoundarys=0,numNeighbors=0,numNull=0;
  for(int i=-1;i<=1;i++){
    for(int j=-1;j<=1;j++){
      for(int k=-1;k<=1;k++){
        if(!i && !j && !k){
          side(i,j,k)=NULL; //this is myself
        }else if(IS_FACE(i,j,k) && NEIGHBOR(i,j,k)>=0){
          //information needs to be passed from one domain to the next
          //but only for direct neighbors, not diagonal ones
          numNeighbors++;
          side(i,j,k)=
          new acousticNeighborProcess(model,dependent,
                                      i,j,k);
        }else if(IS_FACE(i,j,k)){
          //face boundarys are always done
          numBoundarys++;
          //A pressure-free surface
          if(k==-1 && _freeSurfaceMode==SURFACE_PRESS_FREE)
            side(i,j,k) = new boundary(_freeSurfaceMode,i,j,k);
          else {  //any other boundary
            side(i,j,k)=new boundaryCondition(i,j,k);  //does nothing
            //The following values however are used by CPML and sponge boundaries.
            if(i==-1) minXb = fminXb = true;
            else if(i==1) maxXb = fmaxXb = true;
            else if(j==-1) minYb = fminYb = true;
            else if(j==1) maxYb = fmaxYb = true;
            else if(k==-1) minZb = fminZb = true;
            else maxZb = fmaxZb = true;
          }
        }else{
          //edges and corners are only done if they are at the extreme
          // edge of the model. The master boundaryCondition class does
          // nothing at all.
          numNull++;
          side(i,j,k)=new boundaryCondition(i,j,k);
        }
      }
    }
  }
}
