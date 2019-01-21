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
 *  Defines several class functions of the base classes for controlling
 *  boundary conditions.  Mostly these classes hold information concerning
 *  the boundaries and pass the tasks to other routines.  Defines the
 *  static function arrays for calling the correct routines for
 *  extrapolating to the very edge of the model.  Also, some of these
 *  classes control message passing among domains and define
 *  special MPI types to efficiently handle the messages.
 *
 *  Defines the following class functions:
 *  sgfdBoundaries::updateVel
 *  sgfdBoundaries::updateStress
 *  neighborProcess::neighborProcess
 *  neighborProcess::neighborProcess
 *  neighborProcess::startPassVel
 *  neighborProcess::startPassVelXY
 *  neighborProcess::startPassVelZ
 *  neighborProcess::finishPassVel
 *  neighborProcess::finishPassVelZ
 *  neighborProcess::finishPassVelXY
 *  neighborProcess::startPassStress
 *  neighborProcess::finishPassStress
 *  neighborProcess::receiveNeighborVel2
 *  neighborProcess::sendNeighborVel2
 *  neighborProcess::receiveNeighborVel2Z
 *  neighborProcess::sendNeighborVel2Z
 *  neighborProcess::receiveNeighborStress2
 *  neighborProcess::sendNeighborStress2
 *  neighborProcess::sendNeighborVel
 *  neighborProcess::sendNeighborStress
 *  neighborProcess::receiveNeighborStress
 *  neighborProcess::receiveNeighborVel
 *  neighborProcess::buildMPITypes
 *  neighborProcess::setReceiveIndices
 *  neighborProcess::setSendIndices
 *  neighborProcess::setReceiveIndicesMD
 *  neighborProcess::setSendIndicesMD
 *
 *  Defines the following functions:
 *  triCubicCoeff
 *  interpolate
 *
 */

#include "sgfdBoundary.hh"
#include "sgfd.hh"

//the order is important; the boundary conditions must be done before the passing
int sgfdBoundaries::updateVel(int iteration,sgfdDependent* dep){
  //Start the pass
  DEF_PARALLEL(dep->model()->parallelDef());
  // Faces
  for(int i=-1;i<=1;i+=2)
    side(i,0,0)->startPassVel(dep,procI,i);
  for(int j=-1;j<=1;j+=2)
    side(0,j,0)->startPassVel(dep,procJ,j);
  for(int k=-1;k<=1;k+=2)
    side(0,0,k)->startPassVel(dep,procK,k);
  

  //Finish the pass
  // Faces
  for(int i=-1;i<=1;i+=2)
    side(i,0,0)->finishPassVel(dep,procI,i);
  for(int j=-1;j<=1;j+=2)
    side(0,j,0)->finishPassVel(dep,procJ,j);
  for(int k=-1;k<=1;k+=2)
    side(0,0,k)->finishPassVel(dep,procK,k);
  
  return 26;
}
int sgfdBoundaries::updateStress(int iteration,sgfdDependent* dep){
  //Start the pass
  DEF_PARALLEL(dep->model()->parallelDef());
  // Faces
  for(int i=-1;i<=1;i+=2)
    side(i,0,0)->startPassStress(dep,procI,i);
  for(int j=-1;j<=1;j+=2)
    side(0,j,0)->startPassStress(dep,procJ,j);
  for(int k=-1;k<=1;k+=2)
    side(0,0,k)->startPassStress(dep,procK,k);
  
  //Finish the pass
  // Faces
  for(int i=-1;i<=1;i+=2)
    side(i,0,0)->finishPassStress(dep,procI,i);
  for(int j=-1;j<=1;j+=2)
    side(0,j,0)->finishPassStress(dep,procJ,j);
  for(int k=-1;k<=1;k+=2)
    side(0,0,k)->finishPassStress(dep,procK,k);
  
  return 26;
}

//
//now the classes the non-virtual boundaryCondition classes
//

//neighborProcess: internal procedures
neighborProcess::neighborProcess(sgfdModel* model,
                int sideI,int sideJ,int sideK):
boundaryCondition(sideI,sideJ,sideK){
  DEF_PARALLEL(model->parallelDef());
  _tid=BOUNDARY(neighbors,_sideI,_sideJ,_sideK);
#if USE_MPI_SEND
  _splitXY_Z = false;
#endif  //USE_MPI_SEND
}

int neighborProcess::startPassVel(sgfdDependent* dependent,
                 int proc,int dir){
#if USE_MASS_PASS
#if USE_WARNING
#warning "Using mass passing; may not work if default send is blocking"
#endif
  sendNeighborVel2(dependent);
#else
#if USE_WARNING
#warning "Passing neighbors one at a time; slower if default send is nonblocking"
#endif
  if(ISEVEN(proc)){
    //First send from left to right, from top to bottom, and from front to back
    if(dir>0){
      sendNeighborVel2(dependent);
      receiveNeighborVel2(dependent);
    }
  }else{
    if(dir<0){
      receiveNeighborVel2(dependent);
      sendNeighborVel2(dependent);
    }
  }
#endif
  return 6;
}

int neighborProcess::startPassVelXY(sgfdDependent* dependent,
                   int proc,int dir){
#if USE_MASS_PASS
#if USE_WARNING
#warning "Using mass passing; may not work if default send is blocking"
#endif
  sendNeighborVel2(dependent);
#else
#if USE_WARNING
#warning "Passing neighbors one at a time; slower if default send is nonblocking"
#endif
  if(ISEVEN(proc)){
    //First send from left to right, from top to bottom, and from front to back
    if(dir>0){
      sendNeighborVel2(dependent);
      receiveNeighborVel2(dependent);
    }
  }else{
    if(dir<0){
      receiveNeighborVel2(dependent);
      sendNeighborVel2(dependent);
    }
  }
#endif
  return 6;
}

int neighborProcess::startPassVelZ(sgfdDependent* dependent,
                  int proc,int dir){
#if USE_MASS_PASS
#if USE_WARNING
#warning "Using mass passing; may not work if default send is blocking"
#endif
  sendNeighborVel2Z(dependent);
#else
#if USE_WARNING
#warning "Passing neighbors one at a time; slower if default send is nonblocking"
#endif
  if(ISEVEN(proc)){
    //First send from left to right, from top to bottom, and from front to back
    if(dir>0){
      sendNeighborVel2Z(dependent);
      receiveNeighborVel2Z(dependent);
    }
  }else{
    if(dir<0){
      receiveNeighborVel2Z(dependent);
      sendNeighborVel2Z(dependent);
    }
  }
#endif
  return 6;
}

int neighborProcess::finishPassVel(sgfdDependent* dependent,
                  int proc,int dir){
#if USE_MASS_PASS
  receiveNeighborVel2(dependent);
#else
  if(ISEVEN(proc)){
    //First send from left to right, from top to bottom, and from front to back
    if(dir<0){
      sendNeighborVel2(dependent);
      receiveNeighborVel2(dependent);
    }
  }else{
    if(dir>0){
      receiveNeighborVel2(dependent);
      sendNeighborVel2(dependent);
    }
  }
#endif
  return 6;
}

int neighborProcess::finishPassVelZ(sgfdDependent* dependent,
                   int proc,int dir){
#if USE_MASS_PASS
  receiveNeighborVel2Z(dependent);
#else
  if(ISEVEN(proc)){
    //First send from left to right, from top to bottom, and from front to back
    if(dir<0){
      sendNeighborVel2Z(dependent);
      receiveNeighborVel2Z(dependent);
    }
  }else{
    if(dir>0){
      receiveNeighborVel2Z(dependent);
      sendNeighborVel2Z(dependent);
    }
  }
#endif
  return 6;
}

int neighborProcess::finishPassVelXY(sgfdDependent* dependent,
                    int proc,int dir){
#if USE_MASS_PASS
  receiveNeighborVel2(dependent);
#else
  if(ISEVEN(proc)){
    //First send from left to right, from top to bottom, and from front to back
    if(dir<0){
      sendNeighborVel2(dependent);
      receiveNeighborVel2(dependent);
    }
  }else{
    if(dir>0){
      receiveNeighborVel2(dependent);
      sendNeighborVel2(dependent);
    }
  }
#endif
  return 6;
}

int neighborProcess::startPassStress(sgfdDependent* dependent,
                    int proc,int dir){
#if USE_MASS_PASS
  sendNeighborStress2(dependent);
#else
  if(ISEVEN(proc)){
    //First send from left to right, from top to bottom, and from front to back
    if(dir>0){
      sendNeighborStress2(dependent);
      receiveNeighborStress2(dependent);
    }
  }else{
    if(dir<0){
      receiveNeighborStress2(dependent);
      sendNeighborStress2(dependent);
    }
  }
#endif
  return 6;
}

int neighborProcess::finishPassStress(sgfdDependent* dependent,
                     int proc,int dir){
#if USE_MASS_PASS
  receiveNeighborStress2(dependent);
#else
  if(ISEVEN(proc)){
    //First send from left to right, from top to bottom, and from front to back
    if(dir<0){
      sendNeighborStress2(dependent);
      receiveNeighborStress2(dependent);
    }
  }else{
    if(dir>0){
      receiveNeighborStress2(dependent);
      sendNeighborStress2(dependent);
    }
  }
#endif
  return 6;
}

void neighborProcess::receiveNeighborVel2(sgfdDependent* dependent){
  receiveNeighborVel(dependent->modelDef(),
                     dependent->vx(),dependent->vy(),dependent->vz());
}
void neighborProcess::sendNeighborVel2(sgfdDependent* dependent){
  sendNeighborVel(dependent->modelDef(),
                  dependent->vx(),dependent->vy(),dependent->vz());
}

void neighborProcess::receiveNeighborVel2Z(sgfdDependent* dependent){
  receiveNeighborVelZ(dependent->modelDef(),
                      dependent->vz());
}
void neighborProcess::sendNeighborVel2Z(sgfdDependent* dependent){
  sendNeighborVelZ(dependent->modelDef(),
                   dependent->vz());
}

void neighborProcess::receiveNeighborStress2(sgfdDependent* dependent){
  receiveNeighborStress(dependent->modelDef(),
                        dependent->xx(),dependent->yy(),dependent->zz(),
                        dependent->xy(),dependent->xz(),dependent->yz());
}
void neighborProcess::sendNeighborStress2(sgfdDependent* dependent){
  sendNeighborStress(dependent->modelDef(),
                     dependent->xx(),dependent->yy(),dependent->zz(),
                     dependent->xy(),dependent->xz(),dependent->yz());
}

void neighborProcess::receiveNeighborVel(modelDefStruct* modelDef,
                                         float* vx,float* vy,float* vz,
                                         int messageTagInc){
#if USE_MPI_SEND
  getBlocks(_tid,&_vReceiveGroupType,_vReceiveAll);
#else
  DEF_MODEL_SIZE(modelDef);
  
  getMessage(_tid,MESSAGE_SLAVE_TO_SLAVE+messageTagInc,NULL);
#if PASS_NORMAL_VELOCITY_ONLY
#if USE_WARNING
#warning "Passing normal velocities only (Acoustic problem execution)"
#endif
  if(_sideI){
    unpackBlock(_receiveXStartM1,_receiveXStopM1,
                _receiveYStart,_receiveYStop,
                _receiveZStart,_receiveZStop,
                1,NX,NXY,&vx);
  }else if(_sideJ){
    unpackBlock(_receiveXStart,_receiveXStop,
                _receiveYStartM1,_receiveYStopM1,
                _receiveZStart,_receiveZStop,
                1,NX,NXY,&vy);
  }else{
    unpackBlock(_receiveXStart,_receiveXStop,
                _receiveYStart,_receiveYStop,
                _receiveZStartM1,_receiveZStopM1,
                1,NX,NXY,&vz);
  }
#else //#if PASS_NORMAL_VELOCITY_ONLY
  unpackBlock(_receiveXStartM1,_receiveXStopM1,
              _receiveYStart,_receiveYStop,
              _receiveZStart,_receiveZStop,
              1,NX,NXY,&vx);
  unpackBlock(_receiveXStart,_receiveXStop,
              _receiveYStartM1,_receiveYStopM1,
              _receiveZStart,_receiveZStop,
              1,NX,NXY,&vy);
  unpackBlock(_receiveXStart,_receiveXStop,
              _receiveYStart,_receiveYStop,
              _receiveZStartM1,_receiveZStopM1,
              1,NX,NXY,&vz);
#endif //#if PASS_NORMAL_VELOCITY_ONLY
#endif //#if USE_MPI_SEND
}

void neighborProcess::sendNeighborVel(modelDefStruct* modelDef,
                     float* vx,float* vy,float* vz,
                     int messageTagInc){
#if USE_MPI_SEND
#if USE_IMMEDIATE_MPI_SEND
  immediateSendBlocks(&_sendRequest,_tid,
                      &_vSendGroupType,_vSendAll);
#else
  sendBlocks(_tid,&_vSendGroupType,_vSendAll);
#endif //#if USE_IMMEDIATE_MPI_SEND
#else //#if USE_MPI_SEND
  DEF_MODEL_SIZE(modelDef);
  initSend();
#if PASS_NORMAL_VELOCITY_ONLY
  if(_sideI){
    packBlock(_sendXStartM1,_sendXStopM1,
              _sendYStart,_sendYStop,
              _sendZStart,_sendZStop,
              1,NX,NXY,&vx);
  }else if(_sideJ){
    packBlock(_sendXStart,_sendXStop,
              _sendYStartM1,_sendYStopM1,
              _sendZStart,_sendZStop,
              1,NX,NXY,&vy);
  }else{
    packBlock(_sendXStart,_sendXStop,
              _sendYStart,_sendYStop,
              _sendZStartM1,_sendZStopM1,
              1,NX,NXY,&vz);
  }
#else //#if PASS_NORMAL_VELOCITY_ONLY
  packBlock(_sendXStartM1,_sendXStopM1,
            _sendYStart,_sendYStop,
            _sendZStart,_sendZStop,
            1,NX,NXY,&vx);
  packBlock(_sendXStart,_sendXStop,
            _sendYStartM1,_sendYStopM1,
            _sendZStart,_sendZStop,
            1,NX,NXY,&vy);
  packBlock(_sendXStart,_sendXStop,
            _sendYStart,_sendYStop,
            _sendZStartM1,_sendZStopM1,
            1,NX,NXY,&vz);
#endif //#if PASS_NORMAL_VELOCITY_ONLY
#if USE_IMMEDIATE_MPI_SEND
  iSendMessage(_tid,MESSAGE_SLAVE_TO_SLAVE+messageTagInc,NULL);
#else //#if USE_IMMEDIATE_MPI_SEND
  sendMessage(_tid,MESSAGE_SLAVE_TO_SLAVE+messageTagInc,NULL);
#endif //#if USE_IMMEDIATE_MPI_SEND
#endif //#if USE_MPI_SEND
}

void neighborProcess::sendNeighborStress(modelDefStruct* modelDef,
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
  DEF_MODEL_SIZE(modelDef);
  initSend();
  //blocks need to be individually packed to work with the mismatch
  // neighbor class as it is currently working
  packBlock(_sendXStart,_sendXStop,
            _sendYStart,_sendYStop,
            _sendZStart,_sendZStop,
            1,NX,NXY,&xx);
  if(yy && zz){
    packBlock(_sendXStart,_sendXStop,
              _sendYStart,_sendYStop,
              _sendZStart,_sendZStop,
              1,NX,NXY,&yy);
    packBlock(_sendXStart,_sendXStop,
              _sendYStart,_sendYStop,
              _sendZStart,_sendZStop,
              1,NX,NXY,&zz);
    packBlock(_sendXStartM1,_sendXStopM1,
              _sendYStartM1,_sendYStopM1,
              _sendZStart,_sendZStop,
              1,NX,NXY,&xy);
    packBlock(_sendXStartM1,_sendXStopM1,
              _sendYStart,_sendYStop,
              _sendZStartM1,_sendZStopM1,
              1,NX,NXY,&xz);
    packBlock(_sendXStart,_sendXStop,
              _sendYStartM1,_sendYStopM1,
              _sendZStartM1,_sendZStopM1,
              1,NX,NXY,&yz);
  }
  
#if USE_IMMEDIATE_MPI_SEND
  iSendMessage(_tid,MESSAGE_SLAVE_TO_SLAVE+messageTagInc,NULL);
#else //#if USE_IMMEDIATE_MPI_SEND
  sendMessage(_tid,MESSAGE_SLAVE_TO_SLAVE+messageTagInc,NULL);
#endif //#if USE_IMMEDIATE_MPI_SEND
#endif //#if USE_MPI_SEND
}

void neighborProcess::receiveNeighborStress(modelDefStruct* modelDef,
                           float* xx,float* yy,float* zz,
                           float* xy,float* xz,float* yz,
                           int messageTagInc){
#if USE_MPI_SEND
#if USE_IMMEDIATE_MPI_SEND
  getBlocks(_tid,&_sReceiveGroupType,_sReceiveAll);
#else
  getBlocks(_tid,&_sReceiveGroupType,_sReceiveAll);
#endif //#if USE_IMMEDIATE_MPI_SEND
#else
  DEF_MODEL_SIZE(modelDef);
  
  getMessage(_tid,MESSAGE_SLAVE_TO_SLAVE+messageTagInc,NULL);
  //blocks need to be individually packed to work with the mismatch
  // neighbor class as it is currently working
  unpackBlock(_receiveXStart,_receiveXStop,
              _receiveYStart,_receiveYStop,
              _receiveZStart,_receiveZStop,
              1,NX,NXY,&xx);
  if(yy && zz){
    unpackBlock(_receiveXStart,_receiveXStop,
                _receiveYStart,_receiveYStop,
                _receiveZStart,_receiveZStop,
                1,NX,NXY,&yy);
    unpackBlock(_receiveXStart,_receiveXStop,
                _receiveYStart,_receiveYStop,
                _receiveZStart,_receiveZStop,
                1,NX,NXY,&zz);
    unpackBlock(_receiveXStartM1,_receiveXStopM1,
                _receiveYStartM1,_receiveYStopM1,
                _receiveZStart,_receiveZStop,
                1,NX,NXY,&xy);
    unpackBlock(_receiveXStartM1,_receiveXStopM1,
                _receiveYStart,_receiveYStop,
                _receiveZStartM1,_receiveZStopM1,
                1,NX,NXY,&xz);
    unpackBlock(_receiveXStart,_receiveXStop,
                _receiveYStartM1,_receiveYStopM1,
                _receiveZStartM1,_receiveZStopM1,
                1,NX,NXY,&yz);
  }
#endif
}

#if USE_MPI_SEND
void neighborProcess::buildMPITypes(sgfdDependent* dependent){
  DEF_MODEL_SIZE(dependent->modelDef());
  
  //Build locals.
  float *vxV=dependent->vx(),*vyV=dependent->vy(),*vzV=dependent->vz();
  float *xxV=dependent->xx(),*yyV=dependent->yy(),*zzV=dependent->zz();
  float *xyV=dependent->xy(),*xzV=dependent->xz(),*yzV=dependent->yz();
  
  //Receive velocity blocks.
  setBlockType(&_vReceiveBlockTypes[0],&_vReceiveOffset[0],
               NX,NY,NZ,
               _receiveXStartM1,_receiveXStopM1,
               _receiveYStart,_receiveYStop,
               _receiveZStart,_receiveZStop);
  _vReceiveAll[0]=&vxV[_vReceiveOffset[0]];
  setBlockType(&_vReceiveBlockTypes[1],&_vReceiveOffset[1],
               NX,NY,NZ,
               _receiveXStart,_receiveXStop,
               _receiveYStartM1,_receiveYStopM1,
               _receiveZStart,_receiveZStop);
  _vReceiveAll[1]=&vyV[_vReceiveOffset[1]];
  setBlockType(&_vReceiveBlockTypes[2],&_vReceiveOffset[2],
               NX,NY,NZ,
               _receiveXStart,_receiveXStop,
               _receiveYStart,_receiveYStop,
               _receiveZStartM1,_receiveZStopM1);
  _vReceiveAll[2]=&vzV[_vReceiveOffset[2]];
  
  //Send velocity blocks.
  setBlockType(&_vSendBlockTypes[0],&_vSendOffset[0],
               NX,NY,NZ,
               _sendXStartM1,_sendXStopM1,
               _sendYStart,_sendYStop,
               _sendZStart,_sendZStop);
  _vSendAll[0]=&vxV[_vSendOffset[0]];
  setBlockType(&_vSendBlockTypes[1],&_vSendOffset[1],
               NX,NY,NZ,
               _sendXStart,_sendXStop,
               _sendYStartM1,_sendYStopM1,
               _sendZStart,_sendZStop);
  _vSendAll[1]=&vyV[_vSendOffset[1]];
  setBlockType(&_vSendBlockTypes[2],&_vSendOffset[2],
               NX,NY,NZ,
               _sendXStart,_sendXStop,
               _sendYStart,_sendYStop,
               _sendZStartM1,_sendZStopM1);
  _vSendAll[2]=&vzV[_vSendOffset[2]];
  
  //Build groups of the blocks for send and receiving velocity.
#if PASS_NORMAL_VELOCITY_ONLY
#if USE_WARNING
#warning "Passing normal velocities only (Acoustic problem setup)"
#endif
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
#else
  if(_splitXY_Z) {
    setBlocksType(&_vReceiveGroupType,2,
                  _vReceiveAll,_vReceiveBlockTypes);
    setBlocksType(&_vSendGroupType,2,
                  _vSendAll,_vSendBlockTypes);
    setBlocksType(&_vReceiveGroupTypeZ,1,
                  _vReceiveAll+2,_vReceiveBlockTypes+2);
    setBlocksType(&_vSendGroupTypeZ,1,
                  _vSendAll+2,_vSendBlockTypes+2);
  } else {
    setBlocksType(&_vReceiveGroupType,3,
                  _vReceiveAll,_vReceiveBlockTypes);
    setBlocksType(&_vSendGroupType,3,
                  _vSendAll,_vSendBlockTypes);
  }
#endif //#if PASS_NORMAL_VELOCITY_ONLY
  
  //Receive stress blocks.
  setBlockType(&_sReceiveBlockTypes[0],&_sReceiveOffset[0],
               NX,NY,NZ,
               _receiveXStart,_receiveXStop,
               _receiveYStart,_receiveYStop,
               _receiveZStart,_receiveZStop);
  _sReceiveAll[0]=&xxV[_sReceiveOffset[0]];
  if(yyV && zzV && xyV && xzV && yzV){
    setBlockType(&_sReceiveBlockTypes[1],&_sReceiveOffset[1],
                 NX,NY,NZ,
                 _receiveXStart,_receiveXStop,
                 _receiveYStart,_receiveYStop,
                 _receiveZStart,_receiveZStop);
    _sReceiveAll[1]=&yyV[_sReceiveOffset[1]];
    setBlockType(&_sReceiveBlockTypes[2],&_sReceiveOffset[2],
                 NX,NY,NZ,
                 _receiveXStart,_receiveXStop,
                 _receiveYStart,_receiveYStop,
                 _receiveZStart,_receiveZStop);
    _sReceiveAll[2]=&zzV[_sReceiveOffset[2]];
    
    setBlockType(&_sReceiveBlockTypes[3],&_sReceiveOffset[3],
                 NX,NY,NZ,
                 _receiveXStartM1,_receiveXStopM1,
                 _receiveYStartM1,_receiveYStopM1,
                 _receiveZStart,_receiveZStop);
    _sReceiveAll[3]=&xyV[_sReceiveOffset[3]];
    setBlockType(&_sReceiveBlockTypes[4],&_sReceiveOffset[4],
                 NX,NY,NZ,
                 _receiveXStartM1,_receiveXStopM1,
                 _receiveYStart,_receiveYStop,
                 _receiveZStartM1,_receiveZStopM1);
    _sReceiveAll[4]=&xzV[_sReceiveOffset[4]];
    setBlockType(&_sReceiveBlockTypes[5],&_sReceiveOffset[5],
                 NX,NY,NZ,
                 _receiveXStart,_receiveXStop,
                 _receiveYStartM1,_receiveYStopM1,
                 _receiveZStartM1,_receiveZStopM1);
    _sReceiveAll[5]=&yzV[_sReceiveOffset[5]];
  }
  
  //Send stress blocks.
  setBlockType(&_sSendBlockTypes[0],&_sSendOffset[0],
               NX,NY,NZ,
               _sendXStart,_sendXStop,
               _sendYStart,_sendYStop,
               _sendZStart,_sendZStop);
  _sSendAll[0]=&xxV[_sSendOffset[0]];
  if(yyV && zzV && xyV && xzV && yzV){
    setBlockType(&_sSendBlockTypes[1],&_sSendOffset[1],
                 NX,NY,NZ,
                 _sendXStart,_sendXStop,
                 _sendYStart,_sendYStop,
                 _sendZStart,_sendZStop);
    _sSendAll[1]=&yyV[_sSendOffset[1]];
    setBlockType(&_sSendBlockTypes[2],&_sSendOffset[2],
                 NX,NY,NZ,
                 _sendXStart,_sendXStop,
                 _sendYStart,_sendYStop,
                 _sendZStart,_sendZStop);
    _sSendAll[2]=&zzV[_sSendOffset[2]];
    
    setBlockType(&_sSendBlockTypes[3],&_sSendOffset[3],
                 NX,NY,NZ,
                 _sendXStartM1,_sendXStopM1,
                 _sendYStartM1,_sendYStopM1,
                 _sendZStart,_sendZStop);
    _sSendAll[3]=&xyV[_sSendOffset[3]];
    setBlockType(&_sSendBlockTypes[4],&_sSendOffset[4],
                 NX,NY,NZ,
                 _sendXStartM1,_sendXStopM1,
                 _sendYStart,_sendYStop,
                 _sendZStartM1,_sendZStopM1);
    _sSendAll[4]=&xzV[_sSendOffset[4]];
    setBlockType(&_sSendBlockTypes[5],&_sSendOffset[5],
                 NX,NY,NZ,
                 _sendXStart,_sendXStop,
                 _sendYStartM1,_sendYStopM1,
                 _sendZStartM1,_sendZStopM1);
    _sSendAll[5]=&yzV[_sSendOffset[5]];
  }
  //Build groups of the blocks for send and receiving stress.
  setBlocksType(&_sReceiveGroupType,yyV&&zzV?6:1,
                _sReceiveAll,_sReceiveBlockTypes);
  setBlocksType(&_sSendGroupType,yyV&&zzV?6:1,
                _sSendAll,_sSendBlockTypes);
}
#endif //#if USE_MPI_SEND

void neighborProcess::setReceiveIndicesMD(modelDefStruct* modelDef,int widthF){
  DEF_MODEL_SIZE(modelDef);
  if(!_sideI){
    //receive the whole distance
    _receiveXStartM1=_receiveXStart=0;
    _receiveXStopM1=_receiveXStop=NX;
  }else if(_sideI==-1){
    //This is actually a simplification, receive the first two  and uncentered centered nodes.
    _receiveXStart=_receiveXStartM1=0;
    _receiveXStop=_receiveXStopM1=widthF;
  }else{ //sideI==1
    //Receive the last two centered and uncentered nodes.
    _receiveXStart=_receiveXStartM1=NX-widthF;
    _receiveXStop=_receiveXStopM1=NX;
  }

  if(!_sideJ){
    _receiveYStartM1=_receiveYStart=0;
    _receiveYStopM1=_receiveYStop=NY;
  }else if(_sideJ==-1){
    _receiveYStart=_receiveYStartM1=0;
    _receiveYStop=_receiveYStopM1=widthF;
  }else{ //sideJ==1
    _receiveYStart=_receiveYStartM1=NY-widthF;
    _receiveYStop=_receiveYStopM1=NY;
  }

  if(!_sideK){
    _receiveZStartM1=_receiveZStart=0;
    _receiveZStopM1=_receiveZStop=NZ;
  }else if(_sideK==-1){
    _receiveZStart=_receiveZStartM1=0;
    _receiveZStop=_receiveZStopM1=widthF;
  }else{ //sideK==1
    _receiveZStart=_receiveZStartM1=NZ-widthF;
    _receiveZStop=_receiveZStopM1=NZ;
  }
}
void neighborProcess::setSendIndicesMD(modelDefStruct* modelDef,int widthF){
  DEF_MODEL_SIZE(modelDef);
  if(!_sideI){
    //send the whole distance
    _sendXStartM1=_sendXStart=0;
    _sendXStopM1=_sendXStop=NX;
  }else if(_sideI==-1){
    //send the first+2 two centered and non-centerd nodes
    _sendXStart=_sendXStartM1=widthF;
    _sendXStop=_sendXStopM1=2*widthF;
  }else{ //_sideI==1
    //send the last two centered and non-centerednodes
    _sendXStart=_sendXStartM1=NX-2*widthF;
    _sendXStop=_sendXStopM1=NX-widthF;
  }
  if(!_sideJ){
    _sendYStartM1=_sendYStart=0;
    _sendYStopM1=_sendYStop=NY;
  }else if(_sideJ==-1){
    _sendYStart=_sendYStartM1=widthF;
    _sendYStop=_sendYStopM1=2*widthF;
  }else{ //_sideJ==1
    _sendYStart=_sendYStartM1=NY-2*widthF;
    _sendYStop=_sendYStopM1=NY-widthF;
  }
  if(!_sideK){
    _sendZStartM1=_sendZStart=0;
    _sendZStopM1=_sendZStop=NZ;
  }else if(_sideK==-1){
    _sendZStart=_sendZStartM1=widthF;
    _sendZStop=_sendZStopM1=2*widthF;
  }else{ //_sideK==1
    _sendZStart=_sendZStartM1=NZ-2*widthF;
    _sendZStop=_sendZStopM1=NZ-widthF;
  }
}
//boundary class: internal functions

//For efficency purposes we may want to have models running with
// different grid spacing. The mismatchNeighbor class is used if
// the neighbor to this process is runnning on a different grid
inline void triCubicCoeff(int i,int j,int k,
			  float p,float q,float r,
			  int iptr[64],int jptr[64],int kptr[64],float coeff[64]){
  //Compute and store tricubic polynomial interpolator coefficients
  float a[4]={-p*(p-1.0)*(p-2.0)/6.0,
	      (p-1.0)*(p+1.0)*(p-2.0)/2.0,
	      -p*(p+1.0)*(p-2.0)/2.0,
	      p*(p-1.0)*(p+1.0)/6.0};

  float b[4]={-q*(q-1.0)*(q-2.0)/6.0,
	      (q-1.0)*(q+1.0)*(q-2.0)/2.0,
	      -q*(q+1.0)*(q-2.0)/2.0,
	      q*(q-1.0)*(q+1.0)/6.0};

  float c[4]={-r*(r-1.0)*(r-2.0)/6.0,
	      (r-1.0)*(r+1.0)*(r-2.0)/2.0,
	      -r*(r+1.0)*(r-2.0)/2.0,
	      r*(r-1.0)*(r+1.0)/6.0};

  //fill in
  for(int kk=0,kkk=k-1;kk<4;kk++,kkk++){
    for(int jj=0,jjj=j-1;jj<4;jj++,jjj++){
      for(int ii=0,iii=i-1;ii<4;ii++,iii++){
	int index=ii+4*jj+16*kk;
	iptr[index]=iii;
	jptr[index]=jjj;
	kptr[index]=kkk;

	coeff[index]=a[ii]*b[jj]*c[kk];
      }
    }
  }
} 

#define ISZERO(x)(fabs(x)<1e-5)
#define ISONE(x) (fabs(x-1.0)<1e-5)
inline float interpolate(float x,float y,float z,
			 float minX,float dx,float minY,float dy,float minZ,float dz,
			 float *vx,int NX,int NXY,int NY,int NZ){
  //Compute spatial indeces of cell that surrounds receiver.
  int i=(int)floor((x-minX)/dx),j=(int)floor((y-minY)/dy),k=(int)floor((z-minZ)/dz);

  //Compute coordinates of corner of cell.
  float xi=minX+(float)i*dx,yj=minY+(float)j*dy,zk=minZ+(float)k*dz;

  //Compute and store trilinear interpolator coefficients.
  float p=(x-xi)/dx,q=(y-yj)/dy,r=(z-zk)/dz;

  //check for interpolation short circuit
  if(ABS(p)+ABS(q)+ABS(r) < 0.01){
    return VX(i,j,k);
  }

  //do a one step interpolation
  float p1=1.0-p,q1=1.0-q,r1=1.0-r;
  float value=0.0,totalCoeff=0.0;
  for(int kk=0,kkk=k;kk<2;kkk++,kk++){
    for(int jj=0,jjj=j;jj<2;jjj++,jj++){
      for(int ii=0,iii=i;ii<2;iii++,ii++){
	if(ISMID(0,iii,NX-1)&&ISMID(0,jjj,NY-1)&&ISMID(0,kkk,NZ-1)){
	  float coeff=(ii?p:p1)*(jj?q:q1)*(kk?r:r1);
	  value+=coeff*VX(iii,jjj,kkk);
	  totalCoeff+=coeff;
	}
      }
    }
  }

  return totalCoeff>0.0?value/totalCoeff:0.0;
}
