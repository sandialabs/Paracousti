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
 *  fixed_acoustic.cc
 *
 *
 *  This file contains class function definitions for the master
 *  fixed acoustic model.  It also defines a output group type.
 *
 *  Defines the following class functions:
 *  acousticOutputGroup::addExtraOutput
 *  masterFixedMediaModel::masterFixedMediaModel
 *  masterFixedMediaModel::getInitAck
 *  masterFixedMediaModel::getAck
 *
 */

#include "fixed_acoustic.hh"
#include "message_passing.h"
#include "selector.hh"

/*! Output group for the acoustic problem

  This is were the messages and initialization methods are defined that 
  know about the new types of slices and acoustic volume output.
  Note:  This has not been tested recently.
*/
int acousticOutputGroup::addExtraOutput(int& i,int argOffset,
         int argc,char* argv[],
         modelDefStruct* modelDef,
         int runCall){

  //Check for arguments dealing with full wavefield ouput, they must be processed
  // here since the generic extraOuput class does not know what type of variables
  // need to be written.
  if(strlen(argv[i])>=argOffset+1 && argv[i][argOffset]=='w'){
    assert(argc>i+1,"addExtraOutput--adding wavefield output requires time");
    assert(runCall,"addExtraOutput--%s flag only valid at runtime",argv[i]);

    float t=atof(argv[++i]);

    _theOutput->Add(new wavefieldOutput(t,1.0,
          4,"vx","vy","vz","pressure"));
  }else{
    return
      extraOutputGroup::addExtraOutput(i,argOffset,argc,argv,modelDef,runCall);
  }
  return i;
}

//Constructor for the master model class for the fixed media problem.
masterFixedMediaModel::masterFixedMediaModel(int tids[],int mediaType,
                       int surfaceBCType,int gridMultiplier,
                       modelDefStruct* modelDef,parallelDefStruct* parallelDef,
                       char* modelName,int fdCoeffOrder,float* fdCoeffs,
                       int longRead)
:masterSgfdModel(modelDef,parallelDef,modelName,
                 surfaceBCType,gridMultiplier){
  initSend();
  sendMessage(AllProcesses,MESSAGE_INITIALIZE,"ii",
              mediaType,FALSE); //Need to send a pivot type (FALSE) just so the
  //                                        message is the correct length.
  initializeSlaves(tids,
                   modelName,fdCoeffOrder,fdCoeffs,longRead);
  
  getAck(FIXED_MEDIA,tids);
}
  ///This virtual function is called at the end of initializeSlaves to receive an
  /// acknoledgement that the slaves have completed initialization. Redefine here to give some
  /// additional information on the acoustic model.
void masterFixedMediaModel::getInitAck(int tids[]){
  float rhoMin=0.0,rhoMax=0.0,bmMin=0.0,bmMax=0.0;
  for(int i=0;i<NumProcs;i++){
    float localVMin,localVMax,localRMin,localRMax;
    float localBmMin,localBmMax;
    
    getMessage(AllProcesses,MESSAGE_INITIALIZE,"ff ff ff",
   &localVMin,&localVMax,&localRMin,&localRMax,
   &localBmMin,&localBmMax);

    SETMIN(_vMin,localVMin,i);
    SETMAX(_vMax,localVMax,i);

    SETMIN(rhoMin,localRMin,i);
    SETMAX(rhoMax,localRMax,i);

    SETMIN(bmMin,localBmMin,i);
    SETMAX(bmMax,localBmMax,i);
  }
  //and return the global min and max to the slaves as well as printing a diagnostic
  initSend();
  sendMessage(AllProcesses,MESSAGE_INITIALIZE,"ff",_vMin,_vMax);
  tEprintf(Verbose,
     "Model read by slaves; global min vel %.1f; max vel %.1f\n",
     _vMin,_vMax);
  if(rhoMin>0.1){
    tEprintf(Verbose,
       "  Rho %.2f-%.2f; Bulk Modulus %.1f-%.1f; \n",
       rhoMin,rhoMax,
       bmMin,bmMax);
  }else{
    tEprintf(Verbose,
       "  Rho %.3g to %.3g; Bulk Modulus %.3g to %.3g; \n",
       rhoMin,rhoMax,
       bmMin,bmMax);
  }
}

  ///Here is a revised (and non-virtual) acknoledgement method which is used
  /// by the masterSgfdModel constructors.
void masterFixedMediaModel::getAck(int tag,int* tids){
    //Receive and ack so we know we have successfully initialized
  for(int i=0;i<NumProcs;i++){
    int mediaType;
    float localVMin=0.0,localVMax=0.0;
    float localMinVx=0.0,localMinVy=0.0,localMinVz=0.0;
    float localMaxVx=0.0,localMaxVy=0.0,localMaxVz=0.0;
    float localMinMachNum=0.0,localMaxMachNum=0.0;
    getMessage(tids[i],MESSAGE_INITIALIZE,"i ff ff ff ff ff",
   &mediaType,
   &localVMin,&localVMax,
   &localMinMachNum,&localMaxMachNum,
   &localMinVx,&localMaxVx,
   &localMinVy,&localMaxVy,
   &localMinVz,&localMaxVz);
    assert(mediaType==tag,
     "Proc %i, media type set to %i not %i",
     tids[i],mediaType,tag);

    SETMIN(_vMin,localVMin,i);
    SETMAX(_vMax,localVMax,i);
  }

  if(_gridMultiplier!=1){
    DEF_MODEL_SIZE(modelDef());
    DEF_MODEL_LIMITS(modelDef());

    tEprintf(Verbose,"Modified model grid by factor of %i\n",_gridMultiplier);
    tEprintf(Verbose,
       "\tX: start %.1f; dx %.1f; nx %i=>stop %.1f\n",
       minX,dx,NX,minX+dx*(NX-1));
    tEprintf(Verbose,
       "\tY: start %.1f; dy %.1f; ny %i=>stop %.1f\n",
       minY,dy,NY,minY+dy*(NY-1));
    tEprintf(Verbose,
       "\tZ: start %.1f; dz %.1f; nz %i=>stop %.1f\n",
       minZ,dz,NZ,minZ+dz*(NZ-1));
    tEprintf(Verbose,
       "\tT: start %.3fms; dt %.3fms; nt %i=>stop %.3fms\n",
       1000*minT,1000*dt,NT,1000*(minT+dt*(NT-1)));
    tEprintf(Verbose,"\tModel is ~%.2f million nodes\n",
       (float)NXYZ/1e6);
  }

  if(tag==FIXED_MEDIA){
    tEprintf(Verbose,
       "Setup stationary media acoustic model\n");
  }
}

//Now the master fixed media attenuative class functions
masterFixedAttenMediaModel::masterFixedAttenMediaModel(selector* selectors,int tids[],int mediaType,
                           int surfaceBCType,int gridMultiplier,
                           modelDefStruct* modelDef,parallelDefStruct* parallelDef,
                           char* modelName,int fdCoeffOrder,float* fdCoeffs,
                          int longRead)
:masterFixedMediaModel(surfaceBCType,gridMultiplier,modelDef,parallelDef,modelName) {
  _selectors=selectors;
  assert(_selectors->size(),
         "masterAcousticAttenModel--selectors array must have at least 1 member");
  initSend();
  sendMessage(AllProcesses,MESSAGE_INITIALIZE,"ii",
              mediaType,FALSE); //Need to send a pivot type (FALSE) just so the
  //                                        message is the correct length.

  initializeSlaves(tids,
                   modelName,fdCoeffOrder,fdCoeffs,longRead);

  int *temp,*nodesPerSelector;
  assert((temp=(int*)malloc(_selectors->size()*sizeof(int)))!=NULL,
         "masterAnelasticModel--unable to allocate %i ints for temp space",
         _selectors->size());
  assert((nodesPerSelector=(int*)malloc(_selectors->size()*sizeof(int)))!=NULL,
         "masterAnelasticModel--unable to allocate %i ints for nodesPerSelector",
         _selectors->size());
  for(int i=0;i<_selectors->size();nodesPerSelector[i++]=0);
  for(int i=0;i<NumProcs;i++){
    getMessage(tids[i],MESSAGE_INITIALIZE,"I",temp,_selectors->size());
    for(int j=0;j<_selectors->size();j++){
      nodesPerSelector[j]+=temp[j];
    }
  }
  
  for(int i=0;i<_selectors->size();i++)
    tEprintf(Verbose,"Selector #%i <= %i nodes\n",
             i+1,nodesPerSelector[i]);
  
  free(temp);
  free(nodesPerSelector);

  getAck(FIXED_MEDIA,tids);
}

void masterFixedAttenMediaModel::sendInitMessage(char* modelName,int fdCoeffOrder,float* fdCoeffs,int longRead,
                                            parallelDefStruct** slaveParallelDefs,
                                            int tids[],int i,int j,int k,
                                            int nxProc,int nxyProc,
                                            int coworkParadyme) {
  masterFixedMediaModel::sendInitMessage(modelName,fdCoeffOrder,fdCoeffs,longRead,slaveParallelDefs,
                                         tids,i,j,k,nxProc,nxyProc,coworkParadyme);
  int maxMessSize = _selectors->calcMessageSize(0);
  setMessageBuffer(maxMessSize);
  initSend();
  sendMessage(TID(i,j,k),MESSAGE_INITIALIZE,"ii",_selectors->size(),
              maxMessSize);
  for(int ii=0;ii<_selectors->size();++ii) {
    initSend();
    _selectors->packSelectorMessage(ii);
    sendMessage(TID(i,j,k),MESSAGE_INITIALIZE,NULL);
  }
}
