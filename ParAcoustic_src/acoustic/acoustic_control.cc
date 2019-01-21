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
 *  acoustic_control.cc
 *
 *
 *  Define slave model and dependent class functions.  These classes live
 *  on specific domains and handle the actual model specific to that 
 *  domain as well as the dependent variables for that domain.  The
 *  dependent class handles the time step updating for the dependent
 *  variables.
 *
 *  Defines the following class functions:
 *    slaveAcousticModel::slaveAcousticModel
 *    slaveAcousticModel::vp
 *    slaveAcousticModel::rho
 *    slaveAcousticModel::rhox
 *    slaveAcousticModel::rhoy
 *    slaveAcousticModel::rhoz
 *    slaveAcousticModel::resetTimeVector
 *    slaveAcousticReplacementModel::slaveAcousticReplacementModel
 *    slaveAcousticPseudoAttenuateModel::slaveAcousticPseudoAttenuateModel
 *    slaveAcousticDependent::slaveAcousticDependent
 *    slaveAcousticDependent::doCheckpoint
 *    slaveAcousticDependent::readCheckpoint
 *    slaveAcousticDependent::setBoundaryConditions
 *    slaveAcousticDependent::loadSeismogram
 *    slaveAcousticDependent::setInitialConditions
 *    slaveAcousticDependent::advanceVel
 *    slaveAcousticDependent::advanceStress
 *    slaveAcousticDependent::fullFieldOutput
 *    slaveAcousticReplacementDependent::advanceVel
 *    slaveAcousticReplacementDependent::advanceStress
 *    slaveAcousticReplacementDependent::updateAcousticVelReplacement
 *    slaveAcousticReplacementDependent::updateAcousticPressureReplacement
 *    slaveAcousticPseudoAttenuateDependent::updateAcousticVelReplacement
 *    slaveAcousticPseudoAttenuateDependent::updateAcousticPressureReplacement
 *
 */

#include "acoustic_control.hh"
#include "selector.hh"
#include "acousticBoundary.hh"
#include "updateAcousticCPMLFull.h"
#include "updateAcousticMPMLFull.h"
#include "updateAcousticSpongeFull.h"
#include "updateAcousticAttenMPMLFull.h"
#include "xtrautil.hh"
#include "sgfdSources.hh"
#include "sgfdReceiverNetwork.hh"

//The model is read from a cdf file. But the
// model initialization routine receives the message
// from the parent containing the name (as well as this model limits)
slaveAcousticModel::slaveAcousticModel(float rhoRatLim, float cLimit, bool finalizeInit):slaveSgfdModel(FALSE){
  DEF_MODEL_SIZE(modelDef());
  DEF_MODEL_SCALE(modelDef());
  DEF_MODEL_LIMITS(modelDef());

  //Fill in the differentiator coefficients, ignore the holberg stuff for the present.
  float cinner=9.0/8.0,couter=-1.0/24.0;
  float dinner=2.0/3.0,douter=-1.0/12.0;
  _fx[0]=12.0/11.0*dx/dz*cinner;
  _fx[1]=12.0/11.0*dx/dz*couter;
  _fy[0]=12.0/11.0*dy/dz*cinner;
  _fy[1]=12.0/11.0*dy/dz*couter;

  _rx[0]=2.0*scalarSpeed*dt/dx*dinner;
  _rx[1]=2.0*scalarSpeed*dt/dx*douter;
  _ry[0]=2.0*scalarSpeed*dt/dy*dinner;
  _ry[1]=2.0*scalarSpeed*dt/dy*douter;
  _rz[0]=2.0*scalarSpeed*dt/dz*dinner;
  _rz[1]=2.0*scalarSpeed*dt/dz*douter;

  //Allocate space for the material properties (Note, the previous version of the code
  // let the slaveModel allocate the space and then renamed it. The new slaveSgfdModel
  // does not by default allocate any space so it must be done here).
  assert((_bulk=(float*)malloc(NXYZ*sizeof(float)))!=NULL,
   "slaveAcousticModel--unable to allocate %i floats for bulkModulus",
   NXYZ);
  assert((_rho=(float*)malloc(NXYZ*sizeof(float)))!=NULL,
   "slaveAcousticModel--unable to allocate %i floats for rho",
   NXYZ);

  //Define variables related to order-switching
  rhoRatioLimit=rhoRatLim;
  compLimit=cLimit;
  _vxfunc=(unsigned char*)malloc(NXYZ*sizeof(char));
  _vyfunc=(unsigned char*)malloc(NXYZ*sizeof(char));
  _vzfunc=(unsigned char*)malloc(NXYZ*sizeof(char));
  _ssfunc=(unsigned char*)malloc(NXYZ*sizeof(char));
  
  assert(_vxfunc && _vyfunc && _vzfunc && _ssfunc,
         "slaveElasticModel::unable to allocate 4*%i memory for function identifiers\n",
         NXYZ);
  //Now read the values from the model file, the model name is read from the initialization
  // message in the parent initializer.
  readCDFModel(_modelName,NULL,NULL, //Have not implemented a different rho model
   modelDef()->procLim,  // follow elastic_sgfd changes if desired.
   _bulk,NULL,_rho,
   _vMin,_vMax);

  if(finalizeInit) {
    float bMMin=1e50,bMMax=0.0,rhoMin=1e50,rhoMax=0.0;
    convertVels(_vMin,_vMax,rhoMin,rhoMax,bMMin,bMMax);
    
    tEprintf(_doLongRead,
       "[%i]: Read model: velocity range %.2f to %.2f\n",
       messageRank(),_vMin,_vMax);

    //send a completion message
    initSend();
    packMessage("ff ff ff",
    _vMin,_vMax,rhoMin,rhoMax,
    bMMin,bMMax);
    sendMessage(Parent,MESSAGE_INITIALIZE,NULL);
    
    //and receive the global min and max in return
    getMessage(Parent,MESSAGE_INITIALIZE,"ff",&_vMin,&_vMax);
  }
}

void slaveAcousticModel::convertVels(float &vMin, float &vMax, float &rhoMin, float &rhoMax,
                                     float &bMMin, float &bMMax) {
  DEF_MODEL_SIZE(_modelDef);
  DEF_MODEL_SCALE(modelDef());
  
  //And convert the velocities read from the model into the desired parameters
  for(int i=0;i<NXYZ;i++){
    float alfa=_bulk[i];
    float rho=_rho[i];
    _vxfunc[i] = _vyfunc[i] = _vzfunc[i] = _ssfunc[i] = 0;
    
    assert(1e20>alfa && alfa>=0.,
           "slaveAcousticModel--alfa[%i]=%.6g; out of range",
           i,alfa);
    assert(1e20>rho && rho>1e-20,
           "slaveAcousticModel--rho[%i]=%.6g; out of range",
           i,rho);
    
    SETMINMAX(_vMin,_vMax,alfa,i);
    SETMINMAX(rhoMin,rhoMax,rho,i);
    
    //rescale speed and density
    _bulk[i]/=scalarSpeed;
    _rho[i]/=scalarDen;
  }
  
  //Now check to see if any adjacent rho values exceed a certain ratio and use O(2,2) updating for
  //those nodes
  for(int k=0;k<NZ;k++) {
    for(int j=0;j<NY;j++) {
      for(int i=0;i<NX;i++) {
        int index=i+j*NX+k*NXY;
        float comp = compLimit*_rho[index]*_bulk[index];
        if(i>0 && (_rho[index-1]/_rho[index]>rhoRatioLimit || _rho[index-1]*_bulk[index-1]>comp)) {
          _vxfunc[index]=1;
          _vzfunc[index]=_vyfunc[index]=1;
          _ssfunc[index-1]=1;
        }
        if(i<NX-1 && (_rho[index+1]/_rho[index]>rhoRatioLimit || _rho[index+1]*_bulk[index+1]>comp)) {
          if(i>0) _vxfunc[index-1]=1;
          _vzfunc[index]=_vyfunc[index]=1;
          _ssfunc[index+1]=1;
        }
        if(j>0 && (_rho[index-NX]/_rho[index]>rhoRatioLimit || _rho[index-NX]*_bulk[index-NX]>comp)) {
          _vyfunc[index]=1;
          _vzfunc[index]=_vxfunc[index]=1;
          _ssfunc[index-NX]=1;
        }
        if(j<NY-1 && (_rho[index+NX]/_rho[index]>rhoRatioLimit || _rho[index+NX]*_bulk[index+NX]>comp)) {
          if(j>0) _vyfunc[index-NX]=1;
          _vzfunc[index]=_vxfunc[index]=1;
          _ssfunc[index+NX]=1;
        }
        if(k>0 && (_rho[index-NXY]/_rho[index]>rhoRatioLimit || _rho[index-NXY]*_bulk[index-NXY]>comp)) {
          _vzfunc[index]=1;
          _vxfunc[index]=_vyfunc[index]=1;
          _ssfunc[index-NXY]=1;
        }
        if(k<NZ-1 && (_rho[index+NXY]/_rho[index]>rhoRatioLimit || _rho[index+NXY]*_bulk[index+NXY]>comp)) {
          if(k>0) _vzfunc[index-NXY]=1;
          _vxfunc[index]=_vyfunc[index]=1;
          _ssfunc[index+NXY]=1;
        }
      }
    }
  }
  
  int numVx=0, numVy=0, numVz=0, numSS=0;
  //Now convert to divided (more computationally efficient) values.
  for(int i=0;i<NXYZ;i++){
    if(_vzfunc[i]) numVz++;
    if(_vyfunc[i]) numVy++;
    if(_vxfunc[i]) numVx++;
    if(_ssfunc[i]) numSS++;
    //solve for bulk
    _bulk[i]=_rho[i]*SQRNPS(_bulk[i]);
    SETMINMAX(bMMin,bMMax,_bulk[i],i);
    
    //Now convert to divided (more computationally efficient) form.
    _rho[i]/=2.0;
  }
}

//the following are helper routines for retrieving various model parameters
float slaveAcousticModel::vp(int i){
  return modelDef()->scalarSpeed*sqrt(0.5*_bulk[i]/_rho[i]);
}

float slaveAcousticModel::rho(int i){
  return 2.*_rho[i]*modelDef()->scalarDen;
}

//If a new time vector is requested on the command line
int slaveAcousticModel::resetTimeVector(){
  float minT,dt;
  int nt;
  unpackMessage("fif",&minT,&nt,&dt);
  
  //we need to fix some values that were based on the original time step
  _cx[0]*=dt/_modelDef->dt;
  _cx[1]*=dt/_modelDef->dt;
  _cy[0]*=dt/_modelDef->dt;
  _cy[1]*=dt/_modelDef->dt;
  _cz[0]*=dt/_modelDef->dt;
  _cz[1]*=dt/_modelDef->dt;
  
  _rx[0]*=dt/_modelDef->dt;
  _rx[1]*=dt/_modelDef->dt;
  _ry[0]*=dt/_modelDef->dt;
  _ry[1]*=dt/_modelDef->dt;
  _rz[0]*=dt/_modelDef->dt;
  _rz[1]*=dt/_modelDef->dt;
  
  _modelDef->NT=nt;
  _modelDef->dt=dt;
  _modelDef->minT=minT;
  
  return nt;
}

//The constructor for the fixed acoustic attenuative model with standard linear fluid mechanisms
slaveAcousticAttenModel::slaveAcousticAttenModel(float rhoRatLim, float cLimit): slaveAcousticModel(rhoRatLim,cLimit,false) {
  DEF_MODEL_SIZE(_modelDef);
  DEF_MODEL_LIMITS(_modelDef);

  float *vp=_bulk,*rho=_rho;

  //Start the selector.
  _selectors=new selector(NX,NY,NZ,minX,dx,minY,dy,minZ,dz);
  //Add the fields needed for making selections.
  _selectors->addField("C",vp);
  _selectors->addField("Rho",rho);
  
  getMessage(Parent,MESSAGE_INITIALIZE,NULL);
  int nSelector;
  unpackMessage("i",&nSelector);
  
  int maxMessSize;
  unpackMessage("i",&maxMessSize);
  setMessageBuffer(maxMessSize);
  //Now read the selector information from the message.
  for(int ii=0;ii<nSelector;ii++) {
    getMessage(Parent,MESSAGE_INITIALIZE,NULL);
    _selectors->unpackSelectorMessage(ii);
  }

  //Determine the maximum number of relaxation mechanisms, in any of the
  // selectors.
  for(int i=0;i<_selectors->size();i++)
    SETMAX(_numRelaxation,_selectors->intVar(0,i),i);
  
  //Allocate space to count the number of nodes allocated to each selector.
  int *nodesPerSelector;
  assert((nodesPerSelector=(int*)malloc(_selectors->size()*sizeof(int)))!=NULL,
         "masterAnelasticModel--unable to allocate %i ints for nodesPerSelector",
         _selectors->size());
  for(int i=0;i<_selectors->size();nodesPerSelector[i++]=0);
  
  //Get the correct Q for each node.
  _Q=NULL;
  if((_Q=readCDFQindex(_modelName,modelDef()->procLim,_Q))==NULL)
    _Q=_selectors->generateIndexArray();
  else {
    printf("Read Qindex from file\n");
  }
  
  _nMechs=_selectors->generateIntVar(0);
  fprintf(stderr,"nMechs: %d\n",_numRelaxation);
  _nM = -1; //only used if Full updating
  
  _decayRates=_selectors->generateFloatPtrVar(0);
  _ampP=_selectors->generateFloatPtrVar(1);
  
  //Count the number of nodes/selector, this is a good check that the call
  // was done correctly. This also where the velocities are corrected.
  for(int i=0;i<NXYZ;i++){
    assert(_Q[i]>=0 && _Q[i]<_selectors->size(),"Number of selectors (%d) does not match Qindex (%d)",
           _selectors->size(),_Q[i]);
    nodesPerSelector[_Q[i]]++;
    
    vp[i]*= _selectors->floatVar(0,_Q[i]);
  }

  float bMMin=1e50,bMMax=0.0,rhoMin=1e50,rhoMax=0.0;
  convertVels(_vMin,_vMax,rhoMin,rhoMax,bMMin,bMMax);
  
  //Adjust the variables to save on later operations.
  _omegaAmpSum = new float[_selectors->size()];
  for(int i=0;i<_selectors->size();i++){
    _omegaAmpSum[i] = 0.f;
    for(int j=0;j<_selectors->intVar(0,i);j++){
      _selectors->floatPtrVar(0,i)[j]*= dt;
      float om = _selectors->floatPtrVar(0,i)[j];
      _omegaAmpSum[i] += om*_selectors->floatPtrVar(1,i)[j]/(2.f+om);
    }
    _omegaAmpSum[i] -= 1.f;
  }
  
  //swap the order of decayRates and pAmp for full updating
  completeFullSetup();

  tEprintf(_doLongRead,
           "[%i]: Read model: velocity range %.2f to %.2f\n",
           messageRank(),_vMin,_vMax);
  
  //send a completion message
  initSend();
  packMessage("ff ff ff",
              _vMin,_vMax,rhoMin,rhoMax,
              bMMin,bMMax);
  sendMessage(Parent,MESSAGE_INITIALIZE,NULL);
  
  //and receive the global min and max in return
  getMessage(Parent,MESSAGE_INITIALIZE,"ff",&_vMin,&_vMax);

  //Send a final message with the number of nodes per Selector.
  initSend();
  sendMessage(Parent,MESSAGE_INITIALIZE,"I",nodesPerSelector,_selectors->size());
  
  //Clean up.
  free(nodesPerSelector);
}

slaveAcousticAttenModel::~slaveAcousticAttenModel(){
  delete _selectors;
  delete[] _omegaAmpSum;
  free(_Q);
  if(_nM!=-1) {
    for(int i=0;i<_nM;++i) {
      delete[] _decayRates[i];
      delete[] _ampP[i];
    }
    delete[] _decayRates;
    delete[] _ampP;
  }
}

//for the cases that use the Full update we need to swap the order of some of the
//variables and fill the _nM variable
void slaveAcousticAttenModel::completeFullSetup() {
  int ssize = _selectors->size();
  _nM = _nMechs[0];
  for(int i=1;i<ssize;++i) {
    assert(_nMechs[i]==_nM,"All the Q lines must have the same number of mechanisms for Full updating\n");
  }
  float** tmp1 = _decayRates;
  float** tmp2 = _ampP;
  _decayRates = new float*[_nM];
  _ampP = new float*[_nM];
  for(int i=0;i<_nM;++i) {
    _decayRates[i] = new float[ssize];
    _ampP[i] = new float[ssize];
  }
  for(int i=0;i<ssize;++i) { //i is the number of Q lines
    for(int j=0;j<_nM;++j) {
      _decayRates[j][i] = tmp1[i][j];
      _ampP[j][i] = tmp2[i][j];
    }
  }
  free(tmp1); //free only the outermost because the selectors hold the inner pointers
  free(tmp2);
}

//Main dependent constructor.  It handles updates of dependent variables of pressure and velocities
slaveAcousticDependent::slaveAcousticDependent(slaveAcousticModel* model,
     int doAllocate):
  slaveSgfdDependent(model){
  _model=model;
  _sources=NULL;
  _receivers=NULL;
  _boundaries=NULL;

  getMessage(Parent,MESSAGE_INITIALIZE,NULL);

  //Allocate the dependent variables.
  if(doAllocate){
    DEF_MODEL_SIZE(modelDef());
    _vx=(float*)malloc(NXYZ*sizeof(float));
    _vy=(float*)malloc(NXYZ*sizeof(float));
    _vz=(float*)malloc(NXYZ*sizeof(float));
    assert(_vx&&_vy&&_vz,
     "slaveAcousticDependent--unable to allocate 3*%i memory for velocity\n",
     NXYZ);

    assert(( _pressure=(float*)malloc(NXYZ*sizeof(float)))!=NULL,
     "slaveAcousticDependent--unable to allocate %i memory for pressure\n",
     NXYZ);

    initVariables();
  }

  doBarrier();
#if USE_VAMPIR
  VT_classdef("Updates",&_VTUpdateClass_);
  VT_classdef("Boundary_Conditions",&_VTBCClass_);

  VT_funcdef("Vel_Update",_VTUpdateClass_,&_VTVelUpdate_);
  VT_funcdef("Str_Update",_VTUpdateClass_,&_VTStressUpdate_);

  VT_funcdef("Vel_BC",_VTBCClass_,&_VTBCDoVel_);
  VT_funcdef("Str_BC",_VTBCClass_,&_VTBCDoStress_);
#endif
}

///Here is a method that will write a set ofcheckpoint (restart)
/// files to the given directory.
//This has not been tested recently, so it may or may not work
FILE* slaveAcousticDependent::doCheckpoint(int iteration,char* cpDir,int cpID,
         int callDepth,
         FILE* cpFile){
  if(!cpFile)
    cpFile=slaveSgfdDependent::doCheckpoint(iteration,cpDir,cpID,
              callDepth+1);

  //Need to write the dependent variable values out to the file.
  DEF_MODEL_SIZE(modelDef());
  fwrite(vx(),sizeof(float),NXYZ,cpFile);
  fwrite(vy(),sizeof(float),NXYZ,cpFile);
  fwrite(vz(),sizeof(float),NXYZ,cpFile);
  fwrite(P(),sizeof(float),NXYZ,cpFile);

  if(!callDepth){
    fclose(cpFile);
    cpFile=NULL;
  }
  return cpFile;
}
///Here is a method to read a checkpoint file written by doCheckpoint.
//This has not been tested recently, so it may or may not work
FILE* slaveAcousticDependent::readCheckpoint(int& iteration,char* cpDir,int cpID,
           int callDepth,
           FILE* cpFile){
  if(!cpFile)
    cpFile=slaveSgfdDependent::readCheckpoint(iteration,cpDir,cpID,
              callDepth+1);

  //Need to read the dependent variable values out of the file.
  DEF_MODEL_SIZE(modelDef());
  fread(vx(),sizeof(float),NXYZ,cpFile);
  fread(vy(),sizeof(float),NXYZ,cpFile);
  fread(vz(),sizeof(float),NXYZ,cpFile);
  fread(P(),sizeof(float),NXYZ,cpFile);

  if(!callDepth){
    fclose(cpFile);
    cpFile=NULL;
  }
  return cpFile;
}

//Really just a convience method so that the dependent class has more direct access to the boundary condition class values and methods
void slaveAcousticDependent::setBoundaryConditions(acousticSgfdBoundaries* boundaries){
  _boundaries=boundaries;
  //Z-min side.
  if(boundaries->isRealBoundary(0,0,-1) && boundaries->surfaceBC(0,0,-1)){
    _kstart=3;
  }else{
    _kstart=2;
  }
}

//Fill the receiver traces with data
int slaveAcousticDependent::loadSeismogram(int iteration){
  return
  _receivers->fillReceivers(modelDef(),iteration,
                            vx(),vy(),vz(),
                            P(),
                            NULL,NULL,
                            NULL,NULL,NULL);
}

//One can set up non-zero intial conditions for the dependent variables.
//This has nto been tested recently
void slaveAcousticDependent::setInitialConditions(){
  DEF_MODEL_SIZE(modelDef());
  //Read the initial conditions
  char buffer[1024];
  getMessage(Parent,MESSAGE_INITIALIZE,"s",buffer);
  if(!strcmp(buffer,"CHECKPOINT_RESTART")){
    char cpDir[1024];
    int iteration,cpID;
    getMessage(Parent,MESSAGE_GENERAL,"s isi",
   buffer,&iteration,cpDir,&cpID);
    assert(!strcmp(buffer,"READ_CHECKPOINT"),
     "setInitialConditions--message mismatch, should be %s (%s)",
     "READ_CHECKPOINT",buffer);
    readCheckpoint(iteration,cpDir,cpID,0);
    return;
  }

  //Options.
  int negateVelocities;
  unpackMessage("i",&negateVelocities);

  if(buffer[0]){
    assert(readCDFVarBlock(buffer,
         "vx",&_vx,_model->procLim()),
     "setInitialConditions--unable to read vx init cond from %s",
     buffer);
    assert(readCDFVarBlock(buffer,
         "vy",&_vy,_model->procLim()),
     "setInitialConditions--unable to read vy init cond from %s",
     buffer);
    assert(readCDFVarBlock(buffer,
         "vz",&_vz,_model->procLim()),
     "setInitialConditions--unable to read vz init cond from %s",
     buffer);

    assert(readCDFVarBlock(buffer,
         "pressure",&_pressure,_model->procLim()),
     "setInitialConditions--unable to read pressure init cond from %s",
     buffer);
  }

  //I believe this is a means of marching the algorithm back in time if desired, but I have never used it
  if(negateVelocities){
    for(int i=0;i<NXYZ;i++){
_vx[i]*=-1;
_vy[i]*=-1;
_vz[i]*=-1;
    }
  }
  initSend();
  sendMessage(Parent,MESSAGE_INITIALIZE,NULL);
}

//This function takes care of advancing the three components of the
//velocities forward one time step.  It takes care of boundary conditions,
//source insertion.  Either CPML or Sponge boundaries are supported,
//as well as the pressure-free surface
int slaveAcousticDependent::advanceVel(int iteration,int acknowledge){
  //Do the actual update.
#if USE_VAMPIR
  VT_begin(_VTVelUpdate_);
#endif
  if(_boundaries->_usePML) {
    if(_boundaries->_usePML==1 || _boundaries->_usePML==3)  {
      updateVxCPMLBounds(_model->modelDef(),
                         _vx,
                         _pressure,
                         _model->_cx,_model->_cy,_model->_cz,
                         _model->_rho,_model->_vxfunc,
                         _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                         _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                         _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      updateVyCPMLBounds(_model->modelDef(),
                         _vy,
                         _pressure,
                         _model->_cx,_model->_cy,_model->_cz,
                         _model->_rho,_model->_vyfunc,
                         _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                         _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                         _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      updateVzCPMLBounds(_model->modelDef(),
                         _vz,
                         _pressure,
                         _model->_cx,_model->_cy,_model->_cz,
                         _model->_rho,_model->_vzfunc,
                         _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                         _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                         _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      /*if(_boundaries->_freeSurfaceMode==SURFACE_PRESS_FREE && _boundaries->isRealBoundary(0,0,-1))
        updateVzPressFreeCPMLBounds(_model->modelDef(),
                                    _vz,
                                    _pressure,
                                    _model->_cx,_model->_cy,_model->_cz,
                                    _model->_rho,_model->_vzfunc,
                                    _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                    _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);*/
    } else if(_boundaries->_usePML==4) {
      updateVxMPMLBounds(_model->modelDef(),
                         _vx,
                         _pressure,
                         _model->_cx,_model->_cy,_model->_cz,
                         _model->_rho,_model->_vxfunc,
                         _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                         _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                         _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      updateVyMPMLBounds(_model->modelDef(),
                         _vy,
                         _pressure,
                         _model->_cx,_model->_cy,_model->_cz,
                         _model->_rho,_model->_vyfunc,
                         _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                         _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                         _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      updateVzMPMLBounds(_model->modelDef(),
                         _vz,
                         _pressure,
                         _model->_cx,_model->_cy,_model->_cz,
                         _model->_rho,_model->_vzfunc,
                         _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                         _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                         _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      /*if(_boundaries->_freeSurfaceMode==SURFACE_PRESS_FREE && _boundaries->isRealBoundary(0,0,-1))
        updateVzPressFreeMPMLBounds(_model->modelDef(),
                                    _vz,
                                    _pressure,
                                    _model->_cx,_model->_cy,_model->_cz,
                                    _model->_rho,_model->_vzfunc,
                                    _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                    _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);*/
    }
  } else {
    updateVxSpongeBounds(_model->modelDef(),
                       _vx,
                       _pressure,
                       _model->_cx,_model->_cy,_model->_cz,
                       _model->_rho,_model->_vxfunc,
                       _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                       _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                       _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                       _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
    updateVySpongeBounds(_model->modelDef(),
                       _vy,
                       _pressure,
                       _model->_cx,_model->_cy,_model->_cz,
                       _model->_rho,_model->_vyfunc,
                       _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                       _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                       _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                       _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
    updateVzSpongeBounds(_model->modelDef(),
                       _vz,
                       _pressure,
                       _model->_cx,_model->_cy,_model->_cz,
                       _model->_rho,_model->_vzfunc,
                       _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                       _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                       _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                       _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
    /*if(_boundaries->_freeSurfaceMode==SURFACE_PRESS_FREE && _boundaries->isRealBoundary(0,0,-1))
      updateVzPressFreeSpongeBounds(_model->modelDef(),
                                  _vz,
                                  _pressure,
                                  _model->_cx,_model->_cy,_model->_cz,
                                  _model->_rho,_model->_vzfunc,
                                  _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                  _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);*/
  }

  //Add any force or anti-symetric moment sources if not reading
  // time-dependent bc
  if(_sources && iteration<_model->modelDef()->NT)
    _sources->applyAcousticVel(_model->modelDef(),iteration,
       _vx,_vy,_vz,
       _model->_rho);

#if USE_VAMPIR
  VT_end(_VTVelUpdate_);
  VT_begin(_VTBCDoVel_);
#endif
  //this just starts the message passing among domains
  _boundaries->updateVel(iteration,this);

#if USE_VAMPIR
VT_end(_VTBCDoVel_);
#endif
  
  //fill in velocity receivers
  loadSeismogram(iteration);
  if(acknowledge){
    initSend();
    sendMessage(Parent,MESSAGE_UPDATE_VEL,NULL);
  }

  //make sure we are finished with the velocity message passing before
  //proceding
#if USE_MPI_SEND & USE_IMMEDIATE_MPI_SEND
  _boundaries->checkPassComplete();
#endif
  //push vz above the free surface if there is one
  if(_boundaries->_freeSurfaceMode==SURFACE_PRESS_FREE && _boundaries->isRealBoundary(0,0,-1)) {
    if(_boundaries->_usePML==1 || _boundaries->_usePML==3)
      updateVzPressFreeCPMLBounds(_model->modelDef(),
                                  _vz,
                                  _pressure,
                                  _model->_cx,_model->_cy,_model->_cz,
                                  _model->_rho,_model->_vzfunc,
                                  _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                  _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);
    else if(_boundaries->_usePML==4)
      updateVzPressFreeMPMLBounds(_model->modelDef(),
                                  _vz,
                                  _pressure,
                                  _model->_cx,_model->_cy,_model->_cz,
                                  _model->_rho,_model->_vzfunc,
                                  _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                  _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);
    else if(!_boundaries->_usePML)
      updateVzPressFreeSpongeBounds(_model->modelDef(),
                                    _vz,
                                    _pressure,
                                    _model->_cx,_model->_cy,_model->_cz,
                                    _model->_rho,_model->_vzfunc,
                                    _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                    _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);
  }
  
  return iteration;
}
// Perform and update to the stress grid for one time step.  It takes
//care of boundary conditions ((CPML or Sponge) and pressure-free
//surface, if applicable.  It also inserts sources that affect pressure
//directly.
int slaveAcousticDependent::advanceStress(int iteration,int acknowledge){

  //do the actual update
#if USE_VAMPIR
  VT_begin(_VTStressUpdate_);
#endif
  if(_boundaries->_usePML) {
    if(_boundaries->_usePML==1 || _boundaries->_usePML==3) {
      updateAcousticPressureCPML(_model->modelDef(),
                                    _pressure,_vx,_vy,_vz,
                                    _model->_cx,_model->_cy,_model->_cz,
                                    _model->_bulk,_model->_ssfunc,_boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                 _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                                 _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                                 _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      if(_boundaries->_freeSurfaceMode==SURFACE_PRESS_FREE && _boundaries->isRealBoundary(0,0,-1))
        updateAcousticPressurePressFreeCPML(_model->modelDef(),
                                   _pressure,_vx,_vy,_vz,
                                   _model->_cx,_model->_cy,_model->_cz,
                                   _model->_bulk,_model->_ssfunc,_boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                   _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);
    } else if(_boundaries->_usePML==4) {
      updateAcousticPressureMPML(_model->modelDef(),
                                 _pressure,_vx,_vy,_vz,
                                 _model->_cx,_model->_cy,_model->_cz,
                                 _model->_bulk,_model->_ssfunc,_boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                 _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                                 _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                                 _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      if(_boundaries->_freeSurfaceMode==SURFACE_PRESS_FREE && _boundaries->isRealBoundary(0,0,-1))
        updateAcousticPressurePressFreeMPML(_model->modelDef(),
                                            _pressure,_vx,_vy,_vz,
                                            _model->_cx,_model->_cy,_model->_cz,
                                            _model->_bulk,_model->_ssfunc,_boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                            _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);
    }
 } else {
   updateAcousticPressureSponge(_model->modelDef(),
                              _pressure,_vx,_vy,_vz,
                              _model->_cx,_model->_cy,_model->_cz,
                              _model->_bulk,_model->_ssfunc,_boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                              _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                              _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                              _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
   if(_boundaries->_freeSurfaceMode==SURFACE_PRESS_FREE && _boundaries->isRealBoundary(0,0,-1))
     updateAcousticPressurePressFreeSponge(_model->modelDef(),
                                         _pressure,_vx,_vy,_vz,
                                         _model->_cx,_model->_cy,_model->_cz,
                                         _model->_bulk,_model->_ssfunc,_boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);
  }
#if USE_VAMPIR
  VT_end(_VTStressUpdate_);
#endif

  //add any symetric moment and/or traction sources if not reading 
  // time-dependent bc
  if(_sources && iteration<_model->modelDef()->NT)
    _sources->applyPressure(_model->modelDef(),iteration+1,
          _pressure);

  //and update the boundarys and begin message passing
#if USE_VAMPIR
 VT_begin(_VTBCDoStress_);
#endif
  _boundaries->updateStress(iteration,this);
  
#if USE_VAMPIR
 VT_end(_VTBCDoStress_);
#endif
  if(acknowledge){
    initSend();
    sendMessage(Parent,MESSAGE_UPDATE_STRESS,NULL);
  }
  //check and see that we are done passing stresses
#if USE_MPI_SEND & USE_IMMEDIATE_MPI_SEND
  _boundaries->checkPassComplete();
#endif

  return iteration;
}

///Define a real version of fullFieldOutput (required for a real class).
//This has not been tested recently
int slaveAcousticDependent::fullFieldOutput(int recepient,
          char* msgtag){
  //Read the string from the message so we know which field we are sending back.
  char buffer[512];
  if(msgtag){
    strcpy(buffer,msgtag);
  }else{
    unpackMessage("s",buffer);
  }
  initSend();
  sendMessage(recepient,MESSAGE_FULL_FIELD,"P",_model->procLim());

  if(!strcmp(buffer,"vx")){
    return slaveSgfdDependent::fullFieldOutput(recepient,_vx);
  }else if(!strcmp(buffer,"vy")){
    return slaveSgfdDependent::fullFieldOutput(recepient,_vy);
  }else if(!strcmp(buffer,"vz")){
    return slaveSgfdDependent::fullFieldOutput(recepient,_vz);
  }else if(!strcmp(buffer,"pressure")){
    return slaveSgfdDependent::fullFieldOutput(recepient,_pressure);
  }
  return
    assert(FALSE,"slaveAcousticDependent::fullFieldOutput--unknown comp %s",
     buffer);
}

//slave dependent for fixed acoustic attenuative models with standard linear fluid mechs
slaveAcousticAttenDependent::slaveAcousticAttenDependent(slaveAcousticAttenModel* model,
                                                         int doAllocate): slaveAcousticDependent(model) {
  _model = model;
  DEF_MODEL_SIZE(_model->_modelDef);
  
  int numR=_model->_numRelaxation;
  //allocate the memory variables
  assert((_rP=(float**)malloc(numR*sizeof(float*)))!=NULL,
         "slaveAcousticAttenDependent--unable to allocate %i floatPtr for rP",
         numR);
  for(int i=0;i<numR;i++){
    assert((_rP[i]=(float*)malloc(NXYZ*sizeof(float)))!=NULL,
           "slaveAcousticAttenDependent--unable to allocate %i floats for rP[%i]",
           NXYZ,i);
    for(int j=0;j<NXYZ;++j)
      _rP[i][j] = 0.0;
  }
}

//The velocity update for attenuative models.  Only real difference is if the free surface is on
int slaveAcousticAttenDependent::advanceVel(int iteration,int acknowledge) {
  //Do the actual update.
#if USE_VAMPIR
  VT_begin(_VTVelUpdate_);
#endif
  if(_boundaries->_usePML) {
    if(_boundaries->_usePML==1 || _boundaries->_usePML==3)  {
      updateVxCPMLBounds(_model->modelDef(),
                         _vx,
                         _pressure,
                         _model->_cx,_model->_cy,_model->_cz,
                         _model->_rho,_model->_vxfunc,
                         _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                         _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                         _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      updateVyCPMLBounds(_model->modelDef(),
                         _vy,
                         _pressure,
                         _model->_cx,_model->_cy,_model->_cz,
                         _model->_rho,_model->_vyfunc,
                         _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                         _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                         _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      updateVzCPMLBounds(_model->modelDef(),
                         _vz,
                         _pressure,
                         _model->_cx,_model->_cy,_model->_cz,
                         _model->_rho,_model->_vzfunc,
                         _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                         _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                         _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      /*if(_boundaries->_freeSurfaceMode==SURFACE_PRESS_FREE && _boundaries->isRealBoundary(0,0,-1))
       updateVzPressFreeCPMLBounds(_model->modelDef(),
       _vz,
       _pressure,
       _model->_cx,_model->_cy,_model->_cz,
       _model->_rho,_model->_vzfunc,
       _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
       _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);*/
    } else if(_boundaries->_usePML==4) {
      updateVxAttenMPMLBounds(_model->modelDef(),
                         _vx,
                         _pressure,
                         _model->_cx,_model->_cy,_model->_cz,
                         _model->_rho,_model->_vxfunc,
                         _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                         _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                         _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      updateVyAttenMPMLBounds(_model->modelDef(),
                         _vy,
                         _pressure,
                         _model->_cx,_model->_cy,_model->_cz,
                         _model->_rho,_model->_vyfunc,
                         _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                         _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                         _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      updateVzAttenMPMLBounds(_model->modelDef(),
                         _vz,
                         _pressure,
                         _model->_cx,_model->_cy,_model->_cz,
                         _model->_rho,_model->_vzfunc,
                         _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                         _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                         _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      /*if(_boundaries->_freeSurfaceMode==SURFACE_PRESS_FREE && _boundaries->isRealBoundary(0,0,-1))
       updateVzPressFreeMPMLBounds(_model->modelDef(),
       _vz,
       _pressure,
       _model->_cx,_model->_cy,_model->_cz,
       _model->_rho,_model->_vzfunc,
       _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
       _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);*/
    }
  } else {
    updateVxSpongeBounds(_model->modelDef(),
                         _vx,
                         _pressure,
                         _model->_cx,_model->_cy,_model->_cz,
                         _model->_rho,_model->_vxfunc,
                         _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                         _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                         _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
    updateVySpongeBounds(_model->modelDef(),
                         _vy,
                         _pressure,
                         _model->_cx,_model->_cy,_model->_cz,
                         _model->_rho,_model->_vyfunc,
                         _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                         _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                         _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
    updateVzSpongeBounds(_model->modelDef(),
                         _vz,
                         _pressure,
                         _model->_cx,_model->_cy,_model->_cz,
                         _model->_rho,_model->_vzfunc,
                         _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                         _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                         _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                         _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
    /*if(_boundaries->_freeSurfaceMode==SURFACE_PRESS_FREE && _boundaries->isRealBoundary(0,0,-1))
     updateVzPressFreeSpongeBounds(_model->modelDef(),
     _vz,
     _pressure,
     _model->_cx,_model->_cy,_model->_cz,
     _model->_rho,_model->_vzfunc,
     _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
     _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);*/
  }
  
  //Add any force or anti-symetric moment sources if not reading
  // time-dependent bc
  if(_sources && iteration<_model->modelDef()->NT)
    _sources->applyAcousticVel(_model->modelDef(),iteration,
                               _vx,_vy,_vz,
                               _model->_rho);
  
#if USE_VAMPIR
  VT_end(_VTVelUpdate_);
  VT_begin(_VTBCDoVel_);
#endif
  //this just starts the message passing among domains
  _boundaries->updateVel(iteration,this);
  
#if USE_VAMPIR
  VT_end(_VTBCDoVel_);
#endif
  
  //fill in velocity receivers
  loadSeismogram(iteration);
  if(acknowledge){
    initSend();
    sendMessage(Parent,MESSAGE_UPDATE_VEL,NULL);
  }
  
  //make sure we are finished with the velocity message passing before
  //proceding
#if USE_MPI_SEND & USE_IMMEDIATE_MPI_SEND
  _boundaries->checkPassComplete();
#endif
  //push vz above the free surface if there is one
  if(_boundaries->_freeSurfaceMode==SURFACE_PRESS_FREE && _boundaries->isRealBoundary(0,0,-1)) {
    if(_boundaries->_usePML==1 || _boundaries->_usePML==3)
      assert(FALSE,"Not implemented\n");
      /*updateVzPressFreeCPMLBounds(_model->modelDef(),
                                  _vz,
                                  _pressure,
                                  _model->_cx,_model->_cy,_model->_cz,
                                  _model->_rho,_model->_vzfunc,
                                  _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                  _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);*/
    else if(_boundaries->_usePML==4)
      updateVzPressFreeAttenMPMLBounds(_model->modelDef(),
                                  _vz,
                                  _pressure,
                                  _model->_cx,_model->_cy,_model->_cz,
                                  _model->_rho,_model->_vzfunc,
                                  _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                  _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);
    else if(!_boundaries->_usePML)
      assert(FALSE,"Not implemented\n");
      /*updateVzPressFreeSpongeBounds(_model->modelDef(),
                                    _vz,
                                    _pressure,
                                    _model->_cx,_model->_cy,_model->_cz,
                                    _model->_rho,_model->_vzfunc,
                                    _boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                    _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);*/
  }
  
  return iteration;
}

//The pressure update for attenuative models.  Major differences are here
int slaveAcousticAttenDependent::advanceStress(int iteration,int acknowledge) {
  //do the actual update
#if USE_VAMPIR
  VT_begin(_VTStressUpdate_);
#endif
  if(_boundaries->_usePML) {
    if(_boundaries->_usePML==1 || _boundaries->_usePML==3) {
      assert(FALSE,"Not implemented\n");
      updateAcousticPressureCPML(_model->modelDef(),
                                 _pressure,_vx,_vy,_vz,
                                 _model->_cx,_model->_cy,_model->_cz,
                                 _model->_bulk,_model->_ssfunc,_boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                 _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                                 _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                                 _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      if(_boundaries->_freeSurfaceMode==SURFACE_PRESS_FREE && _boundaries->isRealBoundary(0,0,-1))
        updateAcousticPressurePressFreeCPML(_model->modelDef(),
                                            _pressure,_vx,_vy,_vz,
                                            _model->_cx,_model->_cy,_model->_cz,
                                            _model->_bulk,_model->_ssfunc,_boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                            _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);
    } else if(_boundaries->_usePML==4) {
      updateAcousticAttenPressureMPML(_model->modelDef(),
                                      _pressure,_vx,_vy,_vz,_model->_Q,_model->_nM,_model->_decayRates,_model->_ampP,
                                      _rP,
                                 _model->_cx,_model->_cy,_model->_cz,
                                 _model->_bulk,_model->_ssfunc,_boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                 _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                                 _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                                 _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
      if(_boundaries->_freeSurfaceMode==SURFACE_PRESS_FREE && _boundaries->isRealBoundary(0,0,-1))
        updateAcousticAttenPressurePressFreeMPML(_model->modelDef(),
                                            _pressure,_vx,_vy,_vz,
                                            _model->_cx,_model->_cy,_model->_cz,
                                            _model->_bulk,_model->_ssfunc,_boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                            _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);
    }
  } else {
    assert(FALSE,"Not implemented\n");
    updateAcousticPressureSponge(_model->modelDef(),
                                 _pressure,_vx,_vy,_vz,
                                 _model->_cx,_model->_cy,_model->_cz,
                                 _model->_bulk,_model->_ssfunc,_boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                 _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb,
                                 _boundaries->fminXb,_boundaries->fmaxXb,_boundaries->fminYb,
                                 _boundaries->fmaxYb,_boundaries->fminZb,_boundaries->fmaxZb);
    if(_boundaries->_freeSurfaceMode==SURFACE_PRESS_FREE && _boundaries->isRealBoundary(0,0,-1))
      updateAcousticPressurePressFreeSponge(_model->modelDef(),
                                            _pressure,_vx,_vy,_vz,
                                            _model->_cx,_model->_cy,_model->_cz,
                                            _model->_bulk,_model->_ssfunc,_boundaries->minXb,_boundaries->maxXb,_boundaries->minYb,
                                            _boundaries->maxYb,_boundaries->minZb,_boundaries->maxZb);
  }
#if USE_VAMPIR
  VT_end(_VTStressUpdate_);
#endif
  
  //add any symetric moment and/or traction sources if not reading
  // time-dependent bc
  if(_sources && iteration<_model->modelDef()->NT)
    _sources->applyPressure(_model->modelDef(),iteration+1,
                            _pressure);
  
  //and update the boundarys and begin message passing
#if USE_VAMPIR
  VT_begin(_VTBCDoStress_);
#endif
  _boundaries->updateStress(iteration,this);
  
#if USE_VAMPIR
  VT_end(_VTBCDoStress_);
#endif
  if(acknowledge){
    initSend();
    sendMessage(Parent,MESSAGE_UPDATE_STRESS,NULL);
  }
  //check and see that we are done passing stresses
#if USE_MPI_SEND & USE_IMMEDIATE_MPI_SEND
  _boundaries->checkPassComplete();
#endif
  
  return iteration;
}
