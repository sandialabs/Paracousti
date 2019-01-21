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
 *  sgfdSources.cc
 *
 *
 *  Defines class functions used to represent sources.  Classes that represent individual
 *  sources are defined for each source type.  Also, the single class that manages all
 *  sources present in a run is defined as the sourceNetwork class.
 *
 *  Defines the following class functions:
 *  tdbcSource::tdbcSource
 *  tdbcSource::tdbcSource
 *  tdbcSource::~tdbcSource
 *  tdbcSource::packSourceMessage
 *  tdbcSource::applyVel
 *  tdbcSource::applyPressure
 *  tdbcSource::initializeEdgePoints
 *  tdbcSource::initializeUnitPoints
 *  tdbcSource::readNextBuffer
 *  tdbcSource::readNextEdgeBuffer
 *  pointSource::pointSource
 *  pointSource::pointSource
 *  pointSource::pointSource
 *  pointSource::readAuxMessage
 *  pointSource::setActive
 *  pointSource::packSourceMessage
 *  pointSource::writeSource
 *  momentSource::momentSource
 *  momentSource::momentSource
 *  momentSource::momentSource
 *  momentSource::packSourceMessage
 *  momentSource::writeSource
 *  momentSource::applyVel
 *  momentSource::applyStress
 *  momentSource::applyPressure
 *  momentSource::applyPressureMM
 *  momentSource::calcInterpCoeffs
 *  forceSource::forceSource
 *  forceSource::forceSource
 *  forceSource::forceSource
 *  forceSource::packSourceMessage
 *  forceSource::writeSource
 *  forceSource::applyVel
 *  forceSource::applyVelMM
 *  forceSource::calcInterpCoeffs
 *  tractionSource::tractionSource
 *  tractionSource::tractionSource
 *  tractionSource::tractionSource
 *  tractionSource::tractionSource
 *  tractionSource::packSourceMessage
 *  tractionSource::writeSource
 *  tractionSource::applyVel
 *  tractionSource::calcInterpCoeffs
 *  sourceNetwork::sourceNetwork
 *  sourceNetwork::sourceNetwork
 *  sourceNetwork::combine
 *  sourceNetwork::sendSources
 *  sourceNetwork::applyVel
 *  sourceNetwork::applyStress
 *  sourceNetwork::applyAcousticVel
 *  sourceNetwork::applyAcousticVelMM
 *  sourceNetwork::applyPressure
 *  sourceNetwork::applyPressureMM
 *  sourceNetwork::writeSources
 *  sourceNetwork::setDispersionFactor
 *  sourceNetwork::checkDispersion
 *  sourceNetwork::addSources
 *  sourceNetwork::addSources
 *  sourceNetwork::checkWaveletDispersion
 *  sourceNetwork::readWavelet
 *  sourceNetwork::rickerWavelet
 *  sourceNetwork::writeForceSources
 *  sourceNetwork::writeMomentSources
 *  sourceNetwork::writeTractionSources
 *  sourceNetwork::readForceSources
 *  sourceNetwork::readMomentSources
 *  sourceNetwork::readTractionSources
 *  sourceNetwork::cTimesCExp
 *  sourceNetwork::rickerSpec
 *  sourceNetwork::nFreqs
 *  sourceNetwork::sourceSpectrum
 *  sourceNetwork::bandit
 *
 */
//These utility functions are for sources in the SGFD problem

#include "sgfdSources.hh"

#include <math.h>
#include "nstdutil.h"
#include "netcdf.h"
#include "io_procs.h"
#include "message_passing.h"

auto sqr = [](float x){return x*x;};

//
//Here are the distributed source types.
#define TDBC_SOURCE_BUFFER_SIZE 100
tdbcSource::tdbcSource(modelDefStruct* modelDef,char* filename):source(){
  _doDifferences=FALSE;
  strcpy(_infileName,filename);
  
  //Check that the tdbc model size and the modelDef agree.
  int infile=openCDFFile(_infileName,FALSE,NC_NOWRITE);
  
  int dim;
  size_t test;
  DEF_MODEL_SIZE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  //Check NX.
  assert(nc_inq_dimid(infile,"NX",&dim)==NC_NOERR,
         "tdbcSource--unable to get dim NX from %s(%i)",
         _infileName,infile);
  assert(nc_inq_dimlen(infile,dim,&test)==NC_NOERR,
         "tdbcSource--unable to read dim NX from %s(%i,%i)",
         _infileName,infile,dim);
  assert(test==NX,
         "tdbcSource--NX mismatch %i(%s) != %i (model)",
         test,_infileName,NX);
  
  //Check NY.
  assert(nc_inq_dimid(infile,"NY",&dim)==NC_NOERR,
         "tdbcSource--unable to get dim NY from %s(%i)",
         _infileName,infile);
  assert(nc_inq_dimlen(infile,dim,&test)==NC_NOERR,
         "tdbcSource--unable to read dim NY from %s(%i,%i)",
         _infileName,infile,dim);
  assert(test==NY,
         "tdbcSource--NY mismatch %i(%s) != %i (model)",
         test,_infileName,NY);
  
  //Check NZ.
  assert(nc_inq_dimid(infile,"NZ",&dim)==NC_NOERR,
         "tdbcSource--unable to get dim NZ from %s(%i)",
         _infileName,infile);
  assert(nc_inq_dimlen(infile,dim,&test)==NC_NOERR,
         "tdbcSource--unable to read dim NZ from %s(%i,%i)",
         _infileName,infile,dim);
  assert(test==NZ,
         "tdbcSource--NZ mismatch %i(%s) != %i (model)",
         test,_infileName,NZ);
  
  //Check minima.
  float temp[4];
  int var;
  assert(nc_inq_varid(infile,"Minima",&var)==NC_NOERR,
         "tdbcSource--unable to get var Minima from %s(%i)",
         _infileName,infile);
  assert(nc_get_var_float(infile,var,temp)==NC_NOERR,
         "tdbcSource--unable to read var Minima from %s(%i)",
         _infileName,infile,var);
  assert(temp[0]==minX && temp[1]==minY && temp[2]==minZ && temp[3]==minT,
         "tdbcSource--minima mismatch: %.2f, %.2f, %.2f, %.2f != %.2f, %.2f, %.2f, %.2f ",
         temp[0],temp[1],temp[2],temp[3],minX,minY,minZ,minT);
  
  assert(nc_inq_varid(infile,"Increments",&var)==NC_NOERR,
         "tdbcSource--unable to get var Increments from %s(%i)",
         _infileName,infile);
  assert(nc_get_var_float(infile,var,temp)==NC_NOERR,
         "tdbcSource--unable to read var Increments from %s(%i)",
         _infileName,infile,var);
  assert(temp[0]==dx && temp[1]==dy && temp[2]==dz && temp[3]==dt,
         "tdbcSource--increment mismatch: %.2f, %.2f, %.2f, %.2f != %.2f, %.2f, %.2f, %.2f ",
         temp[0],temp[1],temp[2],temp[3],dx,dy,dz,dt);
  
  //Close the file.
  assert(nc_close(infile)==NC_NOERR,
         "tdbcSource--unable to close %s(%i)",
         _infileName,infile);
}

//Read a file corresponding to mapTDBC translation of John Holland's output
// format.
tdbcSource::tdbcSource(modelDefStruct* modelDef):source(){
  _doDifferences=FALSE;
  
  //Read the filename from a message.
  unpackMessage("s",_infileName);
  
  //Open the file for reading.
  _infile=openCDFFile(_infileName,FALSE,NC_NOWRITE);
  
  //Check for the existance of a timeCutoff Variable.
  int tcVar;
  if(nc_inq_varid(_infile,"timeCutoff",&tcVar)!=NC_NOERR){
    _timeCutoffIndex=-1;
  }else{
    assert(nc_get_var_int(_infile,tcVar,&_timeCutoffIndex)==NC_NOERR,
           "initializeUnitPoints--unable to read timeCutoff from %s(%i,%i)",
           _infileName,_infile,tcVar);
  }
  
  //Check for the dimension NumPoints to find out if this is a unit cell or
  // a decomposed type file.
  int numPointsDim;
  size_t numPoints;
  if(nc_inq_dimid(_infile,"NumPoints",&numPointsDim)==NC_NOERR &&
     nc_inq_dimlen(_infile,numPointsDim,&numPoints)==NC_NOERR){
    _nVxPoints=_nVyPoints=_nVzPoints=0;
    initializeUnitPoints(modelDef,numPointsDim,numPoints);
  }else{
    _nPoints=0;
    initializeEdgePoints(modelDef,"Vx",
                         _nVxPoints,_vxVar,
                         _vxPointIndicies,_vxGridIndicies,
                         &_vxBuffer);
    initializeEdgePoints(modelDef,"Vy",
                         _nVyPoints,_vyVar,
                         _vyPointIndicies,_vyGridIndicies,
                         &_vyBuffer);
    initializeEdgePoints(modelDef,"Vz",
                         _nVzPoints,_vzVar,
                         _vzPointIndicies,_vzGridIndicies,
                         &_vzBuffer);
    initializeEdgePoints(modelDef,"Pressure",
                         _nPressPoints,_pressVar,
                         _pressPointIndicies,_pressGridIndicies,
                         &_pressBuffer);
    if(!_nVxPoints && !_nVyPoints && !_nVzPoints && !_nPressPoints){
      assert(nc_close(_infile)==NC_NOERR,
             "initializeUnitPoints--unable to close %s(%i)",
             _infileName,_infile);
      _infile=-1;
    }else{
      readNextEdgeBuffer(modelDef,0);
    }
  }
}

tdbcSource::~tdbcSource(){
  if(_nPoints){
    assert(nc_close(_infile)==NC_NOERR,
           "~tdbcSource--unable to close input file %s(%i)",_infileName,_infile);
    
    free(_gridIndicies);
    free(_gridIndiciesMx);
    free(_gridIndiciesMy);
    free(_gridIndiciesMz);
    
    for(int i=0;i<_nPoints;i++){
      free(_vxBuffer[i]);
      free(_vyBuffer[i]);
      free(_vzBuffer[i]);
    }
    free(_vxBuffer);
    free(_vyBuffer);
    free(_vzBuffer);
  }
}

void tdbcSource::packSourceMessage(modelDefStruct* modelDef,
                                   int startIteration){
  if(!startIteration){
    packMessage("i",type());
    packMessage("s",_infileName);
  }
}
int tdbcSource::applyVel(modelDefStruct* modelDef,int iteration,
                         float *vx,float *vy,float *vz,
                         float* rho){
  if(_timeCutoffIndex>0 && iteration>_timeCutoffIndex){
    return iteration;
  }else if(_timeCutoffIndex>0 && iteration==_timeCutoffIndex){
    tEprintf(Verbose,
             "Applying TDBC for last time before cutoff at %i iterations\n",
             _timeCutoffIndex);
  }
  
  int currBIndex=iteration-_bufferIndex;
  if(_nPoints){
    //Apply the unit cell formulation.
    if(iteration>=_bufferIndex+TDBC_SOURCE_BUFFER_SIZE)
      readNextBuffer(modelDef,_bufferIndex+TDBC_SOURCE_BUFFER_SIZE);
    
    if(_doDifferences){
      assert(FALSE,
             "tdbcSource--difference pseudo-TDBC not implemented for unit cell");
    }else{
      for(int i=0;i<_nPoints;i++){
        int currGIndex=_gridIndicies[i],
        currGIndexMx=_gridIndiciesMx[i],currGIndexMy=_gridIndiciesMy[i],currGIndexMz=_gridIndiciesMz[i];
        if(currGIndex!=-1){
          vx[currGIndex]=_vxBuffer[i][currBIndex]/2.0;
          vy[currGIndex]=_vyBuffer[i][currBIndex]/2.0;
          vz[currGIndex]=_vzBuffer[i][currBIndex]/2.0;
        }
        
        if(currGIndexMx!=-1)
          vx[currGIndexMx]=_vxBuffer[i][currBIndex]/2.0;
        if(currGIndexMy!=-1)
          vy[currGIndexMy]=_vyBuffer[i][currBIndex]/2.0;
        if(currGIndexMz!=-1)
          vz[currGIndexMz]=_vzBuffer[i][currBIndex]/2.0;
      }
    }
  }else if(_nVxPoints || _nVyPoints || _nVzPoints){
    //Apply the decomposed into edge locations form.
    if(iteration>=_bufferIndex+TDBC_SOURCE_BUFFER_SIZE)
      readNextEdgeBuffer(modelDef,_bufferIndex+TDBC_SOURCE_BUFFER_SIZE);
    
    if(_doDifferences){
      for(int i=0;i<_nVxPoints;i++)
        vx[_vxGridIndicies[i]]+=
        (_vxBuffer[i][currBIndex]-
         _vxBuffer[i][currBIndex-1]);
      for(int i=0;i<_nVyPoints;i++)
        vy[_vyGridIndicies[i]]+=
        (_vyBuffer[i][currBIndex]-
         _vyBuffer[i][currBIndex-1]);
      for(int i=0;i<_nVzPoints;i++)
        vz[_vzGridIndicies[i]]+=
        (_vzBuffer[i][currBIndex]-
         _vzBuffer[i][currBIndex-1]);
    }else{
      for(int i=0;i<_nVxPoints;i++)
        vx[_vxGridIndicies[i]]=_vxBuffer[i][currBIndex];
      for(int i=0;i<_nVyPoints;i++)
        vy[_vyGridIndicies[i]]=_vyBuffer[i][currBIndex];
      for(int i=0;i<_nVzPoints;i++)
        vz[_vzGridIndicies[i]]=_vzBuffer[i][currBIndex];
    }
  }
  return TRUE;
}

int tdbcSource::applyPressure(modelDefStruct* modelDef,int iteration,
                              float *pressure){
  if(_timeCutoffIndex>0 && iteration>_timeCutoffIndex){
    return iteration;
  }else if(_timeCutoffIndex>0 && iteration==_timeCutoffIndex){
    tEprintf(Verbose,
             "Applying TDBC for last time before cutoff at %i iterations\n",
             _timeCutoffIndex);
  }
  
  int currBIndex=iteration-_bufferIndex;
  if(_nPoints){
    //Apply the unit cell formulation.
    if(iteration>=_bufferIndex+TDBC_SOURCE_BUFFER_SIZE)
      readNextBuffer(modelDef,_bufferIndex+TDBC_SOURCE_BUFFER_SIZE);
    
    if(_doDifferences){
      assert(FALSE,
             "tdbcSource--difference pseudo-TDBC not implemented for unit cell");
    }else{
      for(int i=0;i<_nPoints;i++){
        int currGIndex=_gridIndicies[i];
        if(currGIndex!=-1){
          pressure[currGIndex]=_pressBuffer[i][currBIndex];
        }
      }
    }
  }else if(_nPressPoints){
    //Apply the decomposed into edge locations form.
    if(iteration>=_bufferIndex+TDBC_SOURCE_BUFFER_SIZE)
      readNextEdgeBuffer(modelDef,_bufferIndex+TDBC_SOURCE_BUFFER_SIZE);
    for(int i=0;i<_nPressPoints;i++)
      pressure[_pressGridIndicies[i]]=0.5*(_pressBuffer[i][currBIndex]+_pressBuffer[i][currBIndex-1]);
  }
  return TRUE;
}

void tdbcSource::initializeEdgePoints(modelDefStruct* modelDef,const char* flag,
                                      int &nPoints,int &var,
                                      intPtr &pointIndicies,intPtr &gridIndicies,
                                      float*** buffer){
  DEF_MODEL_SIZE(modelDef);
  char nameBuffer[1024];
  //Read the dimension and get the number of points.
  int numPointsDim;
  sprintf(nameBuffer,"%sNumPoints",flag);
  assert(nc_inq_dimid(_infile,nameBuffer,&numPointsDim)==NC_NOERR,
         "initializeEdgePoints--unable to get dim %s from %s(%i)",
         nameBuffer,_infileName,_infile);
  size_t numPoints;
  assert(nc_inq_dimlen(_infile,numPointsDim,&numPoints)==NC_NOERR,
         "initializeEdgePoints--unable to read dim %s from %s(%i,%i)",
         nameBuffer,_infileName,_infile,numPointsDim);
  
  //Allocate space for indicies and read.
  int *indicies;
  assert((indicies=(int*)malloc(3*numPoints*sizeof(int)))!=NULL,
         "initializeEdgePoints--unable to allocate %i ints for indicies",
         numPoints);
  int indexVar;
  sprintf(nameBuffer,"%sIndicies",flag);
  assert(nc_inq_varid(_infile,nameBuffer,&indexVar)==NC_NOERR,
         "initializeEdgePoints--unable to get var %s from %s(%i)",
         nameBuffer,_infileName,_infile);
  assert(nc_get_var_int(_infile,indexVar,indicies)==NC_NOERR,
         "initializeEdgePoints--unable to read var %s from %s(%i,%i)",
         nameBuffer,_infileName,_infile,indexVar);
  
  //Get the variable.
  sprintf(nameBuffer,"%s",flag);
  assert(nc_inq_varid(_infile,nameBuffer,&var)==NC_NOERR,
         "initializeEdgePoints--unable to get var %s from %s(%i)",
         nameBuffer,_infileName,_infile);
  
  //Count the number of points that are in bounds for this process.
  nPoints=0;
  for(int i=0;i<numPoints;i++){
    if(ISMID_OC(procLim[0],indicies[0+3*i],procLim[1])&&
       ISMID_OC(procLim[2],indicies[1+3*i],procLim[3])&&
       ISMID_OC(procLim[4],indicies[2+3*i],procLim[5])){
      nPoints++;
    }
  }
  
  if(!nPoints){
    //No points are in bounds.
    gridIndicies=pointIndicies=NULL;
    *buffer=NULL;
  }else{
    //Allocate local indicies and fill in.
    assert((gridIndicies=(int*)malloc(nPoints*sizeof(int)))!=NULL,
           "initializeEdgePoints--unable to allocate %i ints for local indicies",
           nPoints);
    assert((pointIndicies=(int*)malloc(nPoints*sizeof(int)))!=NULL,
           "initializeEdgePoints--unable to allocate %i ints for local indicies",
           nPoints);
    
    for(int i=0,currIndex=0;i<numPoints;i++){
      if(ISMID_OC(procLim[0],indicies[0+3*i],procLim[1]) &&
         ISMID_OC(procLim[2],indicies[1+3*i],procLim[3]) &&
         ISMID_OC(procLim[4],indicies[2+3*i],procLim[5])){
        int index=(indicies[0+3*i]-procLim[0])+
        (indicies[1+3*i]-procLim[2])*NX+
        (indicies[2+3*i]-procLim[4])*NXY;
        gridIndicies[currIndex]=index;
        //Set the point index and increment the currIndex.
        pointIndicies[currIndex++]=i;
      }
    }
    
    //Allocate the buffer.
    assert(((*buffer)=(float**)malloc(nPoints*sizeof(float*)))!=NULL,
           "initializeEdgePoints--unable to allocate %i ptrs for buffer",
           nPoints);
    for(int i=0;i<nPoints;i++){
      assert(((*buffer)[i]=
              (float*)malloc((TDBC_SOURCE_BUFFER_SIZE+1)*sizeof(float)))!=NULL,
             "initializeEdgePoints--unable to allocate %i floats for buffer[%i]",
             TDBC_SOURCE_BUFFER_SIZE,i);
      (*buffer)[i][TDBC_SOURCE_BUFFER_SIZE]=0.0;
    }
  }
}

void tdbcSource::initializeUnitPoints(modelDefStruct* modelDef,
                                      int numPointsDim,size_t numPoints){
  DEF_MODEL_SIZE(modelDef);
  //Allocate space for indicies and read.
  int *indicies;
  assert((indicies=(int*)malloc(3*numPoints*sizeof(int)))!=NULL,
         "initializeUnitPoints--unable to allocate %i ints for indicies",
         numPoints);
  int indexVar;
  assert(nc_inq_varid(_infile,"Indicies",&indexVar)==NC_NOERR,
         "initializeUnitPoints--unable to get var Indicies from %s(%i)",
         _infileName,_infile);
  assert(nc_get_var_int(_infile,indexVar,indicies)==NC_NOERR,
         "initializeUnitPoints--unable to read var Indicies from %s(%i,%i)",
         _infileName,_infile,indexVar);
  
  //Get the variables
  assert(nc_inq_varid(_infile,"Vx",&_vxVar)==NC_NOERR,
         "initializeUnitPoints--unable to get var Vx from %s(%i)",
         _infileName,_infile);
  assert(nc_inq_varid(_infile,"Vy",&_vyVar)==NC_NOERR,
         "initializeUnitPoints--unable to get var Vy from %s(%i)",
         _infileName,_infile);
  assert(nc_inq_varid(_infile,"Vz",&_vzVar)==NC_NOERR,
         "initializeUnitPoints--unable to get var Vz from %s(%i)",
         _infileName,_infile);
  
  //Count the number of points that are in bounds for this process.
  _nPoints=0;
  for(int i=0;i<numPoints;i++){
    if((ISMID_OC(procLim[0],indicies[0+3*i],procLim[1]) ||
        ISMID_OC(procLim[0],indicies[0+3*i]-1,procLim[1])) &&
       (ISMID_OC(procLim[2],indicies[1+3*i],procLim[3]) ||
        ISMID_OC(procLim[2],indicies[1+3*i]-1,procLim[3])) &&
       (ISMID_OC(procLim[4],indicies[2+3*i],procLim[5]) ||
        ISMID_OC(procLim[4],indicies[2+3*i]-1,procLim[5]))){
         _nPoints++;
       }
  }
  
  if(!_nPoints){
    //No points are in bounds.
    _gridIndicies=_gridIndiciesMx=_gridIndiciesMy=_gridIndiciesMz=NULL;
    _vxBuffer=_vyBuffer=_vzBuffer=NULL;
    assert(nc_close(_infile)==NC_NOERR,
           "initializeUnitPoints--unable to close %s(%i)",
           _infileName,_infile);
    _infile=-1;
  }else{
    //Allocate local indicies and fill in.
    assert((_gridIndicies=(int*)malloc(_nPoints*sizeof(int)))!=NULL,
           "initializeUnitPoints--unable to allocate %i ints for local indicies",
           _nPoints);
    assert((_gridIndiciesMx=(int*)malloc(_nPoints*sizeof(int)))!=NULL,
           "initializeUnitPoints--unable to allocate %i ints for local indicies",
           _nPoints);
    assert((_gridIndiciesMy=(int*)malloc(_nPoints*sizeof(int)))!=NULL,
           "initializeUnitPoints--unable to allocate %i ints for local indicies",
           _nPoints);
    assert((_gridIndiciesMz=(int*)malloc(_nPoints*sizeof(int)))!=NULL,
           "initializeUnitPoints--unable to allocate %i ints for local indicies",
           _nPoints);
    
    assert((_pointIndicies=(int*)malloc(_nPoints*sizeof(int)))!=NULL,
           "initializeUnitPoints--unable to allocate %i ints for local indicies",
           _nPoints);
    
    for(int i=0,currIndex=0;i<numPoints;i++){
      if((ISMID_OC(procLim[0],indicies[0+3*i],procLim[1]) ||
          ISMID_OC(procLim[0],indicies[0+3*i]-1,procLim[1])) &&
         (ISMID_OC(procLim[2],indicies[1+3*i],procLim[3]) ||
          ISMID_OC(procLim[2],indicies[1+3*i]-1,procLim[3])) &&
         (ISMID_OC(procLim[4],indicies[2+3*i],procLim[5]) ||
          ISMID_OC(procLim[4],indicies[2+3*i]-1,procLim[5]))){
           //Check the center node.
           if(ISMID_OC(procLim[0],indicies[0+3*i],procLim[1]) &&
              ISMID_OC(procLim[2],indicies[1+3*i],procLim[3]) &&
              ISMID_OC(procLim[4],indicies[2+3*i],procLim[5])){
             int index=(indicies[0+3*i]-procLim[0])+
             (indicies[1+3*i]-procLim[2])*NX+
             (indicies[2+3*i]-procLim[4])*NXY;
             _gridIndicies[currIndex]=index;
           }else{
             _gridIndicies[currIndex]=-1;
           }
           
           //Check the X offset node
           if(ISMID_OC(procLim[0],indicies[0+3*i]-1,procLim[1]) &&
              ISMID_OC(procLim[2],indicies[1+3*i],procLim[3]) &&
              ISMID_OC(procLim[4],indicies[2+3*i],procLim[5])){
             int index=(indicies[0+3*i]-procLim[0]-1)+
             (indicies[1+3*i]-procLim[2])*NX+
             (indicies[2+3*i]-procLim[4])*NXY;
             _gridIndiciesMx[currIndex]=index;
           }else{
             _gridIndiciesMx[currIndex]=-1;
           }
           
           //Check the Y offset node.
           if(ISMID_OC(procLim[0],indicies[0+3*i],procLim[1]) &&
              ISMID_OC(procLim[2],indicies[1+3*i]-1,procLim[3]) &&
              ISMID_OC(procLim[4],indicies[2+3*i],procLim[5])){
             int index=(indicies[0+3*i]-procLim[0])+
             (indicies[1+3*i]-procLim[2]-1)*NX+
             (indicies[2+3*i]-procLim[4])*NXY;
             _gridIndiciesMy[currIndex]=index;
           }else{
             _gridIndiciesMy[currIndex]=-1;
           }
           
           //Check the Z offset node.
           if(ISMID_OC(procLim[0],indicies[0+3*i],procLim[1]) &&
              ISMID_OC(procLim[2],indicies[1+3*i],procLim[3]) &&
              ISMID_OC(procLim[4],indicies[2+3*i]-1,procLim[5])){
             int index=(indicies[0+3*i]-procLim[0])+
             (indicies[1+3*i]-procLim[2])*NX+
             (indicies[2+3*i]-procLim[4]-1)*NXY;
             _gridIndiciesMz[currIndex]=index;
           }else{
             _gridIndiciesMz[currIndex]=-1;
           }
           
           //Set the point index and increment the currIndex.
           _pointIndicies[currIndex++]=i;
         }
    }
    
    //Allocate the buffers.
    assert((_vxBuffer=(float**)malloc(_nPoints*sizeof(float*)))!=NULL,
           "initializeUnitPoints--unable to allocate %i ptrs for vx buffer",
           _nPoints);
    assert((_vyBuffer=(float**)malloc(_nPoints*sizeof(float*)))!=NULL,
           "initializeUnitPoints--unable to allocate %i ptrs for vy buffer",
           _nPoints);
    assert((_vzBuffer=(float**)malloc(_nPoints*sizeof(float*)))!=NULL,
           "initializeUnitPoints--unable to allocate %i ptrs for vz buffer",
           _nPoints);
    
    for(int i=0;i<_nPoints;i++){
      assert((_vxBuffer[i]=(float*)malloc(TDBC_SOURCE_BUFFER_SIZE*sizeof(float)))!=NULL,
             "initializeUnitPoints--unable to allocate %i floats for vx buffer[%i]",
             TDBC_SOURCE_BUFFER_SIZE,i);
      assert((_vyBuffer[i]=(float*)malloc(TDBC_SOURCE_BUFFER_SIZE*sizeof(float)))!=NULL,
             "initializeUnitPoints--unable to allocate %i floats for vy buffer[%i]",
             TDBC_SOURCE_BUFFER_SIZE,i);
      assert((_vzBuffer[i]=(float*)malloc(TDBC_SOURCE_BUFFER_SIZE*sizeof(float)))!=NULL,
             "initializeUnitPoints--unable to allocate %i floats for vz buffer[%i]",
             TDBC_SOURCE_BUFFER_SIZE,i);
    }
  }
  readNextBuffer(modelDef,0);
}

void tdbcSource::readNextBuffer(modelDefStruct* modelDef,int index){
  DEF_MODEL_SCALE(modelDef);
  _bufferIndex=index;
  int readSize=MIN(modelDef->NT-_bufferIndex,TDBC_SOURCE_BUFFER_SIZE);
  for(int i=0;i<_nPoints;i++){
    size_t start[2]={_pointIndicies[i],_bufferIndex},count[2]={1,readSize};
    assert(nc_get_vara_float(_infile,_vxVar,start,count,_vxBuffer[i])==NC_NOERR,
           "readNextBuffer--unable to read Vx[%i,%i,%i,%i] from %s(%i)",
           start[0],start[1],count[0],count[1],_infileName,_infile);
    assert(nc_get_vara_float(_infile,_vyVar,start,count,_vyBuffer[i])==NC_NOERR,
           "readNextBuffer--unable to read Vy[%i,%i,%i,%i] from %i",
           start[0],start[1],count[0],count[1],_infileName,_infile);
    assert(nc_get_vara_float(_infile,_vzVar,start,count,_vzBuffer[i])==NC_NOERR,
           "readNextBuffer--unable to read Vz[%i,%i,%i,%i] from %i",
           start[0],start[1],count[0],count[1],_infileName,_infile);
    for(int j=0;j<TDBC_SOURCE_BUFFER_SIZE;j++){
      _vxBuffer[i][j]/=scalarVel;
      _vyBuffer[i][j]/=scalarVel;
      _vzBuffer[i][j]/=scalarVel;
    }
  }
}
void tdbcSource::readNextEdgeBuffer(modelDefStruct* modelDef,int index){
  _bufferIndex=index;
  float scalar=1.0/modelDef->scalarVel;
  int readSize=MIN(modelDef->NT-_bufferIndex,TDBC_SOURCE_BUFFER_SIZE);
  for(int i=0;i<_nVxPoints;i++){
    size_t start[2]={_vxPointIndicies[i],_bufferIndex},count[2]={1,readSize};
    _vxBuffer[i][0]=_vxBuffer[i][TDBC_SOURCE_BUFFER_SIZE];
    float* currBuffer=&(_vxBuffer[i][1]);
    assert(nc_get_vara_float(_infile,_vxVar,start,count,currBuffer)==NC_NOERR,
           "readNextEdgeBuffer--unable to read Vx[%i,%i,%i,%i] from %s(%i)",
           start[0],start[1],count[0],count[1],_infileName,_infile);
    for(int j=0;j<readSize;currBuffer[j++]*=scalar);
  }
  for(int i=0;i<_nVyPoints;i++){
    size_t start[2]={_vyPointIndicies[i],_bufferIndex},count[2]={1,readSize};
    _vyBuffer[i][0]=_vyBuffer[i][TDBC_SOURCE_BUFFER_SIZE];
    float* currBuffer=&(_vyBuffer[i][1]);
    assert(nc_get_vara_float(_infile,_vyVar,start,count,currBuffer)==NC_NOERR,
           "readNextEdgeBuffer--unable to read Vy[%i,%i,%i,%i] from %s(%i)",
           start[0],start[1],count[0],count[1],_infileName,_infile);
    for(int j=0;j<readSize;currBuffer[j++]*=scalar);
  }
  for(int i=0;i<_nVzPoints;i++){
    size_t start[2]={_vzPointIndicies[i],_bufferIndex},count[2]={1,readSize};
    _vzBuffer[i][0]=_vzBuffer[i][TDBC_SOURCE_BUFFER_SIZE];
    float* currBuffer=&(_vzBuffer[i][1]);
    assert(nc_get_vara_float(_infile,_vzVar,start,count,currBuffer)==NC_NOERR,
           "readNextEdgeBuffer--unable to read Vz[%i,%i,%i,%i] from %s(%i)",
           start[0],start[1],count[0],count[1],_infileName,_infile);
    for(int j=0;j<readSize;currBuffer[j++]*=scalar);
  }
  for(int i=0;i<_nPressPoints;i++){
    size_t start[2]={_pressPointIndicies[i],_bufferIndex},count[2]={1,readSize};
    _pressBuffer[i][0]=_pressBuffer[i][TDBC_SOURCE_BUFFER_SIZE];
    float* currBuffer=&(_pressBuffer[i][1]);
    assert(nc_get_vara_float(_infile,_pressVar,start,count,currBuffer)==NC_NOERR,
           "readNextEdgeBuffer--unable to read Pressure[%i,%i,%i,%i] from %s(%i)",
           start[0],start[1],count[0],count[1],_infileName,_infile);
    for(int j=0;j<readSize;currBuffer[j++]/=modelDef->scalarStress);
  }
  _bufferIndex--;
}

///Parameters passed as arguments.
pointSource::pointSource(modelDefStruct* modelDef,
                         float x,float y,float z,float amp,float* data){
  _x=x;
  _y=y;
  _z=z;
  
  _amp=amp;
  
  if(!data){
    _data=NULL;
  }else{
    DEF_MODEL_SIZE(modelDef);
    assert((_data=(float*)malloc((NT+1)*sizeof(float)))!=NULL,
           "pointSource::pointSource--unable to allocate %i floats for source time function",
           NT);
    for(int i=0;i<NT;i++)
      _data[i]=data[i];
    _data[NT]=_data[NT-1];
  }
  setActive(modelDef);
}

///Parameters read from message.
pointSource::pointSource(modelDefStruct* modelDef){
  DEF_MODEL_SIZE(modelDef);
  assert((_data=(float*)malloc((NT+1)*sizeof(float)))!=NULL,
         "pointSource::pointSource--unable to allocate %i floats for source time function",
         NT);
  
  unpackMessage("ffff",&_x,&_y,&_z,&_amp);
  unpackMessage("F",_data,MIN(MAX_TIME_STEPS_SEND,NT));
  _data[NT]=_data[NT-1];
  
  setActive(modelDef);
}

///Parameters read from CDF file.
pointSource::pointSource(modelDefStruct* modelDef,
                         char* filename,int infile,size_t index,
                         int xVar,int yVar,int zVar,
                         int ampVar,int dataVar){
  if(xVar>=0)
    assert(nc_get_var1_float(infile,xVar,&index,&_x)==NC_NOERR,
           "pointSource::pointSource--unable to read X[%i] from %s(%i,%i)",
           (int)index,filename,infile,xVar);
  if(yVar>=0)
    assert(nc_get_var1_float(infile,yVar,&index,&_y)==NC_NOERR,
           "pointSource::pointSource--unable to read Y[%i] from %s(%i,%i)",
           (int)index,filename,infile,yVar);
  if(zVar>=0)
    assert(nc_get_var1_float(infile,zVar,&index,&_z)==NC_NOERR,
           "pointSource::pointSource--unable to read Z[%i] from %s(%i,%i)",
           (int)index,filename,infile,zVar);
  if(ampVar>=0)
    assert(nc_get_var1_float(infile,ampVar,&index,&_amp)==NC_NOERR,
           "pointSource::pointSource--unable to read amp[%i] from %s(%i,%i)",
           (int)index,filename,infile,ampVar);
  
  if(dataVar < 0){
    _data=NULL;
  }else{
    DEF_MODEL_SIZE(modelDef);
    assert((_data=(float*)malloc((NT+1)*sizeof(float)))!=NULL,
           "pointSource::pointSource--unable to allocate %i floats for source time function",
           NT);
    
    int incVar;
    assert(nc_inq_varid(infile,"increments",&incVar)==NC_NOERR,
           "pointSource::pointSource--unable to get increments  var for %s(%i)",
           filename,infile);
    size_t three=3;
    float modelDt;
    assert(nc_get_var1_float(infile,incVar,&three,&modelDt)==NC_NOERR,
           "pointSource::pointSource--unable to read increment var for %s(%i,%i)",
           filename,infile,incVar);
    if(fabs(modelDt-modelDef->dt)<modelDef->dt/100){
      size_t start[2]={index,0},count[2]={1,NT};
      assert(nc_get_vara_float(infile,dataVar,start,count,_data)==NC_NOERR,
             "pointSource::pointSource--unable to read data[%i,%i,%i,%i] from %s(%i,%i)",
             start[0],start[1],count[0],count[1],filename,infile,dataVar);
      _data[NT]=0.0;
    }else{
      //Model dt is different than run dt, need to do some interpolation/extrapolation
      // to the new time vector.
      float dt=modelDef->dt;
      
      int ntDim;
      assert(nc_inq_dimid(infile,"NT",&ntDim)==NC_NOERR,
             "pointSource::pointSource--unable to get NT dim for %s(%i)",
             filename,infile);
      size_t dataNT;
      assert(nc_inq_dimlen(infile,ntDim,&dataNT)==NC_NOERR,
             "pointSource::pointSource--unable to read NT dim for %s(%i,%i)",
             filename,infile,ntDim);
      
      float *data;
      assert((data=(float*)malloc(dataNT*sizeof(float)))!=NULL,
             "pointSource::pointSource--unable to allocate space to read source data");
      
      size_t start[2]={index,0},count[2]={1,dataNT};
      assert(nc_get_vara_float(infile,dataVar,start,count,data)==NC_NOERR,
             "pointSource::pointSource--unable to read data[%i,%i,%i,%i]  from %s(%i,%i)",
             start[0],start[1],count[0],count[1],filename,infile,dataVar);
      for(int j=0;j<NT;j++){
        float currT=dt*(float)j;
        int index=(int)floor(currT/modelDt);
        if(index>=dataNT){
          _data[j]=data[dataNT-1];
        }else if(index<0){
          _data[j]=data[0];
        }else{
          float h1=currT-modelDt*index;
          float h2=modelDt-h1;
          float value=(h2*data[index]+h1*data[index+1])/modelDt;
          _data[j]=value;
        }
      }
      free(data);
    }
  }
  
  setActive(modelDef);
}

///Additional message reading function. If there are a large number of
/// time-steps we use an additional function to read the rest of the data.
int pointSource::readAuxMessage(int NT,int startIndex,
                                int startTag,int iGetMessg){
  if(iGetMessg)
    getMessage(Parent,startTag+startIndex,NULL);
  unpackMessage("F",
                _data+startIndex,
                MIN(MAX_TIME_STEPS_SEND,NT-startIndex));
  _data[NT]=_data[NT-1];
  
  return (NT+1-startIndex)<=MAX_TIME_STEPS_SEND;
}

int pointSource::setActive(modelDefStruct* modelDef){
  DEF_MODEL_SIZE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  
  //Default must be 1 nodes inside model, this may be changed for the traction source.
  if(!ISMID(minX+dx,_x,minX+dx*(NX-2)) ||
     !ISMID(minY+dy,_y,minY+dy*(NY-2)) ||
     !ISMID(minZ+dz,_z,minZ+dz*(NZ-2))){
    return _active=FALSE;
  }
  return _active=TRUE;
}

//Pack the basic parameters, this should be the same format as is read by the constructor.
void pointSource::packSourceMessage(modelDefStruct* modelDef,
                                    int startIteration){
  if(!startIteration){
    packMessage("i",type());
    packMessage("ffff",_x,_y,_z,_amp);
  }
  packMessage("F",_data+startIteration,
              MIN(MAX_TIME_STEPS_SEND,modelDef->NT-startIteration));
}

//Write the basic parameters to a cdf file, this should be the same format as is read by
// the constructor.
void pointSource::writeSource(int NT,float ampScalar,size_t index,
                              char* filename,int outfile,
                              int xVar,int yVar,int zVar,
                              int ampVar,int dataVar,
                              int* extraVars){
  assert(nc_put_var1_float(outfile,xVar,&index,&_x)==NC_NOERR,
         "pointSource::writeSource--unable to write X[%i] to %s(%i,%i)",
         (int)index,filename,outfile,xVar);
  assert(nc_put_var1_float(outfile,yVar,&index,&_y)==NC_NOERR,
         "pointSource::writeSource--unable to write Y[%i] to %s(%i,%i)",
         (int)index,filename,outfile,yVar);
  assert(nc_put_var1_float(outfile,zVar,&index,&_z)==NC_NOERR,
         "pointSource::writeSource--unable to write Z[%i] to %s(%i,%i)",
         (int)index,filename,outfile,zVar);
  
  float tempAmp=_amp*ampScalar;
  assert(nc_put_var1_float(outfile,ampVar,&index,&tempAmp)==NC_NOERR,
         "pointSource::writeSource--unable to write amp[%i] to %s(%i,%i)",
         (int)index,filename,outfile,ampVar);
  
  assert(_data!=NULL,
         "pointSource::writeSource--attempt to write source at %f,%f,%f with uninitialized data",
         _x,_y,_z);
  size_t start[2]={index,0};
  size_t count[2]={1,NT};
  assert(nc_put_vara_float(outfile,dataVar,start,count,_data)==NC_NOERR,
         "pointSource::writeSource--unable to write data[%i,%i,%i,%i] to %s(%i,%i)",
         start[0],start[1],count[0],count[1],filename,outfile,dataVar);
}

//Constructors
///Parameters passed as arguments.
momentSource::momentSource(modelDefStruct* modelDef,
                           float x,float y,float z,
                           float xxS,float yyS,float zzS,
                           float xyS,float xzS,float yzS,
                           float xyA,float xzA,float yzA,
                           float amp,float* data,
                           float* rho,int useMomentRate):pointSource(modelDef,x,y,z,amp,data){
  _useCubicExtrap=FALSE;
  _useMomentRate=useMomentRate;
  
  _xxS=xxS;
  _yyS=yyS;
  _zzS=zzS;
  _xyS=xyS;
  _xzS=xzS;
  _yzS=yzS;
  _xyA=xyA;
  _xzA=xzA;
  _yzA=yzA;
  
  calcInterpCoeffs(modelDef,rho);
  
  //Scale the amplitude for the model
  DEF_MODEL_SCALE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  _amp/=(-scalarStress*dx*dy*dz);
  
  _iso=-1.0/3.0*_amp*(_xxS+_yyS+_zzS);
}
///Parameters read from message.
momentSource::momentSource(modelDefStruct* modelDef,int startTag,float* rho,
                           int useCubicExtrap):
pointSource(modelDef){
  unpackMessage("i fff fff fff",
                &_useMomentRate,&_xxS,&_yyS,&_zzS,
                &_xyS,&_xzS,&_yzS,
                &_xyA,&_xzA,&_yzA);
  
  _useCubicExtrap=useCubicExtrap;
  calcInterpCoeffs(modelDef,rho);
  
  //Do not scale the amplitude, here, it should have been done when source was read or
  // defined in the sending process.
  _iso=-1.0/3.0*_amp*(_xxS+_yyS+_zzS);
  
  //Check to see if additional data needs to be read because of a large
  // number of time-steps.
  DEF_MODEL_SIZE(modelDef);
  for(int j=MAX_TIME_STEPS_SEND;startTag&&j<NT;j+=MAX_TIME_STEPS_SEND)
    readAuxMessage(NT,j,startTag,TRUE);
}

// Parameters read from CDF file.
momentSource::momentSource(modelDefStruct* modelDef,
                           char* filename,int infile,size_t index,
                           int xVar,int yVar,int zVar,
                           int xxSVar,int yySVar,int zzSVar,
                           int xySVar,int xzSVar,int yzSVar,
                           int xyAVar,int xzAVar,int yzAVar,
                           int ampVar,int dataVar,
                           float* rho):
pointSource(modelDef,filename,infile,index,
            xVar,yVar,zVar,ampVar,dataVar){
  _useCubicExtrap=FALSE;
  _useMomentRate=FALSE;
  
  //Apply the correct scaling to the amplitude variable.
  DEF_MODEL_SCALE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  _amp/=(-scalarStress*dx*dy*dz);
  
  //Read the additional parameters.
  assert(nc_get_var1_float(infile,xxSVar,&index,&_xxS)==NC_NOERR,
         "momentSource::source--unable to read xxS[%i] from %s(%i,%i)",
         (int)index,filename,infile,xxSVar);
  assert(nc_get_var1_float(infile,yySVar,&index,&_yyS)==NC_NOERR,
         "momentSource::source--unable to read yyS[%i] from %s(%i,%i)",
         (int)index,filename,infile,yySVar);
  assert(nc_get_var1_float(infile,zzSVar,&index,&_zzS)==NC_NOERR,
         "momentSource::source--unable to read zzS[%i] from %s(%i,%i)",
         (int)index,filename,infile,zzSVar);
  
  assert(nc_get_var1_float(infile,xySVar,&index,&_xyS)==NC_NOERR,
         "momentSource::source--unable to read xyS[%i] from %s(%i,%i)",
         (int)index,filename,infile,xySVar);
  assert(nc_get_var1_float(infile,xzSVar,&index,&_xzS)==NC_NOERR,
         "momentSource::source--unable to read xzS[%i] from %s(%i,%i)",
         (int)index,filename,infile,xzSVar);
  assert(nc_get_var1_float(infile,yzSVar,&index,&_yzS)==NC_NOERR,
         "momentSource::source--unable to read yzS[%i] from %s(%i,%i)",
         (int)index,filename,infile,yzSVar);
  
  assert(nc_get_var1_float(infile,xyAVar,&index,&_xyA)==NC_NOERR,
         "momentSource::source--unable to read xyA[%i] from %s(%i,%i)",
         (int)index,filename,infile,xyAVar);
  assert(nc_get_var1_float(infile,xzAVar,&index,&_xzA)==NC_NOERR,
         "momentSource::source--unable to read xzA[%i] from %s(%i,%i)",
         (int)index,filename,infile,xzAVar);
  assert(nc_get_var1_float(infile,yzAVar,&index,&_yzA)==NC_NOERR,
         "momentSource::source--unable to read yzA[%i] from %s(%i,%i)",
         (int)index,filename,infile,yzAVar);
  
  calcInterpCoeffs(modelDef,rho);
  _iso=-1.0/3.0*_amp*(_xxS+_yyS+_zzS);
}

void momentSource::packSourceMessage(modelDefStruct* modelDef,
                                     int startIteration){
  pointSource::packSourceMessage(modelDef,startIteration);
  if(!startIteration){
    packMessage("ifffffffff",
                _useMomentRate,_xxS,_yyS,_zzS,
                _xyS,_xzS,_yzS,
                _xyA,_xzA,_yzA);
  }
}
void momentSource::writeSource(int NT,float ampScalar,size_t index,
                               char* filename,int outfile,
                               int xVar,int yVar,int zVar,
                               int ampVar,int dataVar,
                               int* extraVars){
  pointSource::writeSource(NT,ampScalar,index,
                           filename,outfile,
                           xVar,yVar,zVar,ampVar,dataVar,
                           extraVars);
  
  int xxSVar=extraVars[0],yySVar=extraVars[1],zzSVar=extraVars[2];
  int xySVar=extraVars[3],xzSVar=extraVars[4],yzSVar=extraVars[5];
  int xyAVar=extraVars[6],xzAVar=extraVars[7],yzAVar=extraVars[8];
  
  assert(nc_put_var1_float(outfile,xxSVar,&index,&_xxS)==NC_NOERR,
         "momentSource::writeSource--unable to write xxS[%i] to %s(%i,%i)",
         (int)index,filename,outfile,xxSVar);
  assert(nc_put_var1_float(outfile,yySVar,&index,&_yyS)==NC_NOERR,
         "momentSource::writeSource--unable to write yyS[%i] to %s(%i,%i)",
         (int)index,filename,outfile,yySVar);
  assert(nc_put_var1_float(outfile,zzSVar,&index,&_zzS)==NC_NOERR,
         "momentSource::writeSource--unable to write zzS[%i] to %s(%i,%i)",
         (int)index,filename,outfile,zzSVar);
  
  assert(nc_put_var1_float(outfile,xySVar,&index,&_xyS)==NC_NOERR,
         "momentSource::writeSource--unable to write xyS[%i] to %s(%i,%i)",
         (int)index,filename,outfile,xySVar);
  assert(nc_put_var1_float(outfile,xzSVar,&index,&_xzS)==NC_NOERR,
         "momentSource::writeSource--unable to write xzS[%i] to %s(%i,%i)",
         (int)index,filename,outfile,xzSVar);
  assert(nc_put_var1_float(outfile,yzSVar,&index,&_yzS)==NC_NOERR,
         "momentSource::writeSource--unable to write yzS[%i] to %s(%i,%i)",
         (int)index,filename,outfile,yzSVar);
  
  assert(nc_put_var1_float(outfile,xyAVar,&index,&_xyA)==NC_NOERR,
         "momentSource::writeSource--unable to write xyA[%i] to %s(%i,%i)",
         (int)index,filename,outfile,xyAVar);
  assert(nc_put_var1_float(outfile,xzAVar,&index,&_xzA)==NC_NOERR,
         "momentSource::writeSource--unable to write xzA[%i] to %s(%i,%i)",
         (int)index,filename,outfile,xzAVar);
  assert(nc_put_var1_float(outfile,yzAVar,&index,&_yzA)==NC_NOERR,
         "momentSource::writeSource--unable to write yzA[%i] to %s(%i,%i)",
         (int)index,filename,outfile,yzAVar);
}

//In the elastic/anelastic case this source is applied to velocities and stresses.
int momentSource::applyVel(modelDefStruct* modelDef,int iteration,
                           float *vx,float *vy,float *vz,
                           float* rho){
  static double pMomRate = 0.;
  if(!_isAnisotropic) return FALSE;
  
  float sx=modelDef->dt/modelDef->dx*modelDef->scalarSpeed;
  float sy=modelDef->dt/modelDef->dy*modelDef->scalarSpeed;
  float sz=modelDef->dt/modelDef->dz*modelDef->scalarSpeed;
  DEF_MODEL_SIZE(modelDef);
  float factr;
  //for moment rate sources, the STF must be integrated.  Using cumulative Simpson's Rule from Blake 1971
  if(_useMomentRate) {  //dt has already been included in amp in addSources in sourceNetwork
    if(iteration==0) {  //assume 0 for times < 0
      factr=_amp/3.*_data[iteration]*1.25;
      pMomRate = factr;
    } else if(iteration>NT-2) {  //assume all t>tmax equal to final point in STF
      if(iteration==NT-1) {
        if(iteration%2)
          factr=pMomRate+_amp/3.*(1.25*_data[iteration-1]+1.75*_data[iteration]);
        else
          factr=pMomRate+_amp/3.*(-.25*_data[iteration-2]+2.*_data[iteration-1]+1.25*_data[iteration]);
      } else if(iteration==NT) {
        if(iteration%2)
          factr=pMomRate+_amp*_data[iteration-1];
        else
          factr=pMomRate+_amp/3.*(-.25*_data[iteration-2]+3.25*_data[iteration-1]);
      } else {
        factr=pMomRate+_amp*_data[NT-1];
      }
      pMomRate = factr;
    } else {
      if(iteration%2)
        factr=pMomRate+_amp/3.*(1.25*_data[iteration-1]+2.*_data[iteration]-.25*_data[iteration+1]);
      else
        factr=pMomRate+_amp/3.*(-.25*_data[iteration-2]+2.*_data[iteration-1]+1.25*_data[iteration]);
      pMomRate = factr;
    }
  } else
    factr=_amp*_data[iteration];
  
  //Check for near 0 value.
  if(fabs(factr)<1e-30) return FALSE;
  
  //Update vx and vy components.
  if(fabs(_xyA)>1e-30){
    for(int n=0;n<8;n++){
      float value=factr*_xyA*_xy_vInterp.coeff[n];
      int i=_xy_vInterp.iptr[n];
      int j=_xy_vInterp.jptr[n];
      int k=_xy_vInterp.kptr[n];
      int index = i+j*NX+k*NXY;
      
      if(ISMID(0,i,NX-1) &&
         ISMID(0,k,NZ-1)){
        if(ISMID(0,j  ,NY-1)){
          VX(i,j,  k)+=  sy*value/(rho[index]+rho[index+1]);
        }
        if(ISMID(0,j+1,NY-1)){
          VX(i,j+1,k)+= -sy*value/(rho[index+NX]+rho[index+NX+1]);
        }
      }
      if(ISMID(0,j,NY-1) &&
         ISMID(0,k,NZ-1)){
        if(ISMID(0,i  ,NX-1)){
          VY(i,  j,k)+= -sx*value/(rho[index]+rho[index+NX]);
        }
        if(ISMID(0,i+1,NX-1)){
          VY(i+1,j,k)+=  sx*value/(rho[index+1]+rho[index+1+NX]);
        }
      }
    }
  }
  //Update vy and vz components.
  if(fabs(_yzA)>1e-30){
    for(int n=0;n<8;n++){
      float value=factr*_yzA*_yz_vInterp.coeff[n];
      int i=_yz_vInterp.iptr[n];
      int j=_yz_vInterp.jptr[n];
      int k=_yz_vInterp.kptr[n];
      int index = i+j*NX+k*NXY;
      
      if(ISMID(0,i,NX-1) &&
         ISMID(0,j,NY-1)){
        if(ISMID(0,k  ,NZ-1))
          VY(i,j,k  )+=  sz*value/(rho[index]+rho[index+NX]);
        if(ISMID(0,k+1,NZ-1))
          VY(i,j,k+1)+= -sz*value/(rho[index+NXY]+rho[index+NXY+NX]);
      }
      if(ISMID(0,i,NX-1) &&
         ISMID(0,k,NZ-1)){
        if(ISMID(0,j,NY-1))
          VZ(i,j,  k)+= -sy*value/(rho[index]+rho[index+NXY]);
        if(ISMID(0,j+1,NY-1))
          VZ(i,j+1,k)+=  sy*value/(rho[index+NX]+rho[index+NX+NXY]);
      }
    }
  }
  //Update vx and vz components.
  if(fabs(_xzA)>1e-30){
    for(int n=0;n<8;n++){
      float value=factr*_xzA*_xz_vInterp.coeff[n];
      int i=_xz_vInterp.iptr[n];
      int j=_xz_vInterp.jptr[n];
      int k=_xz_vInterp.kptr[n];
      int index = i+j*NX+k*NXY;
      
      if(ISMID(0,i,NX-1) &&
         ISMID(0,j,NY-1)){
        if(ISMID(0,k  ,NZ-1))
          VX(i,j,k  )+=  sz*value/(rho[index]+rho[index+1]);
        if(ISMID(0,k+1,NZ-1))
          VX(i,j,k+1)+= -sz*value/(rho[index+NXY]+rho[index+NXY+1]);
      }
      if(ISMID(0,j,NY-1) &&
         ISMID(0,k,NZ-1)){
        if(ISMID(0,i  ,NX-1))
          VZ(i,  j,k)+= -sx*value/(rho[index]+rho[index+NXY]);
        if(ISMID(0,i+1,NX-1))
          VZ(i+1,j,k)+=  sx*value/(rho[index+1]+rho[index+1+NXY]);
      }
    }
  }
  return iteration;
}

int momentSource::applyStress(modelDefStruct* modelDef,int iteration,
                              float *xx,float* yy,float *zz,
                              float *xy,float* xz,float *yz){
  if(!_isIsotropic) return FALSE;
  
  float value;
  if(_useMomentRate) {  //dt has already been included in amp in addSources in sourceNetwork
    if(iteration)
      value = _amp*0.5*(_data[iteration]+_data[iteration-1]);
    else
      value = _amp*0.5*_data[iteration];
  } else {
    if(iteration){
      value=_amp*(_data[iteration]-_data[iteration-1]);
    }else{
      value=_amp*_data[iteration];
    }
  }
  
  //Check for near 0 value.
  if(fabs(value)<1e-30) return FALSE;
  
  //xx-component.
  trilinExtrap(modelDef,
               xx,&_xxInterp,_xxS*value);
  
  //yy-component.
  trilinExtrap(modelDef,
               yy,&_xxInterp,_yyS*value);
  
  //zz-component.
  trilinExtrap(modelDef,
               zz,&_xxInterp,_zzS*value);
  
  //xy-component.
  trilinExtrap(modelDef,
               xy,&_xyInterp,_xyS*value);
  
  //yz-component.
  trilinExtrap(modelDef,
               yz,&_yzInterp,_yzS*value);
  
  //xz-component.
  trilinExtrap(modelDef,
               xz,&_xzInterp,_xzS*value);
  
  return iteration;
}

int momentSource::applyPressure(modelDefStruct* modelDef,int iteration,
                                float *pressure){
  if(!_isIsotropic) return FALSE;
  
  float value;
  if(_useMomentRate) {
    if(iteration){
      value=_iso*0.5*(_data[iteration]+_data[iteration-1]);
    }else{
      value=_iso*0.5*_data[iteration];
    }
  } else {
    if(iteration){
      value=_iso*(_data[iteration]-_data[iteration-1]);
    }else{
      value=_iso*_data[iteration];
    }
  }
  
  //Check for near 0 value.
  if(fabs(value)<1e-30) return FALSE;
  
  //Only need to do the pressure update with the sum.
  if(_useCubicExtrap){
    triCubicExtrap(modelDef,&_cubicXxInterp,
                   pressure,value);
  }else{
    trilinExtrap(modelDef,
                 pressure,&_xxInterp,value);
  }
  return iteration;
}

int momentSource::applyPressureMM(modelDefStruct* modelDef,int iteration,
                                  float *pressure){
  if(!_isIsotropic) return FALSE;
  
  float value;
  if(_useMomentRate) {
    if(iteration>0){
      value=_iso*_data[iteration-1];
    }else{
      value=_iso*0.5*_data[iteration];
    }
  } else {
    if(iteration>1){
      value=_iso*(_data[iteration]-_data[iteration-2]);
    }else{
      value=_iso*_data[iteration];
    }
  }
  
  //Check for near 0 value.
  if(fabs(value)<1e-30) return FALSE;
  
  //Only need to do the pressure update with the sum.
  if(_useCubicExtrap){
    triCubicExtrap(modelDef,&_cubicXxInterp,
                   pressure,value);
  }else{
    trilinExtrap(modelDef,
                 pressure,&_xxInterp,value);
  }
  return iteration;
}

void momentSource::calcInterpCoeffs(modelDefStruct* modelDef,float* rho){
  //check if source has isotropic and anisotropic components.
  _isAnisotropic=(fabs(_xyA)+fabs(_xzA)+fabs(_yzA))>0.01;
  _isIsotropic=(fabs(_xxS)+fabs(_yyS)+fabs(_zzS)+fabs(_xyS)+fabs(_xzS)+fabs(_yzS))>0.01;
  
  DEF_MODEL_LIMITS(modelDef);
  
  //some useful numbers
  float xmin_dx=minX+0.5*dx;
  float ymin_dy=minY+0.5*dy;
  float zmin_dz=minZ+0.5*dz;
  
  //Anti-symetric Components; extrapolate to the velocity vector
  trilinCoeff(_x,_y,_z,
              xmin_dx,dx,ymin_dy,dy,minZ,dz,
              &(_xy_vInterp));
  trilinCoeff(_x,_y,_z,
              minX,dx,ymin_dy,dy,zmin_dz,dz,
              &(_yz_vInterp));
  trilinCoeff(_x,_y,_z,
              xmin_dx,dx,minY,dy,zmin_dz,dz,
              &(_xz_vInterp));
  
  //Symetric Components; extrapolate to the stress tensor
  // diagonal component coeffs
  if(_useCubicExtrap)
    triCubicCoeff(_x,_y,_z,
                  minX,dx,minY,dy,minZ,dz,
                  &_cubicXxInterp);
  trilinCoeff(_x,_y,_z,
              minX,dx,minY,dy,minZ,dz,
              &(_xxInterp));
  
  //Off diagional components.
  trilinCoeff(_x,_y,_z,
              xmin_dx,dx,ymin_dy,dy,minZ,dz,
              &(_xyInterp));
  trilinCoeff(_x,_y,_z,
              minX,dx,ymin_dy,dy,zmin_dz,dz,
              &(_yzInterp));
  trilinCoeff(_x,_y,_z,
              xmin_dx,dx,minY,dy,zmin_dz,dz,
              &(_xzInterp));
}

///Values passed.
forceSource::forceSource(modelDefStruct* modelDef,
                         float x,float y,float z,
                         float ax,float ay,float az,
                         float amp,float* data,
                         float* rho):pointSource(modelDef,x,y,z,amp,data){
  _ax=ax;
  _ay=ay;
  _az=az;
  
  calcInterpCoeffs(modelDef,rho);
  
  //Also scale the source amplitude to the non-dimensionalizing factors and the grid size.
  DEF_MODEL_LIMITS(modelDef);
  DEF_MODEL_SCALE(modelDef);
  _amp*= dt/(scalarVel*scalarDen*dx*dy*dz);
}

///Parameters read from message.
forceSource::forceSource(modelDefStruct* modelDef,int startTag,float* rho):
pointSource(modelDef){
  unpackMessage("fff",&_ax,&_ay,&_az);
  
  calcInterpCoeffs(modelDef,rho);
  
  //Do not scale the amplitude, here, it should have been done when source was read or
  // defined in the sending process.
  
  //Check to see if additional data needs to be read because of a large
  // number of time-steps.
  DEF_MODEL_SIZE(modelDef);
  for(int j=MAX_TIME_STEPS_SEND;startTag&&j<NT;j+=MAX_TIME_STEPS_SEND)
    readAuxMessage(NT,j,startTag,TRUE);
}

///Parameters read from CDF file.
forceSource::forceSource(modelDefStruct* modelDef,
                         char* filename,int infile,size_t index,
                         int xVar,int yVar,int zVar,
                         int axVar,int ayVar,int azVar,
                         int ampVar,int dataVar,
                         float* rho):
pointSource(modelDef,filename,infile,index,
            xVar,yVar,zVar,ampVar,dataVar){
  //Apply the correct scaling to the amplitude variable.
  DEF_MODEL_SCALE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  _amp*= dt/(scalarVel*scalarDen*dx*dy*dz); //Apply the correct scaling.
  
  //Read the additional parameters.
  assert(nc_get_var1_float(infile,axVar,&index,&_ax)==NC_NOERR,
         "forceSource::source--unable to read aX[%i] from %s(%i,%i)",
         (int)index,filename,infile,axVar);
  assert(nc_get_var1_float(infile,ayVar,&index,&_ay)==NC_NOERR,
         "forceSource::source--unable to read aY[%i] from %s(%i,%i)",
         (int)index,filename,infile,ayVar);
  assert(nc_get_var1_float(infile,azVar,&index,&_az)==NC_NOERR,
         "forceSource::source--unable to read aZ[%i] from %s(%i,%i)",
         (int)index,filename,infile,azVar);
  
  calcInterpCoeffs(modelDef,rho);
}

void forceSource::packSourceMessage(modelDefStruct* modelDef,
                                    int startIteration){
  pointSource::packSourceMessage(modelDef,startIteration);
  if(!startIteration)
    packMessage("fff",_ax,_ay,_az);
}
void forceSource::writeSource(int NT,float ampScalar,size_t index,
                              char* filename,int outfile,
                              int xVar,int yVar,int zVar,
                              int ampVar,int dataVar,
                              int* extraVars){
  pointSource::writeSource(NT,ampScalar,index,
                           filename,outfile,
                           xVar,yVar,zVar,ampVar,dataVar,
                           extraVars);
  
  int axVar=extraVars[0],ayVar=extraVars[1],azVar=extraVars[2];
  assert(nc_put_var1_float(outfile,axVar,&index,&_ax)==NC_NOERR,
         "forceSource::writeSource--unable to write ax[%i] to %s(%i,%i)",
         (int)index,filename,outfile,axVar);
  assert(nc_put_var1_float(outfile,ayVar,&index,&_ay)==NC_NOERR,
         "forceSource::writeSource--unable to write ay[%i] to %s(%i,%i)",
         (int)index,filename,outfile,ayVar);
  assert(nc_put_var1_float(outfile,azVar,&index,&_az)==NC_NOERR,
         "forceSource::writeSource--unable to write az[%i] to %s(%i,%i)",
         (int)index,filename,outfile,azVar);
}

//Need to apply this source to velocities only.
int forceSource::applyVel(modelDefStruct* modelDef,int iteration,
                          float *vx,float *vy,float *vz,
                          float* rho){
  float factr=_amp*_data[iteration];  //rhox, y and z are actually buoyancies and are applied below
  //Check for near 0 value.
  if(fabs(factr)<1e-30) return FALSE;
  
  
  DEF_MODEL_SIZE(modelDef);
  //Extrapolate components of the body force vector to the grid points
  // that store the components of the particle velocity vector.
  if(fabs(_ax)>1e-30) {
    for(int n=0;n<8;++n) {
      float value = factr*_ax*_vxInterp.coeff[n];
      int i=_vxInterp.iptr[n];
      int j=_vxInterp.jptr[n];
      int k=_vxInterp.kptr[n];
      if(ISMID(0,i,NX-1) && ISMID(0,j,NY-1) && ISMID(0,k,NZ-1) && _rhox[n]<1e12)
        VX(i,j,k) += value*_rhox[n];
    }
  }
  
  if(fabs(_ay)>1e-30) {
    for(int n=0;n<8;++n) {
      float value = factr*_ay*_vyInterp.coeff[n];
      int i=_vyInterp.iptr[n];
      int j=_vyInterp.jptr[n];
      int k=_vyInterp.kptr[n];
      if(ISMID(0,i,NX-1) && ISMID(0,j,NY-1) && ISMID(0,k,NZ-1) && _rhoy[n]<1e12)
        VY(i,j,k) += value*_rhoy[n];
    }
  }
  
  if(fabs(_az)>1e-30) {
    for(int n=0;n<8;++n) {
      float value = factr*_az*_vzInterp.coeff[n];
      int i=_vzInterp.iptr[n];
      int j=_vzInterp.jptr[n];
      int k=_vzInterp.kptr[n];
      if(ISMID(0,i,NX-1) && ISMID(0,j,NY-1) && ISMID(0,k,NZ-1) && _rhoz[n]<1e12)
        VZ(i,j,k) += value*_rhoz[n];
    }
  }
  
  return iteration;
}
int forceSource::applyVelMM(modelDefStruct* modelDef,int iteration,
                            float *vx,float *vy,float *vz,
                            float* rho){
  float factr;
  if(iteration>0)
    factr=_amp*(_data[iteration]+_data[iteration-1]);
  else
    factr=_amp*_data[iteration];
  //Check for near 0 value.
  if(fabs(factr)<1e-30) return FALSE;
  
  
  DEF_MODEL_SIZE(modelDef);
  //Extrapolate components of the body force vector to the grid points
  // that store the components of the particle velocity vector.
  if(fabs(_ax)>1e-30) {
    for(int n=0;n<8;++n) {
      float value = factr*_ax*_vxInterp.coeff[n];
      int i=_vxInterp.iptr[n];
      int j=_vxInterp.jptr[n];
      int k=_vxInterp.kptr[n];
      if(ISMID(0,i,NX-1) && ISMID(0,j,NY-1) && ISMID(0,k,NZ-1) && _rhox[n]<1e12)
        VX(i,j,k) += value*_rhox[n];
    }
  }
  
  if(fabs(_ay)>1e-30) {
    for(int n=0;n<8;++n) {
      float value = factr*_ay*_vyInterp.coeff[n];
      int i=_vyInterp.iptr[n];
      int j=_vyInterp.jptr[n];
      int k=_vyInterp.kptr[n];
      if(ISMID(0,i,NX-1) && ISMID(0,j,NY-1) && ISMID(0,k,NZ-1) && _rhoy[n]<1e12)
        VY(i,j,k) += value*_rhoy[n];
    }
  }
  
  if(fabs(_az)>1e-30) {
    for(int n=0;n<8;++n) {
      float value = factr*_az*_vzInterp.coeff[n];
      int i=_vzInterp.iptr[n];
      int j=_vzInterp.jptr[n];
      int k=_vzInterp.kptr[n];
      if(ISMID(0,i,NX-1) && ISMID(0,j,NY-1) && ISMID(0,k,NZ-1) && _rhoz[n]<1e12)
        VZ(i,j,k) += value*_rhoz[n];
    }
  }
  
  return iteration;
}

void forceSource::calcInterpCoeffs(modelDefStruct* modelDef,float* rho){
  DEF_MODEL_SIZE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  
  //some useful numbers
  float xmin_dx=minX+0.5*dx;
  float ymin_dy=minY+0.5*dy;
  float zmin_dz=minZ+0.5*dz;
  
  //Extrapolate x-component of the body force vector to the grid points
  // that store the x-component of the particle velocity vector.
  trilinCoeff(_x,_y,_z,
              xmin_dx,dx,minY,dy,minZ,dz,
              &(_vxInterp));
  
  //Extrapolate y-component of the body force vector to the grid points
  // that store the y-component of the particle velocity vector.
  trilinCoeff(_x,_y,_z,
              minX,dx,ymin_dy,dy,minZ,dz,
              &(_vyInterp));
  
  //Extrapolate z-component of the body force vector to the grid points
  // that store the z-component of the particle velocity vector.
  trilinCoeff(_x,_y,_z,
              minX,dx,minY,dy,zmin_dz,dz,
              &(_vzInterp));
  
  //interpolate the bouyency value for this source location
  if(!rho){
    for(int ii=0;ii<8;++ii) {
      _rhox[ii] = _rhoy[ii] = _rhoz[ii] = 1e13;
    }
  }else{
    for(int ii=0;ii<8;++ii) {
      int i=_vxInterp.iptr[ii];
      int j=_vxInterp.jptr[ii];
      int k=_vxInterp.kptr[ii];
      if(ISMID(0,i,NX-2) && ISMID(0,j,NY-1) && ISMID(0,k,NZ-1)) _rhox[ii] = 1.0/(rho[i+j*NX+k*NXY]+rho[i+1+j*NX+k*NXY]);
      else _rhox[ii] = 1e13;
      i=_vyInterp.iptr[ii];
      j=_vyInterp.jptr[ii];
      k=_vyInterp.kptr[ii];
      if(ISMID(0,i,NX-1) && ISMID(0,j,NY-2) && ISMID(0,k,NZ-1)) _rhoy[ii] = 1.0/(rho[i+j*NX+k*NXY]+rho[i+(j+1)*NX+k*NXY]);
      else _rhoy[ii] = 1e13;
      i=_vzInterp.iptr[ii];
      j=_vzInterp.jptr[ii];
      k=_vzInterp.kptr[ii];
      if(ISMID(0,i,NX-1) && ISMID(0,j,NY-1) && ISMID(0,k,NZ-2)) _rhoz[ii] = 1.0/(rho[i+j*NX+k*NXY]+rho[i+j*NX+(k+1)*NXY]);
      else _rhoz[ii] = 1e13;
    }
  }
}

//Constructors
// Parameters passed as arguments.
tractionSource::tractionSource(modelDefStruct* modelDef,
                               float x,float y,float z,
                               float ax,float ay,float az,
                               float amp,float* data,
                               float rho):pointSource(modelDef,x,y,z,amp,data){
  _ax=ax;
  _ay=ay;
  _az=az;
  
  calcInterpCoeffs(modelDef,rho);
  
  //Also scale the source amplitude to the non-dimensionalizing factors and the grid size.
  DEF_MODEL_LIMITS(modelDef);
  DEF_MODEL_SCALE(modelDef);
  _amp*= dt/(scalarVel*scalarDen*dx*dy*dz);
}
tractionSource::tractionSource(modelDefStruct* modelDef,
                               float x,float y,float z,
                               float ax,float ay,float az,
                               float amp,float* data,
                               float* rho):pointSource(modelDef,x,y,z,amp,data){
  _ax=ax;
  _ay=ay;
  _az=az;
  
  calcInterpCoeffs(modelDef,rho);
  
  //Also scale the source amplitude to the non-dimensionalizing factors and the grid size.
  DEF_MODEL_LIMITS(modelDef);
  DEF_MODEL_SCALE(modelDef);
  _amp*= dt/(scalarVel*scalarDen*dx*dy*dz);
}

// Parameters read from message.
tractionSource::tractionSource(modelDefStruct* modelDef,int startTag,float* rho):
pointSource(modelDef){
  unpackMessage("fff",&_ax,&_ay,&_az);
  
  calcInterpCoeffs(modelDef,rho);
  
  //Do not scale the amplitude, here, it should have been done when source was read or
  // defined in the sending process.
  
  //Check to see if additional data needs to be read because of a large
  // number of time-steps.
  DEF_MODEL_SIZE(modelDef);
  for(int j=MAX_TIME_STEPS_SEND;startTag&&j<NT;j+=MAX_TIME_STEPS_SEND)
    readAuxMessage(NT,j,startTag,TRUE);
}

// Parameters read from CDF file.
tractionSource::tractionSource(modelDefStruct* modelDef,
                               char* filename,int infile,size_t index,
                               int xVar,int yVar,int zVar,
                               int axVar,int ayVar,int azVar,
                               int ampVar,int dataVar,
                               float* rho):
pointSource(modelDef,filename,infile,index,
            xVar,yVar,zVar,ampVar,dataVar){
  //Apply the correct scaling to the amplitude variable.
  DEF_MODEL_SCALE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  _amp*= dt/(scalarVel*scalarDen*dx*dy*dz); //Apply the correct scaling.
  
  //Read the additional parameters.
  assert(nc_get_var1_float(infile,axVar,&index,&_ax)==NC_NOERR,
         "tractionSource::source--unable to read aX[%i] from %s(%i,%i)",
         (int)index,filename,infile,axVar);
  assert(nc_get_var1_float(infile,ayVar,&index,&_ay)==NC_NOERR,
         "tractionSource::source--unable to read aY[%i] from %s(%i,%i)",
         (int)index,filename,infile,ayVar);
  assert(nc_get_var1_float(infile,azVar,&index,&_az)==NC_NOERR,
         "tractionSource::source--unable to read aZ[%i] from %s(%i,%i)",
         (int)index,filename,infile,azVar);
  
  calcInterpCoeffs(modelDef,rho);
}

void tractionSource::packSourceMessage(modelDefStruct* modelDef,
                                       int startIteration){
  pointSource::packSourceMessage(modelDef,startIteration);
  if(!startIteration)
    packMessage("fff",_ax,_ay,_az);
}
void tractionSource::writeSource(int NT,float ampScalar,size_t index,
                                 char* filename,int outfile,
                                 int xVar,int yVar,int zVar,
                                 int ampVar,int dataVar,
                                 int* extraVars){
  pointSource::writeSource(NT,ampScalar,index,
                           filename,outfile,
                           xVar,yVar,zVar,ampVar,dataVar,
                           extraVars);
  
  int axVar=extraVars[0],ayVar=extraVars[1],azVar=extraVars[2];
  assert(nc_put_var1_float(outfile,axVar,&index,&_ax)==NC_NOERR,
         "tractionSource::writeSource--unable to write ax[%i] to %s(%i,%i)",
         (int)index,filename,outfile,axVar);
  assert(nc_put_var1_float(outfile,ayVar,&index,&_ay)==NC_NOERR,
         "tractionSource::writeSource--unable to write ay[%i] to %s(%i,%i)",
         (int)index,filename,outfile,ayVar);
  assert(nc_put_var1_float(outfile,azVar,&index,&_az)==NC_NOERR,
         "tractionSource::writeSource--unable to write az[%i] to %s(%i,%i)",
         (int)index,filename,outfile,azVar);
}

//Need to apply this source to velocities only.
int tractionSource::applyVel(modelDefStruct* modelDef,int iteration,
                             float *vx,float *vy,float *vz,
                             float* rho){
  //     float value=_amp*_data[iteration]*2.0*_boy;
  float value=-_amp*(_data[iteration]-_data[iteration-1])/(2.0*_rho);
  //Check for near 0 value.
  if(fabs(value)<1e-30) return FALSE;
  
  
  //Extrapolate components of the body force vector to the grid points
  // that store the components of the particle velocity vector.
  trilinExtrap(modelDef,
               vx,&_vxInterp,
               _ax*value);
  
  trilinExtrap(modelDef,
               vy,&_vyInterp,
               _ay*value);
  
  trilinExtrap(modelDef,
               vz,&_vzInterp,
               _az*value);
  
  return iteration;
}

void tractionSource::calcInterpCoeffs(modelDefStruct* modelDef,float* rho){
  DEF_MODEL_LIMITS(modelDef);
  
  //some useful numbers
  float xmin_dx=minX+0.5*dx;
  float ymin_dy=minY+0.5*dy;
  float zmin_dz=minZ+0.5*dz;
  
  //first interpolate the boyency value for this source location
  if(!rho){
    _rho=-1;
  }else{
    interpStruct interpCoeff;
    trilinCoeff(_x,_y,_z,
                minX,dx,minY,dy,minZ,dz,
                &interpCoeff);
    _rho=trilinInterp(modelDef,&interpCoeff,rho);
  }
  
  //Extrapolate x-component of the body force vector to the grid points
  // that store the x-component of the particle velocity vector.
  trilinCoeff(_x,_y,minZ,
              xmin_dx,dx,minY,dy,minZ,dz,
              &(_vxInterp));
  
  //Extrapolate y-component of the body force vector to the grid points
  // that store the y-component of the particle velocity vector.
  trilinCoeff(_x,_y,minZ,
              minX,dx,ymin_dy,dy,minZ,dz,
              &(_vyInterp));
  
  //Extrapolate z-component of the body force vector to the grid points
  // that store the z-component of the particle velocity vector.
  trilinCoeff(_x,_y,zmin_dz,
              minX,dx,minY,dy,zmin_dz,dz,
              &(_vzInterp));
}

// Sources passed in message.
sourceNetwork::sourceNetwork(modelDefStruct* modelDef,float *rho,
                             int useCubicExtrap){
  _useModelSources=TRUE;
  _useMomentRate = false;
  _bandCheckPCent=DEFAULT_CHECK_BANDWIDTH_PCENT;
  _wavelet=NULL;
  _nfs=_nms=_nts=_ntdbc=0;
  //Read the number of sources.
  int numSources;
  unpackMessage("i",&numSources);
  
  setMessageBuffer(numSources*((20+modelDef->NT)*sizeof(float)));
  
  _sources=new sourceArray(0,numSources);
  for(int i=0;i<numSources;i++){
    //New 8/18/05: send/receive each source in a seperate message
    getMessage(Parent,MESSAGE_SET_SOURCES+1+i,NULL);
    
    //For each source, read the type and then call the constructor corresponding to
    // that source type.
    int type;
    unpackMessage("i",&type);
    switch(type){
      case FORCE_SOURCE:
        _sources->Add(new forceSource(modelDef,MESSAGE_SET_SOURCES+1+i,rho));
        _nfs++;
        break;
      case MOMENT_SOURCE:
        _sources->Add(new momentSource(modelDef,
                                       MESSAGE_SET_SOURCES+1+i,rho,
                                       useCubicExtrap));
        _nms++;
        break;
      case TRACTION_SOURCE:
        _sources->Add(new tractionSource(modelDef,
                                         MESSAGE_SET_SOURCES+1+i,rho));
        _nts++;
        break;
        
      case TDBC_SOURCE:
        _sources->Add(new tdbcSource(modelDef));
        _ntdbc++;
        break;
        
        //Traction source not yet implemented
      default:
        assert(FALSE,"sourceNetwork--unable create new source[%i] of type %i",
               i,type);
    }
  }
}

// Sources read from cdf file
sourceNetwork::sourceNetwork(modelDefStruct* modelDef,char* filename,
                             int useModelSources,float *rho){
  _useModelSources=useModelSources;
  _bandCheckPCent=DEFAULT_CHECK_BANDWIDTH_PCENT;
  _nfs=_nms=_nts=_ntdbc=0;
  _wavelet=NULL;
  _useMomentRate = false;
  _sources=new sourceArray;
  
  if(_useModelSources){
    //Open the file for input.
    int infile=openCDFFile(filename,FALSE,NC_NOWRITE);
    
    readForceSources(modelDef,filename,infile,rho);
    readMomentSources(modelDef,filename,infile,rho);
    readTractionSources(modelDef,filename,infile,rho);
    assert(nc_close(infile)==NC_NOERR,
           "sourceNetwork--unable to close file %s",
           filename);
  }
}

int sourceNetwork::combine(sourceNetwork* extraSources){
  if(extraSources){
    for(int i=0;i<extraSources->size();i++){
      _sources->Add((*extraSources->_sources)[i]);
      switch((*extraSources->_sources)[i]->type()){
        case FORCE_SOURCE:
          _nfs++;
          break;
        case MOMENT_SOURCE:
          _nms++;
          break;
        case TDBC_SOURCE:
          _ntdbc++;
          break;
        case TRACTION_SOURCE:
          _nts++;
          break;
          
        default:
          assert(FALSE,"combine--extraSources[%i], unknown or unsupported type %i",
                 i,(*extraSources->_sources)[i]->type());
      }
    }
  }
  return size();
}

int sourceNetwork::sendSources(modelDefStruct* modelDef,
                               int target){
  DEF_MODEL_SIZE(modelDef);
  
  //Prepare to send info to the slave processes
  initSend();
  packMessage("i",size());
  ((target!=AllProcesses)?iSendMessage:sendMessage)(target,MESSAGE_SET_SOURCES,NULL);
  
  for(int i=0;i<size();i++){
    for(int j=0,count=0;j<NT;count++,j+=MAX_TIME_STEPS_SEND){
      initSend();
      (*_sources)[i]->packSourceMessage(modelDef,j);
      int tag=MESSAGE_SET_SOURCES+1+i+j;
      ((target!=AllProcesses)?iSendMessage:sendMessage)
      (target,tag,NULL);
    }
  }
  
  //receive response from slaves
  if(target==AllProcesses){
    for(int i=0;i<NumProcs;i++)
      getMessage(AllProcesses,MESSAGE_SET_SOURCES,NULL);
  }
  
  return size();
}

int sourceNetwork::applyVel(modelDefStruct* modelDef,int iteration,
                            float *vx,float *vy,float *vz,
                            float* rho){
  int numUsed=0;
  for(int i=0;i<size();i++){
    source* curr=(*_sources)[i];
    numUsed+=curr->applyVel(modelDef,iteration,vx,vy,vz,rho);
  }
  return numUsed;
}
int sourceNetwork::applyStress(modelDefStruct* modelDef,int iteration,
                               float *xx,float* yy,float *zz,
                               float *xy,float* xz,float *yz){
  int numUsed=0;
  for(int i=0;i<size();i++){
    source* curr=(*_sources)[i];
    numUsed+=curr->applyStress(modelDef,iteration,xx,yy,zz,xy,xz,yz);
  }
  return numUsed;
}

int sourceNetwork::applyAcousticVel(modelDefStruct* modelDef,int iteration,
                                    float *vx,float *vy,float *vz,
                                    float* rho){
  int numUsed=0;
  for(int i=0;i<size();i++){
    source* currSource=(*_sources)[i];
    numUsed+=currSource->
    applyAcousticVel(modelDef,iteration,vx,vy,vz,rho);
  }
  return numUsed;
}
int sourceNetwork::applyAcousticVelMM(modelDefStruct* modelDef,int iteration,
                                      float *vx,float *vy,float *vz,
                                      float* rho){
  int numUsed=0;
  for(int i=0;i<size();i++){
    source* currSource=(*_sources)[i];
    numUsed+=currSource->
    applyAcousticVelMM(modelDef,iteration,vx,vy,vz,rho);
  }
  return numUsed;
}
int sourceNetwork::applyPressure(modelDefStruct* modelDef,int iteration,
                                 float *pressure){
  int numUsed=0;
  for(int i=0;i<size();i++){
    source* currSource=(*_sources)[i];
    numUsed+=currSource->applyPressure(modelDef,iteration,pressure);
  }
  return numUsed;
}
int sourceNetwork::applyPressureMM(modelDefStruct* modelDef,int iteration,
                                   float *pressure){
  int numUsed=0;
  for(int i=0;i<size();i++){
    source* currSource=(*_sources)[i];
    numUsed+=currSource->applyPressureMM(modelDef,iteration,pressure);
  }
  return numUsed;
}

//Bigger routines.
// Write all the current sources to a cdf file.
int sourceNetwork::writeSources(modelDefStruct* modelDef,char* filename,int createNew){
  int outfile,ntDim;
  if(!createNew){
    //Open it to write, note this assumes the file already exists and has no sources already.
    outfile=openCDFFile(filename,FALSE,NC_WRITE);
    //get the time dimesion id.
    if(nc_inq_dimid(outfile,"NT",&ntDim)!=NC_NOERR){
      assert(nc_redef(outfile)==NC_NOERR,
             "writeSources--unable to move %s(%i) into define mode",
             filename,outfile);
      assert(nc_def_dim(outfile,"NT",(size_t)modelDef->NT,&ntDim)==NC_NOERR,
             "writeSources--unable to create NT in file %s(%i)",
             filename,outfile);
      assert(nc_enddef(outfile)==NC_NOERR,
             "writeSources--unable to move %s(%i) out of define mode",
             filename,outfile);
    }
  }else{
    outfile=openCDFFile(filename,TRUE,NC_WRITE);
    assert(nc_def_dim(outfile,"NT",(size_t)modelDef->NT,&ntDim)==NC_NOERR,
           "writeSources--unable to create NT in file %s(%i)",
           filename,outfile);
    assert(nc_enddef(outfile)==NC_NOERR,
           "writeSources--unable to move %s(%i) out of define mode",
           filename,outfile);
  }
  
  //Write the three types of sources
  writeForceSources(modelDef,filename,outfile,ntDim);
  writeMomentSources(modelDef,filename,outfile,ntDim);
  writeTractionSources(modelDef,filename,outfile,ntDim);
  
  assert(nc_close(outfile)==NC_NOERR,
         "writeSources--unable to close file %s",
         filename);
  return _nfs+_nms+_nts;
}

///Set the source differentiation factor.
int sourceNetwork::setDispersionFactor(int i,int doPrint){
  int factor=0;
  switch(type(i)){
    case FORCE_SOURCE:
      tEprintf(doPrint," Source #%i: Force\n",i+1);
      factor=1;
      break;
    case MOMENT_SOURCE:
      tEprintf(doPrint," Source #%i: Moment \n",i+1);
      factor=2;
      break;
    case TDBC_SOURCE:
      tEprintf(doPrint," Source #%i: TDBC (No Dispersion Check)\n",i+1);
      factor=-1;
      break;
      
    default:
      fprintf(stderr,
              " WARNING: checkDispersion--source %i; unknown factor type %i\n",
              i,type(i));
  }
  return factor;
}

///Perform dispersion analysis and check CFL criteria.
int sourceNetwork::checkDispersion(modelDefStruct* modelDef,
                                   float vMin,float vMax,
                                   int doPrint,float *wavelengths){
  DEF_MODEL_SIZE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  int nf=nFreqs(NT);
  float df=dFreq(nf,dt);
  float* wrk=NULL; //variable to hold the spectra
  for(int i=0;i<size();i++){
    int factor=setDispersionFactor(i,_bandCheckPCent>0.0);
    if(factor>=0){
      sourceSpectrum(data(i),factor,
                     NT,dt,
                     nf,df,wrk);
      float curr=
      bandit(_bandCheckPCent,
             dx,dy,dz,df,nf,vMin,wrk,_bandCheckPCent>0.0);
      if(wavelengths)
        SETMINMAX(wavelengths[0],wavelengths[1],curr,i);
    }
  }
  free(wrk);
  
  //Report CFL stability condition criteria.
  tEprintf(doPrint," X-dimension CFL stability ratio = %.3f\n",2.0*vMax*dt/dx);
  tEprintf(doPrint," Y-dimension CFL stability ratio = %.3f\n",2.0*vMax*dt/dy);
  tEprintf(doPrint," Z-dimension CFL stability ratio = %.3f\n",2.0*vMax*dt/dz);
  
  return size();
}

int sourceNetwork::addSources(sourceNetwork* extraSources){
  //First set any parameters that might be needed.
  _bandCheckPCent=extraSources->_bandCheckPCent;
  _useMomentRate = extraSources->_useMomentRate;
  
  //Now add the sources.
  for(int i=0;i<extraSources->size();i++){
    _sources->Add(extraSources->s(i));
    if(extraSources->s(i)->type()==FORCE_SOURCE){
      _nfs++;
    }else if(extraSources->s(i)->type()==MOMENT_SOURCE){
      _nms++;
    } else if(extraSources->s(i)->type()==TDBC_SOURCE){
      ;
    }else{
      assert(FALSE,"Source type %i (S%i) not yet implemented.",
             extraSources->s(i)->type(),i);
    }
  }
  
  //Manually zero the array size in extra sources so the elements
  // are not deleted.
  extraSources->_sources->EmptyArray();
  
  return size();
}
int sourceNetwork::addSources(int& i,int argOffset,
                              int argc,char* argv[],
                              modelDefStruct* modelDef){
  DEF_MODEL_SIZE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  
  static float angleMult=DEG_TO_RAD;
  
  if(strlen(argv[i])<argOffset+1){
    assert(argc>i+1,"addSource--defining wavelet requires 1 argument");
    readWavelet(NT,argv[++i]);
  }else{
    switch(argv[i][argOffset]){
      case 'M':
        _useModelSources=!_useModelSources;
        tEprintf(Verbose,"%s sources from model file\n",
                 _useModelSources?"Using":"Not using");
        break;
        
      case 'F':
        assert(argc>i+1,"addSource--source file requires 1 argument");
      {
        char* filename=argv[++i];
        tEprintf(Verbose,"Adding sources from\n\t\"%s\"\n",filename);
        
        //Open the file for input.
        int infile=openCDFFile(filename,FALSE,NC_NOWRITE);
        
        readForceSources(modelDef,filename,infile,(float*)NULL);
        readMomentSources(modelDef,filename,infile,(float*)NULL);
        readTractionSources(modelDef,filename,infile,(float*)NULL);
        if(nfs()+nms()+nts()>100){
          tEprintf(Verbose," Now %i source; turned off bandwidth check\n",
                   nfs()+nms()+nts());
          _bandCheckPCent=-1;
        }
        
        //Close the file.
        assert(nc_close(infile)==NC_NOERR,
               "addSources--unable to close %s(%i)",filename,infile);
      }
        break;
        
      case 'T':
        assert(argc>i+1,"addSource--TDBC source requires 1 argument");
        tEprintf(Verbose,"Adding TDBC source from\n\t\"%s\"\n",argv[++i]);
      {
        tdbcSource* newSource=new tdbcSource(modelDef,argv[i]);
        _sources->Add(newSource);
      }
        break;
        
      case 'R':
        assert(argc>i+1,"addSource--ramp function wavelet requires shift");
      {
        int shift;
        if(strlen(argv[i])>argOffset+1 && argv[i][argOffset+1]=='t'){
          float otime=atof(argv[++i]);
          shift=(int)floor(otime/dt);
        }else{
          shift=atoi(argv[++i]);
        }
        if(!_wavelet)
          assert((_wavelet=(float*)malloc(NT*sizeof(float)))!=NULL,
                 "addSource--unable to allocate %i floats for wavelet",NT);
        for(int iii=0;iii<NT;++iii) _wavelet[iii]=iii<shift?0.0:(iii-shift)*dt;
        
        tEprintf(Verbose,
                 "Ramp function source; origin index %i (t %.3fms)\n",
                 shift,minT+shift*dt);
      }
        break;
      case 'g':
        assert(argc>i+1,"addSource--gaussian function wavelet requires central frequency");
      {
        float freq = atof(argv[++i]);
        if(!_wavelet)
          assert((_wavelet=(float*)malloc(NT*sizeof(float)))!=NULL,
                 "addSource--unable to allocate %i floats for wavelet",NT);
        buildGaussianWavelet(freq,NT,dt);
        tEprintf(Verbose,"Gaussian function source with central frequency %g\n",freq);
      }
        break;
      case 'H':
        assert(argc>i+1,"addSource--heavyside function wavelet requires shift");
      {
        int shift;
        if(strlen(argv[i])>argOffset+1 && argv[i][argOffset+1]=='t'){
          float otime=atof(argv[++i]);
          shift=(int)floor(otime/dt);
        }else{
          shift=atoi(argv[++i]);
        }
        if(!_wavelet)
          assert((_wavelet=(float*)malloc(NT*sizeof(float)))!=NULL,
                 "addSource--unable to allocate %i floats for wavelet",NT);
        for(int iii=0;iii<NT;++iii) _wavelet[iii]=iii<shift?0.0:1.0;
        
        tEprintf(Verbose,
                 "Heavyside function source; origin index %i (t %.3fms)\n",
                 shift,minT+shift*dt);
      }
        break;
        
      case 'D':
        assert(argc>i+1,"addSource--delta function wavelet requires shift");
      {
        int shift;
        if(strlen(argv[i])>argOffset+1 && argv[i][argOffset+1]=='t'){
          float otime=atof(argv[++i]);
          shift=(int)floor(otime/dt);
        }else{
          shift=atoi(argv[++i]);
        }
        if(!_wavelet)
          assert((_wavelet=(float*)malloc(NT*sizeof(float)))!=NULL,
                 "addSource--unable to allocate %i floats for wavelet",NT);
        for(int iii=0;iii<NT;_wavelet[iii++]=0.0);
        
        if(strlen(argv[i-1])>argOffset+1 && argv[i-1][argOffset+1]=='1'){
          _wavelet[shift]=1.0;
        }else{
          _wavelet[shift]=1.0/dt;
        }
        tEprintf(Verbose,
                 "Delta function source; origin index %i (t %.3fms), amp %.2f\n",
                 shift,minT+shift*dt,_wavelet[shift]);
      }
        break;
        
      case 'b':
        if(strlen(argv[i])>argOffset+1 &&
           (argv[i][argOffset+1]=='n' || argv[i][argOffset+1]=='N')){
          _bandCheckPCent=-1;
          tEprintf(Verbose,"No bandwidth check on sources\n");
        }else{
          assert(argc>i+1,"addSource--changing bandwidth check percentage requires aregument");
          _bandCheckPCent=atof(argv[++i]);
          tEprintf(Verbose,"Set bandwidth check percentage to %.1f%%\n",
                   _bandCheckPCent);
        }
        break;
      case 'N':
      {
        char* callFlag=argv[i];
        
        //Time reverse the wavelet and shift by the given time amount.
        assert(argc>i+1,"addSource--wavelet reversal requires shift");
        int shift=atoi(argv[++i]);
        int mid=(int)floor((float)NT/2.0);
        for(int ii=0;ii<mid;ii++){
          float temp=_wavelet[ii];
          _wavelet[ii]=_wavelet[NT-1-ii];
          _wavelet[NT-1-ii]=temp;
        }
        if(shift>0){
          for(int ii=0;ii<NT;ii++){
            _wavelet[ii]=((shift+ii)<NT)?_wavelet[ii+shift]:0.0;
          }
        }else if(shift<0){
          for(int ii=NT-1;ii>0;ii--){
            _wavelet[ii]=((ii+shift)>0)?_wavelet[ii+shift]:0.0;
          }
        }
        
        if(strlen(callFlag)>3 && callFlag[3]=='n'){
          for(int ii=0;ii<NT;ii++)
            _wavelet[ii]*=-1.0;
          tEprintf(Verbose,
                   "Time reversed and negated wavelet; shifted by %i iterations\n",
                   shift);
        }else{
          tEprintf(Verbose,
                   "Time reversed wavelet; shifted by %i iterations (%.2fms)\n",
                   shift,1000.0*shift*dt);
        }
      }
        break;
        
      case 'w':
      {
        assert(argc>i+1,"addSource--defining wavelet requires 1 argument");
        readWavelet(NT,argv[++i]);
      }
        break;
        
      case 'r':
        assert(argc>i+1,
               "addSource--Ricker wavelet definition requires peak frequency");
      {
        char* callFlag=argv[i];
        float freq=atof(argv[++i]);
        int derivative=0;
        if(argc>i+1 && !strcmp(argv[i+1],"0")){
          i++;
          derivative=0;
        }else if(argc>i+1 && !strcmp(argv[i+1],"1")){
          i++;
          derivative=1;
        }else if(argc>i+1 && !strcmp(argv[i+1],"2")){
          i++;
          derivative=2;
        }
        
        if(_wavelet) free(_wavelet);
        if(strlen(callFlag)<argOffset+1 || callFlag[argOffset+1]!='z'){
          rickerWavelet(modelDef->NT,modelDef->dt,freq,derivative);
        }else{
          rickerWavelet(modelDef->NT,modelDef->dt,freq,derivative,
                        TRUE,-modelDef->minT);
        }
        tEprintf(Verbose,
                 "Built Ricker wavelet: frequency %.2f; normalized %i derivative\n",
                 freq,derivative);
      }
        break;
        
      case 'u':
        if(strlen(argv[i])<argOffset+1 || argv[i][argOffset+1] == 'r'){
          angleMult=1.0;
          tEprintf(Verbose,"Source Angles in Degrees (Default)\n");
        }else if(argv[i][argOffset+1] == 'd'){
          angleMult=DEG_TO_RAD;
          tEprintf(Verbose,"Source Angles in Radians\n");
        }else{
          assert(FALSE,"addReceivers--u option is [d|r] only; arg is %s",
                 argv[i]);
        }
        break;
      case 'f':
        assert(argc>i+4,"addSource--force source requires loc (3) and amp");
      {
        char *currArg=argv[i];
        
        //read source location and orientation
        float x=atof(argv[++i]);
        float y=atof(argv[++i]);
        float z=atof(argv[++i]);
        
        float amp=atof(argv[++i]);
        
        //Check for along axis direction in flag.
        float ax,ay,az;
        if(strlen(currArg)>3 && currArg[3]=='x'){
          az=ay=0.0;
          ax=1.0;
        }else if(strlen(currArg)>3 && currArg[3]=='y'){
          ax=az=0.0;
          ay=1.0;
        }else if(strlen(currArg)>3 && currArg[3]=='z'){
          ax=ay=0.0;
          az=1.0;
        }else{
          assert(argc>i+2,"addSource--free source direction requires 2 angles");
          float theta=atof(argv[++i])*angleMult;
          float phi=atof(argv[++i])*angleMult;
          
          ax=sin(theta)*cos(phi);
          ay=sin(theta)*sin(phi);
          az=cos(theta);
          
          //Check for near along axis.
          if(fabs(ax-1.0)<1e-3){
            ax=1.0;
            ay=az=0.0;
          }else if(fabs(ax+1.0)<1e-3){
            ax=-1.0;
            ay=az=0.0;
          }else if(fabs(ay-1.0)<1e-3){
            ay=1.0;
            ax=az=0.0;
          }else if(fabs(ay+1.0)<1e-3){
            ay=-1.0;
            ax=az=0.0;
          }else if(fabs(az-1.0)<1e-3){
            az=1.0;
            ax=ay=0.0;
          }else if(fabs(az+1.0)<1e-3){
            az=-1.0;
            ax=ay=0.0;
          }
        }
        assert(_wavelet!=NULL,
               "addSource--adding force source requires wavelet be initialized");
        _sources->Add(new forceSource(modelDef,
                                      x,y,z,ax,ay,az,
                                      amp,_wavelet,(float*)NULL));
        _nfs++;
        
        tEprintf(Verbose,
                 "Added force source %i (%.1f,%.1f,%.1f) (%.1f,%.1f,%.1f,%.1f)\n",
                 size(),x,y,z,amp,ax,ay,az);
      }
        break;
      case 'e':
        assert(argc>i+4,"addSource--explosion requires 4 arguments");
      {
        float x=atof(argv[++i]);
        float y=atof(argv[++i]);
        float z=atof(argv[++i]);
        float amp=atof(argv[++i]);
        if(_useMomentRate) amp*=dt;
        
        float xxS=1.0;
        
        assert(_wavelet!=NULL,
               "addSource--adding explosion requires wavelet be initialized");
        _sources->Add(new momentSource(modelDef,
                                       x,y,z,
                                       xxS,xxS,xxS,
                                       0.0,0.0,0.0,
                                       0.0,0.0,0.0,
                                       amp,_wavelet,NULL,_useMomentRate));
        _nms++;
        
        tEprintf(Verbose,
                 "Added explosive moment source %i (%.1f,%.1f,%.1f,%.3f)\n",
                 size(),x,y,z,amp);
      }
        break;
      case 'd':
        //Definitions from Gubbins; pp198-201
        assert(argc>i+7,"addSource--double couple requires 7 arguments");
      {
        //allocate the source
        float x=atof(argv[++i]);
        float y=atof(argv[++i]);
        float z=atof(argv[++i]);
        float amp=atof(argv[++i]);
        if(_useMomentRate) amp*=dt;
        
        //Read the angles.
        // Note the strike is measured from x-direction north.
        double strike=atof(argv[++i])*angleMult;
        double dip=atof(argv[++i])*angleMult;
        double rake=atof(argv[++i])*angleMult;
        
        //calculate the components of the fault normal vector
        double normalVec[3]={-sin(dip)*sin(strike),sin(dip)*cos(strike),-cos(dip)};
        tEprintf(Verbose,
                 "Double Couple vectors: normal (%.4f,%.4f,%.4f)\n",
                 normalVec[0],normalVec[1],normalVec[2]);
        
        double slipVec[3]={0.0,0.0,0.0};
        if(fabs(normalVec[0])+fabs(normalVec[1])<0.01){
          //Fault plane is in the z plane, strike is undefined, calculate rake from north.
          slipVec[0]=cos(rake);
          slipVec[1]=sin(rake);
        }else{
          //Start with a north pointing vector. And rotate into the correct rake direction
          slipVec[0]=cos(rake);
          slipVec[2]=sin(rake);
          
          //Rotate by the dip angle.
          slipVec[1] =slipVec[2]*cos(dip);
          slipVec[2]*=sin(dip);
          
          //Rotate by the strike.
          double temp=cos(strike)*slipVec[0]-sin(strike)*slipVec[1];
          slipVec[1]= sin(strike)*slipVec[0]+cos(strike)*slipVec[1];
          slipVec[0]=temp;
        }
        
        tEprintf(Verbose,
                 "                       slip   (%.4f,%.4f,%.4f)\n",
                 slipVec[0],slipVec[1],slipVec[2]);
        assert(fabs(slipVec[0]*normalVec[0]+
                    slipVec[1]*normalVec[1]+
                    slipVec[2]*normalVec[2])<0.01,
               "addSource--Fault Source; Normal and slip not orthogonal (dot product %f)",
               slipVec[0]*normalVec[0]+slipVec[1]*normalVec[1]+slipVec[2]*normalVec[2]);
        
        //Fill in the moment components
        // Aki & Richards 3.21 (p52)
        assert(_wavelet!=NULL,
               "addSource--adding force source requires wavelet be initialized");
        _sources->Add(new momentSource(modelDef,x,y,z,
                                       2*normalVec[0]*slipVec[0],
                                       2*normalVec[1]*slipVec[1],
                                       2*normalVec[2]*slipVec[2],
                                       normalVec[0]*slipVec[1]+normalVec[1]*slipVec[0],
                                       normalVec[0]*slipVec[2]+normalVec[2]*slipVec[0],
                                       normalVec[1]*slipVec[2]+normalVec[2]*slipVec[1],
                                       0.0,0.0,0.0,
                                       amp,_wavelet,NULL,_useMomentRate));
        _nms++;
        
        // Trace should be zero (vectors are supposed to be orthogonal)
        assert(fabs(xxS(size()-1)+yyS(size()-1)+zzS(size()-1)) < 0.01,
               "addSource--Fault Source; Non-zero trace in moment tensor");
        
        
        tEprintf(Verbose,
                 "Added moment source %i (%.1f,%.1f,%.1f,%.3f)\n",
                 size(),
                 x,y,z,amp);
        tEprintf(Verbose,"                      (%6.1f,%6.1f,%6.1f)\n",
                 xxS(size()-1),xyS(size()-1)+xyA(size()-1),xzS(size()-1)+xzA(size()-1));
        tEprintf(Verbose,"                      (%6.1f,%6.1f,%6.1f)\n",
                 xyS(size()-1)+xyA(size()-1),yyS(size()-1),yzS(size()-1)-yzA(size()-1));
        tEprintf(Verbose,"                      (%6.1f,%6.1f,%6.1f)\n",
                 xzS(size()-1)-xzA(size()-1),yzS(size()-1)-yzA(size()-1),zzS(size()-1));
      }
        break;
        
      case 'm':
        if(strlen(argv[i])>3 && argv[i][3]=='r') {
          tEprintf(Verbose,"All moment waveforms will be assumed moment rate waveforms\n");
          _useMomentRate = true;
        } else {
          assert(argc>i+13,"addSource--moment source requires 13 arguments");
          {
            float x=atof(argv[++i]);
            float y=atof(argv[++i]);
            float z=atof(argv[++i]);
            float amp=atof(argv[++i]);
            if(_useMomentRate) amp*=dt;
            
            //Read the moment components
            float axx=atof(argv[++i]);
            float axy=atof(argv[++i]);
            float axz=atof(argv[++i]);
            float ayx=atof(argv[++i]);
            float ayy=atof(argv[++i]);
            float ayz=atof(argv[++i]);
            float azx=atof(argv[++i]);
            float azy=atof(argv[++i]);
            float azz=atof(argv[++i]);
            
            assert(_wavelet!=NULL,
                   "addSource--adding force source requires wavelet be initialized");
            _sources->Add(new momentSource(modelDef,x,y,z,
                                           axx,ayy,azz,
                                           (axy+ayx)/2.0,(axz+azx)/2.0,(ayz+azy)/2.0,
                                           (axy-ayx)/2.0,(axz-azx)/2.0,(ayz-azy)/2.0,
                                           amp,_wavelet,NULL,_useMomentRate));
            _nms++;
            
            tEprintf(Verbose,
                     "Added moment source %i (%.1f,%.1f,%.1f,%.3f)\n",
                     size(),x,y,z,amp);
            tEprintf(Verbose,"                               (%6.1f,%6.1f,%6.1f)\n",
                     xxS(size()-1),xyS(size()-1)+xyA(size()-1),xzS(size()-1)+xzA(size()-1));
            tEprintf(Verbose,"                               (%6.1f,%6.1f,%6.1f)\n",
                     xyS(size()-1)-xyA(size()-1),yyS(size()-1),yzS(size()-1)+yzA(size()-1));
            tEprintf(Verbose,"                               (%6.1f,%6.1f,%6.1f)\n",
                     xzS(size()-1)-xzA(size()-1),yzS(size()-1)-yzA(size()-1),zzS(size()-1));
          }
        }
        break;
        
      default:
        assert(FALSE,"addSource--unable to process argument %s (%c) ",
               argv[i],argv[i][argOffset]);
    }
  }
  return size();
}

void sourceNetwork::checkWaveletDispersion(float dx,float dy,float dz,
                                           int nt,float dt,
                                           float minVel,int derivative){
  float *fullSpec=NULL;
  
  //Calculate the spectrum in the full format and run bandit.
  float df;
  int nf;
  sourceSpectrum(wavelet(),derivative,nt,dt,nf,df,fullSpec);
  bandit(_bandCheckPCent,
         dx,dy,dz,df,nf,minVel,fullSpec);
  free(fullSpec);
}

float* sourceNetwork::readWavelet(int NT,char* filename){
  if(!_wavelet){
    assert((_wavelet=(float*)malloc(NT*sizeof(float)))!=NULL,
           "readWavelet--unable to allocate %i floats for wavelet",NT);
  }
  
  FILE* infile=fopen(filename,"r");
  assert(infile!=NULL,
         "sourceTimeFunc--unable to open source file \"%s\"\n",
         filename);
  for(int i=0;i<NT;i++){
    float time,value;
    if(fscanf(infile,"%g %g\n",&time,&value)==EOF){
      for(int j=i;j<NT;_wavelet[j++]=0.0);
      break;
    }
    _wavelet[i]=value;
  }
  fclose(infile);
  
  tEprintf(Verbose,
           "Read wavelet from file \"%s\"\n",
           filename);
  return _wavelet;
}

float* sourceNetwork::rickerWavelet(int nt,float dt,float fmode,
                                    int normalizationDerivative,
                                    int setPeakTime,float peakTime){
  //Determine the number of frequencies required for this length of time series.
  double exponent=ceil(log((float)nt)/log((float)2));
  int nsamps=2*(int)ceil(pow(2.0,exponent));
  int nf=nsamps/2;
  
  //Determine the frequency step
  double df=1.0/((float)nsamps*dt);
  double tmin=1.25/fmode;
  
  //Get the spectrum.
  double *spect;
  assert((spect=(double*)malloc((1+2*nsamps)*sizeof(double)))!=NULL,
         "ricker--unable to allocate 2*%i floats for spectrum",
         nsamps);
  rickerSpec(nf,df,fmode,tmin,spect);
  
  if(Verbose>1){
    FILE* out=fopen("ricker_spectrum.txt","w");
    for(int j=0;j<nf;j++)
      fprintf(out,"%f\t%f\n",spect[2*j+0],spect[2*j+1]);
    fclose(out);
  }
  
  //Transform into the time domain. Remember that numerical recipies subroutines are
  // indexed from 1, not 0.
  NUMREC_dfour1(spect-1,
                nsamps,-1);
  
  if(Verbose>1){
    FILE* out=fopen("ricker_raw_td.txt","w");
    for(int j=0;j<nt;j++)
      fprintf(out,"%f\n",spect[2*j+0]);
    fclose(out);
  }
  
  //Get the scale factor to set the max amplitude to 1.
  int maxValIndex;
  double maxVal=0.0;
  for(int i=0;i<nt;i++){
    switch(normalizationDerivative){
      case 0:
        SETMAXINDEX(maxVal,spect[2*i],maxValIndex,i);
        break;
      case 1:
        if(i>0){
          double diff=spect[2*i]-spect[2*(i-1)];
          SETMAXINDEX(maxVal,diff,maxValIndex,i-1);
        }
        break;
      case 2:
        if(i>2){
          double diff=(spect[2*i]-spect[2*(i-1)])-(spect[2*(i-2)]-spect[2*(i-3)]);
          SETMAXINDEX(maxVal,diff,maxValIndex,i-3);
        }
        break;
        
      default:
        assert(FALSE,"ricker--normalizationDerivatives of 0, 1, or 2 only supported not %i",
               normalizationDerivative);
    }
  }
  
  float *data;
  assert((data=(float*)malloc(nt*sizeof(float)))!=NULL,
         "ricker--unable to allocate %i floats for data",
         nt);
  for(int i=0;i<nt;data[i++]=0.0);
  
  if(!setPeakTime){
    for(int i=0;i<nt;i++)
      data[i]=spect[2*i]/maxVal;
  }else{
    //Calculate the index offset
    int offset=maxValIndex-(int)floor(peakTime/dt);
    if(offset==0){
      for(int i=0;i<nt;i++)
        data[i]=spect[2*i]/maxVal;
    }else if(offset>0){
      for(int i=0;i<nt-offset;i++)
        data[i]=spect[2*(i+offset)]/maxVal;
    }else{
      for(int i=0;i<nt+offset;i++)
        data[i-offset]=spect[2*i]/maxVal;
    }
  }
  free(spect);
  
  if(Verbose>1){
    FILE* out=fopen("ricker_final_td.txt","w");
    for(int j=0;j<nt;j++)
      fprintf(out,"%f\n",data[j]);
    fclose(out);
  }
  
  return _wavelet=data;
}

void sourceNetwork::buildGaussianWavelet(float freq,int nt,float dt) {
  float tmin = -1.2/freq;
  if(nt<ceil(2.4/freq/dt)) {
    fprintf(stderr,"Error: nt (%d) is too small for given frequency (min: %d)\n",nt,static_cast<int>(ceil(2.4/freq/dt)));
  }
  for(int i=0;i<nt;++i) {
    float t = tmin+i*dt;
    float alpha = freq*t*PI;
    _wavelet[i] = exp(-sqr(alpha));
  }
}

//
//Helper routines to write sources to CDF file.

///Write force sources to a NetCDF file.
int sourceNetwork::writeForceSources(modelDefStruct* modelDef,char* filename,
                                     int outfile,int ntDim){
  if(!_nfs) return FALSE;
  
  //Define the required variables.
  DEF_MODEL_SIZE(modelDef);
  
  nc_redef(outfile);
  int nfsDim;
  assert(nc_def_dim(outfile,"numFSources",_nfs,&nfsDim)==NC_NOERR,
         "writeForceSources--unable to define dimension numFSources[%i] in %s(%i)",
         _nfs,filename,outfile);
  
  int fsXsVar,fsYsVar,fsZsVar,fsSampVar;
  assert(nc_def_var(outfile,"fsourcesXs",NC_FLOAT,1,&nfsDim,&fsXsVar)==NC_NOERR,
         "writeForceSources--unable to define var fsourcesXs in %s(%i,%i)",
         filename,outfile,nfsDim);
  assert(nc_def_var(outfile,"fsourcesYs",NC_FLOAT,1,&nfsDim,&fsYsVar)==NC_NOERR,
         "writeForceSources--unable to define var fsourcesYs in %s(%i,%i)",
         filename,outfile,nfsDim);
  assert(nc_def_var(outfile,"fsourcesZs",NC_FLOAT,1,&nfsDim,&fsZsVar)==NC_NOERR,
         "writeForceSources--unable to define var fsourcesZs in %s(%i,%i)",
         filename,outfile,nfsDim);
  assert(nc_def_var(outfile,"fsourcesSamp",NC_FLOAT,1,&nfsDim,&fsSampVar)==NC_NOERR,
         "writeForceSources--unable to define var fsourcesSamp in %s(%i,%i)",
         filename,outfile,nfsDim);
  
  int extraVars[3];
  assert(nc_def_var(outfile,"fsourcesAx",NC_FLOAT,1,&nfsDim,&extraVars[0])==NC_NOERR,
         "writeForceSources--unable to define var fsourcesAx in %s(%i,%i)",
         filename,outfile,nfsDim);
  assert(nc_def_var(outfile,"fsourcesAy",NC_FLOAT,1,&nfsDim,&extraVars[1])==NC_NOERR,
         "writeForceSources--unable to define var fsourcesAy in %s(%i,%i)",
         filename,outfile,nfsDim);
  assert(nc_def_var(outfile,"fsourcesAz",NC_FLOAT,1,&nfsDim,&extraVars[2])==NC_NOERR,
         "writeForceSources--unable to define var fsourcesAz in %s(%i,%i)",
         filename,outfile,nfsDim);
  
  int fsDims[2]={nfsDim,ntDim};
  int fsDataVar;
  assert(nc_def_var(outfile,"fsourcesData",NC_FLOAT,2,fsDims,&fsDataVar)==NC_NOERR,
         "writeForceSources--unable to define var fsourcesData in %s(%i,%i,%i)",
         filename,outfile,nfsDim,ntDim);
  assert(nc_enddef(outfile)==NC_NOERR,
         "writeForceSources--unable to take %s(%i) out of define mode",
         filename,outfile);
  
  //Now do the actual writing of source data.
  int numWritten=0;
  DEF_MODEL_SCALE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  float forceAmpScalar=1.0/(dt/(scalarVel*scalarDen*dx*dy*dz));
  for(int i=0;i<size();i++){
    if((*_sources)[i]->type() == FORCE_SOURCE){
      (*_sources)[i]->writeSource(NT,forceAmpScalar,i,filename,outfile,
                                  fsXsVar,fsYsVar,fsZsVar,fsSampVar,fsDataVar,
                                  extraVars);
      numWritten++;
    }
  }
  assert(numWritten==_nfs,
         "writeForceSources-- _nfs is %i, only wrote (found) %i force sources",
         _nfs,numWritten);
  
  //KLUDGE: this should no be required but I have been getting dropouts in the sourceTime
  // function for certian values of NT.
  assert(nc_sync(outfile)==NC_NOERR,
         "writeForceSources--unable to sync %s(%i)",
         filename,outfile);
  return _nfs;
}

///Write moment sources to a NetCDF file.
int sourceNetwork::writeMomentSources(modelDefStruct* modelDef,char* filename,
                                      int outfile,int ntDim){
  DEF_MODEL_SIZE(modelDef);
  if(!_nms)
    return FALSE;
  
  assert(nc_redef(outfile)==NC_NOERR,
         "writeMomentSources--unable to put %s(%i) into define mode",
         filename,outfile);
  int nmsDim;
  assert(nc_def_dim(outfile,"numMSources",_nms,&nmsDim)==NC_NOERR,
         "writeMomentSources--unable to define dimension numMSources[%i] in %s(%i)",
         _nms,filename,outfile);
  
  int msXsVar,msYsVar,msZsVar,msSampVar;
  assert(nc_def_var(outfile,"mSourcesXs",NC_FLOAT,1,&nmsDim,&msXsVar)==NC_NOERR,
         "writeMomentSources--unable to define var mSourcesXs in %s(%i,%i)",
         filename,outfile,nmsDim);
  assert(nc_def_var(outfile,"mSourcesYs",NC_FLOAT,1,&nmsDim,&msYsVar)==NC_NOERR,
         "writeMomentSources--unable to define var mSourcesYs in %s(%i,%i)",
         filename,outfile,nmsDim);
  assert(nc_def_var(outfile,"mSourcesZs",NC_FLOAT,1,&nmsDim,&msZsVar)==NC_NOERR,
         "writeMomentSources--unable to define var mSourcesZs in %s(%i,%i)",
         filename,outfile,nmsDim);
  assert(nc_def_var(outfile,"mSourcesSamp",NC_FLOAT,1,&nmsDim,&msSampVar)==NC_NOERR,
         "writeMomentSources--unable to define var mSourcesSamp in %s(%i,%i)",
         filename,outfile,nmsDim);
  
  int extraVars[9];
  assert(nc_def_var(outfile,"mSourcesXxS",NC_FLOAT,1,&nmsDim,&extraVars[0])==NC_NOERR,
         "writeMomentSources--unable to define var mSourcesXxS in %s(%i,%i)",
         filename,outfile,nmsDim);
  assert(nc_def_var(outfile,"mSourcesYyS",NC_FLOAT,1,&nmsDim,&extraVars[1])==NC_NOERR,
         "writeMomentSources--unable to define var mSourcesYyS in %s(%i,%i)",
         filename,outfile,nmsDim);
  assert(nc_def_var(outfile,"mSourcesZzS",NC_FLOAT,1,&nmsDim,&extraVars[2])==NC_NOERR,
         "writeMomentSources--unable to define var mSourcesZzS in %s(%i,%i)",
         filename,outfile,nmsDim);
  assert(nc_def_var(outfile,"mSourcesXyS",NC_FLOAT,1,&nmsDim,&extraVars[3])==NC_NOERR,
         "writeMomentSources--unable to define var mSourcesXyS in %s(%i,%i)",
         filename,outfile,nmsDim);
  assert(nc_def_var(outfile,"mSourcesXzS",NC_FLOAT,1,&nmsDim,&extraVars[4])==NC_NOERR,
         "writeMomentSources--unable to define var mSourcesXzS in %s(%i,%i)",
         filename,outfile,nmsDim);
  assert(nc_def_var(outfile,"mSourcesYzS",NC_FLOAT,1,&nmsDim,&extraVars[5])==NC_NOERR,
         "writeMomentSources--unable to define var mSourcesYzS in %s(%i,%i)",
         filename,outfile,nmsDim);
  assert(nc_def_var(outfile,"mSourcesXyA",NC_FLOAT,1,&nmsDim,&extraVars[6])==NC_NOERR,
         "writeMomentSources--unable to define var mSourcesXyA in %s(%i,%i)",
         filename,outfile,nmsDim);
  assert(nc_def_var(outfile,"mSourcesXzA",NC_FLOAT,1,&nmsDim,&extraVars[7])==NC_NOERR,
         "writeMomentSources--unable to define var mSourcesXzA in %s(%i,%i)",
         filename,outfile,nmsDim);
  assert(nc_def_var(outfile,"mSourcesYzA",NC_FLOAT,1,&nmsDim,&extraVars[8])==NC_NOERR,
         "writeMomentSources--unable to define var mSourcesYzA in %s(%i,%i)",
         filename,outfile,nmsDim);
  
  int msDims[2]={nmsDim,ntDim};
  int msDataVar;
  assert(nc_def_var(outfile,"mSourcesData",NC_FLOAT,2,msDims,&msDataVar)==NC_NOERR,
         "writeMomentSources--unable to define mSourcesData(%i,%i) in %s(%i)",
         msDims[0],msDims[1],filename,outfile);
  assert(nc_enddef(outfile)==NC_NOERR,
         "writeMomentSources--unable to take %s(%i) out of define mode",
         filename,outfile);
  
  int numWritten=0;
  DEF_MODEL_SCALE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  float momentAmpScalar= -scalarStress*dx*dy*dz;
  for(int i=0;i<size();i++){
    if((*_sources)[i]->type() == MOMENT_SOURCE){
      (*_sources)[i]->writeSource(NT,momentAmpScalar,i,filename,outfile,
                                  msXsVar,msYsVar,msZsVar,msSampVar,msDataVar,
                                  extraVars);
      numWritten++;
    }
  }
  assert(numWritten==_nms,
         "writeMomentSources-- _nms is %i, only wrote (found) %i force sources",
         _nms,numWritten);
  
  //KLUDGE: this should no be required but I have been getting dropouts in the sourceTime
  // function for certian values of NT.
  assert(nc_sync(outfile)==NC_NOERR,
         "writeMomentSources--unable to sync %s(%i)",
         filename,outfile);
  return _nms;
}

///Write traction sources to a NetCDF file.
int sourceNetwork::writeTractionSources(modelDefStruct* modelDef,char* filename,
                                        int outfile,int ntDim){
  if(!_nts) return FALSE;
  
  nc_redef(outfile);
  int ntsDim;
  assert(nc_def_dim(outfile,"numTSources",_nts,&ntsDim)==NC_NOERR,
         "writeTractionSources--unable to define dimension numTSources[%i] in %s(%i)",
         _nts,filename,outfile);
  
  int tsXsVar,tsYsVar,tsZsVar,tsSampVar;
  assert(nc_def_var(outfile,"tSourcesXs",NC_FLOAT,1,&ntsDim,&tsXsVar)==NC_NOERR,
         "writeTractionSources--unable to define var tSourcesXs in %s(%i,%i)",
         filename,outfile,ntsDim);
  assert(nc_def_var(outfile,"tSourcesYs",NC_FLOAT,1,&ntsDim,&tsYsVar)==NC_NOERR,
         "writeTractionSources--unable to define var tSourcesYs in %s(%i,%i)",
         filename,outfile,ntsDim);
  assert(nc_def_var(outfile,"tSourcesZs",NC_FLOAT,1,&ntsDim,&tsZsVar)==NC_NOERR,
         "writeTractionSources--unable to define var tSourcesZs in %s(%i,%i)",
         filename,outfile,ntsDim);
  assert(nc_def_var(outfile,"tSourcesSamp",NC_FLOAT,1,&ntsDim,&tsSampVar)==NC_NOERR,
         "writeTractionSources--unable to define var tSourcesSamp in %s(%i,%i)",
         filename,outfile,ntsDim);
  
  int extraVars[3];
  assert(nc_def_var(outfile,"tSourcesAx",NC_FLOAT,1,&ntsDim,&extraVars[0])==NC_NOERR,
         "writeTractionSources--unable to define var tSourcesAx in %s(%i,%i)",
         filename,outfile,ntsDim);
  assert(nc_def_var(outfile,"tSourcesAy",NC_FLOAT,1,&ntsDim,&extraVars[1])==NC_NOERR,
         "writeTractionSources--unable to define var tSourcesAy in %s(%i,%i)",
         filename,outfile,ntsDim);
  assert(nc_def_var(outfile,"tSourcesAz",NC_FLOAT,1,&ntsDim,&extraVars[2])==NC_NOERR,
         "writeTractionSources--unable to define var tSourcesAz in %s(%i,%i)",
         filename,outfile,ntsDim);
  
  int tsDims[2]={ntsDim,ntDim};
  int tsDataVar;
  assert(nc_def_var(outfile,"tSourcesData",NC_FLOAT,2,tsDims,&tsDataVar)==NC_NOERR,
         "writeTractionSources--unable to define var tSourcesData in %s(%i,%i,%i)",
         filename,outfile,ntsDim,ntDim);
  assert(nc_enddef(outfile)==NC_NOERR,
         "writeTractionSources--unable to take %s(%i) out of define mode",
         filename,outfile);
  
  //Now do the actual writing of source data.
  int numWritten=0;
  for(int i=0;i<size();i++){
    if((*_sources)[i]->type() == TRACTION_SOURCE){
      numWritten++;
    }
  }
  assert(numWritten==_nts,
         "writeTractionSources-- _nts is %i, only wrote (found) %i force sources",
         _nts,numWritten);
  
  //KLUDGE: this should no be required but I have been getting dropouts in the sourceTime
  // function for certian values of NT.
  assert(nc_sync(outfile)==NC_NOERR,
         "writeTractionSources--unable to sync %s(%i)",
         filename,outfile);
  return _nts;
}

//Helper routines to read sources from CDF file.
///Read force sources.
int sourceNetwork::readForceSources(modelDefStruct* modelDef,
                                    char* filename,int infile,
                                    float* rho){
  int nfsDim;
  if(nc_inq_dimid(infile,"numFSources",&nfsDim)!=NC_NOERR)
    return _nfs;
  
  //Read the dimension to find the actual number of sources.
  size_t nfs;
  assert(nc_inq_dimlen(infile,nfsDim,&nfs)==NC_NOERR,
         "sourceNetwork--unable read dim numFSources in %s(%i,%i)",
         filename,infile,nfsDim);
  _nfs+=nfs;
  
  //Get the required variable id's.
  int fsXsVar,fsYsVar,fsZsVar;
  assert(nc_inq_varid(infile,"fsourcesXs",&fsXsVar)==NC_NOERR,
         "sourceNetwork--unable to open variable fsourcesXs in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"fsourcesYs",&fsYsVar)==NC_NOERR,
         "sourceNetwork--unable to open variable fsourcesYs in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"fsourcesZs",&fsZsVar)==NC_NOERR,
         "sourceNetwork--unable to open variable fsourcesZs in file %s(%i)",
         filename,infile);
  
  int fsSampVar,fsAxVar,fsAyVar,fsAzVar;
  assert(nc_inq_varid(infile,"fsourcesSamp",&fsSampVar)==NC_NOERR,
         "sourceNetwork--unable to open variable fsourcesSamp in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"fsourcesAx",&fsAxVar)==NC_NOERR,
         "sourceNetwork--unable to open variable fsourcesAx in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"fsourcesAy",&fsAyVar)==NC_NOERR,
         "sourceNetwork--unable to open variable fsourcesAy in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"fsourcesAz",&fsAzVar)==NC_NOERR,
         "sourceNetwork--unable to open variable fsourcesAz in file %s(%i)",
         filename,infile);
  
  int fsDataVar;
  assert(nc_inq_varid(infile,"fsourcesData",&fsDataVar)==NC_NOERR,
         "sourceNetwork--unable to open variable fsourcesData in file %s(%i)",
         filename,infile);
  
  //Create a new source for each.
  for(int i=0;i<_nfs;i++)
    _sources->Add(new forceSource(modelDef,filename,infile,i,
                                  fsXsVar,fsYsVar,fsZsVar,
                                  fsAxVar,fsAyVar,fsAzVar,
                                  fsSampVar,fsDataVar,
                                  rho));
  
  return _nfs;
}

///Method to read any moment sources form the CDF file.
int sourceNetwork::readMomentSources(modelDefStruct* modelDef,
                                     char* filename,int infile,
                                     float* rho){
  int nmsDim;
  if(nc_inq_dimid(infile,"numMSources",&nmsDim)!=NC_NOERR)
    return _nms=0;
  
  //Read the dimension to find the actual number of sources.
  assert(nc_inq_dimlen(infile,nmsDim,&_nms)==NC_NOERR,
         "sourceNetwork--unable read dim numMSources in %s(%i,%i)",
         filename,infile,nmsDim);
  
  //Get the required variable id's.
  int msXsVar,msYsVar,msZsVar;
  assert(nc_inq_varid(infile,"mSourcesXs",&msXsVar)==NC_NOERR,
         "sourceNetwork--unable to open variable mSourcesXs in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"mSourcesYs",&msYsVar)==NC_NOERR,
         "sourceNetwork--unable to open variable mSourcesYs in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"mSourcesZs",&msZsVar)==NC_NOERR,
         "sourceNetwork--unable to open variable mSourcesZs in file %s(%i)",
         filename,infile);
  
  int msSampVar;
  assert(nc_inq_varid(infile,"mSourcesSamp",&msSampVar)==NC_NOERR,
         "sourceNetwork--unable to open variable mSourcesSamp in file %s(%i)",
         filename,infile);
  
  int msXxSVar,msYySVar,msZzSVar,msXySVar,msXzSVar,msYzSVar;
  assert(nc_inq_varid(infile,"mSourcesXxS",&msXxSVar)==NC_NOERR,
         "sourceNetwork--unable to open variable mSourcesXxS in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"mSourcesYyS",&msYySVar)==NC_NOERR,
         "sourceNetwork--unable to open variable mSourcesYyS in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"mSourcesZzS",&msZzSVar)==NC_NOERR,
         "sourceNetwork--unable to open variable mSourcesZzS in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"mSourcesXyS",&msXySVar)==NC_NOERR,
         "sourceNetwork--unable to open variable mSourcesXyS in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"mSourcesXzS",&msXzSVar)==NC_NOERR,
         "sourceNetwork--unable to open variable mSourcesXzS in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"mSourcesYzS",&msYzSVar)==NC_NOERR,
         "sourceNetwork--unable to open variable mSourcesYzS in file %s(%i)",
         filename,infile);
  
  int msXyAVar,msXzAVar,msYzAVar;
  assert(nc_inq_varid(infile,"mSourcesXyA",&msXyAVar)==NC_NOERR,
         "sourceNetwork--unable to open variable mSourcesXyA in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"mSourcesXzA",&msXzAVar)==NC_NOERR,
         "sourceNetwork--unable to open variable mSourcesXzA in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"mSourcesYzA",&msYzAVar)==NC_NOERR,
         "sourceNetwork--unable to open variable mSourcesYzA in file %s(%i)",
         filename,infile);
  
  int msDataVar;
  assert(nc_inq_varid(infile,"mSourcesData",&msDataVar)==NC_NOERR,
         "sourceNetwork--unable to open variable mSourcesData in file %s(%i)",
         filename,infile);
  
  //Create a new source for each.
  for(int i=0;i<_nms;i++)
    _sources->Add(new momentSource(modelDef,filename,infile,i,
                                   msXsVar,msYsVar,msZsVar,
                                   msXxSVar,msYySVar,msZzSVar,msXySVar,msXzSVar,msYzSVar,
                                   msXyAVar,msXzAVar,msYzAVar,
                                   msSampVar,msDataVar,
                                   rho));
  
  return _nms;
}

///Read traction sources from a NetCDF file.
int sourceNetwork::readTractionSources(modelDefStruct* modelDef,
                                       char* filename,int infile,
                                       float* rho){
  int ntsDim;
  if(nc_inq_varid(infile,"numTSources",&ntsDim)!=NC_NOERR)
    return _nts=0;
  
  assert(nc_inq_dimlen(infile,ntsDim,&_nts),
         "sourceNetwork--unable read dim numTSources in %s(%i,%i)",
         filename,infile,ntsDim);
  
  int tsXsVar,tsYsVar;
  assert(nc_inq_varid(infile,"tSourcesXs",&tsXsVar)==NC_NOERR,
         "sourceNetwork--unable to open variable tSourcesXs in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"tSourcesYs",&tsYsVar)==NC_NOERR,
         "sourceNetwork--unable to open variable tSourcesYs in file %s(%i)",
         filename,infile);
  
  int tsSampVar,tsAxVar,tsAyVar,tsAzVar;
  assert(nc_inq_varid(infile,"tSourcesSamp",&tsSampVar)==NC_NOERR,
         "sourceNetwork--unable to open variable tSourcesSamp in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"tSourcesAx",&tsAxVar)==NC_NOERR,
         "sourceNetwork--unable to open variable tSourcesAx in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"tSourcesAy",&tsAyVar)==NC_NOERR,
         "sourceNetwork--unable to open variable tSourcesAy in file %s(%i)",
         filename,infile);
  assert(nc_inq_varid(infile,"tSourcesAz",&tsAzVar)==NC_NOERR,
         "sourceNetwork--unable to open variable tSourcesAz in file %s(%i)",
         filename,infile);
  
  int tsDataVar;
  assert(nc_inq_varid(infile,"tSourcesData",&tsDataVar)==NC_NOERR,
         "sourceNetwork--unable to open variable tSourcesData in file %s(%i)",
         filename,infile);
  
  for(size_t i=0;i<_nts;i++);
  
  return _nts;
}

//
//Direct wavelet generation subroutines. 
void sourceNetwork::cTimesCExp(double* inOut,double r,double c){
  //Determine the length and angle of the complex represented by inOut
  double length1=sqrt(inOut[0]*inOut[0]+inOut[1]*inOut[1]);
  double theta1=atan2(inOut[1],inOut[0]);
  
  inOut[0]=length1*exp(r)*cos(theta1+c);
  inOut[1]=length1*exp(r)*sin(theta1+c);
} 
double* sourceNetwork::rickerSpec(int nf,double df,double fmode,double tmin,
                                  double *spect){
  //All the complex parts are 0.
  for(int i=0;i<4*nf;spect[i++]=0.0);
  
  //set the real parts
  double pi2=2*PI;
  
  for(int i=0;i<=nf;i++){
    double f=df*i;
    double arg=(f/fmode);
    arg*=arg;
    
    if(arg>70.0){
      spect[2*i]=0.0;
    }else{
      spect[2*i]=arg*exp(-arg);
    }
  }
  
  //Now account for the non-zero origin time
  for(int i=0;i<nf;i++)
    cTimesCExp(&spect[2*i],0.0,pi2*df*i*tmin);
  
  //Fill in the negative frequencys with the complex conjugate values.
  for(int i=1;i<nf;i++){
    spect[2*nf+2*i]=   spect[2*nf-2*i];
    spect[2*nf+2*i+1]=-spect[2*nf-2*i+1];
  }
  
  return spect;
}

//
//And some more internal functions to calculate the dispersion limits for a given source.
int sourceNetwork::nFreqs(int nt){
  double exponent=ceil(log((float)nt)/log((float)2));
  int nsamps=(int)ceil(pow((double)2.0,exponent));
  return (int)floor((float)nsamps/2);
}
//calculate the fourier aspec of the data
// returns wrk array and allocates it as required
float* sourceNetwork::sourceSpectrum(float* data,int scaleFactor,
                                     int nt,float dt,
                                     int& nf,float& df,
                                     floatPtr& wrk){
  //fill in the parameters and allocate space if required
  nf=nFreqs(nt);
  df=dFreq(nf,dt);
  int workSize=2*(2*nf+1);
  if(!wrk)
    assert((wrk=(float*)malloc(workSize*sizeof(float)))!=NULL,
           "sourceSpectrum--unable to allocate work array (%i)",
           2*(2*nf+1));
  
  //Initialize work array.
  for(int it=0;it<workSize;wrk[it++]=0.0);
  //Retrieve force source waveform; load into every other word of
  // the real-valued work array.
  for(int it=0;it<nt;it++)
    wrk[1+2*it]=data[it];
  
  //Compute FFT of source waveform.
  NUMREC_four1(wrk,2*nf,-1);
  
  //Calculate amplitude spectrum.  Multiply spectrum by factor power of
  // frequency f to account for the far-field propagation effect. 
  for(int it=0;it<=nf;it++){
    float f=it*df;
    wrk[it]=sqrt(wrk[1+2*it]*wrk[1+2*it]+wrk[2+2*it]*wrk[2+2*it]);
    
    if(scaleFactor==1){
      wrk[it]*=f;
    }else if(scaleFactor==2){
      wrk[it]*=f*f;
    }else if(scaleFactor>2){
      for(int j=0;j<scaleFactor;j++)
        wrk[it]*=f;
    }
  }
  
  return wrk;
}
//Subroutine BANDIT determines the bandwidth of a frequency spectrum
// between amplitude levels that are a specified percentage of the 
// maximum amplitude.
float sourceNetwork::bandit(float pcent,
                            float dx,float dy,float dz,
                            float df,int nf,float vmin,
                            float *data,int doPrint){
  int if_lo=1;
  int if_hi=nf-2;
  
  //Find mode of amplitude spectrum.
  float ampmax=0.0;
  int if_mode=0;
  for(int i=0;i<nf;i++)
    SETMAXINDEX(ampmax,data[i],if_mode,i);
  
  //Compute minimum nonzero value of amplitude spectrum.
  float ampmin=(pcent/100.0)*ampmax;
  
  //Determine upper cutoff frequency.
  for(int i=nf-1;i>=0;i--){
    if (data[i]>ampmin){
      if_hi=i+1;
      break;
    }
  }
  
  //Determine lower cutoff frequency.
  for(int i=0;i<=nf;i++){
    if (data[i]>ampmin){
      if_lo=i-1;
      break;
    }
  }
  
  //Calculate low-cut, high-cut, and mode frequencies.
  float flo  =if_lo*df;
  float fhi  =if_hi*df;
  float fmode=if_mode*df;
  
  //Write diagnostic information to screen.
  tEprintf(doPrint,"  Peak frequency                     = %7.2f Hz\n",
           fmode);
  tEprintf(doPrint,"  Effective bandwidth (%5.1f%% level) = %7.2f to %7.2f Hz\n",
           pcent,flo,fhi);
  
  //Calculate minimum propagating wavelength.
  float lamdamin=vmin/fhi;
  
  //Report ratios of the minimum wavelength to grid intervals.
  tEprintf(doPrint,
           "  Min dx/wavelength=%.4f\n  Min dy/wavelength=%.4f\n  Min dz/wavelength=%.4f\n",
           lamdamin/dx,lamdamin/dy,lamdamin/dz);
  return MIN3(lamdamin/dx,lamdamin/dy,lamdamin/dz);
} 
