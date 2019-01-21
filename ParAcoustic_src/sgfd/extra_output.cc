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
 *  extra_output.cc
 *
 *
 *  Define several class functions that control slice and volume
 *  output for wavefield visualization.
 *
 *  Defines the following functions (mostly class):
 *  time (needed if _fixTimeOnMac is defined)
 *  sliceOutput::doOutput
 *  sliceOutput::doOutput
 *  sliceOutput::writeCDFSlice
 *  wavefieldOutput::wavefieldOutput
 *  wavefieldOutput::doOutput
 *  wavefieldOutput::doOutput
 *  extraOutputGroup::extraOutputGroup
 *  extraOutputGroup::doCheckpoint
 *  extraOutputGroup::readCheckpoint
 *  extraOutputGroup::correctGlobalOrigin
 *  extraOutputGroup::initOutput
 *  extraOutputGroup::doOutput
 *  extraOutputGroup::doOutput
 *  extraOutputGroup::writeOutput
 *  extraOutputGroup::readExtraOutput
 *  extraOutputGroup::addExtraOutput
 *  extraOutputGroup::addWavefieldOutput
 *  extraOutputGroup::doWavefieldOutput
 *  extraOutputGroup::initializeSliceVariables
 *  extraOutputGroup::sliceComp2Num
 *  extraOutputGroup::slicePlane2Num
 *
 */
/*These classes produce extra output.
 
 Extra output consists of (slices, full wavefield, etc???) for the
 SGFD problem. Slice output is the only one of these that has been
 extensivly tested.  The full waveform output creates large files but
 can provide some really cool visualizations.
 */

#include "extra_output.hh"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "netcdf.h"
#include "nstdutil.hh"
#include "nstdutil.h"
#include "xtrautil.hh"
#include "io_procs.h"

#include "message_passing.h"
#include "sgfdSources.hh"

#include "sgfd.hh"

#ifdef _serial_elasti_hh_
int sliceOutput::doOutput(char *outfileName,serialDependent* forwardProblem,
                          const char* varName,float* work,
                          int **sliceIndices){
  //find the correct index to write to
  int index=sliceIndices[_plane][_comp]++;
  
  forwardProblem->getSliceData(work,_plane,_comp,_coord);
  writeCDFSlice(outfileName,index,varName,
                _plane,_comp,_t,_coord,
                work,
                forwardProblem->modelDef()->procLim);
  return index;
}
#endif
int sliceOutput::doOutput(char *outfileName,masterSgfdDependent* forwardProblem,
                          const char* varName,float* work,
                          int **sliceIndices){
  //find the correct index to write to
  int index=sliceIndices[_plane][_comp]++;
  
  //Ask the slaves to send this slice.
  int numWriting=0;
  //Send a message to all procs
  initSend();
  sendMessage(AllProcesses,MESSAGE_SLICE,"iif",
              _plane,_comp,_coord);
  
  if(forwardProblem->doIHaveData()){
    //The caller also has data that must be added to the slice.
    float *sliceData;
    int currProcLim[6];
    int inLimits=forwardProblem->slicer(_plane,_comp,_coord,
                                        sliceData,currProcLim);
    if(inLimits){
      numWriting++;
      //Write the result
      writeCDFSlice(outfileName,index,varName,
                    _plane,_comp,_t,_coord,
                    sliceData,
                    currProcLim);
    }
  }
  
  for(int ii=0;ii<NumProcs;ii++){
    //receive the result
    int inLimits;
    initSend();
    sendMessage(Tids[ii],MESSAGE_SLICE,NULL);
    getMessage(Tids[ii],MESSAGE_SLICE,"i",
               &inLimits);
    
    if(inLimits){
      int currProcLim[6];
      unpackMessage("iiiiii",
                    &currProcLim[0],&currProcLim[1],
                    &currProcLim[2],&currProcLim[3],
                    &currProcLim[4],&currProcLim[5]);
      
      numWriting++;
      //calculate the result length
      int resultLength;
      if(_plane==SLICE_YZ_PLANE){
        resultLength=
        (currProcLim[3]-currProcLim[2])*
        (currProcLim[5]-currProcLim[4]);
      }else if(_plane==SLICE_XZ_PLANE){
        resultLength=
        (currProcLim[1]-currProcLim[0])*
        (currProcLim[5]-currProcLim[4]);
      }else if(_plane==SLICE_XY_PLANE){
        resultLength=
        (currProcLim[1]-currProcLim[0])*
        (currProcLim[3]-currProcLim[2]);
      }
      assert(resultLength>0,
             "Received bad procLim: %i-%i, %i-%i, %i-%i from proc %i",
             currProcLim[0],currProcLim[1],
             currProcLim[2],currProcLim[3],
             currProcLim[4],currProcLim[5],
             ii);
      
      //unpack and write the result
      unpackMessage("F",work,resultLength);
      writeCDFSlice(outfileName,index,varName,
                    _plane,_comp,_t,_coord,
                    work,
                    currProcLim);
    }
  }
  assert(numWriting,
         "sliceOutput::doOutput--no procs wrote to slice %s[%i], %.3f",
         varName,index,_coord);
  return index;
}
int sliceOutput::writeCDFSlice(char* outfileName,int sliceNum,const char* varName,
                               int plane,int comp,float time,float position,
                               float* data,
                               int procLim[6]){
  //Make sure filename ends in .cdf
  char buffer[512];
  if(!strcmp(outfileName+strlen(outfileName)-4,".cdf")){
    sprintf(buffer,"%s",outfileName);
  }else{
    sprintf(buffer,"%s.cdf",outfileName);
  }
  
  //open the file
  int outFile;
  assert(nc_open(buffer,NC_WRITE,&outFile)==NC_NOERR,
         "writeCDFSlice--unable to open file %s\n",buffer);
  
  //Set the start and count and variable for the current slice
  // Also get the variables to for the current slice type count
  // and location
  size_t start[3]={sliceNum},count[3]={1};
  int currSliceTypeVar;
  int currSliceTimeVar,currSlicePosVar;
  
  //Identify the plane, use a switch here since until I actually add some new
  // type of plane?.
  switch(plane){
    case SLICE_YZ_PLANE:
      start[2]=procLim[2];
      start[1]=procLim[4];
      count[2]=procLim[3]-procLim[2];
      count[1]=procLim[5]-procLim[4];
      break;
    case SLICE_XZ_PLANE:
      start[2]=procLim[0];
      start[1]=procLim[4];
      count[2]=procLim[1]-procLim[0];
      count[1]=procLim[5]-procLim[4];
      break;
    case SLICE_XY_PLANE:
      start[2]=procLim[0];
      start[1]=procLim[2];
      count[2]=procLim[1]-procLim[0];
      count[1]=procLim[3]-procLim[2];
      break;
    default:
      assert(FALSE,"writeCDFSlice--unknown plane %i for slice %i",
             plane,sliceNum);
  }
  //Identify the component, we can just index to the comp name and generate the
  // required variable names on the fly.
  char timeVarName[128],posVarName[128];
  sprintf(timeVarName,"%sTime",varName);
  sprintf(posVarName,"%sPos",varName);
  assert(nc_inq_varid(outFile,varName,&currSliceTypeVar)==NC_NOERR,
         "writeCDFSlice--unable to open variable %s in file %s",
         varName,buffer);
  assert(nc_inq_varid(outFile,timeVarName,&currSliceTimeVar)==NC_NOERR,
         "writeCDFSlice--unable to open variable %s in file %s",
         timeVarName,buffer);
  assert(nc_inq_varid(outFile,posVarName,&currSlicePosVar)==NC_NOERR,
         "writeCDFSlice--unable to open variable %s in file %s",
         posVarName,buffer);
  
  //now write the data
  assert(nc_put_vara_float(outFile,currSliceTypeVar,start,count,data)==NC_NOERR,
         "writeCDFSlice--unable to write (%i,%i,%i) block to %s at (%i,%i,%i)",
         count[0],count[1],count[2],
         buffer,
         start[0],start[1],start[2]);
  
  assert(nc_put_var1_float(outFile,currSliceTimeVar,start,&time)==NC_NOERR,
         "writeCDFSlice--unable to write time to %s(%i,%i) at index %i",
         buffer,outFile,currSliceTimeVar,start[0]);
  assert(nc_put_var1_float(outFile,currSlicePosVar,start,&position)==NC_NOERR,
         "writeCDFSlice--unable to write position to %s(%i,%i) at index %i",
         buffer,outFile,currSlicePosVar,start[0]);
  
  //close the file
  assert(nc_close(outFile)!=-1,
         "writeCDFSlice--nc_close of %s unsuccessful",buffer);
  return count[0]*count[1]*count[2];
}

wavefieldOutput::wavefieldOutput(float t,float fieldScalar,
                                 int nVars,...)
:extraOutput(t){
  _outfileName=NULL;
  _fieldScalar=fieldScalar;
  
  _numVars=nVars;
  
  va_list args;
  va_start(args,nVars);
  for(int i=0;i<nVars;i++){
    assert((_varNames[i]=(char*)malloc(512*sizeof(char)))!=NULL,
           "wavefieldOutput--unable to allocate %i chars for variable name %i",
           512,i);
    char* buffer=va_arg(args,char*);
    strcpy(_varNames[i],buffer);
  }
  va_end(args);
  
  if(_numVars){
    tEprintf(Verbose,"New wavefield output for %i variables:\n\t",_numVars);
    for(int i=0;i<_numVars;i++)
      tEprintf(Verbose,"%s ",_varNames[i]);
    tEprintf(Verbose,"\n");
  }
}

#ifdef _serial_elasti_hh_
int wavefieldOutput::doOutput(char *outfileName,serialDependent* forwardProblem,
                              char* varName,float* work,
                              int **sliceIndices){
  DEF_MODEL_SIZE(forwardProblem->modelDef());
  DEF_MODEL_LIMITS(forwardProblem->modelDef());
  //Make an array of the required variables
  if(!_numVars)
    return; //nothing to do, return
  
  //Write a cdf header
  writeCDFHeader(outfileName,
                 NX,NY,NZ,NT,NULL,
                 dx,minX,dy,minY,dz,minZ,dt,minT,
                 _numVars,_varNames);
  
  //write the current iteration and time
  int iteration=(int)floor(t()/dt);
  float time=dt*iteration;
  writeCDFVarFloat(outfileName,"currTime",time);
  
  //Now open my local copy of the file.
  int outFile=openCDFFile(outfileName,FALSE,NC_WRITE);
  
  //Write the variables to the cdf file maintain space for the result.
  for(int i=0;i<_numVars;i++){
    int varID;
    assert(nc_inq_varid(outFile,_varNames[i],&varID)==NC_NOERR,
           "doOutput--unable to get var %s from %s(%i)",
           _varNames[i],outfileName,outFile);
    
    float *var;
    if(!strcmp("vx",_varNames[i])){
      var=forwardProblem->_vx;
    }else if(!strcmp("vy",_varNames[i])){
      var=forwardProblem->_vy;
    }else if(!strcmp("vz",_varNames[i])){
      var=forwardProblem->_vz;
    }else if(!strcmp("xx",_varNames[i])){
      var=forwardProblem->_xx;
    }else if(!strcmp("yy",_varNames[i])){
      var=forwardProblem->_yy;
    }else if(!strcmp("zz",_varNames[i])){
      var=forwardProblem->_zz;
    }else if(!strcmp("xy",_varNames[i])){
      var=forwardProblem->_xy;
    }else if(!strcmp("xz",_varNames[i])){
      var=forwardProblem->_xz;
    }else if(!strcmp("yz",_varNames[i])){
      var=forwardProblem->_yz;
    }else{
      assert(FALSE,"doOutput--unknown var name %s[%i]",
             _varNames[i],i);
    }
    
    assert(nc_put_var_float(outFile,varID,var)==NC_NOERR,
           "doOutput--unable to write %s %s(%i,%i)",
           _varNames[i],
           outfileName,outFile,varID);
  }
  
  assert(nc_close(outFile)==NC_NOERR,
         "doOutput--unable to close %s(%i)",outfileName,outFile);
  
  tEprintf(Verbose,"  Wrote %i wavefield variables:\n\t",
           _numVars);
  for(int i=0;i<_numVars;i++)
    tEprintf(Verbose,"%s ",_varNames[i]);
  tEprintf(Verbose,"\n\tto file: %s\n",outfileName);
  
  return TRUE;
}
#endif
int wavefieldOutput::doOutput(char *outfileName,masterSgfdDependent* forwardProblem,
                              const char* varName,float* work,
                              int **sliceIndices){
  DEF_MODEL_SIZE(forwardProblem->modelDef());
  DEF_MODEL_LIMITS(forwardProblem->modelDef());
  
  if(_outfileName)
    outfileName=_outfileName;
  
  //Make an array of the required variables
  if(!_numVars)
    return FALSE; //nothing to do, return
  
  //Write a cdf header.
  writeCDFHeader(outfileName,
                 NX,NY,NZ,NT,NULL,
                 dx,minX,dy,minY,dz,minZ,dt,minT,
                 _numVars,_varNames);
  
  //Write the current iteration and time.
  int iteration=(int)floor(t()/dt);
  float time=dt*iteration;
  writeCDFVarFloat(outfileName,"currTime",time);
  
  //Now open my local copy of the file.
  int outFile=openCDFFile(outfileName,FALSE,NC_WRITE);
  
  //Write the variables to the cdf file maintain space for the result.
  for(int i=0;i<_numVars;i++){
    int varID;
    assert(nc_inq_varid(outFile,_varNames[i],&varID)==NC_NOERR,
           "doOutput--unable to get var %s from %s(%i)",
           _varNames[i],outfileName,outFile);
    
    for(int j=0;j<NumProcs;j++){
      //Send a message to the current process to send this variable. Note the new
      // format for the send, each layer will be packed individually.
      initSend();
      sendMessage(Tids[j],MESSAGE_FULL_FIELD,
                  "s",_varNames[i]);
      
      //Receive and write the results
      int currProcLim[6];
      getMessage(Tids[j],MESSAGE_FULL_FIELD,"P",currProcLim);
      
      int currNXY=
      (currProcLim[1]-currProcLim[0])*
      (currProcLim[3]-currProcLim[2]);
      int startZ=currProcLim[4],stopZ=currProcLim[5];
      
      //Set the unchanging part of the start and count.
      size_t start[3],count[3];
      start[1]=currProcLim[2];
      start[2]=currProcLim[0];
      
      count[0]=1;
      count[1]=currProcLim[3]-currProcLim[2];
      count[2]=currProcLim[1]-currProcLim[0];
      int totalCount=count[0]*count[1]*count[2];
      
      for(int k=startZ;k<stopZ;k++){
        //Unpack and write the result ONE LAYER AT A TIME.
        getMessage(Tids[j],MESSAGE_FULL_FIELD,"F",work,currNXY);
        
        for(int iii=0;iii<totalCount;work[iii++]*=_fieldScalar);
        
        start[0]=k;
        assert(nc_put_vara_float(outFile,varID,start,count,work)==NC_NOERR,
               "doOutput--unable to write %s[%i,%i,%i;%i,%i,%i] %s(%i,%i)",
               _varNames[i],
               start[0],start[1],start[2],count[0],count[1],count[2],
               outfileName,outFile,varID);
      }
    }
  }
  
  assert(nc_close(outFile)==NC_NOERR,
         "doOutput--unable to close %s(%i)",outfileName,outFile);
  
  tEprintf(Verbose,"  Wrote %i wavefield variables:\n\t",
           _numVars);
  for(int i=0;i<_numVars;i++)
    tEprintf(Verbose,"%s ",_varNames[i]);
  tEprintf(Verbose,"\n\tto file: %s\n",outfileName);
  return TRUE;
}

/*! \brief Here is container class that holds all of the extra output.
 
 This is the class that the program will allocate and use to check if it is time
 to create some extra output.
 */
//Initializer.
extraOutputGroup::extraOutputGroup(modelDefStruct* modelDef){
  //Initialize the array.
  _theOutput=new extraOutputArray;
  _nextOutput=0;
  _useModelExtraOutput=TRUE;
  
  //Set default names.
  strcpy(_sliceOutputName,"slice.cdf");
  strcpy(_wavefieldOutputName,"wavefield.cdf");
  _maxNumPerSliceFile = 1000000000LL;
  
  //Set default slice types.
  _sliceCompNames=new stringArray;
  addNewSliceComp("");
  addNewSliceComp("Vx");
  addNewSliceComp("Vy");
  addNewSliceComp("Vz");
  addNewSliceComp("Pressure");
  
  addNewSliceComp("Vp");
  addNewSliceComp("Vs");
  addNewSliceComp("Rho");
  
  addNewSliceComp("Rx");
  addNewSliceComp("Ry");
  addNewSliceComp("Rz");
  
  //Set default slice planes.
  _slicePlaneNames=new stringArray;
  addNewSlicePlane("");
  addNewSlicePlane("yz");
  addNewSlicePlane("xz");
  addNewSlicePlane("xy");
}

///Method to write restart information required to correctly restart the extra
/// output.
FILE* extraOutputGroup::doCheckpoint(int iteration,char* cpDir,int cpID,
                                     int callDepth,
                                     FILE* cpFile){
  //Since this is the base the file has not been open yet,
  // open the binary file that will hold all the info need for
  // the restart.
  char filename[1024];
  sprintf(filename,"%s/checkpoint_extra_output_%i.bin",cpDir,cpID);
  cpFile=fopen(filename,"wb");
  assert(cpFile!=NULL,
         "extraOutputGroup::doCheckpoint--unable to open file\n\t%s\n",
         filename);
  
  fwrite(&_nextOutput,sizeof(int),1,cpFile);
  
  for(int j=1;j<=numPlanes();j++)
    fwrite(_sliceIndices[j],sizeof(int),numComps(),cpFile);
  
  if(!callDepth){
    fclose(cpFile);
    cpFile=NULL;
  }
  return cpFile;
}
///Here is a method to read a checkpoint file written by doCheckpoint.
FILE* extraOutputGroup::readCheckpoint(masterSgfdModel* model,
                                       int& iteration,char* cpDir,int cpID,
                                       int callDepth,
                                       FILE* cpFile){
  //Since this is the base the file has not been open yet,
  // open the binary file that will hold all the info need for
  // the restart.
  char filename[1024];
  sprintf(filename,"%s/checkpoint_extra_output_%i.bin",cpDir,cpID);
  cpFile=fopen(filename,"rb");
  assert(cpFile!=NULL,
         "extraOutputGroup::readCheckpoint--unable to open file\n\t%s\n",
         filename);
  
  
  initializeSliceVariables(model,TRUE);
  
  fread(&_nextOutput,sizeof(int),1,cpFile);
  
  for(int j=1;j<=numPlanes();j++)
    fread(_sliceIndices[j],sizeof(int),numComps(),cpFile);
  
  if(!callDepth){
    fclose(cpFile);
    cpFile=NULL;
  }
  return cpFile;
}

//Initialize output files (actully only the slice file needs to be initialized).
int extraOutputGroup::initOutput(masterSgfdModel* model,
                                 sourceNetwork* sources,int allocateStorage,
                                 int fileInit){
  modelDefStruct* modelDef = model->modelDef();
  DEF_MODEL_SIZE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  
  //Initialize the counters for indexing each slice type.
  initializeSliceVariables(model,allocateStorage);
  
  if(!fileInit)
    return FALSE;
  
  if(Verbose>1){
    tEprintf(Verbose,"Extra Output\n");
    for(int i=0;i<_theOutput->size();i++){
      extraOutput* curr=(*_theOutput)[i];
      switch(curr->type()){
        case SLICE_OUTPUT_TYPE:
        {
          sliceOutput* currSlice=(sliceOutput*)curr;
          tEprintf(Verbose,
                   " #%3i: slice, %s%s: t %.2fms\n",
                   i,compText(currSlice->comp()),
                   planeText(currSlice->plane()),
                   1000.0*currSlice->t());
        }
          break;
        case WAVEFIELD_OUTPUT_TYPE:
          tEprintf(Verbose,
                   " #%3i: wavefield: t %.2fms\n",
                   i,1000.0*curr->t());
          break;
        default:
          assert(FALSE,"initOutput--unknown output type %i for output %i",
                 curr->type(),i);
      }
    }
  }
  
  //Loop through all the output points and count different types of slices,
  // there are 12 possible types.
  //      size_t numSliceTypes[numPlanes()][numComps()];
  size_t **numSliceTypes;
  size_t **maxNumSliceTypes;
  assert((numSliceTypes=(size_t**)malloc(numPlanes()*sizeof(size_t*)))!=NULL,
         "initOutput--unable to allocate %i size** for slicetypes",
         numPlanes());
  assert((maxNumSliceTypes=(size_t**)malloc(numPlanes()*sizeof(size_t*)))!=NULL,
         "initOutput--unable to allocate %i size** for slicetypes",
         numPlanes());
  for(int i=0;i<numPlanes();i++) {
    assert((numSliceTypes[i]=(size_t*)malloc(numComps()*sizeof(size_t)))!=
           NULL,
           "initOutput--unable to allocate %i size_t for slicetypes[%i]",
           numComps(),i);
    assert((maxNumSliceTypes[i]=(size_t*)malloc(numComps()*sizeof(size_t)))!=
           NULL,
           "initOutput--unable to allocate %i size_t for slicetypes[%i]",
           numComps(),i);
  }
  for(int i=0;i<numPlanes();i++){
    for(int j=0;j<numComps();j++){
      numSliceTypes[i][j]=0;
      maxNumSliceTypes[i][j]=0;
    }
  }
  //determine the total number of slices requested for each type
  int totalSlices=0;
  for(int i=0;i<_theOutput->size();i++){
    extraOutput* curr=(*_theOutput)[i];
    if(curr->type()!=SLICE_OUTPUT_TYPE) continue;
    
    //numSliceTypes[curr->plane()][curr->comp()]++; //this is needed if C1 below is uncommented
    totalSlices++;
  }
  //determine the number of output files needed for the requested slices
  int numFiles = 1;
  for(int i=0;i<numPlanes();i++){
    long long planeSize;
    switch (i) {
      case SLICE_YZ_PLANE:
        planeSize=NY*NZ;
        break;
      case SLICE_XZ_PLANE:
        planeSize=NX*NZ;
        break;
      case SLICE_XY_PLANE:
        planeSize=NX*NY;
        break;
      default:
        continue;
    }
    int maxSlices = static_cast<int>(_maxNumPerSliceFile/planeSize);
    for(int j=0;j<numComps();j++){
      maxNumSliceTypes[i][j] = maxSlices;
    }
  }
  for(int ii=0;ii<_theOutput->size();ii++){
    extraOutput* curr=(*_theOutput)[ii];
    if(curr->type()!=SLICE_OUTPUT_TYPE) continue;
    
    if(numSliceTypes[curr->plane()][curr->comp()]+1>maxNumSliceTypes[curr->plane()][curr->comp()]) {
      ++numFiles;
      for(int i=0;i<numPlanes();i++){
        for(int j=0;j<numComps();j++){
          numSliceTypes[i][j]=0;
        }
      }
    }
    numSliceTypes[curr->plane()][curr->comp()]++;
  }
  
  if(!totalSlices)
    return 0;
  int totalSlices0=totalSlices;
  for(int i=0;i<numPlanes();i++){
    for(int j=0;j<numComps();j++){
      numSliceTypes[i][j]=0;
    }
  }
  
  //There are slices so write the output file.
  float *data=(float*)malloc(MAX3(NX,NY,NZ)*sizeof(float));
  char sliceOutputName0[1024];
  char ending[256] = "";
  int beginI = 0;
  _sliceFileIndices = new int[numFiles];
  totalSlices=0;
  if(!strcmp(_sliceOutputName+strlen(_sliceOutputName)-4,".cdf")) {
    strncpy(sliceOutputName0,_sliceOutputName,strlen(_sliceOutputName)-4);
    sliceOutputName0[strlen(_sliceOutputName)-4]='\0';
  } else
    strcpy(sliceOutputName0,_sliceOutputName);
  for(int ii=0;ii<numFiles;++ii) {
    for(int i=beginI;i<_theOutput->size();++i,++beginI){
      extraOutput* curr=(*_theOutput)[i];
      if(curr->type()!=SLICE_OUTPUT_TYPE) continue;
      if((numSliceTypes[curr->plane()][curr->comp()]+1)>maxNumSliceTypes[curr->plane()][curr->comp()])
        break;
      totalSlices++;
      numSliceTypes[curr->plane()][curr->comp()]++;
    }
    _sliceFileIndices[ii] = beginI;
    if(numFiles>1)
      sprintf(ending,"_%d",ii);
    sprintf(_sliceOutputName,"%s%s.cdf",sliceOutputName0,ending);
    //Create the file. Define standard dimensions, coordinates, and attributes.
    int outFile=openCDFFile(_sliceOutputName,TRUE,NC_WRITE);
    //Global attributes
    char buffer[1024];
    sprintf(buffer,"parallel_elasti slice file");
    nc_put_att_text(outFile,NC_GLOBAL,"title",strlen(buffer)+1,buffer);
    nc_put_att_text(outFile,NC_GLOBAL,"history",strlen(CommandLine)+1,CommandLine);
    
    //Define dimensions and variables for coordinate vectors.
    int nxDim,nyDim,nzDim;
    assert(nc_def_dim(outFile,"NX",(size_t)NX,&nxDim)==NC_NOERR,
           "initOutput--unable to define dim NX in %s(%i)",
           _sliceOutputName,outFile);
    assert(nc_def_dim(outFile,"NY",(size_t)NY,&nyDim)==NC_NOERR,
           "initOutput--unable to define dim NY in %s(%i)",
           _sliceOutputName,outFile);
    assert(nc_def_dim(outFile,"NZ",(size_t)NZ,&nzDim)==NC_NOERR,
           "initOutput--unable to define dim NZ in %s(%i)",
           _sliceOutputName,outFile);
    
    int xVar,yVar,zVar;
    assert(nc_def_var(outFile,"x",NC_FLOAT,1,&nxDim,&xVar)==NC_NOERR,
           "initOutput--unable to define var x in %s(%i,%i)",
           _sliceOutputName,outFile,nxDim);
    assert(nc_def_var(outFile,"y",NC_FLOAT,1,&nyDim,&yVar)==NC_NOERR,
           "initOutput--unable to define var y in %s(%i,%i)",
           _sliceOutputName,outFile,nyDim);
    assert(nc_def_var(outFile,"z",NC_FLOAT,1,&nzDim,&zVar)==NC_NOERR,
           "initOutput--unable to define var z in %s(%i,%i)",
           _sliceOutputName,outFile,nzDim);
    
    //define a dimension and variable for each used type of slice
    for(int i=1;i<numPlanes();i++){
      int xDim,yDim;
      switch(i){
        case SLICE_YZ_PLANE:
          xDim=nyDim;
          yDim=nzDim;
          break;
        case SLICE_XZ_PLANE:
          xDim=nxDim;
          yDim=nzDim;
          break;
        case SLICE_XY_PLANE:
          xDim=nxDim;
          yDim=nyDim;
          break;
      }
      
      for(int j=1;j<numComps();j++){
        if(!numSliceTypes[i][j]) continue;
        
        //Generate the names for the dimension and variables.
        char dimName[64],tName[64],pName[64],varName[64];
        sprintf(varName,"%s%s",
                planeText(i),compText(j));
        sprintf(dimName,"%sDim",varName);
        sprintf(tName,"%sTime",varName);
        sprintf(pName,"%sPos",varName);
        
        //Create the
        int dim,var,tVar,pVar;
        assert(nc_def_dim(outFile,dimName,
                          numSliceTypes[i][j],&dim)==NC_NOERR,
               "initOutput--unable to define dim %s[%i] in %s(%i)",
               dimName,numSliceTypes[i][j],_sliceOutputName,outFile);
        
        //Create variables for the slice times and positions
        assert(nc_def_var(outFile,tName,NC_FLOAT,1,&dim,&tVar)==NC_NOERR,
               "initOutput--unable to define var %s in %s(%i)",
               tName,_sliceOutputName,outFile);
        assert(nc_def_var(outFile,pName,NC_FLOAT,1,&dim,&pVar)==NC_NOERR,
               "initOutput--unable to define var %s in %s(%i)",
               pName,_sliceOutputName,outFile);
        
        //Create the variable for the data
        int dims[3]={dim,yDim,xDim};
        assert(nc_def_var(outFile,varName,NC_FLOAT,3,dims,&var)==NC_NOERR,
               "initOutput--unable to define var %s in %s(%i)",
               varName,_sliceOutputName,outFile);
      }
    }
    
    //Now fill in the data variables.
    assert(nc_enddef(outFile)==NC_NOERR,
           "initOutput--unable to take %s(%i) out of define mode",
           _sliceOutputName,outFile);
    size_t start=0;
    
    for(int i=0;i<NX;++i) data[i]=minX+i*dx;
    size_t count=NX;
    assert(nc_put_vara_float(outFile,xVar,&start,&count,data)==NC_NOERR,
           "initOutput--unable to fill var x in %s(%i,%i)",
           _sliceOutputName,outFile,xVar);
    
    for(int i=0;i<NY;++i) data[i]=minY+i*dy;
    count=NY;
    assert(nc_put_vara_float(outFile,yVar,&start,&count,data)==NC_NOERR,
           "initOutput--unable to fill var y in %s(%i,%i)",
           _sliceOutputName,outFile,yVar);
    
    for(int i=0;i<NZ;++i) data[i]=minZ+i*dz;
    count=NZ;
    assert(nc_put_vara_float(outFile,zVar,&start,&count,data)==NC_NOERR,
           "initOutput--unable to fill var z in %s(%i,%i)",
           _sliceOutputName,outFile,zVar);
    
    //close the file
    assert(nc_close(outFile)==NC_NOERR,
           "initOutput--unable to close %s(%i)",
           _sliceOutputName,outFile);
    
    if(sources &&  sources->size() && sources->_bandCheckPCent>0)
      sources->writeSources(modelDef,_sliceOutputName);
    
    tEprintf(Verbose,"Wrote header for %i slices to\n\t%s\n",
             totalSlices,_sliceOutputName);
    for(int i=0;i<numPlanes();i++){
      for(int j=0;j<numComps();j++){
        numSliceTypes[i][j]=0;
      }
    }
    totalSlices=0;
  }
  strcpy(_sliceOutputName,sliceOutputName0);
  _currFileIndex = 0;
  _numSliceFiles = numFiles;
  for(int i=0;i<numPlanes();i++) {
    free(numSliceTypes[i]);
    free(maxNumSliceTypes[i]);
  }
  free(numSliceTypes);
  free(maxNumSliceTypes);
  free(data);
  return totalSlices0;
}

//Check if extra output is required at the current iteration.
#ifdef _serial_elasti_hh_
int extraOutputGroup::doOutput(serialDependent* forwardProblem,
                               float t_vel,float t_press){
  while(_nextOutput<_theOutput->size()&&
        (*_theOutput)[_nextOutput]->t()<=t_press){
    extraOutput* curr=(*_theOutput)[_nextOutput++];
    char* filename;
    switch(curr->type()){
      case SLICE_OUTPUT_TYPE:
        filename=_sliceOutputName;
        break;
      case WAVEFIELD_OUTPUT_TYPE:
        filename=_wavefieldOutputName;
        break;
      default:
        assert(FALSE,"doOutput--unable to process output %i: type %i",
               _nextOutput-1,curr->type());
    }
    curr->doOutput(filename,forwardProblem,
                   _work,_sliceIndices);
  }
  return _nextOutput;
}
#endif
int extraOutputGroup::doOutput(masterSgfdDependent* forwardProblem,
                               float t_vel,float t_press){
  while(_nextOutput<_theOutput->size()&&
        (*_theOutput)[_nextOutput]->t()<=t_press){
    extraOutput* curr=(*_theOutput)[_nextOutput++];
    char filename[1024],ending[128]="",varName[128]="";
    switch(curr->type()){
      case SLICE_OUTPUT_TYPE:
        if(_nextOutput>_sliceFileIndices[_currFileIndex]) {
          ++_currFileIndex;
          for(int j=1;j<=numPlanes();j++){
            for(int i=0;i<=numComps();i++){
              _sliceIndicesGlobal[j][i]+=_sliceIndices[j][i];
              _sliceIndices[j][i]=0;
            }
          }
        }
        if(_numSliceFiles>1)
          sprintf(ending,"_%d",_currFileIndex);
        sprintf(filename,"%s%s.cdf",_sliceOutputName,ending);
        sprintf(varName,"%s%s",
                planeText(curr->plane()),compText(curr->comp()));
        break;
      case WAVEFIELD_OUTPUT_TYPE:
        strcpy(filename,_wavefieldOutputName);
        sprintf(varName,"%s",compText(curr->comp()));
        break;
      default:
        assert(FALSE,"doOutput--unable to process output %i: type %i",
               _nextOutput-1,curr->type());
    }
    
    int num=
    curr->doOutput(filename,forwardProblem,varName,
                   _work,_sliceIndices);
    
    if(curr->type()==SLICE_OUTPUT_TYPE)
      tEprintf(Verbose,
               "  Generated output for slice %s%s %i\n",
               planeText(curr->plane()),compText(curr->comp()),num+_sliceIndicesGlobal[curr->plane()][curr->comp()]);
    
  }
  return _nextOutput;
}

//Write the extra output to a file for later reading, this should match the
// format that is read by readExtraOutput.
int extraOutputGroup::writeOutput(char* fileName){
  //Only slices are written, count the number of slices.
  int nsl=0;
  for(int i=0;_theOutput && i<_theOutput->size();i++){
    if((*_theOutput)[i]->type()==SLICE_OUTPUT_TYPE)
      nsl++;
  }
  if(!nsl) return FALSE;
  
  //Open the file for writing.
  int outFile=openCDFFile(fileName,FALSE,NC_WRITE);
  
  //Define a dimension for the number of slices.
  assert(nc_redef(outFile)==NC_NOERR,
         "writeOutput--unable to put %s(%i) into define mode",
         outFile,fileName);
  int nSlicesDim;
  assert(nc_def_dim(outFile,"numSlices",nsl,&nSlicesDim)==NC_NOERR,
         "writeOutput--unable to define dim numSlices[%i] in %s(%i)",
         nsl,outFile,fileName);
  
  //Define the required variables.
  int sliceTimeVar,sliceCompVar,slicePlaneVar,sliceCoordVar;
  assert(nc_def_var(outFile,"sliceTime",NC_FLOAT,1,&nSlicesDim,&sliceTimeVar)==NC_NOERR,
         "writeOutput--unable to define var sliceTime in %s(%i,%i)",
         nsl,outFile,fileName,nSlicesDim);
  assert(nc_def_var(outFile,"sliceComp",NC_INT,1,&nSlicesDim,&sliceCompVar)==NC_NOERR,
         "writeOutput--unable to define var sliceComp in %s(%i,%i)",
         nsl,outFile,fileName,nSlicesDim);
  assert(nc_def_var(outFile,"slicePlane",NC_INT,1,&nSlicesDim,&slicePlaneVar)==NC_NOERR,
         "writeOutput--unable to define var slicePlane in %s(%i,%i)",
         nsl,outFile,fileName,nSlicesDim);
  assert(nc_def_var(outFile,"sliceCoord",NC_FLOAT,1,&nSlicesDim,&sliceCoordVar)==NC_NOERR,
         "writeOutput--unable to define var sliceCoord in %s(%i,%i)",
         nsl,outFile,fileName,nSlicesDim);
  assert(nc_enddef(outFile)==NC_NOERR,
         "writeOutput--unable to take %s(%i) out of define mode",
         outFile,fileName);
  
  size_t index=0;
  for(int i=0;i<_theOutput->size();i++){
    extraOutput* curr=(*_theOutput)[i];
    if(curr->type() != SLICE_OUTPUT_TYPE) continue;
    
    float temp=curr->outputTime();
    assert(nc_put_var1_float(outFile,sliceTimeVar,&index,&temp)==NC_NOERR,
           "writeOutput--unable to write var sliceTime[%i] in %s(%i,%i)",
           index,outFile,fileName,sliceTimeVar);
    int temp2=curr->comp();
    assert(nc_put_var1_int(outFile,sliceCompVar,&index,&temp2)==NC_NOERR,
           "writeOutput--unable to write var sliceComp[%i] in %s(%i,%i)",
           index,outFile,fileName,sliceCompVar);
    temp2=curr->plane();
    assert(nc_put_var1_int(outFile,slicePlaneVar,&index,&temp2)==NC_NOERR,
           "writeOutput--unable to write var slicePlane[%i] in %s(%i,%i)",
           index,outFile,fileName,slicePlaneVar);
    temp=curr->coord();
    assert(nc_put_var1_float(outFile,sliceCoordVar,&index,&temp)==NC_NOERR,
           "writeOutput--unable to write var sliceCoord[%i] in %s(%i,%i)",
           index,outFile,fileName,sliceCoordVar);
    
    index++;
  }
  
  assert(nc_close(outFile)==NC_NOERR,
         "writeOutput--unable to close %s(%i)",
         outFile,fileName);
  return nsl;
}

//Read extra output from a cdf file (model definition, this is the old way
// of defining slices while generating the model).
int extraOutputGroup::readExtraOutput(char* infileName){
  int inFile=openCDFFile(infileName,FALSE,NC_NOWRITE);
  
  //Check for a numSlices dimension, if it exists add slices.
  int nSlicesDim;
  if(nc_inq_dimid(inFile,"numSlices",&nSlicesDim)==NC_NOERR){
    size_t nSlices;
    nc_inq_dimlen(inFile,nSlicesDim,&nSlices);
    
    tEprintf(Verbose,"Loading %i slices from\n\t%s\n",
             nSlices,infileName);
    
    int sliceTimeVar,sliceCompVar,slicePlaneVar,sliceCoordVar;
    assert(nc_inq_varid(inFile,"sliceTime",&sliceTimeVar)==NC_NOERR,
           "readExtraOutput--unable to open var sliceTime in %s(%i)",
           infileName,inFile);
    assert(nc_inq_varid(inFile,"sliceComp",&sliceCompVar)==NC_NOERR,
           "readExtraOutput--unable to open var sliceComp in %s(%i)",
           infileName,inFile);
    assert(nc_inq_varid(inFile,"slicePlane",&slicePlaneVar)==NC_NOERR,
           "readExtraOutput--unable to open var slicePlane in %s(%i)",
           infileName,inFile);
    assert(nc_inq_varid(inFile,"sliceCoord",&sliceCoordVar)==NC_NOERR,
           "readExtraOutput--unable to open var sliceCoord in %s(%i)",
           infileName,inFile);
    
    for(size_t i=0;i<nSlices;i++){
      float t,coord;
      int comp,plane;
      assert(nc_get_var1_float(inFile,sliceTimeVar,&i,&t)==NC_NOERR,
             "readExtraOuput--unable to read var sliceTime from %s(%i,%i)",
             infileName,inFile,sliceTimeVar);
      assert(nc_get_var1_int(inFile,sliceCompVar,&i,&comp)==NC_NOERR,
             "readExtraOuput--unable to read var sliceComp from %s(%i,%i)",
             infileName,inFile,sliceCompVar);
      assert(nc_get_var1_int(inFile,slicePlaneVar,&i,&plane)==NC_NOERR,
             "readExtraOuput--unable to read var slicePlane from %s(%i,%i)",
             infileName,inFile,slicePlaneVar);
      assert(nc_get_var1_float(inFile,sliceCoordVar,&i,&coord)==NC_NOERR,
             "readExtraOuput--unable to read var sliceCoord from %s(%i,%i)",
             infileName,inFile,sliceCoordVar);
      _theOutput->Add(new sliceOutput(t,comp,plane,coord));
    }
  }
  assert(nc_close(inFile)==NC_NOERR,
         "readExtraOutput--unable to close %s(%i)",
         infileName,inFile);
  return _theOutput->size();
}

///Parse requests for additional output from the command line, note that some
/// requests can not be honored unless this is a run call and not a generate
/// model call.
int extraOutputGroup::addExtraOutput(int& i,int argOffset,
                                     int argc,char* argv[],
                                     modelDefStruct* modelDef,
                                     int runCall){
  DEF_MODEL_SIZE(modelDef);
  DEF_MODEL_LIMITS(modelDef);

  if(strlen(argv[i])<argOffset+1 || argv[i][argOffset]=='o'){
    assert(runCall,"addExtraOutput--%s flag only valid at runtime",argv[i]);
    assert(argc>i+1,"addExtraOutput--%s flag requires filename",argv[i]);
    if(strlen(argv[i])<argOffset+2){
      strcpy(_sliceOutputName,argv[++i]);
      tEprintf(Verbose,"Slice output to \"%s\"\n",
               _sliceOutputName);
    }else{
      strcpy(_wavefieldOutputName,argv[++i]);
      tEprintf(Verbose,"Wavefield output to \"%s\"\n",
               _wavefieldOutputName);
    }
  }else{
    switch(argv[i][argOffset]){
      case 'M': //Use extra output defined in the model file.
        _useModelExtraOutput=!_useModelExtraOutput;
        tEprintf(Verbose,"%s slice output from model file\n",
                 _useModelExtraOutput?"Using":"Not using");
        break;
      case 's': //Add a single slice.
        assert(argc>i+4,"addExtraOutput--adding slice requires 4 arguments");
      {
        float t=atof(argv[++i]);
        int comp=sliceComp2Num(argv[++i]);
        int plane=slicePlane2Num(argv[++i]);
        float coord=atof(argv[++i]);
        makeNewSlice(t,comp,plane,coord,TRUE);
      }
        break;
      case 'n': //Add n slices evenly distributed in time.
        assert(argc>i+4,"addExtraOutput--adding n slices requires 4 arguments");
      {
        int nToAdd=atoi(argv[++i]);
        
        int comp=sliceComp2Num(argv[++i]);
        int plane=slicePlane2Num(argv[++i]);
        float coord=atof(argv[++i]);
        
        float startTime=minT,stopTime=minT+(NT-1)*dt;
        for(int j=0;j<nToAdd;j++){
          float t=floor((stopTime-startTime)*(j+1)/nToAdd/dt);
          t=startTime+(t+sliceCompTimeRasterOffset(comp))*dt;
          
          makeNewSlice(t,comp,plane,coord);
        }
        tEprintf(Verbose,
                 "Added %i slices (%s.%s,%.1f) from time %.4g to %.4g\n",
                 nToAdd,
                 planeText(plane),compText(comp),
                 coord,startTime,stopTime);
      }
        break;
      case 't': //Add slices evenly distributed in time defined by Matlab style time vector.
        assert(argc>i+4,"addExtraOutput--adding slices requires 4 arguments");
      {
        float startTime, dtSlice;
        int nToAdd;
        setVectorValues(argc,argv,i,startTime,dtSlice,nToAdd);
        if(startTime<minT) {  //catch errors
          int dj = static_cast<int>(ceil((minT-startTime)/dtSlice));
          startTime += dj*dtSlice;
          nToAdd -= dj;
        }
        int comp=sliceComp2Num(argv[++i]);
        int plane=slicePlane2Num(argv[++i]);
        float coord=atof(argv[++i]);
        
        float stopTime=startTime+(nToAdd-1)*dtSlice;
        if(stopTime>minT+(NT-1)*dt) {  //catch errors
          stopTime = minT+(NT-1)*dt;
          nToAdd = static_cast<int>(floor((stopTime-startTime)/dtSlice))+1;
        }
        for(int j=0;j<nToAdd;j++){
          float t=floor((startTime+j*dtSlice-minT)/dt);
          t=minT+(t+sliceCompTimeRasterOffset(comp))*dt;
          makeNewSlice(t,comp,plane,coord);
        }
        tEprintf(Verbose,
                 "Added %i slices (%s.%s,%.1f) from time %.4g to %.4g\n",
                 nToAdd,
                 planeText(plane),compText(comp),
                 coord,startTime,stopTime);
      }
        break;
      case 'f': //set the max number of points per slice in a file
        assert(argc>i+1,"AddExtraOutput--setting max number of slice points per file requires 1 argument");
      {
        _maxNumPerSliceFile = atoll(argv[++i]);
        tEprintf(Verbose,"Set max number of slice points per file to %lld\n",_maxNumPerSliceFile);
      }
        break;
      case 'A': //For testing add a slice at every time step.
        assert(argc>i+3,"addExtraOutput--adding debug slices requires 3 arguments");
      {
        int comp=sliceComp2Num(argv[++i]);
        int plane=slicePlane2Num(argv[++i]);
        float coord=atof(argv[++i]);
        
        for(int j=0;j<NT;j++){
          float t=minT+(j+sliceCompTimeRasterOffset(comp))*dt;
          
          makeNewSlice(t,comp,plane,coord);
        }
        tEprintf(Verbose,
                 "Added %i slices (%s.%s,%.1f) at every time step\n",
                 NT,
                 planeText(plane),compText(comp),
                 coord);
      }
        break;
        
      default:
        assert(FALSE,"Unknown ExtraOuput flag %s",argv[0]);
    }
  }
  return _theOutput->size();
}

void extraOutputGroup::addWavefieldOutput(char* filename,
                                          float t,float fieldScalar,
                                          int nVars,...){
  wavefieldOutput* newOutput=
  new wavefieldOutput(t,fieldScalar,
                      0);
  newOutput->setFileName(filename);
  _theOutput->Add(newOutput);
  
  va_list args;
  va_start(args,nVars);
  for(int i=0;i<nVars;i++){
    char* buffer=va_arg(args,char*);
    newOutput->addVar(buffer);
  }
  va_end(args);
  sortTimes();
}
void extraOutputGroup::doWavefieldOutput(char* filename,
                                         masterSgfdDependent* forwardProblem,
                                         float t,float fieldScalar,
                                         int nVars,...){
  wavefieldOutput* newOutput=
  new wavefieldOutput(t,fieldScalar,
                      0);
  newOutput->setFileName(filename);
  
  va_list args;
  va_start(args,nVars);
  for(int i=0;i<nVars;i++){
    char* buffer=va_arg(args,char*);
    newOutput->addVar(buffer);
  }
  va_end(args);
  
  newOutput->doOutput(filename,forwardProblem,"",_work,_sliceIndices);
  delete newOutput;
}

void extraOutputGroup::initializeSliceVariables(masterSgfdModel* model,int allocateStorage){
  modelDefStruct* modelDef = model->modelDef();
  if(allocateStorage){
    DEF_MODEL_SIZE(modelDef);
    //Initialize the indices array.
    assert((_sliceIndices=(int**)malloc((numPlanes()+1)*sizeof(int*)))!=NULL,
           "initializeSliceVariables--unable to allocate space for sliceIndices");
    _sliceIndices[0]=NULL;
    assert((_sliceIndicesGlobal=(int**)malloc((numPlanes()+1)*sizeof(int*)))!=NULL,
           "initializeSliceVariables--unable to allocate space for sliceIndices");
    _sliceIndicesGlobal[0]=NULL;
    for(int j=1;j<=numPlanes();j++){
      assert((_sliceIndices[j]=(int*)malloc((numComps()+1)*sizeof(int)))!=NULL,
             "initializeSliceVariables--unable to allocate space for sliceIndices");
      assert((_sliceIndicesGlobal[j]=(int*)malloc((numComps()+1)*sizeof(int)))!=NULL,
             "initializeSliceVariables--unable to allocate space for sliceIndices");
      for(int i=0;i<=numComps();i++){
        _sliceIndices[j][i]=0;
        _sliceIndicesGlobal[j][i]=0;
      }
    }
    
    //Allocate space for the slices.
    //int numToAllocate=MAX(NX,NY)*MAX(NY,NZ);
    int numToAllocate = MAX(model->_maxSlaveNX,model->_maxSlaveNY)*MAX(MIN(model->_maxSlaveNX,model->_maxSlaveNY),model->_maxSlaveNZ);
    fprintf(stderr,"maxSlave: %d %d %d : %d\n",model->_maxSlaveNX,model->_maxSlaveNY,model->_maxSlaveNZ,numToAllocate);
    assert((_work=(float*)malloc(numToAllocate*sizeof(float)))!=NULL,
           "initializeSliceVariables--unable to allocate %ix%i block for result",
           MAX(NX,NY),MAX(NY,NZ));
    for(int i=0;i<numToAllocate;_work[i++]=0.0);
    setMessageBuffer(numToAllocate*sizeof(float)+100);  //technically 28+sliceDim1*sliceDim2*sizeof(float)
  }
  
  //Here is where the output is sorted.
  sortTimes();
}

//Make these virtual so derived classes that need to define new slice types can
// do so.
int extraOutputGroup::sliceComp2Num(char* compName){
  if(isdigit(*compName))
    return atoi(compName);
  
  char buffer1[128],buffer2[128];
  for(int i=1;i<numComps();i++){
    if(!strncmp(lower2(compName,buffer1),lower2(compText(i),buffer2),
                MIN(strlen(compName),strlen(compText(i)))))
      return i;
  }
  assert(FALSE,
         "extraOutputGroup::sliceComp2Num--unable to match slice component %s",
         compName);
  return FALSE; //never executed
}
int extraOutputGroup::slicePlane2Num(char* planeName){
  if(isdigit(*planeName))
    return atoi(planeName);
  
  char buffer1[128],buffer2[128];
  for(int i=1;i<numPlanes();i++){
    if(!strcmp(lower2(planeName,buffer1),lower2(planeText(i),buffer2)))
      return i;
  }
  assert(FALSE,"slicePlane2Num--unable to match slice plane %s",planeName);
  return FALSE; //never executed
}

