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
 *  movingSgfdSources.cc
 *
 *
 *  Define several class functions that control experimental moving
 *  sources.  Even stationary sources can call these classes because
 *  they simply override the standard stationary sources.
 *
 *  Defines the following class functions:
 *  movingSource::movingSource
 *  movingSource::movingSource
 *  movingSource::readAuxMessage
 *  movingSource::packSourceMessage
 *  movingForceSource::movingForceSource
 *  movingForceSource::applyVel
 *  movingForceSource::applyAcousticVel
 *  movingMomentSource::movingMomentSource
 *  movingMomentSource::applyVel
 *  movingMomentSource::applyAcousticVel
 *  movingMomentSource::applyStress
 *  movingMomentSource::applyPressure
 *  movingSourceNetwork::movingSourceNetwork
 *  movingSourceNetwork::movingSourceNetwork
 *  movingSourceNetwork::readForceSources
 *  movingSourceNetwork::readMomentSources
 *  movingSourceNetwork::setDispersionFactor
 *
 */

#include "movingSgfdSources.hh"
#include "netcdf.h"
#include "io_procs.h"
#include "message_passing.h"

///Here is constructor that reads the values from a message.
movingSource::movingSource(modelDefStruct* modelDef,int startTag,float *rho):
pointSource(modelDef){
  DEF_MODEL_SIZE(modelDef);
  _rhoField=rho;
  
  //Allocate space for the location variables and read the time variable
  // locations from the message.
  assert((_xMoving=(float*)malloc((NT+1)*sizeof(float)))!=NULL,
         "movingSource--unable to allocate %i floats for x-loc time function",
         NT);
  assert((_yMoving=(float*)malloc((NT+1)*sizeof(float)))!=NULL,
         "movingSource--unable to allocate %i floats for y-loc time function",
         NT);
  assert((_zMoving=(float*)malloc((NT+1)*sizeof(float)))!=NULL,
         "movingSource--unable to allocate %i floats for z-loc time function",
         NT);
  
  unpackMessage("FFF",
                _xMoving,MIN(MAX_TIME_STEPS_SEND,NT),
                _yMoving,MIN(MAX_TIME_STEPS_SEND,NT),
                _zMoving,MIN(MAX_TIME_STEPS_SEND,NT));
  
  _xMoving[NT]=_xMoving[NT-1];
  _yMoving[NT]=_yMoving[NT-1];
  _zMoving[NT]=_zMoving[NT-1];
}

///Constructor that reads parameters read from CDF file.
movingSource::movingSource(modelDefStruct* modelDef,
                           char* filename,int infile,size_t index,
                           int xVar,int yVar,int zVar,
                           int ampVar,int dataVar):
pointSource(modelDef,filename,infile,index,-1,-1,-1,ampVar,dataVar){
  DEF_MODEL_SIZE(modelDef);
  //Allocate space for the location variables.
  assert((_xMoving=(float*)malloc((NT+1)*sizeof(float)))!=NULL,
         "movingSource--unable to allocate %i floats for x-loc time function",
         NT);
  assert((_yMoving=(float*)malloc((NT+1)*sizeof(float)))!=NULL,
         "movingSource--unable to allocate %i floats for y-loc time function",
         NT);
  assert((_zMoving=(float*)malloc((NT+1)*sizeof(float)))!=NULL,
         "movingSource--unable to allocate %i floats for z-loc time function",
         NT);
  
  //Read the time-varing locations from the file.
  size_t start[2]={index,0},count[2]={1,NT};
  assert(nc_get_vara_float(infile,xVar,start,count,_xMoving)==NC_NOERR,
         "movingSource--unable to read X[%i] from %s(%i,%i)",
         (int)index,filename,infile,xVar);
  assert(nc_get_vara_float(infile,yVar,start,count,_yMoving)==NC_NOERR,
         "movingSource--unable to read Y[%i] from %s(%i,%i)",
         (int)index,filename,infile,yVar);
  assert(nc_get_vara_float(infile,zVar,start,count,_zMoving)==NC_NOERR,
         "movingSource--unable to read Z[%i] from %s(%i,%i)",
         (int)index,filename,infile,zVar);
  
  //Set the zero-time locations to the first point.
  _x=_xMoving[0];
  _y=_yMoving[0];
  _z=_zMoving[0];
  
  setActive(modelDef);
}

///Additional message reading function. If there are a large number of
/// time-steps we use an additional function to read the rest of the data.
int movingSource::readAuxMessage(int NT,int startIndex,
                                 int startTag,int iGetMessg){
  if(iGetMessg)
    getMessage(Parent,startTag+startIndex,NULL);
  unpackMessage("FFF",
                _xMoving+startIndex,MIN(MAX_TIME_STEPS_SEND,NT-startIndex),
                _yMoving+startIndex,MIN(MAX_TIME_STEPS_SEND,NT-startIndex),
                _zMoving+startIndex,MIN(MAX_TIME_STEPS_SEND,NT-startIndex));
  _xMoving[NT]=_xMoving[NT-1];
  _yMoving[NT]=_yMoving[NT-1];
  _zMoving[NT]=_zMoving[NT-1];
  
  return (NT+1-startIndex)<=MAX_TIME_STEPS_SEND;
}

///Pack the correct source message. This is just type time variable location,
/// we assume the other parameters are packed by other inhereited classes.
void movingSource::packSourceMessage(modelDefStruct* modelDef,
                                     int startIteration){
  int NT=modelDef->NT;
  packMessage("FFF",
              _xMoving+startIteration,MIN(MAX_TIME_STEPS_SEND,NT),
              _yMoving+startIteration,MIN(MAX_TIME_STEPS_SEND,NT),
              _zMoving+startIteration,MIN(MAX_TIME_STEPS_SEND,NT));
}

///Parameters read from a message.
movingForceSource::movingForceSource(modelDefStruct* modelDef,int startTag,
                                     float* rho):
pointSource(modelDef),
forceSource(modelDef,FALSE,rho),
movingSource(modelDef,FALSE,rho){
  calcInterpCoeffs(modelDef,rho);
  
  //Check to see if additional data needs to be read because of a large
  // number of time-steps.
  DEF_MODEL_SIZE(modelDef);
  int call=1;
  for(int j=MAX_TIME_STEPS_SEND;startTag&&j<NT;call++,j+=MAX_TIME_STEPS_SEND){
    pointSource::readAuxMessage(NT,j,startTag,TRUE);
    movingSource::readAuxMessage(NT,j,FALSE,FALSE);
  }
  
  //Optional very verbose output.
  int veryVerbose=FALSE;
  if(veryVerbose){
    tEprintf(veryVerbose,
             "Read Moving Source Parameters from Message w/%i segments\nInitial Values\n",
             call);
    tEprintf(veryVerbose,
             "\tLoc (%0.2f, %0.2f, %0.2f) Dir (%0.2f, %0.2f, %0.2f) Scale (%0.2f, %0.2f)\n",
             _x,_y,_z,_ax,_ay,_az,_amp,_rhox[0]);
    tEprintf(veryVerbose,"Time Varying Values\n");
    for(int i=0;i<10;i++)
      tEprintf(veryVerbose,"\t%.2f\t%.2f\t%.2f\t%.2f\n",
               _xMoving[i],_yMoving[i],_zMoving[i],_data[i]);
    tEprintf(veryVerbose,"\t.....................................\n");
    for(int i=NT-10;i<NT;i++)
      tEprintf(veryVerbose,"\t%.2f\t%.2f\t%.2f\t%.2f\n",
               _xMoving[i],_yMoving[i],_zMoving[i],_data[i]);
    
    DEF_MODEL_LIMITS(modelDef);
    FILE* out=fopen("mfsource.txt","w");
    for(int i=0;i<NT;i++)
      fprintf(out,"%.5g\t%.5g\t%.5g\t%.5g\t%.5g\n",
              minT+dt*i,_xMoving[i],_yMoving[i],_zMoving[i],_data[i]);
    fclose(out);
  }
}

///The applyVel method is easy, just call movingSource::applyVel
/// and then forceSource::applyVel.
int movingForceSource::applyVel(modelDefStruct* modelDef,int iteration,
                                float *vx,float *vy,float *vz,
                                float* rho){
  movingSource::applyVel(modelDef,iteration,
                         vx,vy,vz,rho);
  if(!setActive(modelDef)){
    return 0;
  }
  
  calcInterpCoeffs(modelDef,_rhoField);
  return
  forceSource::applyVel(modelDef,iteration,
                        vx,vy,vz,rho);
}
///The applyAcousticVel method is easy, just call movingSource::applyAcousticVel
/// and then forceSource::applyAcousticVel.
int movingForceSource::applyAcousticVel(modelDefStruct* modelDef,int iteration,
                                        float *vx,float *vy,float *vz,
                                        float* rho){
  movingSource::applyAcousticVel(modelDef,iteration,
                                 vx,vy,vz,rho);
  if(!setActive(modelDef)){
    return 0;
  }
  
  calcInterpCoeffs(modelDef,_rhoField);
  return
  forceSource::applyAcousticVel(modelDef,iteration,
                                vx,vy,vz,rho);
}

///Parameters read from a message.
movingMomentSource::movingMomentSource(modelDefStruct* modelDef,int startTag,float* rho):
pointSource(modelDef),
momentSource(modelDef,FALSE,rho),
movingSource(modelDef,FALSE,rho){
  calcInterpCoeffs(modelDef,rho);
  
  //Check to see if additional data needs to be read because of a large
  // number of time-steps.
  DEF_MODEL_SIZE(modelDef);
  for(int j=MAX_TIME_STEPS_SEND;startTag&&j<NT;j+=MAX_TIME_STEPS_SEND){
    pointSource::readAuxMessage(NT,j,startTag,TRUE);
    movingSource::readAuxMessage(NT,j,FALSE,FALSE);
  }
  
  //Optional very verbose output can be set in a debugger.
  int veryVerbose=FALSE;
  if(veryVerbose){
    tEprintf(veryVerbose,"Read Moving Source Parameters from Message\nInitial Values\n");
    tEprintf(veryVerbose,"\tLoc (%0.2f, %0.2f, %0.2f) Scale (%0.2f)\n",
             _x,_y,_z,_amp);
    tEprintf(veryVerbose,"Time Varying Values\n");
    
    for(int i=0;i<10;i++)
      tEprintf(veryVerbose,"\t%.2f\t%.2f\t%.2f\t%.2f\n",
               _xMoving[i],_yMoving[i],_zMoving[i],_data[i]);
    tEprintf(veryVerbose,"\t.....................................\n");
    for(int i=NT-10;i<NT;i++)
      tEprintf(veryVerbose,"\t%.2f\t%.2f\t%.2f\t%.2f\n",
               _xMoving[i],_yMoving[i],_zMoving[i],_data[i]);
  }
}

///The applyVel method is easy, just call movingSource::applyVel
/// and then momentSource::applyVel.
int movingMomentSource::applyVel(modelDefStruct* modelDef,int iteration,
                                 float *vx,float *vy,float *vz,
                                 float* rho){
  movingSource::applyVel(modelDef,iteration,
                         vx,vy,vz,rho);
  if(!setActive(modelDef)){
    return 0;
  }
  
  calcInterpCoeffs(modelDef,_rhoField);
  return
  momentSource::applyVel(modelDef,iteration,
                         vx,vy,vz,rho);
}
///The applyAcousticVel method is easy, just call movingSource::applyAcousticVel
/// and then momentSource::applyAcousticVel.
int movingMomentSource::applyAcousticVel(modelDefStruct* modelDef,int iteration,
                                         float *vx,float *vy,float *vz,
                                         float* rho){
  movingSource::applyAcousticVel(modelDef,iteration,
                                 vx,vy,vz,rho);
  if(!setActive(modelDef)){
    return 0;
  }
  
  calcInterpCoeffs(modelDef,_rhoField);
  return
  momentSource::applyAcousticVel(modelDef,iteration,
                                 vx,vy,vz,rho);
}
int movingMomentSource::applyStress(modelDefStruct* modelDef,int iteration,
                                    float *xx,float* yy,float *zz,
                                    float *xy,float* xz,float *yz){
  movingSource::applyStress(modelDef,iteration,
                            xx,yy,zz,xy,xz,yz);
  if(!setActive(modelDef)){
    return 0;
  }
  
  calcInterpCoeffs(modelDef,_rhoField);
  return
  momentSource::applyStress(modelDef,iteration,
                            xx,yy,zz,xy,xz,yz);
}
int movingMomentSource::applyPressure(modelDefStruct* modelDef,int iteration,
                                      float *pressure){
  movingSource::applyPressure(modelDef,iteration,pressure);
  if(!setActive(modelDef)){
    return 0;
  }
  
  calcInterpCoeffs(modelDef,_rhoField);
  return
  momentSource::applyPressure(modelDef,iteration,pressure);
}


///Sources passed in message.
movingSourceNetwork::movingSourceNetwork(modelDefStruct* modelDef,float *rho,
                                         int useCubicExtrap):sourceNetwork(){
  _bandCheckPCent=DEFAULT_CHECK_BANDWIDTH_PCENT;
  _wavelet=NULL;
  _nfs=_nms=_nts=_ntdbc=0;
  //Read the number of sources.
  int numSources;
  unpackMessage("i",&numSources);
  
  setMessageBuffer(numSources*((20+MIN(modelDef->NT,MAX_TIME_STEPS_SEND))*sizeof(float)));
  
  _sources=new sourceArray(0,numSources);
  for(int i=0;i<numSources;i++){
    //For each source, read the type and then call the constructor corresponding to
    // that source type.
    getMessage(Parent,MESSAGE_SET_SOURCES+1+i,NULL);
    int type;
    unpackMessage("i",&type);
    switch(type){
      case FORCE_SOURCE:
        _sources->Add(new forceSource(modelDef,MESSAGE_SET_SOURCES+1+i,rho));
        _nfs++;
        break;
      case MOMENT_SOURCE:
        _sources->Add(new momentSource(modelDef,MESSAGE_SET_SOURCES+1+i,
                                       rho,useCubicExtrap));
        _nms++;
        break;
        
      case MOVING_FORCE_SOURCE:
        _sources->Add(new movingForceSource(modelDef,MESSAGE_SET_SOURCES+1+i,rho));
        _nfs++;
        break;
      case MOVING_MOMENT_SOURCE:
        _sources->Add(new movingMomentSource(modelDef,MESSAGE_SET_SOURCES+1+i,
                                             rho));
        _nms++;
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

///Sources read from cdf file
movingSourceNetwork::movingSourceNetwork(modelDefStruct* modelDef,char* filename,
                                         int useModelSources,float *rho):sourceNetwork(){
  _useModelSources=useModelSources;
  _bandCheckPCent=DEFAULT_CHECK_BANDWIDTH_PCENT;
  _nfs=_nms=_nts=_ntdbc=0;
  _wavelet=NULL;
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

///Read force sources with special treatment for the moving ones.
int movingSourceNetwork::readForceSources(modelDefStruct* modelDef,
                                          char* filename,int infile,
                                          float* rho){
  //Strip off any moving sources.
  int nMovingFsDim;
  if(nc_inq_dimid(infile,"numMFSources",&nMovingFsDim)==NC_NOERR){
    //Read the dimension to find the actual number of sources.
    size_t nmfs;
    assert(nc_inq_dimlen(infile,nMovingFsDim,&nmfs)==NC_NOERR,
           "sourceNetwork--unable read dim numMFSources in %s(%i,%i)",
           filename,infile,nMovingFsDim);
    _nfs+=nmfs;
    
    //Get the required variable id's.
    int fsXsVar,fsYsVar,fsZsVar;
    assert(nc_inq_varid(infile,"mfsourcesXs",&fsXsVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mfsourcesXs in file %s(%i)",
           filename,infile);
    assert(nc_inq_varid(infile,"mfsourcesYs",&fsYsVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mfsourcesYs in file %s(%i)",
           filename,infile);
    assert(nc_inq_varid(infile,"mfsourcesZs",&fsZsVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mfsourcesZs in file %s(%i)",
           filename,infile);
    
    int fsSampVar,fsAxVar,fsAyVar,fsAzVar;
    assert(nc_inq_varid(infile,"mfsourcesSamp",&fsSampVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mfsourcesSamp in file %s(%i)",
           filename,infile);
    assert(nc_inq_varid(infile,"mfsourcesAx",&fsAxVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mfsourcesAx in file %s(%i)",
           filename,infile);
    assert(nc_inq_varid(infile,"mfsourcesAy",&fsAyVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mfsourcesAy in file %s(%i)",
           filename,infile);
    assert(nc_inq_varid(infile,"mfsourcesAz",&fsAzVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mfsourcesAz in file %s(%i)",
           filename,infile);
    
    int fsDataVar;
    assert(nc_inq_varid(infile,"mfsourcesData",&fsDataVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mfsourcesData in file %s(%i)",
           filename,infile);
    
    //Create a new source for each.
    for(unsigned int i=0;i<nmfs;i++)
      _sources->Add(new movingForceSource(modelDef,filename,infile,i,
                                          fsXsVar,fsYsVar,fsZsVar,
                                          fsAxVar,fsAyVar,fsAzVar,
                                          fsSampVar,fsDataVar,
                                          rho));
  }
  
  //Now call the parent class to get any non-moving force sources.
  return sourceNetwork:: readForceSources(modelDef,filename,infile,rho);
}

///Modification of the readMomentSources method. First strip off any
/// moving-moment sources, then call the sourceNetwork method.
int movingSourceNetwork::readMomentSources(modelDefStruct* modelDef,
                                           char* filename,int infile,
                                           float* rho){
  //Strip off any moving-moment sources.
  int nmmsDim;
  if(nc_inq_dimid(infile,"numMMSources",&nmmsDim)==NC_NOERR){
    //Read the dimension to find the actual number of sources.
    size_t nmms;
    assert(nc_inq_dimlen(infile,nmmsDim,&nmms)==NC_NOERR,
           "sourceNetwork--unable read dim numMMSources in %s(%i,%i)",
           filename,infile,nmmsDim);
    _nms+=nmms;
    
    //Get the required variable id's.
    int msXsVar,msYsVar,msZsVar;
    assert(nc_inq_varid(infile,"mMSourcesXs",&msXsVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mMSourcesXs in file %s(%i)",
           filename,infile);
    assert(nc_inq_varid(infile,"mMSourcesYs",&msYsVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mMSourcesYs in file %s(%i)",
           filename,infile);
    assert(nc_inq_varid(infile,"mMSourcesZs",&msZsVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mMSourcesZs in file %s(%i)",
           filename,infile);
    
    int msSampVar;
    assert(nc_inq_varid(infile,"mMSourcesSamp",&msSampVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mMSourcesSamp in file %s(%i)",
           filename,infile);
    
    int msXxSVar,msYySVar,msZzSVar,msXySVar,msXzSVar,msYzSVar;
    assert(nc_inq_varid(infile,"mMSourcesXxS",&msXxSVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mMSourcesXxS in file %s(%i)",
           filename,infile);
    assert(nc_inq_varid(infile,"mMSourcesYyS",&msYySVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mMSourcesYyS in file %s(%i)",
           filename,infile);
    assert(nc_inq_varid(infile,"mMSourcesZzS",&msZzSVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mMSourcesZzS in file %s(%i)",
           filename,infile);
    assert(nc_inq_varid(infile,"mMSourcesXyS",&msXySVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mMSourcesXyS in file %s(%i)",
           filename,infile);
    assert(nc_inq_varid(infile,"mMSourcesXzS",&msXzSVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mMSourcesXzS in file %s(%i)",
           filename,infile);
    assert(nc_inq_varid(infile,"mMSourcesYzS",&msYzSVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mMSourcesYzS in file %s(%i)",
           filename,infile);
    
    int msXyAVar,msXzAVar,msYzAVar;
    assert(nc_inq_varid(infile,"mMSourcesXyA",&msXyAVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mMSourcesXyA in file %s(%i)",
           filename,infile);
    assert(nc_inq_varid(infile,"mMSourcesXzA",&msXzAVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mMSourcesXzA in file %s(%i)",
           filename,infile);
    assert(nc_inq_varid(infile,"mMSourcesYzA",&msYzAVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mMSourcesYzA in file %s(%i)",
           filename,infile);
    
    int msDataVar;
    assert(nc_inq_varid(infile,"mMSourcesData",&msDataVar)==NC_NOERR,
           "sourceNetwork--unable to open variable mMSourcesData in file %s(%i)",
           filename,infile);
    
    //Create a new source for each.
    for(unsigned int i=0;i<nmms;i++)
      _sources->Add(new movingMomentSource(modelDef,filename,infile,i,
                                           msXsVar,msYsVar,msZsVar,
                                           msXxSVar,msYySVar,msZzSVar,
                                           msXySVar,msXzSVar,msYzSVar,
                                           msXyAVar,msXzAVar,msYzAVar,
                                           msSampVar,msDataVar,
                                           rho));
  }
  
  //Now call the parent class to get any non-moving moment sources.
  return sourceNetwork:: readMomentSources(modelDef,filename,infile,rho);
}

///Pick of the new types and set the factor. The dispersion check should really
/// acount for the Doppler shift but I will ignore that for now.
int movingSourceNetwork::setDispersionFactor(int i,int doPrint){
  if(type(i)==MOVING_FORCE_SOURCE){
    tEprintf(doPrint," Source #%i: Moving Force\n",i+1);
    return 1;
  }else if(type(i)==MOVING_MOMENT_SOURCE){
    tEprintf(doPrint," Source #%i: Moving Moment\n",i+1);
    return 2;
  }
  return sourceNetwork::setDispersionFactor(i,doPrint);
}

