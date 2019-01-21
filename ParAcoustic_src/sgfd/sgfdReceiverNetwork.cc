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
 *  sgfdReceiverNetwork.cc
 *
 *
 *  Defines class functions used in the classes that hold all the receivers
 *  for both master and slave.   These receiverNetwork classes allow a single
 *  class to handle all receiver functions as a group and are basically containers of
 *  receivers.
 *
 *  Defines the following class functions:
 *  receiverNetwork::receiverNetwork
 *  receiverNetwork::~receiverNetwork
 *  receiverNetwork::combine
 *  receiverNetwork::work
 *  receiverNetwork::allocateData
 *  receiverNetwork::addVelocityReceiver
 *  receiverNetwork::fillReceivers
 *  masterReceiverNetwork::masterReceiverNetwork
 *  masterReceiverNetwork::checkActive
 *  masterReceiverNetwork::addReceivers
 *  masterReceiverNetwork::addReceivers
 *  masterReceiverNetwork::initializeNetwork
 *  masterReceiverNetwork::writeTraces
 *  masterReceiverNetwork::getReceiverData
 *  masterReceiverNetwork::getReceiverGridData
 *  masterReceiverNetwork::writeSingleHeader
 *  masterReceiverNetwork::writeSingleLocations
 *  masterReceiverNetwork::rewriteSingleLocations
 *  masterReceiverNetwork::writeSingleTraces
 *  masterReceiverNetwork::generateSplitFilename
 *  masterReceiverNetwork::writeCompSplitHeader
 *  masterReceiverNetwork::writeCompSplitLocations
 *  masterReceiverNetwork::writeCompSplitTraces
 *  masterReceiverNetwork::writeCompSplitTraces
 *  masterReceiverNetwork::writeIndexSplitHeader
 *  masterReceiverNetwork::writeIndexSplitLocations
 *  masterReceiverNetwork::writeIndexSplitTraces
 *  masterReceiverNetwork::writeIndexSplitTraces
 *  masterReceiverNetwork::writeSingleGrid
 *  masterReceiverNetwork::add3CReceiver
 *  masterReceiverNetwork::add4CReceiver
 *  slaveReceiverNetwork::slaveReceiverNetwork
 *  slaveReceiverNetwork::doCheckpoint
 *  slaveReceiverNetwork::readCheckpoint
 *  slaveReceiverNetwork::initializeNetwork
 *
 */

#include "sgfdReceiverNetwork.hh"
#include <time.h>
#include "netcdf.h"
#include "io_procs.h"
#include "message_passing.h"

///Need a null constructor for cases where receivers are specified entirely on
/// the command line
receiverNetwork::receiverNetwork(){
  _receivers=new receiverArray;
  _receiverGrids=new receiverGridArray;
  _wrk=NULL;
  
  strcpy(_traceOutputName,"trace.cdf");
}

receiverNetwork::~receiverNetwork(){
  if(_receivers){
    for(int i=0;i<_receivers->size();i++){
      //receiver* currR=(*_receivers)[i];
      // 	delete currR; //KLUDGE: getting an error when freeing the receivers.
    }
    delete _receivers;
    _receivers=NULL;
  }
  if(_receiverGrids){
    for(int i=0;i<_receiverGrids->size();i++){
      //receiverGrid* currR=(*_receiverGrids)[i];
      // 	delete currR; //KLUDGE: getting an error when freeing the receivers.
    }
    delete _receiverGrids;
    _receiverGrids=NULL;
  }
  
  //     if(_wrk) //KLUDGE: getting an error when freeing the receivers.
  //       free(_wrk);
}

///Method to combine receivers (some read from file; others
/// defined on the command line).
int receiverNetwork::combine(receiverNetwork* extraReceivers){
  if(extraReceivers){
    for(int i=0;i<extraReceivers->size();i++){
      _receivers->Add((*extraReceivers->_receivers)[i]);
    }
  }
  return size();
}

float* receiverNetwork::work(modelDefStruct* modelDef){
  if(!_wrk)
    assert((_wrk=(float*)malloc((modelDef->NT+2)*sizeof(float)))!=NULL,
           "receiverNetwork::work--unable to allocate %i floats",
           modelDef->NT+2);
  return _wrk;
}

//May want to allocate data for all of the receivers.
int receiverNetwork::allocateData(modelDefStruct* modelDef){
  int numAllocated=0;
  if(_receivers){
    for(int ir=0;ir<_receivers->size();ir++){
      receiver* currR=(*_receivers)[ir];
      if(currR->active()){
        numAllocated++;
        currR->allocateData(modelDef);
      }
    }
  }
  return numAllocated;
}

int receiverNetwork::addVelocityReceiver(modelDefStruct* modelDef,
                                         float xr,float yr,float zr,
                                         float bx,float by,float bz,
                                         float ramp,int integrate,
                                         int print,int useCubicInterp,int surfaceReceiver){
  velocityReceiver *currReceiver;
  _receivers->Add(currReceiver=new velocityReceiver(modelDef,FALSE,
                                                    xr,yr,zr,bx,by,bz,ramp,useCubicInterp,surfaceReceiver));
  if(integrate>0)
    currReceiver->setDisplacement();
  else if(integrate<0)
    currReceiver->setAcceleration();
  
  tEprintf(print && Verbose,
           "Added velocity receiver %i (%.1f,%.1f,%.1f) (%.1f,%.1f,%.1f,%.1f)\n",
           _receivers->size(),xr,yr,zr,ramp,bx,by,bz);
  return size();
}

//Fill the data for this network of receivers. Moved from slaveReceiverNetwork
// to here for serial implementation.
int receiverNetwork::fillReceivers(modelDefStruct* modelDef,int iteration,
                                   float* vx,float* vy,float* vz,
                                   float* xx,float* yy,float* zz,
                                   float* xy,float* xz,float* yz){
  if(_receivers){
    for(int ir=0;ir<_receivers->size();ir++){
      receiver* currR=(*_receivers)[ir];
      if(currR->active())
        currR->fill(modelDef,iteration,
                    vx,vy,vz,
                    xx,yy,zz,
                    xy,xz,yz);
    }
  }
  if(_receiverGrids){
    for(int ir=0;ir<_receiverGrids->size();ir++){
      receiver* currR=(*_receiverGrids)[ir];
      if(currR->active())
        currR->fill(modelDef,iteration,
                    vx,vy,vz,
                    xx,yy,zz,
                    xy,xz,yz);
    }
  }
  return size();
}

#define REC_SPLIT_BY_COMP 1

///A new receiverNetwork can be read directly from a cdf file. This subsumes the previous function
/// readCDFReceivers.
masterReceiverNetwork::masterReceiverNetwork(modelDefStruct* modelDef,const char* fileName,int freeSurface,
                                             int useModelReceivers,int timeSubsample):
receiverNetwork(){
  _useModelReceivers=useModelReceivers;
  _splitFiles=FALSE;
  
  _receiverProcs=NULL;
  _timeSubsample=timeSubsample;
  
  _rTypeFlag=NULL;
  if(!_useModelReceivers){
    _nVxR=_nVyR=_nVzR=_nPrR=0;
  }else{
    //open cdf file for input
    int inFile=openCDFFile(fileName,FALSE,NC_NOWRITE);
    
    //read receivers.
    int nReceiversDim;
    if(nc_inq_dimid(inFile,"numReceivers",&nReceiversDim)==NC_NOERR){
      size_t nReceivers;
      assert(nc_inq_dimlen(inFile,nReceiversDim,&nReceivers)==NC_NOERR,
             "masterReceiverNetworks--unable to read dimension numReceivers from %s(%i,%i)",
             fileName,inFile,nReceiversDim);
      
      int receiverRtypeVar,receiverRampVar;
      assert(nc_inq_varid(inFile,"receiverType",&receiverRtypeVar)==NC_NOERR,
             "masterReceiverNetworks--unable to open variable receiverRtype in file %s",
             fileName);
      assert(nc_inq_varid(inFile,"receiverAmp",&receiverRampVar)==NC_NOERR,
             "masterReceiverNetworks--unable to open variable receiverAmp in file %s",
             fileName);
      
      int receiverXrVar,receiverYrVar,receiverZrVar;
      assert(nc_inq_varid(inFile,"receiverX",&receiverXrVar)==NC_NOERR,
             "masterReceiverNetworks--unable to open variable receiverX in file %s",
             fileName);
      assert(nc_inq_varid(inFile,"receiverY",&receiverYrVar)==NC_NOERR,
             "masterReceiverNetworks--unable to open variable receiverY in file %s",
             fileName);
      assert(nc_inq_varid(inFile,"receiverZ",&receiverZrVar)==NC_NOERR,
             "masterReceiverNetworks--unable to open variable receiverZ in file %s",
             fileName);
      
      int receiverBxVar,receiverByVar,receiverBzVar;
      assert(nc_inq_varid(inFile,"receiverBx",&receiverBxVar)==NC_NOERR,
             "masterReceiverNetworks--unable to open variable receiverBx in file %s",
             fileName);
      assert(nc_inq_varid(inFile,"receiverBy",&receiverByVar)==NC_NOERR,
             "masterReceiverNetworks--unable to open variable receiverBy in file %s",
             fileName);
      assert(nc_inq_varid(inFile,"receiverBz",&receiverBzVar)==NC_NOERR,
             "masterReceiverNetworks--unable to open variable receiverBz in file %s",
             fileName);
      
      int receiverIntegrateVar;
      assert(nc_inq_varid(inFile,"receiverIntegrate",&receiverIntegrateVar)==NC_NOERR,
             "masterReceiverNetworks--unable to open variable receiverIntegrate in file %s",
             fileName);
      
      for(size_t i=0;i<nReceivers;i++){
        int rtype;
        nc_get_var1_int(inFile,receiverRtypeVar,&i,&rtype);
        
        //get the location and amplitude
        float x,y,z,amp;
        nc_get_var1_float(inFile,receiverXrVar,&i,&x);
        nc_get_var1_float(inFile,receiverYrVar,&i,&y);
        nc_get_var1_float(inFile,receiverZrVar,&i,&z);
        nc_get_var1_float(inFile,receiverRampVar,&i,&amp);
        
        if(rtype == VELOCITY_RECEIVER){
          //get the direction
          float bx,by,bz;
          nc_get_var1_float(inFile,receiverBxVar,&i,&bx);
          nc_get_var1_float(inFile,receiverByVar,&i,&by);
          nc_get_var1_float(inFile,receiverBzVar,&i,&bz);
          
          velocityReceiver* currReceiver;
          _receivers->Add(currReceiver=new velocityReceiver(modelDef,FALSE,
                                                            x,y,z,bx,by,bz,amp));
          
          //determine if this is a displacment or acceleration receiver
          int integrate;
          assert(nc_get_var1_int(inFile,receiverIntegrateVar,&i,&integrate)==NC_NOERR,
                 "masterReceiverNetworks--unable to read receiverIntegrate[%i] from %s(%i,%i)",
                 i,fileName,inFile,receiverIntegrateVar);
          if(integrate>0){
            currReceiver->setDisplacement();
          }else if(integrate<0){
            currReceiver->setAcceleration();
          }
        }else if(rtype == PRESSURE_RECEIVER){
          _receivers->Add(new pressureReceiver(modelDef,FALSE,x,y,z,amp));
        }else{
          assert(FALSE,"masterReceiverNetworks--unknown receiver type %i for receiver %i in %s",
                 rtype,i,fileName);
        }
      }
    }
    
    assert(nc_close(inFile)==NC_NOERR,
           "masterReceiverNetwork--unable to close %s(%i)",
           fileName,inFile);
  }
}

///Method to make sure that all receivers in the domain are
/// active. Failure of this check is a fatal error.
int masterReceiverNetwork::checkActive(int fail){
  int numActive=0;
  for(int i=0;i<_receivers->size();i++){
    int currActive=_receiverProcs && _receiverProcs[i]>=0;
    if(currActive){
      numActive++;
    }else{
      assert(!fail,
             "checkActive--receiver[%i] (%.2f, %.2f, %.2f) is NOT active",
             i,x(i),y(i),z(i));
    }
  }
  
  return numActive;
}

//Can add receivers to an existing receiverNetwork by processing command line arguments. This replaces
// addReceivers. The argument argOffset is the index of the first character after the one that
// tells the calling process to call this function.
int masterReceiverNetwork::addReceivers(int& i,int argOffset,
                                        int argc,char* argv[],
                                        modelDefStruct* modelDef){
  static int integrate=FALSE;
  static float angleMult=DEG_TO_RAD;
  static int useCubicInterp = TRUE;
  static int surfaceReceiver = FALSE;
  
  if(strlen(argv[i])<argOffset+1 || argv[i][argOffset]=='1'){
    assert(argc>i+5,"addReceivers--requires at least 5 arguments to add basic receiver");
    
    //read the receiver type
    int rtype=receiverType(argv[++i]);
    // all types have a location and amplitude
    float xr=atof(argv[++i]);
    float yr=atof(argv[++i]);
    float zr=atof(argv[++i]);
    float ramp=atof(argv[++i]);
    
    switch(rtype){
      case PRESSURE_RECEIVER:
        _receivers->Add(new pressureReceiver(modelDef,FALSE,xr,yr,zr,ramp,useCubicInterp));
        tEprintf(Verbose,
                 "Added pressure receiver %i (%.1f,%.1f,%.1f) (%.1f)\n",
                 _receivers->size(),xr,yr,zr,ramp);
        break;
      case VELOCITY_RECEIVER:
        //need to read the angles and generate the direction cosines
        assert(argc>i+2,
               "addReceivers--requires 2 arguments to specify direction of velocity receiver");
      {
        float theta=atof(argv[++i])*angleMult;
        float phi=atof(argv[++i])*angleMult;
        float bx=sin(theta)*cos(phi),
        by=sin(theta)*sin(phi),
        bz=cos(theta);
        addVelocityReceiver(modelDef,xr,yr,zr,bx,by,bz,ramp,integrate,TRUE,useCubicInterp,surfaceReceiver);
      }
        break;
      case THREE_C_RECEIVER:
        add3CReceiver(modelDef,xr,yr,zr,ramp,integrate,TRUE,useCubicInterp,surfaceReceiver);
        break;
      case FOUR_C_RECEIVER:
        add4CReceiver(modelDef,xr,yr,zr,ramp,integrate,TRUE,useCubicInterp,surfaceReceiver);
        break;
      case VX_RECEIVER:
        addVelocityReceiver(modelDef,xr,yr,zr,1.0,0.0,0.0,1.0,integrate,TRUE,useCubicInterp,surfaceReceiver);
        break;
      case VY_RECEIVER:
        addVelocityReceiver(modelDef,xr,yr,zr,0.0,1.0,0.0,1.0,integrate,TRUE,useCubicInterp,surfaceReceiver);
        break;
      case VZ_RECEIVER:
        addVelocityReceiver(modelDef,xr,yr,zr,0.0,0.0,1.0,1.0,integrate,TRUE,useCubicInterp,surfaceReceiver);
        break;
        
      default:
        assert(FALSE,"Receiver %i; type %i (%s) is unknown",
               _receivers->size()+1,rtype,ReceiverTypeNames[rtype]);
    }
  }else{
    switch(argv[i][argOffset]){
      case 'M':
        _useModelReceivers=!_useModelReceivers;
        tEprintf(Verbose,"%s receivers from model file\n",
                 _useModelReceivers?"Using":"Not using");
        break;
      case 'o':
        assert(argc>i+1,"addReceivers--1 argument required to set output file");
        strcpy(_traceOutputName,argv[++i]);
        tEprintf(Verbose,"Set trace output file to\n\t%s\n",_traceOutputName);
        break;
      case 'L':
        if(strlen(argv[i])<argOffset+1 || argv[i][argOffset+1] == 'c'){
          _splitFiles=REC_SPLIT_BY_COMP;
          tEprintf(Verbose,
                   "Splitting trace output by component to reduce file size\n");
        }else if(argv[i][argOffset+1] == 'n'){
          assert(argc>i+1,"addReceivers--1 argument required to split traces by index");
          _splitFiles=atoi(argv[++i]);
          tEprintf(Verbose,
                   "Splitting trace output with %i traces/file\n",
                   _splitFiles);
        }
        break;
      case 'v':
        integrate=FALSE;
        tEprintf(Verbose,"Subsequent velocity receivers record velocity\n");
        break;
      case 'u':
        if(strlen(argv[i])<argOffset+1 || argv[i][argOffset+1] == 'r'){
          angleMult=1.0;
          tEprintf(Verbose,"Receiver Angles in Degrees (Default)\n");
        }else if(argv[i][argOffset+1] == 'd'){
          angleMult=DEG_TO_RAD;
          tEprintf(Verbose,"Receiver Angles in Degrees (Default)\n");
        }else{
          assert(FALSE,"addReceivers--u option is [d|r] only; arg is %s",
                 argv[i]);
        }
        break;
      case 'd':
        integrate=1;
        tEprintf(Verbose,"Subsequent velocity receivers record displacement\n");
        break;
      case 'n':
        integrate=0;
        tEprintf(Verbose,"Subsequent velocity receivers record velocity\n");
        break;
      case 'a':
        integrate=-1;
        tEprintf(Verbose,"Subsequent velocity receivers record acceleration\n");
        break;
      case 's':
        assert(argc>i+1,
               "addReceivers--1 argument required to add set TimeSubsample");
        _timeSubsample=atoi(argv[++i]);
        tEprintf(Verbose,"Set ReceiverSubsample to %i\n",_timeSubsample);
        break;
      case 'c':
        useCubicInterp = TRUE;
        tEprintf(Verbose,"Subsequent velocity and pressure receivers use cubic interpolation\n");
        break;
      case 'l':
        useCubicInterp = FALSE;
        tEprintf(Verbose,"Subsequent velocity and pressure receivers use trilinear interpolation\n");
        break;
      case 'S':
        surfaceReceiver = TRUE;
        tEprintf(Verbose,"Subsequent velocity receivers will assume in rock (not air) near interface\n");
        break;
      case 'i':
        surfaceReceiver = FALSE;
        tEprintf(Verbose,"Subsequent velocity receivers will interpolate as smoothly varying properties\n");
        break;
        
      case '3':
      {
        assert(argc>i+4,"addReceivers--4 arguments required to add 3C receiver");
        // all types have a location and amplitude
        float xr=atof(argv[++i]);
        float yr=atof(argv[++i]);
        float zr=atof(argv[++i]);
        float ramp=atof(argv[++i]);
        add3CReceiver(modelDef,xr,yr,zr,ramp,integrate,TRUE,useCubicInterp,surfaceReceiver);
      }
        break;
      case '4':
      {
        assert(argc>i+4,"addReceivers--4 arguments required to add 4C receiver");
        
        // all types have a location and amplitude
        float xr=atof(argv[++i]);
        float yr=atof(argv[++i]);
        float zr=atof(argv[++i]);
        float ramp=atof(argv[++i]);
        add4CReceiver(modelDef,xr,yr,zr,ramp,integrate,TRUE,useCubicInterp,surfaceReceiver);
      }
        break;
      case 'G': //Add a 3D receiverGrid.
      {
        assert(argc>i+4,
               "addReceivers--adding receiver grid requires at least 4 arguments (type,x,y,z)");
        char* callFlag=argv[i];
        int rawData=FALSE;
        char rawDataDir[512]="";
        if(strlen(callFlag)>argOffset+1 && callFlag[argOffset+1]=='r'){
          rawData=1;
          strcpy(rawDataDir,argv[++i]);
          assert(rawDataDir[0]='/',"Raw grid directory must be fully qualified (%s).",
                 rawDataDir);
          argOffset++;
        }else if(strlen(callFlag)>argOffset+1 && callFlag[argOffset+1]=='x'){
          rawData=2;
          strcpy(rawDataDir,argv[++i]);
          assert(rawDataDir[0]='/',"Raw grid directory must be fully qualified (%s).",
                 rawDataDir);
          argOffset++;
        }
        
        int rtype=receiverType(argv[++i]);
        
        float ix,dx;
        int nx;
        setVectorValues(argc,argv,i,ix,dx,nx);
        
        float iy,dy;
        int ny;
        setVectorValues(argc,argv,i,iy,dy,ny);
        
        float iz,dz;
        int nz;
        setVectorValues(argc,argv,i,iz,dz,nz);
        
        switch(rtype){
          case PRESSURE_RECEIVER:
            _receiverGrids->Add(new
                                pressureReceiverGrid(modelDef,FALSE,1.0,
                                                     ix,dx,nx,iy,dy,ny,iz,dz,nz,
                                                     _timeSubsample,rawData,rawDataDir));
            break;
          case VX_RECEIVER:
            _receiverGrids->Add(new
                                vxReceiverGrid(modelDef,FALSE,1.0,
                                               ix,dx,nx,iy,dy,ny,iz,dz,nz,
                                               _timeSubsample,rawData,rawDataDir));
            break;
          case VY_RECEIVER:
            _receiverGrids->Add(new
                                vyReceiverGrid(modelDef,FALSE,1.0,
                                               ix,dx,nx,iy,dy,ny,iz,dz,nz,
                                               _timeSubsample,rawData,rawDataDir));
            break;
          case VZ_RECEIVER:
            _receiverGrids->Add(new
                                vzReceiverGrid(modelDef,FALSE,1.0,
                                               ix,dx,nx,iy,dy,ny,iz,dz,nz,
                                               _timeSubsample,rawData,rawDataDir));
            break;
            
          case FOUR_C_RECEIVER:
            _receiverGrids->Add(new
                                pressureReceiverGrid(modelDef,FALSE,1.0,
                                                     ix,dx,nx,iy,dy,ny,iz,dz,nz,
                                                     _timeSubsample,rawData,rawDataDir));
            //Fall through to get the velocity components.
          case THREE_C_RECEIVER:
            _receiverGrids->Add(new
                                vxReceiverGrid(modelDef,FALSE,1.0,
                                               ix,dx,nx,iy,dy,ny,iz,dz,nz,
                                               _timeSubsample,rawData,rawDataDir));
            _receiverGrids->Add(new
                                vyReceiverGrid(modelDef,FALSE,1.0,
                                               ix,dx,nx,iy,dy,ny,iz,dz,nz,
                                               _timeSubsample,rawData,rawDataDir));
            _receiverGrids->Add(new
                                vzReceiverGrid(modelDef,FALSE,1.0,
                                               ix,dx,nx,iy,dy,ny,iz,dz,nz,
                                               _timeSubsample,rawData,rawDataDir));
            break;
            
          default:
            assert(FALSE,
                   "addReceivers--receiver grid only defined for P, Vx, Vy, or Vz receivers not %i",
                   rtype);
        }
        
        tEprintf(Verbose,"Added %ix%ix%i (%i) %s %s receiver-grid\n",
                 nx,ny,nz,nx*ny*nz,
                 rawData?"Raw Binary":"NetCDF",ReceiverTypeNames[rtype]);
        tEprintf(Verbose,"  Grid %i spans: X %.2f to %.2f, Y %.2f to %.2f, Z %.2f to %.2f\n",
                 _receiverGrids->size(),
                 ix,ix+dx*(nx-1),iy,iy+dy*(ny-1),iz,iz+dz*(nz-1));
        tEprintf(Verbose&&rawData==2,"  Creating on the fly Cross-correlation with\n\t%s\n",
                 rawDataDir);
      }
        break;
        
      case 'g': //Add a 3D grid of receivers. Note a string is a special case nx=nz=1.
      {
        assert(argc>i+4,
               "addReceivers--adding receiver grid requires at least 4 arguments (type,x,y,z)");
        int rtype=receiverType(argv[++i]);
        
        float ix,dx;
        int nx;
        setVectorValues(argc,argv,i,ix,dx,nx);
        
        float iy,dy;
        int ny;
        setVectorValues(argc,argv,i,iy,dy,ny);
        
        float iz,dz;
        int nz;
        setVectorValues(argc,argv,i,iz,dz,nz);
        
        for(int k=0;k<nz;k++){
          float z=iz+dz*k;
          for(int j=0;j<ny;j++){
            float y=iy+dy*j;
            for(int i=0;i<nx;i++){
              float x=ix+dx*i;
              
              switch(rtype){
                case PRESSURE_RECEIVER:
                  _receivers->Add(new pressureReceiver(modelDef,FALSE,x,y,z,1.0,useCubicInterp));
                  break;
                case THREE_C_RECEIVER:
                  add3CReceiver(modelDef,x,y,z,1.0,integrate,FALSE,useCubicInterp,surfaceReceiver);
                  break;
                case FOUR_C_RECEIVER:
                  add4CReceiver(modelDef,x,y,z,1.0,integrate,FALSE,useCubicInterp,surfaceReceiver);
                  break;
                  
                case VX_RECEIVER:
                  addVelocityReceiver(modelDef,x,y,z,1.0,0.0,0.0,1.0,integrate,FALSE,useCubicInterp,surfaceReceiver);
                  break;
                case VY_RECEIVER:
                  addVelocityReceiver(modelDef,x,y,z,0.0,1.0,0.0,1.0,integrate,FALSE,useCubicInterp,surfaceReceiver);
                  break;
                case VZ_RECEIVER:
                  addVelocityReceiver(modelDef,x,y,z,0.0,0.0,1.0,1.0,integrate,FALSE,useCubicInterp,surfaceReceiver);
                  break;
                  
                default:
                  assert(FALSE,
                         "addReceivers--grid only works with P, 3, 4, I , Vx, Vy, or Vz receivers not %i",
                         rtype);
              }
            }
          }
        }
        tEprintf(Verbose,"Added %ix%ix%i (%i) grid of %s receivers\n",
                 nx,ny,nz,nx*ny*nz,
                 ReceiverTypeNames[rtype]);
        tEprintf(Verbose,"  Grid spans: X %.2f to %.2f, Y %.2f to %.2f, Z %.2f to %.2f\n",
                 ix,ix+dx*(nx-1),iy,iy+dy*(ny-1),iz,iz+dz*(nz-1));
      }
        break;
      case 'f':
      {
        assert(argc>i+2,
               "addReceivers--adding receivers from file requires at 2 arguments (type,filename)");
        int lineCount=4;
        if(strlen(argv[i])>argOffset+1 && argv[i][argOffset+1]=='3'){
          lineCount=3;
        }
        
        int rtype=receiverType(argv[++i]);
        FILE* inFile=fopen(argv[++i],"r");
        
        assert(inFile!=NULL,
               "processArgs--unable to open receiver file %s (arg %i)",
               argv[i],i);
        
        for(;;){
          char buffer[1024],testChar=getc(inFile);
          if(testChar==EOF)
            break;
          
          if(testChar=='#'){
            //Discard rest of line and continue.
            fgets(buffer,1024,inFile);
            continue;
          }else{
            //Push back this char.
            ungetc(testChar,inFile);
          }
          
          float x,y,z,amp=1;
          int numRead;
          if(lineCount==3){
            numRead=fscanf(inFile,"%f %f %f",&x,&y,&z);
          }else{
            numRead=fscanf(inFile,"%f %f %f %f",&x,&y,&z,&amp);
          }
          assert(numRead==lineCount,
                 "Error reading receiver %i; read %i of %i elements from %s; ",
                 _receivers->size(),numRead,lineCount,argv[i]);
          
          switch(rtype){
            case VELOCITY_RECEIVER:
            {
              float bx,by,bz;
              fscanf(inFile,"%f %f %f",&bx,&by,&bz);
              addVelocityReceiver(modelDef,x,y,z,bx,by,bz,amp,integrate,FALSE,useCubicInterp,surfaceReceiver);
            }
              break;
            case PRESSURE_RECEIVER:
              _receivers->Add(new pressureReceiver(modelDef,FALSE,x,y,z,amp,useCubicInterp));
              break;
            case THREE_C_RECEIVER:
              add3CReceiver(modelDef,x,y,z,amp,integrate,FALSE,useCubicInterp,surfaceReceiver);
              break;
            case FOUR_C_RECEIVER:
              add4CReceiver(modelDef,x,y,z,amp,integrate,FALSE,useCubicInterp,surfaceReceiver);
              break;
              
            case VX_RECEIVER:
              addVelocityReceiver(modelDef,x,y,z,1.0,0.0,0.0,1.0,integrate,FALSE,useCubicInterp,surfaceReceiver);
              break;
            case VY_RECEIVER:
              addVelocityReceiver(modelDef,x,y,z,0.0,1.0,0.0,1.0,integrate,FALSE,useCubicInterp,surfaceReceiver);
              break;
            case VZ_RECEIVER:
              addVelocityReceiver(modelDef,x,y,z,0.0,0.0,1.0,1.0,integrate,FALSE,useCubicInterp,surfaceReceiver);
              break;
              
            default:
              assert(FALSE,
                     "addReceivers--grid only works with P, 3, 4, I , Vx, Vy, or Vz receivers not %i",
                     rtype);
          }
          fgets(buffer,1024,inFile);
        }
        fclose(inFile);
        
        tEprintf(Verbose," Added %s receivers from file; now %i\n",
                 ReceiverTypeNames[rtype],_receivers->size());
      }
        break;
      default:
        assert(FALSE,"addReceivers--unknown receiver flag \"%s\"",argv[i]);
    }
  }
  
  return _receivers->size();
}

int masterReceiverNetwork::addReceivers(masterReceiverNetwork* extraReceivers){
  if(!extraReceivers) return 0;
  
  strcpy(_traceOutputName,extraReceivers->_traceOutputName);
  _timeSubsample=extraReceivers->_timeSubsample;
  _splitFiles=extraReceivers->_splitFiles;
  
  int numRegularAdded=0;
  for(int i=0;i<extraReceivers->_receivers->size();i++,numRegularAdded++)
    _receivers->Add((*extraReceivers->_receivers)[i]);
  tEprintf(Verbose && numRegularAdded,
           "Added %i receivers to existing receiverNetwork\n",
           numRegularAdded);
  
  int numGridsAdded=0;
  for(int i=0;i<extraReceivers->_receiverGrids->size();i++,numGridsAdded++)
    _receiverGrids->Add((*extraReceivers->_receiverGrids)[i]);
  tEprintf(Verbose && numGridsAdded,
           "Added %i receiver-grids to existing receiverNetwork\n",
           numGridsAdded);
  
  //Manually zero the size of the arrays in extraReceivers so the elements are safe.
  extraReceivers->_receivers->EmptyArray();
  extraReceivers->_receiverGrids->EmptyArray();
  
  return numRegularAdded;
}

//Send the data for these receivers back and forth between master and slave process.
int masterReceiverNetwork::initializeNetwork(modelDefStruct* modelDef,int allocateData,
                                             int target){
  //Make sure the buffer is big enough, this has been causing problems for
  // small models with lots of receivers.
  DEF_MODEL_SIZE(modelDef);
  //uncommented setMessageBuffer is based on an actual count.  It will be less than NT+2 if
  //traces are output every n time steps, but n would have to be passed in and it can be
  //approx 1.5*n time steps between outputs at the end of a run (see mainLoopWork)
  //setMessageBuffer(4*(MAX_R_PER_SEND+10)*(NT+10));
  
  //Start the initialization.
  int numReceivers=_receivers->size();
  setMessageBuffer(MAX(MIN(MAX_R_PER_SEND,numReceivers)*44+100,(NT+2)*sizeof(float)+100));
  initSend();
  packMessage("ii",numReceivers,_timeSubsample);
  
  for(int i=0,currPacked=0;i<numReceivers;i++){
    (*_receivers)[i]->packInit();
    if(++currPacked>=MAX_R_PER_SEND){
      (target!=AllProcesses?iSendMessage:sendMessage)(target,MESSAGE_SEND_RECEIVER,NULL);
      initSend();
      currPacked=0;
    }
  }
  int totalReceivers=_receivers->size();
  (target!=AllProcesses?iSendMessage:sendMessage)(target,MESSAGE_SEND_RECEIVER,NULL);
  
  //Wait for an acknowlegement that all receivers have been received. This should also contain
  // an array of receivers that are active for this slave process.
  if(!_receivers->size()){
    for(int i=0;i<NumProcs;i++){
      initSend();
      sendMessage(Tids[i],MESSAGE_SEND_RECEIVER,NULL);
      getMessage(Tids[i],MESSAGE_SEND_RECEIVER,NULL);
    }
  }else{
    if(target==AllProcesses){
      if(!_receiverProcs){
        assert((_receiverProcs=(int*)malloc(numReceivers*sizeof(int)))!=NULL,
               "masterReceiverNetwork::initializeNetwork--unable to allocate %i ints for receiverProcs",
               numReceivers);
        for(int i=0;i<numReceivers;_receiverProcs[i++]=-1);
      }
      
      int *temp;
      assert((temp=(int*)malloc(numReceivers*sizeof(int)))!=NULL,
             "masterReceiverNetwork::initializeNetwork--unable to allocate %i ints for temp receiverProcs",
             numReceivers);
      for(int i=0;i<NumProcs;i++){
        initSend();
        sendMessage(Tids[i],MESSAGE_SEND_RECEIVER,NULL);
        for(int j=0;j<numReceivers;j+=MAX_R_PER_SEND){
          int *ptr=temp+j,n=MIN(MAX_R_PER_SEND,numReceivers-j);
          getMessage(Tids[i],MESSAGE_SEND_RECEIVER,"I",
                     ptr,n);
        }
        for(int j=0;j<numReceivers;j++){
          if(_receiverProcs[j]<0 && temp[j])
            _receiverProcs[j]=Tids[i];
        }
      }
      free(temp);
    }
  }
  
  //Now send the receiver grids.
  int n=_receiverGrids->size();
  initSend();
  sendMessage(target,MESSAGE_SEND_RECEIVER_GRID,"i",n);
  for(int i=0;i<_receiverGrids->size();i++){
    initSend();
    (*_receiverGrids)[i]->packInit();
    sendMessage(target,MESSAGE_SEND_RECEIVER_GRID+i+1,NULL);
    for(int j=0;j<NumProcs;j++){
      getMessage(AllProcesses,MESSAGE_SEND_RECEIVER_GRID+i+1,NULL);
    }
  }
  
  
  return totalReceivers;
}

int masterReceiverNetwork::writeTraces(modelDefStruct* modelDef,int startI,int stopI,
                                       int useLocalData){
  if(_receiverGrids && _receiverGrids->size()){
    tEprintf(Verbose,"Writing %i receiver-grids\n",
             _receiverGrids->size());
    writeSingleGrid(modelDef);
  }
  
  if(_splitFiles==1){
    return writeCompSplitTraces(modelDef,startI,stopI,useLocalData);
  }else if(_splitFiles){
    return writeIndexSplitTraces(modelDef,startI,stopI,useLocalData);
  }
  return writeSingleTraces(modelDef,startI,stopI,useLocalData);
}

//Get the data for a specific receiver.
int masterReceiverNetwork::getReceiverData(int index,int startI,int stopI,float* data){
  initSend();
  sendMessage(_receiverProcs[index],MESSAGE_SEND_RECEIVER,"iii",
              index,startI,stopI);
  
  int firstI=(int)((float)(startI)/_timeSubsample);
  int nSamps=(int)((float)(stopI)/_timeSubsample)-firstI+1;
  getMessage(_receiverProcs[index],MESSAGE_SEND_RECEIVER,"F",
             data,nSamps);
  return stopI-startI;
}

///Get data for a specific receiver-grid from a specific process.
int masterReceiverNetwork::getReceiverGridData(int index,int tid,int isXcorr,
                                               int &n,floatPtr& data,
                                               size_t start[4],size_t count[4]){
  initSend();
  sendMessage(tid,MESSAGE_SEND_RECEIVER_GRID,"i",index);
  
  int nt,ix,nx,iy,ny,iz,nz;
  getMessage(tid,MESSAGE_SEND_RECEIVER_GRID,"i ii ii ii",
             &nt,&ix,&nx,&iy,&ny,&iz,&nz);
  start[0]=iz;start[1]=iy;start[2]=ix;start[3]=0;
  count[0]=nz;count[1]=ny;count[2]=nx;count[3]=nt;
  
  int nxyz=nx*ny*nz,na;
  if(isXcorr){
    na=nxyz;
  }else{
    na=nxyz*nt;
  }
  if(na>n){
    n=na;
    if(!data){
      assert((data=(float*)malloc(na*sizeof(float)))!=NULL,
             "getReceiverGridData--unable to allocate data for a %ix%ix%ix%i grid\n",
             nx,ny,nz,nt);
    }else{
      assert((data=(float*)realloc(data,na*sizeof(float)))!=NULL,
             "getReceiverGridData--unable to reallocate data for a %ix%ix%ix%i grid\n",
             nx,ny,nz,nt);
    }
  }
  if(nxyz){
    if(isXcorr){
      int test;
      unpackMessage("Fi",data,nxyz,&test);
      assert(test==index,
             "getReceiverGridData--receiver grid index mismatch %i!=%i.",
             test,index);
    }else{
      int currSend=1;
      getMessage(tid,MESSAGE_SEND_RECEIVER_GRID+currSend++,NULL);
      for(int i=0,currPacked=0;i<nxyz;i++){
        int rn;
        unpackMessage("iF",&rn,data+i*nt,nt);
        assert(rn==i,
               "getReceiverGridData--receiver number mismatch at receiver %i, tid %i",
               i,tid);
        if(++currPacked==MAX_R_PER_SEND && i+1<nxyz){
          currPacked=0;
          getMessage(tid,MESSAGE_SEND_RECEIVER_GRID+currSend++,NULL);
        }
      }
    }
  }
  return nxyz*nt;
}

int masterReceiverNetwork::writeSingleHeader(modelDefStruct* modelDef,int defData){
  if(!_receivers->size()) return 0;
  DEF_MODEL_SIZE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  
  //Create the file with a standard header.
  writeCDFHeader(_traceOutputName,
                 NX,NY,NZ,NT,NULL,
                 dx,minX,dy,minY,dz,minZ,dt,minT,
                 0,NULL);
  
  //Fill in the receiver locations.
  writeSingleLocations(modelDef,defData,-1);
  tEprintf(Verbose,"Wrote header for %i receivers to\n\t%s\n",
           _receivers->size(),_traceOutputName);
  
  return _receivers->size();
}

int masterReceiverNetwork::writeSingleLocations(modelDefStruct* modelDef,
                                                int defData,int outFile){
  int iOpenedFile=FALSE;
  if(outFile<0){
    iOpenedFile=TRUE;
    outFile=openCDFFile(_traceOutputName,FALSE,NC_WRITE);
  }
  
  //get required dimension
  int ntDim;
  assert(nc_inq_dimid(outFile,"NT",&ntDim)==NC_NOERR,
         "writeSingleLocations--unable to open dimension NT from file %s",
         _traceOutputName);
  
  //Define a dimension for the decimation factor and another for
  // the actual number of samples in each trace.
  nc_redef(outFile);
  int receiverSkipDim,nSamplesDim;
  assert(nc_def_dim(outFile,"receiverDecimate",_timeSubsample,
                    &receiverSkipDim)==NC_NOERR,
         "writeSingleLocations--unable to create dimension receiverDecimate in %s",
         _traceOutputName);
  assert(nc_def_dim(outFile,"nSamples",nSamples(modelDef),
                    &nSamplesDim)==NC_NOERR,
         "writeSingleLocations--unable to create dimension nSamples in %s",
         _traceOutputName);
  
  //Define a dimension and variables for each type
  int receiversDim;
  assert(nc_def_dim(outFile,"numReceivers",_receivers->size(),&receiversDim)==NC_NOERR,
         "writeSingleLocations--unable to create dimensioon numReceivers in %s",
         _traceOutputName);
  
  int typeVar,ampVar;
  assert(nc_def_var(outFile,"receiverType",NC_FLOAT,1,
                    &receiversDim,&typeVar)==NC_NOERR,
         "writeSingleLocations--unable to define variable receiverType in file %s(%i,%i)",
         _traceOutputName,outFile,receiversDim);
  assert(nc_def_var(outFile,"receiverAmp",NC_FLOAT,1,
                    &receiversDim,&ampVar)==NC_NOERR,
         "writeSingleLocations--unable to define variable receiverAmp in file %s(%i,%i)",
         _traceOutputName,outFile,receiversDim);
  
  int xVar,yVar,zVar;
  int bxVar,byVar,bzVar;
  assert(nc_def_var(outFile,"receiverX",NC_FLOAT,1,
                    &receiversDim,&xVar)==NC_NOERR,
         "writeSingleLocations--unable to create variable receiverX in %s(%i,%i)",
         _traceOutputName,outFile,receiversDim);
  assert(nc_def_var(outFile,"receiverY",NC_FLOAT,1,
                    &receiversDim,&yVar)==NC_NOERR,
         "writeSingleLocations--unable to create variable receiverY in %s(%i,%i)",
         _traceOutputName);
  assert(nc_def_var(outFile,"receiverZ",NC_FLOAT,1,
                    &receiversDim,&zVar)==NC_NOERR,
         "writeSingleLocations--unable to create variable receiverZ in %s(%i,%i)",
         _traceOutputName,outFile,receiversDim);
  
  assert(nc_def_var(outFile,"receiverBx",NC_FLOAT,1,
                    &receiversDim,&bxVar)==NC_NOERR,
         "writeSingleLocations--unable to create variable receiverBx in %s(%i,%i)",
         _traceOutputName,outFile,receiversDim);
  assert(nc_def_var(outFile,"receiverBy",NC_FLOAT,1,
                    &receiversDim,&byVar)==NC_NOERR,
         "writeSingleLocations--unable to create variable receiverBy in %s(%i,%i)",
         _traceOutputName,outFile,receiversDim);
  assert(nc_def_var(outFile,"receiverBz",NC_FLOAT,1,
                    &receiversDim,&bzVar)==NC_NOERR,
         "writeSingleLocations--unable to create variable receiverBz in %s(%i,%i)",
         _traceOutputName,outFile,receiversDim);
  
  int integrateVar;
  assert(nc_def_var(outFile,"receiverIntegrate",NC_INT,1,
                    &receiversDim,&integrateVar)==NC_NOERR,
         "writeSingleLocations--unable to create variable receiverIntegrate in %s(%i)",
         _traceOutputName,outFile);
  
  if(defData){
    int receiverDims[2]={receiversDim,nSamplesDim};
    int dataVar;
    nc_enddef(outFile);
    nc_redef(outFile);
    assert(nc_def_var(outFile,"receiverData",NC_FLOAT,2,
                      receiverDims,&dataVar)==NC_NOERR,
           "writeSingleLocations--unable to create variable receiverData in %s(%i,%i,%i)",
           _traceOutputName,outFile,receiverDims[0],receiverDims[1]);
  }
  
  //Fill the receiver variables.
  nc_enddef(outFile);
  
  float *tempflt;
  assert((tempflt=(float*)malloc(_receivers->size()*sizeof(float)))!=NULL,
         "writeSingleLocations--unable to allocate %i floats",
         _receivers->size());
  
  //X locations
  for(int i=0;i<_receivers->size();i++)
    tempflt[i]=(*_receivers)[i]->_x;
  assert(nc_put_var_float(outFile,xVar,tempflt)==NC_NOERR,
         "writeSingleLocations--unable to write variable xr in %s(%i,%i)",
         _traceOutputName,outFile,xVar);
  
  //Y locations
  for(int i=0;i<_receivers->size();i++)
    tempflt[i]=(*_receivers)[i]->_y;
  assert(nc_put_var_float(outFile,yVar,tempflt)==NC_NOERR,
         "writeSingleLocations--unable to write variable yr in %s(%i,%i)",
         _traceOutputName,outFile,yVar);
  
  //Z locations
  for(int i=0;i<_receivers->size();i++)
    tempflt[i]=(*_receivers)[i]->_z;
  assert(nc_put_var_float(outFile,zVar,tempflt)==NC_NOERR,
         "writeSingleLocations--unable to write variable zr in %s(%i,%i)",
         _traceOutputName,outFile,zVar);
  
  //X orientations
  for(int i=0;i<_receivers->size();i++)
    tempflt[i]=(*_receivers)[i]->_bx;
  assert(nc_put_var_float(outFile,bxVar,tempflt)==NC_NOERR,
         "writeSingleLocations--unable to write variable bx in %s(%i,%i)",
         _traceOutputName,outFile,xVar);
  
  //Y orientations
  for(int i=0;i<_receivers->size();i++)
    tempflt[i]=(*_receivers)[i]->_by;
  assert(nc_put_var_float(outFile,byVar,tempflt)==NC_NOERR,
         "writeSingleLocations--unable to write variable by in %s(%i,%i)",
         _traceOutputName,outFile,yVar);
  
  //Z orientations
  for(int i=0;i<_receivers->size();i++)
    tempflt[i]=(*_receivers)[i]->_bz;
  assert(nc_put_var_float(outFile,bzVar,tempflt)==NC_NOERR,
         "writeSingleLocations--unable to write variable bz in %s(%i,%i)",
         _traceOutputName,outFile,zVar);
  
  //Amplitude factor
  for(int i=0;i<_receivers->size();i++)
    tempflt[i]=(*_receivers)[i]->_amp;
  assert(nc_put_var_float(outFile,ampVar,tempflt)==NC_NOERR,
         "writeSingleLocations--unable to write variable amp in %s(%i,%i)",
         _traceOutputName,outFile,ampVar);
  
  //And the integer variables.
  free(tempflt);
  int *tempint;
  assert((tempint=(int*)malloc(_receivers->size()*sizeof(int)))!=NULL,
         "writeSingleLocations--unable to allocate %i ints",
         _receivers->size());
  
  //Type
  for(int i=0;i<_receivers->size();i++)
    tempint[i]=(*_receivers)[i]->type();
  assert(nc_put_var_int(outFile,typeVar,tempint)==NC_NOERR,
         "writeSingleLocations--unable to write variable type in %s(%i,%i)",
         _traceOutputName,outFile,typeVar);
  
  //Integrate
  for(int i=0;i<_receivers->size();i++)
    tempint[i]=(*_receivers)[i]->integrate();
  assert(nc_put_var_int(outFile,integrateVar,tempint)==NC_NOERR,
         "writeSingleLocations--unable to write variable integrate in %s(%i,%i)",
         _traceOutputName,outFile,integrateVar);
  
  free(tempint);
  
  if(iOpenedFile)
    assert(nc_close(outFile)==NC_NOERR,
           "writeSingleLocations--unable to close file %s(%i)",
           _traceOutputName,outFile);
  
  return size();
}
int masterReceiverNetwork::rewriteSingleLocations(){
  int outFile=openCDFFile(_traceOutputName,FALSE,NC_WRITE);
  
  int xVar,yVar,zVar;
  assert(nc_inq_varid(outFile,"receiverX",&xVar)==NC_NOERR,
         "rewriteSingleLocations--unable to get variable receiverX in %s(%i)",
         _traceOutputName,outFile);
  assert(nc_inq_varid(outFile,"receiverY",&yVar)==NC_NOERR,
         "rewriteSingleLocations--unable to get variable receiverY in %s(%i)",
         _traceOutputName,outFile);
  assert(nc_inq_varid(outFile,"receiverZ",&zVar)==NC_NOERR,
         "rewriteSingleLocations--unable to get variable receiverZ in %s(%i)",
         _traceOutputName,outFile);
  
  for(size_t i=0;i<_receivers->size();i++){
    assert(nc_put_var1_float(outFile,xVar,&i,
                             &(*_receivers)[i]->_x)==NC_NOERR,
           "writeSingleLocations--unable to write variable xr in %s(%i,%i)",
           _traceOutputName,outFile,xVar);
    assert(nc_put_var1_float(outFile,yVar,&i,
                             &(*_receivers)[i]->_y)==NC_NOERR,
           "writeSingleLocations--unable to write variable yr in %s(%i,%i)",
           _traceOutputName,outFile,yVar);
    assert(nc_put_var1_float(outFile,zVar,&i,
                             &(*_receivers)[i]->_z)==NC_NOERR,
           "writeSingleLocations--unable to write variable zr in %s(%i,%i)",
           _traceOutputName,outFile,zVar);
  }
  
  assert(nc_close(outFile)==NC_NOERR,
         "Unable to close file %s(%i)",
         _traceOutputName,outFile);
  
  return size();
}

int masterReceiverNetwork::writeSingleTraces(modelDefStruct* modelDef,int startI,int stopI,
                                             int useLocalData){
  if(!_receivers->size()) return 0;
  
  //Calculate the sub-sampled trace indicies.
  int startSample=(int)floor((float)(startI)/_timeSubsample);
  int stopSample=MIN(nSamples(modelDef)-1,
                     (int)floor((float)(stopI)/_timeSubsample));
  int nSamps=stopSample-startSample+1;
  
  tEprintf(Verbose,
           "Writing trace output to file\n\t %s\n",_traceOutputName);
  if(_timeSubsample==1){
    tEprintf(Verbose,
             "\tIteration %i to %i\n",
             startI,stopI);
  }else{
    tEprintf(Verbose,
             "\tIteration %i to %i; Subsampled indicies %i to %i\n",
             startI,stopI,startSample,stopSample);
  }
  
  //open the file and get the data variable
  int dataVar,outFile=openCDFFile(_traceOutputName,FALSE,NC_WRITE);
  assert(nc_inq_varid(outFile,"receiverData",&dataVar)==NC_NOERR,
         "writeSingleTraces--unable to open var data in %s(%i)",
         _traceOutputName,outFile);
  
  //Now fill in the data variables.
  for(size_t i=0;i<_receivers->size();i++){
    receiver* currR=(*_receivers)[i];
    
    //where to write this data
    size_t start[2]={i,startSample};
    size_t count[2]={1,nSamps};
    
    if(useLocalData || _receiverProcs[i]==Parent){
      currR->scaleData(modelDef,startI,stopI,currR->_data,
                       work(modelDef));
    }else{
      getReceiverData(i,startI,stopI,work(modelDef));
    }
    assert(nc_put_vara_float(outFile,dataVar,
                             start,count,work(modelDef))==NC_NOERR,
           "writeSingleTraces--unable to write data[%i,%i,%i,%i] to %s(%i,%i)",
           start[0],start[1],count[0],count[1],
           _traceOutputName,outFile,dataVar);
    
    if(i && !(i%2000)){
      tEprintf(Verbose,
               "\tCompleted receiver %i\n",
               i);
    }
  }
  
  //close the file and exit
  assert(nc_close(outFile)!=-1,
         "writeSingleTraces--nc_close of %s unsuccessful",
         _traceOutputName);
  return _receivers->size()*nSamps;
}

char* masterReceiverNetwork::generateSplitFilename(char* target,const char* tag){
  //Generate the correct file name for this split file.
  strcpy(target,_traceOutputName);
  
  // Rremove trailing .cdf if present, make sure we get the last one if there are
  //  multiple .cdf's in the filename.
  char *endPtr=target;
  for(;;){
    char* tEndPtr;
    if(!(tEndPtr=strstr(endPtr+1,".cdf"))){
      break; //no more substrings
    }else{
      //found substring
      endPtr=tEndPtr;
    }
  }
  //check for no substrings
  if(endPtr==target)
    endPtr+=strlen(target);
  
  //add correct suffix
  sprintf(endPtr,"_%s.cdf",tag);
  
  return target;
}

//Methods to split files by component.
int masterReceiverNetwork::writeCompSplitHeader(modelDefStruct* modelDef,int defData){
  if(!_receivers->size()) return 0;
  DEF_MODEL_SIZE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  
  //Determine the number of each type of receiver.
  _nVxR=0;
  _nVyR=0;
  _nVzR=0;
  _nPrR=0;
  assert((_rTypeFlag=(int*)malloc(_receivers->size()*sizeof(int)))!=NULL,
         "writeCompSplitHeader--unable to allocate %i ints for rTypeFlag",
         _receivers->size());
  
  for(int i=0;i<_receivers->size();i++){
    receiver* curr=(*_receivers)[i];
    if(curr->type()==VELOCITY_RECEIVER){
      //Determine the dominant type.
      if(curr->_bx>0.5){
        _nVxR++;
        _rTypeFlag[i]=1;
      }else if(curr->_by>0.5){
        _nVyR++;
        _rTypeFlag[i]=2;
      }else if(curr->_bz>0.5){
        _nVzR++;
        _rTypeFlag[i]=3;
      }else{
        assert(FALSE,
               "writeCompSplitHeader--receiver %i: type is velocity can not group %.1f, %.1f, %.1f",
               i,curr->_bx,curr->_by,curr->_bz);
      }
    }else if(curr->type()==PRESSURE_RECEIVER){
      _nPrR++;
      _rTypeFlag[i]=4;
    }else{
      assert(FALSE,
             "writeCompSplitHeader--receiver %i: type %i not supported for split writing",
             i,curr->type());
    }
  }
  tEprintf(Verbose,
           "Writing Split Header for %i receivers: Vx %i, Vy %i, Vz %i, Pressure %i\n",
           _receivers->size(),_nVxR,_nVyR,_nVzR,_nPrR);
  
  char buffer[1024];
  if(_nVxR){
    //Generate the correct file name for this split file.
    generateSplitFilename(buffer,"Vx");
    
    //Create the file with a standard header.
    writeCDFHeader(buffer,
                   NX,NY,NZ,NT,NULL,
                   dx,minX,dy,minY,dz,minZ,dt,minT,
                   0,NULL);
    
    //Fill in the receiver locations.
    int count=writeCompSplitLocations(buffer,_nVxR,1,_rTypeFlag,
                                      modelDef,defData,-1);
    tEprintf(Verbose," Wrote header for %i receivers to\n   %s\n",
             count,buffer);
  }
  if(_nVyR){
    //Generate the correct file name for this split file.
    generateSplitFilename(buffer,"Vy");
    
    //Create the file with a standard header.
    writeCDFHeader(buffer,
                   NX,NY,NZ,NT,NULL,
                   dx,minX,dy,minY,dz,minZ,dt,minT,
                   0,NULL);
    
    //Fill in the receiver locations.
    int count=writeCompSplitLocations(buffer,_nVyR,2,_rTypeFlag,
                                      modelDef,defData,-1);
    tEprintf(Verbose," Wrote header for %i receivers to\n   %s\n",
             count,buffer);
  }
  if(_nVzR){
    //Generate the correct file name for this split file.
    generateSplitFilename(buffer,"Vz");
    
    //Create the file with a standard header.
    writeCDFHeader(buffer,
                   NX,NY,NZ,NT,NULL,
                   dx,minX,dy,minY,dz,minZ,dt,minT,
                   0,NULL);
    
    //Fill in the receiver locations.
    int count=writeCompSplitLocations(buffer,_nVzR,3,_rTypeFlag,
                                      modelDef,defData,-1);
    tEprintf(Verbose," Wrote header for %i receivers to\n   %s\n",
             count,buffer);
  }
  if(_nPrR){
    //Generate the correct file name for this split file.
    generateSplitFilename(buffer,"Pressure");
    
    //Create the file with a standard header.
    writeCDFHeader(buffer,
                   NX,NY,NZ,NT,NULL,
                   dx,minX,dy,minY,dz,minZ,dt,minT,
                   0,NULL);
    
    //Fill in the receiver locations.
    int count=writeCompSplitLocations(buffer,_nPrR,4,_rTypeFlag,
                                      modelDef,defData,-1);
    tEprintf(Verbose," Wrote header for %i receivers to\n   %s\n",
             count,buffer);
  }
  
  return _receivers->size();
}
int masterReceiverNetwork::writeCompSplitLocations(char* filename,int count,
                                                   int tag,int *tagArray,
                                                   modelDefStruct* modelDef,
                                                   int defData,int outFile){
  int iOpenedFile=FALSE;
  if(outFile<0){
    iOpenedFile=TRUE;
    outFile=openCDFFile(filename,FALSE,NC_WRITE);
  }
  
  //get required dimension
  int ntDim;
  assert(nc_inq_dimid(outFile,"NT",&ntDim)==NC_NOERR,
         "writeCompSplitLocations--unable to open dimension NT from file %s",
         filename);
  
  //Define a dimension for the decimation factor and another for
  // the actual number of samples in each trace.
  nc_redef(outFile);
  int receiverSkipDim,nSamplesDim;
  assert(nc_def_dim(outFile,"receiverDecimate",_timeSubsample,
                    &receiverSkipDim)==NC_NOERR,
         "writeCompSplitLocations--unable to create dimension receiverDecimate in %s",
         filename);
  assert(nc_def_dim(outFile,"nSamples",nSamples(modelDef),
                    &nSamplesDim)==NC_NOERR,
         "writeCompSplitLocations--unable to create dimension nSamples in %s",
         filename);
  
  //Define a dimension and variables for each type
  int receiversDim;
  assert(nc_def_dim(outFile,"numReceivers",count,&receiversDim)==NC_NOERR,
         "writeCompSplitLocations--unable to create dimension numReceivers in %s",
         filename);
  
  int typeVar,ampVar;
  assert(nc_def_var(outFile,"receiverType",NC_FLOAT,1,
                    &receiversDim,&typeVar)==NC_NOERR,
         "writeCompSplitLocations--unable to define variable receiverType in file %s(%i,%i)",
         filename,outFile,receiversDim);
  assert(nc_def_var(outFile,"receiverAmp",NC_FLOAT,1,
                    &receiversDim,&ampVar)==NC_NOERR,
         "writeCompSplitLocations--unable to define variable receiverAmp in file %s(%i,%i)",
         filename,outFile,receiversDim);
  
  int xVar,yVar,zVar;
  int bxVar,byVar,bzVar;
  assert(nc_def_var(outFile,"receiverX",NC_FLOAT,1,
                    &receiversDim,&xVar)==NC_NOERR,
         "writeCompSplitLocations--unable to create variable receiverX in %s(%i,%i)",
         filename,outFile,receiversDim);
  assert(nc_def_var(outFile,"receiverY",NC_FLOAT,1,
                    &receiversDim,&yVar)==NC_NOERR,
         "writeCompSplitLocations--unable to create variable receiverY in %s(%i,%i)",
         filename);
  assert(nc_def_var(outFile,"receiverZ",NC_FLOAT,1,
                    &receiversDim,&zVar)==NC_NOERR,
         "writeCompSplitLocations--unable to create variable receiverZ in %s(%i,%i)",
         filename,outFile,receiversDim);
  
  assert(nc_def_var(outFile,"receiverBx",NC_FLOAT,1,
                    &receiversDim,&bxVar)==NC_NOERR,
         "writeCompSplitLocations--unable to create variable receiverBx in %s(%i,%i)",
         filename,outFile,receiversDim);
  assert(nc_def_var(outFile,"receiverBy",NC_FLOAT,1,
                    &receiversDim,&byVar)==NC_NOERR,
         "writeCompSplitLocations--unable to create variable receiverBy in %s(%i,%i)",
         filename,outFile,receiversDim);
  assert(nc_def_var(outFile,"receiverBz",NC_FLOAT,1,
                    &receiversDim,&bzVar)==NC_NOERR,
         "writeCompSplitLocations--unable to create variable receiverBz in %s(%i,%i)",
         filename,outFile,receiversDim);
  
  int integrateVar;
  assert(nc_def_var(outFile,"receiverIntegrate",NC_INT,1,
                    &receiversDim,&integrateVar)==NC_NOERR,
         "writeCompSplitLocations--unable to create variable receiverIntegrate in %s(%i)",
         filename,outFile);
  
  if(defData){
    int receiverDims[2]={receiversDim,nSamplesDim};
    int dataVar;
    assert(nc_def_var(outFile,"receiverData",NC_FLOAT,2,
                      receiverDims,&dataVar)==NC_NOERR,
           "writeCompSplitLocations--unable to create variable receiverData in %s",
           filename);
  }
  
  //Fill the receiver variables.
  nc_enddef(outFile);
  
  float *tempflt;
  assert((tempflt=(float*)malloc(count*sizeof(float)))!=NULL,
         "writeCompSplitLocations--unable to allocate %i floats",
         count);
  
  //X locations
  for(int i=0,currI=0;i<_receivers->size();i++){
    if(tag==tagArray[i]){
      tempflt[currI++]=(*_receivers)[i]->_x;
    }
  }
  assert(nc_put_var_float(outFile,xVar,tempflt)==NC_NOERR,
         "writeCompSplitLocations--unable to write variable xr in %s(%i,%i)",
         filename,outFile,xVar);
  
  //Y locations
  for(int i=0,currI=0;i<_receivers->size();i++){
    if(tag==tagArray[i]){
      tempflt[currI++]=(*_receivers)[i]->_y;
    }
  }
  assert(nc_put_var_float(outFile,yVar,tempflt)==NC_NOERR,
         "writeCompSplitLocations--unable to write variable yr in %s(%i,%i)",
         filename,outFile,yVar);
  
  //Z locations
  for(int i=0,currI=0;i<_receivers->size();i++){
    if(tag==tagArray[i]){
      tempflt[currI++]=(*_receivers)[i]->_z;
    }
  }
  assert(nc_put_var_float(outFile,zVar,tempflt)==NC_NOERR,
         "writeCompSplitLocations--unable to write variable zr in %s(%i,%i)",
         filename,outFile,zVar);
  
  //X orientations
  for(int i=0,currI=0;i<_receivers->size();i++){
    if(tag==tagArray[i]){
      tempflt[currI++]=(*_receivers)[i]->_bx;
    }
  }
  assert(nc_put_var_float(outFile,bxVar,tempflt)==NC_NOERR,
         "writeCompSplitLocations--unable to write variable bx in %s(%i,%i)",
         filename,outFile,xVar);
  
  //Y orientations
  for(int i=0,currI=0;i<_receivers->size();i++){
    if(tag==tagArray[i]){
      tempflt[currI++]=(*_receivers)[i]->_by;
    }
  }
  assert(nc_put_var_float(outFile,byVar,tempflt)==NC_NOERR,
         "writeCompSplitLocations--unable to write variable by in %s(%i,%i)",
         filename,outFile,yVar);
  
  //Z orientations
  for(int i=0,currI=0;i<_receivers->size();i++){
    if(tag==tagArray[i]){
      tempflt[currI++]=(*_receivers)[i]->_bz;
    }
  }
  assert(nc_put_var_float(outFile,bzVar,tempflt)==NC_NOERR,
         "writeCompSplitLocations--unable to write variable bz in %s(%i,%i)",
         filename,outFile,zVar);
  
  //Amplitude factor
  for(int i=0,currI=0;i<_receivers->size();i++){
    if(tag==tagArray[i]){
      tempflt[currI++]=(*_receivers)[i]->_amp;
    }
  }
  assert(nc_put_var_float(outFile,ampVar,tempflt)==NC_NOERR,
         "writeCompSplitLocations--unable to write variable amp in %s(%i,%i)",
         filename,outFile,ampVar);
  
  //And the integer variables.
  free(tempflt);
  int *tempint;
  assert((tempint=(int*)malloc(count*sizeof(int)))!=NULL,
         "writeCompSplitLocations--unable to allocate %i ints",
         count);
  
  //Type
  for(int i=0,currI=0;i<_receivers->size();i++){
    if(tag==tagArray[i]){
      tempint[currI++]=(*_receivers)[i]->type();
    }
  }
  assert(nc_put_var_int(outFile,typeVar,tempint)==NC_NOERR,
         "writeCompSplitLocations--unable to write variable type in %s(%i,%i)",
         filename,outFile,typeVar);
  
  //Integrate
  for(int i=0,currI=0;i<_receivers->size();i++){
    if(tag==tagArray[i]){
      tempint[currI++]=(*_receivers)[i]->integrate();
    }
  }
  assert(nc_put_var_int(outFile,integrateVar,tempint)==NC_NOERR,
         "writeCompSplitLocations--unable to write variable integrate in %s(%i,%i)",
         filename,outFile,integrateVar);
  
  free(tempint);
  
  if(iOpenedFile)
    assert(nc_close(outFile)==NC_NOERR,
           "writeCompSplitLocations--unable to close file %s(%i)",
           filename,outFile);
  
  return count;
}
int masterReceiverNetwork::writeCompSplitTraces(modelDefStruct* modelDef,int iteration,int nSamps,
                                                int useLocalData){
  char buffer[1024];
  tEprintf(Verbose,
           " Writing component split trace output (%i receivers to iteration %i)\n",
           _receivers->size(),iteration);
  if(_nVxR){
    //Generate the correct file name for this split file.
    generateSplitFilename(buffer,"Vx");
    writeCompSplitTraces(buffer,_nVxR,
                         1,_rTypeFlag,
                         modelDef,iteration,nSamps,useLocalData);
  }
  if(_nVyR){
    //Generate the correct file name for this split file.
    generateSplitFilename(buffer,"Vy");
    writeCompSplitTraces(buffer,_nVyR,
                         2,_rTypeFlag,
                         modelDef,iteration,nSamps,useLocalData);
  }
  if(_nVzR){
    //Generate the correct file name for this split file.
    generateSplitFilename(buffer,"Vz");
    writeCompSplitTraces(buffer,_nVzR,
                         3,_rTypeFlag,
                         modelDef,iteration,nSamps,useLocalData);
  }
  if(_nPrR){
    //Generate the correct file name for this split file.
    generateSplitFilename(buffer,"Pressure");
    writeCompSplitTraces(buffer,_nPrR,
                         4,_rTypeFlag,
                         modelDef,iteration,nSamps,useLocalData);
  }
  return _receivers->size();
}
int masterReceiverNetwork::writeCompSplitTraces(char* filename,int count,
                                                int tag,int *tagArray,
                                                modelDefStruct* modelDef,
                                                int startI,int stopI,
                                                int useLocalData){
  tEprintf(Verbose,
           " %i receivers in file: %s\n",
           count,filename);
  
  //open the file and get the data variable
  int dataVar,outFile=openCDFFile(filename,FALSE,NC_WRITE);
  assert(nc_inq_varid(outFile,"receiverData",&dataVar)==NC_NOERR,
         "writeCompSplitTraces--unable to open var data in %s(%i)",
         filename,outFile);
  
  //Now fill in the data variables
  //Calculate the sub-sampled trace indicies.
  int startSample=(int)floor((float)(startI)/_timeSubsample);
  int stopSample=MIN(nSamples(modelDef)-1,
                     (int)floor((float)(stopI)/_timeSubsample));
  int nSamps=stopSample-startSample+1;
  
  for(size_t i=0,currI=0;i<_receivers->size();i++){
    if(tagArray[i]==tag){
      receiver* currR=(*_receivers)[i];
      
      //where to write this data
      size_t start[2]={currI++,startSample};
      size_t count[2]={1,nSamps};
      
      if(useLocalData){
        currR->scaleData(modelDef,startI,stopI,currR->_data,work(modelDef));
      }else{
        getReceiverData(i,startI,stopI,work(modelDef));
      }
      assert(nc_put_vara_float(outFile,dataVar,
                               start,count,work(modelDef))==NC_NOERR,
             "writeCompSplitTraces--unable to write data[%i,%,%i,%i] to %s(%i,%i)",
             start[0],start[1],count[0],count[1],
             filename,outFile,dataVar);
      
      if(i && !(i%2000)){
        tEprintf(Verbose,
                 "\tCompleted receiver %i\n",
                 i);
      }
    }
  }
  
  //close the file and exit
  assert(nc_close(outFile)!=-1,
         "writeCompSplitTraces--nc_close of %s unsuccessful",
         filename);
  return count*nSamps;
}

//Methods to split files by index.
int masterReceiverNetwork::writeIndexSplitHeader(modelDefStruct* modelDef,int defData){
  if(!_receivers->size()) return 0;
  DEF_MODEL_SIZE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  
  int nFiles=(int)ceil((float)_receivers->size()/(float)_splitFiles);
  tEprintf(Verbose,
           "Writing Split Header for %i receivers with %i/file, %i files\n",
           _receivers->size(),_splitFiles,nFiles);
  
  char buffer[1024],tag[10];
  for(int i=0;i<nFiles;i++){
    sprintf(tag,"Seg%02i",i+1);
    generateSplitFilename(buffer,tag);
    
    //Determine the start and end receiver index for this segment.
    int start=i*_splitFiles,end=(i+1)*_splitFiles-1;
    end=MIN(end,_receivers->size()-1);
    
    //Create the file with a standard header.
    writeCDFHeader(buffer,
                   NX,NY,NZ,NT,NULL,
                   dx,minX,dy,minY,dz,minZ,dt,minT,
                   0,NULL);
    
    //Fill in the receiver locations.
    writeIndexSplitLocations(buffer,start,end,
                             modelDef,defData,-1);
    tEprintf(Verbose," Wrote header for receivers %i to %i (%i) to file\n   %s\n",
             start,end,end-start+1,buffer);
  }
  
  return _receivers->size();
}
int masterReceiverNetwork::writeIndexSplitLocations(char* filename,
                                                    int start,int end,
                                                    modelDefStruct* modelDef,
                                                    int defData,int outFile){
  int iOpenedFile=FALSE;
  if(outFile<0){
    iOpenedFile=TRUE;
    outFile=openCDFFile(filename,FALSE,NC_WRITE);
  }
  
  //get required dimension
  int ntDim;
  assert(nc_inq_dimid(outFile,"NT",&ntDim)==NC_NOERR,
         "writeIndexSplitLocations--unable to open dimension NT from file %s",
         filename);
  
  //Define a dimension for the decimation factor and another for
  // the actual number of samples in each trace.
  nc_redef(outFile);
  int receiverSkipDim,nSamplesDim;
  assert(nc_def_dim(outFile,"receiverDecimate",_timeSubsample,
                    &receiverSkipDim)==NC_NOERR,
         "writeIndexSplitLocations--unable to create dimension receiverDecimate in %s",
         filename);
  assert(nc_def_dim(outFile,"nSamples",nSamples(modelDef),
                    &nSamplesDim)==NC_NOERR,
         "writeIndexSplitLocations--unable to create dimension nSamples in %s",
         filename);
  
  //Define two variables for the start and end index of the receivers in this file.
  int zero=0;
  int startIndexVar,endIndexVar;
  assert(nc_def_var(outFile,"startIndex",NC_INT,0,&zero,&startIndexVar)==NC_NOERR,
         "writeIndexSplitLocations--unable to define variable startIndex in file %s(%i)",
         filename,outFile);
  assert(nc_def_var(outFile,"endIndex",NC_INT,0,&zero,&endIndexVar)==NC_NOERR,
         "writeIndexSplitLocations--unable to define variable endIndex in file %s(%i)",
         filename,outFile);
  
  //Define a dimension and variables for each type
  int count=end-start+1;
  int receiversDim;
  assert(nc_def_dim(outFile,"numReceivers",count,&receiversDim)==NC_NOERR,
         "writeIndexSplitLocations--unable to create dimension numReceivers in %s",
         filename);
  
  int typeVar,ampVar;
  assert(nc_def_var(outFile,"receiverType",NC_FLOAT,1,
                    &receiversDim,&typeVar)==NC_NOERR,
         "writeIndexSplitLocations--unable to define variable receiverType in file %s(%i,%i)",
         filename,outFile,receiversDim);
  assert(nc_def_var(outFile,"receiverAmp",NC_FLOAT,1,
                    &receiversDim,&ampVar)==NC_NOERR,
         "writeIndexSplitLocations--unable to define variable receiverAmp in file %s(%i,%i)",
         filename,outFile,receiversDim);
  
  int xVar,yVar,zVar;
  int bxVar,byVar,bzVar;
  assert(nc_def_var(outFile,"receiverX",NC_FLOAT,1,
                    &receiversDim,&xVar)==NC_NOERR,
         "writeIndexSplitLocations--unable to create variable receiverX in %s(%i,%i)",
         filename,outFile,receiversDim);
  assert(nc_def_var(outFile,"receiverY",NC_FLOAT,1,
                    &receiversDim,&yVar)==NC_NOERR,
         "writeIndexSplitLocations--unable to create variable receiverY in %s(%i,%i)",
         filename);
  assert(nc_def_var(outFile,"receiverZ",NC_FLOAT,1,
                    &receiversDim,&zVar)==NC_NOERR,
         "writeIndexSplitLocations--unable to create variable receiverZ in %s(%i,%i)",
         filename,outFile,receiversDim);
  
  assert(nc_def_var(outFile,"receiverBx",NC_FLOAT,1,
                    &receiversDim,&bxVar)==NC_NOERR,
         "writeIndexSplitLocations--unable to create variable receiverBx in %s(%i,%i)",
         filename,outFile,receiversDim);
  assert(nc_def_var(outFile,"receiverBy",NC_FLOAT,1,
                    &receiversDim,&byVar)==NC_NOERR,
         "writeIndexSplitLocations--unable to create variable receiverBy in %s(%i,%i)",
         filename,outFile,receiversDim);
  assert(nc_def_var(outFile,"receiverBz",NC_FLOAT,1,
                    &receiversDim,&bzVar)==NC_NOERR,
         "writeIndexSplitLocations--unable to create variable receiverBz in %s(%i,%i)",
         filename,outFile,receiversDim);
  
  int integrateVar;
  assert(nc_def_var(outFile,"receiverIntegrate",NC_INT,1,
                    &receiversDim,&integrateVar)==NC_NOERR,
         "writeIndexSplitLocations--unable to create variable receiverIntegrate in %s(%i)",
         filename,outFile);
  
  if(defData){
    int receiverDims[2]={receiversDim,nSamplesDim};
    int dataVar;
    assert(nc_def_var(outFile,"receiverData",NC_FLOAT,2,
                      receiverDims,&dataVar)==NC_NOERR,
           "writeIndexSplitLocations--unable to create variable receiverData in %s",
           filename);
  }
  
  //Fill the receiver variables.
  nc_enddef(outFile);
  
  //The starting and ending indicies of receivers in this file.
  assert(nc_put_var_int(outFile,startIndexVar,&start)==NC_NOERR,
         "writeIndexSplitLocations--unable to write variable startIndex in %s(%i,%i)",
         filename,outFile,startIndexVar);
  assert(nc_put_var_int(outFile,endIndexVar,&end)==NC_NOERR,
         "writeIndexSplitLocations--unable to write variable endIndex in %s(%i,%i)",
         filename,outFile,endIndexVar);
  
  //Use a temp to fill in the locations and flags, this makes for a much faster
  // write.
  float *tempflt;
  assert((tempflt=(float*)malloc(count*sizeof(float)))!=NULL,
         "writeIndexSplitLocations--unable to allocate %i floats",
         count);
  
  //X locations
  for(int i=start;i<=end;i++){
    tempflt[i-start]=(*_receivers)[i]->_x;
  }
  assert(nc_put_var_float(outFile,xVar,tempflt)==NC_NOERR,
         "writeIndexSplitLocations--unable to write variable xr in %s(%i,%i)",
         filename,outFile,xVar);
  
  //Y locations
  for(int i=start;i<=end;i++){
    tempflt[i-start]=(*_receivers)[i]->_y;
  }
  assert(nc_put_var_float(outFile,yVar,tempflt)==NC_NOERR,
         "writeIndexSplitLocations--unable to write variable yr in %s(%i,%i)",
         filename,outFile,yVar);
  
  //Z locations
  for(int i=start;i<=end;i++){
    tempflt[i-start]=(*_receivers)[i]->_z;
  }
  assert(nc_put_var_float(outFile,zVar,tempflt)==NC_NOERR,
         "writeIndexSplitLocations--unable to write variable zr in %s(%i,%i)",
         filename,outFile,zVar);
  
  //X orientations
  for(int i=start;i<=end;i++){
    tempflt[i-start]=(*_receivers)[i]->_bx;
  }
  assert(nc_put_var_float(outFile,bxVar,tempflt)==NC_NOERR,
         "writeIndexSplitLocations--unable to write variable bx in %s(%i,%i)",
         filename,outFile,xVar);
  
  //Y orientations
  for(int i=start;i<=end;i++){
    tempflt[i-start]=(*_receivers)[i]->_by;
  }
  assert(nc_put_var_float(outFile,byVar,tempflt)==NC_NOERR,
         "writeIndexSplitLocations--unable to write variable by in %s(%i,%i)",
         filename,outFile,yVar);
  
  //Z orientations
  for(int i=start;i<=end;i++){
    tempflt[i-start]=(*_receivers)[i]->_bz;
  }
  assert(nc_put_var_float(outFile,bzVar,tempflt)==NC_NOERR,
         "writeIndexSplitLocations--unable to write variable bz in %s(%i,%i)",
         filename,outFile,zVar);
  
  //Amplitude factor
  for(int i=start;i<=end;i++){
    tempflt[i-start]=(*_receivers)[i]->_amp;
  }
  assert(nc_put_var_float(outFile,ampVar,tempflt)==NC_NOERR,
         "writeIndexSplitLocations--unable to write variable amp in %s(%i,%i)",
         filename,outFile,ampVar);
  
  //And the integer variables.
  free(tempflt);
  int *tempint;
  assert((tempint=(int*)malloc(count*sizeof(int)))!=NULL,
         "writeIndexSplitLocations--unable to allocate %i ints",
         count);
  
  //Type
  for(int i=start;i<=end;i++){
    tempint[i-start]=(*_receivers)[i]->type();
  }
  assert(nc_put_var_int(outFile,typeVar,tempint)==NC_NOERR,
         "writeIndexSplitLocations--unable to write variable type in %s(%i,%i)",
         filename,outFile,typeVar);
  
  //Integrate
  for(int i=start;i<=end;i++){
    tempint[i-start]=(*_receivers)[i]->integrate();
  }
  assert(nc_put_var_int(outFile,integrateVar,tempint)==NC_NOERR,
         "writeIndexSplitLocations--unable to write variable integrate in %s(%i,%i)",
         filename,outFile,integrateVar);
  
  free(tempint);
  
  if(iOpenedFile)
    assert(nc_close(outFile)==NC_NOERR,
           "writeIndexSplitLocations--unable to close file %s(%i)",
           filename,outFile);
  
  return count;
}
int masterReceiverNetwork::writeIndexSplitTraces(modelDefStruct* modelDef,int iteration,int nSamps,
                                                 int useLocalData){
  char buffer[1024],tag[10];
  int nFiles=(int)ceil((float)_receivers->size()/(float)_splitFiles);
  tEprintf(Verbose,
           " Writing trace output to file for %i receivers to iteration %i\n",
           _receivers->size(),iteration);
  
  for(int i=0;i<nFiles;i++){
    sprintf(tag,"Seg%02i",i+1);
    generateSplitFilename(buffer,tag);
    
    //Determine the start and end receiver index for this segment.
    int start=i*_splitFiles,end=(i+1)*_splitFiles-1;
    end=MIN(end,_receivers->size()-1);
    
    writeIndexSplitTraces(buffer,start,end,
                          modelDef,iteration,nSamps,useLocalData);
  }
  return _receivers->size();
}
int masterReceiverNetwork::writeIndexSplitTraces(char* filename,
                                                 int rStartI,int rEndI,
                                                 modelDefStruct* modelDef,
                                                 int startI,int stopI,
                                                 int useLocalData){
  int count=rEndI-rStartI+1;
  tEprintf(Verbose,
           "  %i receivers in file: %s \n",
           count,filename);
  
  //open the file and get the data variable
  int dataVar,outFile=openCDFFile(filename,FALSE,NC_WRITE);
  assert(nc_inq_varid(outFile,"receiverData",&dataVar)==NC_NOERR,
         "writeIndexSplitTraces--unable to open var data in %s(%i)",
         filename,outFile);
  
  //Now fill in the data variables
  //Calculate the sub-sampled trace indicies.
  int startSample=(int)floor((float)(startI)/_timeSubsample);
  int stopSample=MIN(nSamples(modelDef)-1,
                     (int)floor((float)(stopI)/_timeSubsample));
  int nSamps=stopSample-startSample+1;
  
  for(size_t i=rStartI;i<=rEndI;i++){
    receiver* currR=(*_receivers)[i];
    
    //where to write this data
    size_t start[2]={i-rStartI,startSample};
    size_t count[2]={1,nSamps};
    
    if(useLocalData){
      currR->scaleData(modelDef,startI,stopI,currR->_data,work(modelDef));
    }else{
      getReceiverData(i,startI,stopI,work(modelDef));
    }
    assert(nc_put_vara_float(outFile,dataVar,
                             start,count,work(modelDef))==NC_NOERR,
           "writeIndexSplitTraces--unable to write data[%i,%,%i,%i] to %s(%i,%i)",
           start[0],start[1],count[0],count[1],
           filename,outFile,dataVar);
  }
  
  //close the file and exit
  assert(nc_close(outFile)!=-1,
         "writeIndexSplitTraces--nc_close of %s unsuccessful",
         filename);
  return count*nSamps;
}

///Here is the method for writing out receiver grids.
int masterReceiverNetwork::writeSingleGrid(modelDefStruct* modelDef){
  DEF_MODEL_SIZE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  
  float* data=NULL;
  int ndata=0;
  for(int ii=0;ii<_receiverGrids->size();ii++){
    receiverGrid* currGrid=(*_receiverGrids)[ii];
    tEprintf(Verbose,"  Writing receiver-grid %i\n",ii);
    
    if(currGrid->isRaw()==1){
      initSend();
      sendMessage(AllProcesses,MESSAGE_SEND_RECEIVER_GRID,"ii",ii,ii);
      
      for(int jj=0;jj<NumProcs;jj++){
        int checkIndex;
        getMessage(Tids[jj],MESSAGE_SEND_RECEIVER_GRID,"i",&checkIndex);
        assert(checkIndex==ii,
               "writeSingleGrid--index mismatch for raw output %i!=%i; Proc %i",
               checkIndex,ii,Tids[jj]);
      }
      tEprintf(Verbose,"    Slaves ouput raw grid directly.\n");
      continue;
    }
    
    //Create a modified version of the trace filename.
    char buffer[1024];
    strcpy(buffer,_traceOutputName);
    char* tag=strstr(buffer,".cdf");
    if(tag){
      *tag='\0';
    }else{
      tag=buffer+strlen(buffer);
    }
    
    if(currGrid->isRaw()==2){
      sprintf(tag,"%s%iXCorr.cdf",currGrid->typeName(),ii);
    }else{
      sprintf(tag,"%s%i.cdf",currGrid->typeName(),ii);
    }
    
    //Create the file with a standard header.
    writeCDFHeader(buffer,
                   NX,NY,NZ,NT,NULL,
                   dx,minX,dy,minY,dz,minZ,dt,minT,
                   0,NULL);
    
    //Fill in the receiver locations.
    int outFile=openCDFFile(buffer,FALSE,NC_WRITE);
    
    //get required dimension
    int ntDim;
    assert(nc_inq_dimid(outFile,"NT",&ntDim)==NC_NOERR,
           "writeSingleLocations--unable to open dimension NT from file %s",
           _traceOutputName);
    
    //Define a dimension for the decimation factor and another for
    // the actual number of samples in each trace.
    nc_redef(outFile);
    int receiverSkipDim,rtDim;
    assert(nc_def_dim(outFile,"receiverDecimate",_timeSubsample,
                      &receiverSkipDim)==NC_NOERR,
           "writeSingleLocations--unable to create dimension receiverDecimate in %s",
           _traceOutputName);
    assert(nc_def_dim(outFile,"RT",currGrid->_nt,
                      &rtDim)==NC_NOERR,
           "writeSingleLocations--unable to create dimension nSamples in %s",
           _traceOutputName);
    
    //Define dimensions for the number of receivers in each direction.
    int rxDim,ryDim,rzDim;
    assert(nc_def_dim(outFile,"RX",(*_receiverGrids)[ii]->_nx,&rxDim)==NC_NOERR,
           "writeSingleGrid--unable to create dimension RXDim in %s(%i)",
           buffer,outFile);
    assert(nc_def_dim(outFile,"RY",(*_receiverGrids)[ii]->_ny,&ryDim)==NC_NOERR,
           "writeSingleGrid--unable to create dimension RYDim in %s(%i)",
           buffer,outFile);
    assert(nc_def_dim(outFile,"RZ",(*_receiverGrids)[ii]->_nz,&rzDim)==NC_NOERR,
           "writeSingleGrid--unable to create dimension RZDim in %s(%i)",
           buffer,outFile);
    
    //Define variables.
    int rxVar,ryVar,rzVar,rdVar;
    assert(nc_def_var(outFile,"rx",NC_FLOAT,1,&rxDim,&rxVar)==NC_NOERR,
           "writeSingleGrid--unable to create variable rx in %s(%i,%i)",
           buffer,outFile,rxDim);
    assert(nc_def_var(outFile,"ry",NC_FLOAT,1,&ryDim,&ryVar)==NC_NOERR,
           "writeSingleGrid--unable to create variable ry in %s(%i,%i)",
           buffer,outFile,ryDim);
    assert(nc_def_var(outFile,"rz",NC_FLOAT,1,&rzDim,&rzVar)==NC_NOERR,
           "writeSingleGrid--unable to create variable rz in %s(%i,%i)",
           buffer,outFile,rzDim);
    
    
    int rDims[4]={rzDim,ryDim,rxDim,rtDim};
    assert(nc_def_var(outFile,
                      currGrid->isRaw()==2?"receiverXCorr":"receiverData",
                      NC_FLOAT,currGrid->isRaw()==2?3:4,
                      rDims,&rdVar)==NC_NOERR,
           "writeSingleGrid--unable to create variable receiverData in %s(%i,%i,%i,%i,%i)",
           buffer,outFile,rzDim,ryDim,rxDim,rtDim);
    
    //Fill the receiver variables.
    assert(nc_enddef(outFile)==NC_NOERR,
           "writeSingleGrid--unable to take %s(%i) out of define mode\n",
           buffer,outFile);
    
    float *pos=(float*)malloc(MAX3((*_receiverGrids)[ii]->_nx,
                                   (*_receiverGrids)[ii]->_ny,
                                   (*_receiverGrids)[ii]->_nz)*sizeof(float));
    for(int i=0;i<(*_receiverGrids)[ii]->_nx;i++)
      pos[i]=(*_receiverGrids)[ii]->_x+dx*i*(*_receiverGrids)[ii]->_dx;
    assert(nc_put_var_float(outFile,rxVar,pos)==NC_NOERR,
           "writeSingleGrid--unable to write variable rx in %s(%i,%i)",
           buffer,outFile,rxVar);
    
    for(int i=0;i<(*_receiverGrids)[ii]->_ny;i++)
      pos[i]=(*_receiverGrids)[ii]->_y+dy*i*(*_receiverGrids)[ii]->_dy;
    assert(nc_put_var_float(outFile,ryVar,pos)==NC_NOERR,
           "writeSingleGrid--unable to write variable ry in %s(%i,%i)",
           buffer,outFile,ryVar);
    
    for(int i=0;i<(*_receiverGrids)[ii]->_nz;i++)
      pos[i]=(*_receiverGrids)[ii]->_z+dz*i*(*_receiverGrids)[ii]->_dz;
    assert(nc_put_var_float(outFile,rzVar,pos)==NC_NOERR,
           "writeSingleGrid--unable to write variable rz in %s(%i,%i)",
           buffer,outFile,rzVar);
    free(pos);
    
    //And write the receiver data.
    size_t start[4]={0,0,0,0},count[4]={0,0,0,0};
    for(int jj=0;jj<NumProcs;jj++){
      getReceiverGridData(ii,Tids[jj],currGrid->isRaw()==2,
                          ndata,data,start,count);
      if(currGrid->isRaw()==2){
        if(!(count[0]*count[1]*count[2])){
          tEprintf(Verbose>1," Proc %i: No receivers\n",jj);
        }else{
          tEprintf(Verbose>1," Proc %i: writing [%i:%i,%i:%i,%i:%i]\n",
                   jj,
                   start[0],start[0]+count[0]-1,
                   start[1],start[1]+count[1]-1,
                   start[2],start[2]+count[2]-1);
          assert(nc_put_vara_float(outFile,rdVar,
                                   start,count,data)==NC_NOERR,
                 "writeSingleGrid--unable to write grid%i[%i:%i,%i:%i,%i:%i] to %s(%i,%i); proc %i",
                 ii,
                 start[0],count[0],start[1],count[1],start[2],count[2],
                 buffer,outFile,rdVar,jj);
        }
      }else{
        assert(count[3]==currGrid->_nt,
               "writeSingleGrid--nt mismatch %i in master, %i from slave.",
               currGrid->_nt,count[3]);
        if(!(count[0]*count[1]*count[2]*count[3])){
          tEprintf(Verbose>1," Proc %i: No receivers\n",jj);
        }else{
          tEprintf(Verbose>1," Proc %i: writing [%i:%i,%i:%i,%i:%i,%i:%i]\n",
                   jj,
                   start[0],start[0]+count[0]-1,
                   start[1],start[1]+count[1]-1,
                   start[2],start[2]+count[2]-1,
                   start[3],start[3]+count[3]-1);
          assert(nc_put_vara_float(outFile,rdVar,
                                   start,count,data)==NC_NOERR,
                 "writeSingleGrid--unable to write grid%i[%i:%i,%i:%i,%i:%i,%i:%i] to %s(%i,%i); proc %i",
                 ii,
                 start[0],count[0],start[1],count[1],start[2],count[2],start[3],count[3],
                 buffer,outFile,rdVar,jj);
        }
      }
    }
    
    //Close the file.
    assert(nc_close(outFile)==NC_NOERR,
           "writeSingleGrid--unable to close file %s(%i)",
           buffer,outFile);
  }
  if(data)free(data);
  return _receivers->size();
}

//Convience subroutines for adding receivers
int masterReceiverNetwork::add3CReceiver(modelDefStruct* modelDef,
                                         float xr,float yr,float zr,float ramp,int integrate,
                                         int print,int useCubicInterp,int surfaceReceiver){
  velocityReceiver *currZReceiver,*currXReceiver,*currYReceiver;
  _receivers->Add(currZReceiver=new velocityReceiver(modelDef,FALSE,
                                                     xr,yr,zr,1.0,0.0,0.0,ramp,useCubicInterp,surfaceReceiver));
  _receivers->Add(currXReceiver=new velocityReceiver(modelDef,FALSE,
                                                     xr,yr,zr,0.0,1.0,0.0,ramp,useCubicInterp,surfaceReceiver));
  _receivers->Add(currYReceiver=new velocityReceiver(modelDef,FALSE,
                                                     xr,yr,zr,0.0,0.0,1.0,ramp,useCubicInterp,surfaceReceiver));
  
  const char* compName[]={"acceleration","velocity","displacement"};
  if(integrate>0){
    currZReceiver->setDisplacement();
    currXReceiver->setDisplacement();
    currYReceiver->setDisplacement();
  }else if(integrate<0){
    currZReceiver->setAcceleration();
    currXReceiver->setAcceleration();
    currYReceiver->setAcceleration();
  }
  tEprintf(print && Verbose,
           "Added 3-comp %s receiver %i-%i (%.1f,%.1f,%.1f) (%.1f)\n",
           compName[integrate+1],
           _receivers->size()-2,_receivers->size(),
           xr,yr,zr,ramp);
  
  return _receivers->size();
}
int masterReceiverNetwork::add4CReceiver(modelDefStruct* modelDef,
                                         float xr,float yr,float zr,float ramp,int integrate,
                                         int print,int useCubicInterp,int surfaceReceiver){
  velocityReceiver *currZReceiver,*currXReceiver,*currYReceiver;
  _receivers->Add(currZReceiver=new velocityReceiver(modelDef,FALSE,
                                                     xr,yr,zr,1.0,0.0,0.0,ramp,useCubicInterp,surfaceReceiver));
  _receivers->Add(currXReceiver=new velocityReceiver(modelDef,FALSE,
                                                     xr,yr,zr,0.0,1.0,0.0,ramp,useCubicInterp,surfaceReceiver));
  _receivers->Add(currYReceiver=new velocityReceiver(modelDef,FALSE,
                                                     xr,yr,zr,0.0,0.0,1.0,ramp,useCubicInterp,surfaceReceiver));
  
  char compName[20];
  if(!integrate){
    strcpy(compName,"velocity");
  }else if(integrate>0){
    strcpy(compName,"displacement");
    currZReceiver->setDisplacement();
    currXReceiver->setDisplacement();
    currYReceiver->setDisplacement();
  }else if(integrate<0){
    strcpy(compName,"acceleration");
    currZReceiver->setAcceleration();
    currXReceiver->setAcceleration();
    currYReceiver->setAcceleration();
  }
  
  _receivers->Add(new pressureReceiver(modelDef,FALSE,
                                       xr,yr,zr,ramp,useCubicInterp));
  tEprintf(print && Verbose,
           "Added 4-comp receiver %s-pressure %i-%i (%.1f,%.1f,%.1f) (%.1f)\n",
           compName,
           _receivers->size()-3,_receivers->size(),
           xr,yr,zr,ramp);
  
  return _receivers->size();
}

slaveReceiverNetwork::slaveReceiverNetwork(modelDefStruct* modelDef,slaveReceiverNetwork* theSource){
  _timeSubsample=theSource->_timeSubsample;
  for(int i=0;i<theSource->size();i++){
    receiver* curr=(*theSource->_receivers)[i],*currCopy;
    switch(curr->type()){
      case VELOCITY_RECEIVER:
        currCopy=new velocityReceiver(modelDef,TRUE,
                                      curr->_x,curr->_y,curr->_z,
                                      curr->_bx,curr->_by,curr->_bz,curr->_amp);
        break;
      case PRESSURE_RECEIVER:
        currCopy=new pressureReceiver(modelDef,TRUE,
                                      curr->_x,curr->_y,curr->_z,curr->_amp);
        break;
        
      default:
        assert(FALSE,"slaveReceiverNetwork::operator(=)--unknow type for theSource receiver %i, %i",
               i,curr->type());
    }
    _receivers->Add(currCopy);
  }
}

///Method for checkpoint writing of receiver data.
FILE* slaveReceiverNetwork::doCheckpoint(int iteration,char* cpDir,int cpID,
                                         int callDepth,
                                         FILE* cpFile){
  //Since this is the base the file has not been open yet,
  // open the binary file that will hold all the info need for
  // the restart.
  char filename[1024];
  sprintf(filename,"%s/checkpoint_receivers_%i.bin",cpDir,cpID);
  cpFile=fopen(filename,"wb");
  assert(cpFile!=NULL,
         "slaveReceiverNetwork::doCheckpoint--unable to open file\n\t%s\n",
         filename);
  
  //Save the current data for each active receiver.
  for(int i=0;i<size();i++){
    if((*_receivers)[i]->active()){
      (*_receivers)[i]->doCheckpoint(iteration,cpFile);
    }
  }
  
  if(!callDepth){
    fclose(cpFile);
    cpFile=NULL;
  }
  return cpFile;
}

///Method for checkpoint reading of receiver data.
FILE* slaveReceiverNetwork::readCheckpoint(int iteration,char* cpDir,int cpID,
                                           int callDepth,
                                           FILE* cpFile){
  //Since this is the base the file has not been open yet,
  // open the binary file that will hold all the info need for
  // the restart.
  char filename[1024];
  sprintf(filename,"%s/checkpoint_receivers_%i.bin",cpDir,cpID);
  cpFile=fopen(filename,"rb");
  assert(cpFile!=NULL,
         "slaveReceiverNetwork::doCheckpoint--unable to open file\n\t%s\n",
         filename);
  
  //Read the saved data for each active receiver.
  for(int i=0;i<size();i++){
    if((*_receivers)[i]->active()){
      (*_receivers)[i]->readCheckpoint(iteration,cpFile);
    }
  }
  
  if(!callDepth){
    fclose(cpFile);
    cpFile=NULL;
  }
  return cpFile;
}

///
///Method to initialize a network using information that is
/// sent from the master process.
int slaveReceiverNetwork::initializeNetwork(modelDefStruct* modelDef,int allocateData,
                                            int useCubicInterp){
  int numReceivers;
  getMessage(Parent,MESSAGE_SEND_RECEIVER,"ii",
             &numReceivers,&_timeSubsample);
  
  int *activeReceivers;
  if(!numReceivers){
    _receivers=NULL;
    activeReceivers=NULL;
  }else{
    _receivers=new receiverArray(0,numReceivers);
    assert((activeReceivers=(int*)malloc(numReceivers*sizeof(int)))!=NULL,
           "slaveReceiverNetwork::initializeNetwork--uable to allocate %i ints for activeReceivers",
           numReceivers);
    
    for(int i=0,currRead=0;i<numReceivers;i++){
      int type;
      unpackMessage("i",&type);
      switch(type){
        case VELOCITY_RECEIVER:
          _receivers->Add(new velocityReceiver(modelDef,allocateData));
          break;
        case PRESSURE_RECEIVER:
          _receivers->Add(new pressureReceiver(modelDef,allocateData));
          break;
          
        default:
          assert(FALSE,"receiverNetwork(modelDefStruct* modelDef,int theSource)--unknown receiver type %i; %i",
                 type,i);
      }
      if(++currRead>=MAX_R_PER_SEND){
        currRead=0;
        getMessage(Parent,MESSAGE_SEND_RECEIVER,NULL);
      }
      activeReceivers[i]=(*_receivers)[i]->active();
    }
  }
  
  //Send an acknowlegement that all receivers have been received.
  getMessage(Parent,MESSAGE_SEND_RECEIVER,NULL);
  if(!numReceivers){
    sendMessage(Parent,MESSAGE_SEND_RECEIVER,NULL);
  }else{
    for(int i=0;i<numReceivers;i+=MAX_R_PER_SEND){
      initSend();
      sendMessage(Parent,MESSAGE_SEND_RECEIVER,"I",
                  &(activeReceivers[i]),MIN(MAX_R_PER_SEND,numReceivers-i));
    }
    free(activeReceivers);
  }
  
  //Now check for receiverGrids.
  getMessage(Parent,MESSAGE_SEND_RECEIVER_GRID,"i",&numReceivers);
  for(int i=0;i<numReceivers;i++){
    int type;
    getMessage(Parent,MESSAGE_SEND_RECEIVER_GRID+i+1,"i",&type);
    switch(type){
      case PRESSURE_RECEIVER_GRID:
        _receiverGrids->Add(new pressureReceiverGrid(modelDef,i,TRUE));
        break;
      case VX_RECEIVER_GRID:
        _receiverGrids->Add(new vxReceiverGrid(modelDef,i,TRUE));
        break;
      case VY_RECEIVER_GRID:
        _receiverGrids->Add(new vyReceiverGrid(modelDef,i,TRUE));
        break;
      case VZ_RECEIVER_GRID:
        _receiverGrids->Add(new vzReceiverGrid(modelDef,i,TRUE));
        break;
        
      default:
        assert(FALSE,"initializeNetwork--bad grid type %i\n",type);
    }
    
    //Make sure the buffer is big enough, this has been causing problems for
    // small models with lots of receivers.
    DEF_MODEL_SIZE(modelDef);
    setMessageBuffer(4*(MAX_R_PER_SEND+10)*(NT+10));
    
    //Send the final receiver acknoledgment.
    initSend();
    sendMessage(Parent,MESSAGE_SEND_RECEIVER_GRID+i+1,NULL);
  }
  
  return size();
}
