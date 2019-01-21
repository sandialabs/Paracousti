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
 *  sgfdReceiverNetwork.hh
 *
 *
 *  Declares classes that hold all the receivers for the master and slave processes.  These
 *  are basically containers that serve as a single class to handle all receiver functions
 *  as a group.
 *
 *  Declares the following classes:
 *  receiverNetwork
 *  masterReceiverNetwork
 *  slaveReceiverNetwork
 *
 *  Defines the following macros:
 *  THREE_C_RECEIVER
 *  FOUR_C_RECEIVER
 *  VX_RECEIVER
 *  VY_RECEIVER
 *  VZ_RECEIVER
 *  REC_SPLIT_BY_COMP
 *
 */

#ifndef _sgfdReceiverNetwork_hh_
#define _sgfdReceiverNetwork_hh_

#include "sgfdReceivers.hh"
#include "sgfdReceiverGroups.hh"

#define THREE_C_RECEIVER 3
#define FOUR_C_RECEIVER  4

#define VX_RECEIVER 6
#define VY_RECEIVER 7
#define VZ_RECEIVER 8

/*
 // Global Variables
 */
extern int NumReceiverTypeNames;
extern const char* ReceiverTypeNames[]; //{"","Velocity","Pressure","3C","4C",
//                                   "Vx","Vy","Vz"}

//
///Also want a container class to do the work on a group of receivers. This allows
/// the incorporation of receivers with vastly different properties (such as a
/// stress receiver that records 9 data fields).
/// Break this class into a master and slave version. This is done to
/// improve the handleing of a large number of receivers. Make the original
/// receiverNetwork a pure virtual class.
class receiverNetwork{
public:
  receiverArray* _receivers;
  receiverGridArray* _receiverGrids;
  
  char _traceOutputName[1024];
  int _timeSubsample;   //may not want to write every time.
  
  float* _wrk;
  
  ///Need a null constructor for cases where receivers are specified entirely on
  /// the command line
  receiverNetwork();
  
  virtual ~receiverNetwork();
  
  ///Operator to index receivers using array notation, don't know
  /// if this ever used or not but it might be nice to have.
  receiver* operator[](int index){return (*_receivers)[index];}
  
  ///Pure virtual ethod for checkpoint writing of receiver data.
  virtual FILE* doCheckpoint(int iteration,char* cpDir,int cpID,
                             int callDepth,
                             FILE* cpFile=NULL)=0;
  
  ///Pure virtual method for checkpoint reading of receiver data.
  virtual FILE* readCheckpoint(int iteration,char* cpDir,int cpID,
                               int callDepth,
                               FILE* cpFile=NULL)=0;
  
  
  ///Method to combine receivers (some read from file; others
  /// defined on the command line).
  int combine(receiverNetwork* extraReceivers);
  
  //Here are some methods to help with subsampled (decimated) output.
  int timeSubsample(){return _timeSubsample;}
  int subsampleIndex(int i){
    return i*_timeSubsample;
  }
  int nSamples(modelDefStruct* modelDef){
    return 1+(int)floor((float)(modelDef->NT-1)/_timeSubsample);
  }
  
  char* setOutputName(char* name){
    strcpy(_traceOutputName,name);
    return _traceOutputName;
  }
  int size(){return _receivers?_receivers->size():0;}
  int active(int index){return (*_receivers)[index]->active();}
  
  float* work(modelDefStruct* modelDef);
  
  //May want to allocate data for all of the receivers.
  int allocateData(modelDefStruct* modelDef);
  
  //in some cases we want direct access to the individual receiver parameters. The locations
  // might need to be changes (for a coordinate tranformation)
  float& x(int index){return (*_receivers)[index]->_x;}
  float& y(int index){return (*_receivers)[index]->_y;}
  float& z(int index){return (*_receivers)[index]->_z;}
  
  float& bx(int index){return (*_receivers)[index]->_bx;}
  float& by(int index){return (*_receivers)[index]->_by;}
  float& bz(int index){return (*_receivers)[index]->_bz;}
  float& amp(int index){return (*_receivers)[index]->_amp;}
  
  float* data(int index){return (*_receivers)[index]->_data;}
  float indexData(int index,int i){
    return (*_receivers)[index]->_data[i*_timeSubsample];
  }
  
  int type(int index){return (*_receivers)[index]->type();}
  int integrate(int index){return (*_receivers)[index]->integrate();}
  
  //Send the data for these receivers back and forth between master and slave process.
  virtual int initializeNetwork(modelDefStruct* modelDef,int allocateData,
                                int useCubicInterp=TRUE)=0;
  
  int addVelocityReceiver(modelDefStruct* modelDef,
                          float xr,float yr,float zr,
                          float bx,float by,float bz,
                          float ramp,int integrate,
                          int print=TRUE,int useCubicInterp=TRUE,int surfaceReceiver=FALSE);
  
  //Fill the data for this network of receivers. Moved from slaveReceiverNetwork
  // to here for serial implementation.
  virtual int fillReceivers(modelDefStruct* modelDef,int iteration,
                            float* vx,float* vy,float* vz,
                            float* xx,float* yy,float* zz,
                            float* xy,float* xz,float* yz);
};

#define REC_SPLIT_BY_COMP 1
///Here is the class for the receiver network that is used by the
/// master process. This does the IO into a combined file.
class masterReceiverNetwork:public receiverNetwork{
public:
  int* _receiverProcs; //this is an array the same size as _receivers with indicies for
  // which process owns each receiver.
  
  //Some flags to control aspects of the output specific to the master processes.
  int _splitFiles;
  int _nVxR,_nVyR,_nVzR,_nPrR,*_rTypeFlag;
  
  int _useModelReceivers; //-RM flag
  
  /// Need a null constructor for cases where receivers are specified entirely on the command line
  masterReceiverNetwork():receiverNetwork(){
    _receiverProcs=NULL;
    _useModelReceivers=TRUE;
    
    _splitFiles=FALSE;
    _timeSubsample=1;
    _nVxR=_nVyR=_nVzR=_nPrR=0;
    _rTypeFlag=NULL;
  }
  
  ///A new receiverNetwork can be read directly from a cdf file. This subsumes the previous function
  /// readCDFReceivers.
  masterReceiverNetwork(modelDefStruct* modelDef,const char* fileName,int freeSurface,
                        int useModelReceivers=TRUE,int timeSubsample=1);
  
  virtual ~masterReceiverNetwork(){
    if(_receiverProcs)
      free(_receiverProcs);
    if(_rTypeFlag)
      free(_rTypeFlag);
  }
  
  ///Method for checkpoint writing of receiver data.
  virtual FILE* doCheckpoint(int iteration,char* cpDir,int cpID,
                             int callDepth,
                             FILE* cpFile=NULL){
    //As far as I can figure out at this time the master does
    // not need to do anything. The slaves store the data and
    // will need to do all the checkpointing.
    
    return NULL;
  }
  
  ///Method for checkpoint reading of receiver data.
  virtual FILE* readCheckpoint(int iteration,char* cpDir,int cpID,
                               int callDepth,
                               FILE* cpFile=NULL){
    //As far as I can figure out at this time the master does
    // not need to do anything. The slaves store the data and
    // will need to do all the checkpointing.
    
    return NULL;
  }
  
  
  ///Method to make sure that all receivers in the domain are
  /// active. Failure of this check is a fatal error.
  int checkActive(int fail);
  
  ///Return the slave process id that contains a given receiver.
  int proc(int index){return _receiverProcs[index];}
  
  //Can add receivers to an existing receiverNetwork by processing command line arguments. This replaces
  // addReceivers. The argument argOffset is the index of the first character after the one that
  // tells the calling process to call this function.
  virtual int addReceivers(int& i,int argOffset,
                           int argc,char* argv[],
                           modelDefStruct* modelDef=NULL);
  
  virtual int addReceivers(masterReceiverNetwork* extraReceivers);
  
  //Send the data for these receivers back and forth between master and slave process.
  virtual int initializeNetwork(modelDefStruct* modelDef,int allocateData,
                                int target=AllProcesses);
  
  
  //Also want to be able to generate netcdf output to a trace file. Write the header and
  // the traces seperatly
  int writeHeader(modelDefStruct* modelDef,int defData){
    if(_splitFiles==1){
      return writeCompSplitHeader(modelDef,defData);
    }else if(_splitFiles){
      return writeIndexSplitHeader(modelDef,defData);
    }
    return writeSingleHeader(modelDef,defData);
  }
  int writeLocations(modelDefStruct* modelDef,
                     int defData,int outFile=-1){
    if(_splitFiles)
      assert(FALSE,
             "writeLocations--can not be called directly with split files");
    return writeSingleLocations(modelDef,defData,outFile);
  }
  int rewriteLocations(){
    if(_splitFiles)
      assert(FALSE,
             "rewriteLocations--can not be called directly with split files");
    return rewriteSingleLocations();
  }
  virtual int writeTraces(modelDefStruct* modelDef,int startI,int stopI,
                          int useLocalData=FALSE);
  
  //Get the data for a specific receiver.
  int getReceiverData(int index,int startI,int stopI,float* data);
  
  ///Get data for a specific receiver-grid from a specific process.
  int getReceiverGridData(int index,int tid,int isXcorr,
                          int &n,floatPtr& data,
                          size_t start[4],size_t count[4]);
  
protected:
  int receiverType(char* argv){
    if(isdigit(argv[0]))
      return atoi(argv);
    
    for(int i=1;i<=NumReceiverTypeNames;i++){
      if(!strcmp(argv,ReceiverTypeNames[i]))
        return i;
    }
    
    assert(FALSE,"rtype--unknown receiver type %s",argv);
    return FALSE; //never executed
  }
  
  int writeSingleHeader(modelDefStruct* modelDef,int defData);
  
  int writeSingleLocations(modelDefStruct* modelDef,
                           int defData,int outFile);
  int rewriteSingleLocations();
  
  int writeSingleTraces(modelDefStruct* modelDef,int startI,int stopI,
                        int useLocalData);
  
  char* generateSplitFilename(char* target,const char* tag);
  
  //Methods to split files by component.
  int writeCompSplitHeader(modelDefStruct* modelDef,int defData);
  int writeCompSplitLocations(char* filename,int count,
                              int tag,int *tagArray,
                              modelDefStruct* modelDef,
                              int defData,int outFile);
  int writeCompSplitTraces(modelDefStruct* modelDef,int iteration,int nSamps,
                           int useLocalData);
  int writeCompSplitTraces(char* filename,int count,
                           int tag,int *tagArray,
                           modelDefStruct* modelDef,
                           int startI,int stopI,
                           int useLocalData);
  
  //Methods to split files by index.
  int writeIndexSplitHeader(modelDefStruct* modelDef,int defData);
  int writeIndexSplitLocations(char* filename,
                               int start,int end,
                               modelDefStruct* modelDef,
                               int defData,int outFile);
  int writeIndexSplitTraces(modelDefStruct* modelDef,int iteration,int nSamps,
                            int useLocalData);
  int writeIndexSplitTraces(char* filename,
                            int rStartI,int rEndI,
                            modelDefStruct* modelDef,
                            int startI,int stopI,
                            int useLocalData);
  
  ///Here is the method for writing out receiver grids.
  int writeSingleGrid(modelDefStruct* modelDef);
  
  //Convience subroutines for adding receivers
  int add3CReceiver(modelDefStruct* modelDef,
                    float xr,float yr,float zr,float ramp,int integrate,
                    int print=TRUE,int useCubicInterp=TRUE,int surfaceReceiver=FALSE);
  int add4CReceiver(modelDefStruct* modelDef,
                    float xr,float yr,float zr,float ramp,int integrate,
                    int print=TRUE,int useCubicInterp=TRUE,int surfaceReceiver=FALSE);
};
typedef masterReceiverNetwork* masterReceiverNetworkPtr;

///Here is the class for the slave receiver networks. These do
/// the actual work of offloading and interpolating the data.
class slaveReceiverNetwork:public receiverNetwork{
public:
  slaveReceiverNetwork(modelDefStruct* modelDef,
                       int allocateData=TRUE,int mesgFromParent=TRUE,
                       int useCubicInterp=TRUE){
    if(mesgFromParent){
      initializeNetwork(modelDef,allocateData,useCubicInterp);
    }else{
      _receivers=NULL;
      _timeSubsample=1;
    }
  }
  
  slaveReceiverNetwork(modelDefStruct* modelDef,slaveReceiverNetwork* theSource);
  
  ///Method for checkpoint writing of receiver data.
  virtual FILE* doCheckpoint(int iteration,char* cpDir,int cpID,
                             int callDepth,
                             FILE* cpFile=NULL);
  
  ///Method for checkpoint reading of receiver data.
  virtual FILE* readCheckpoint(int iteration,char* cpDir,int cpID,
                               int callDepth,
                               FILE* cpFile=NULL);
  
  ///
  ///Method to initialize a network using information that is
  /// sent from the master process.
  virtual int initializeNetwork(modelDefStruct* modelDef,int allocateData,
                                int useCubicInterp=TRUE);
};

#endif
