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
 *  paracoustic.cc
 *
 *
 *  Main entry point.
 *
 *  Defines the following function:
 *  main
 *
 */

#include <time.h>

#include "nstdutil.hh"

#include "message_passing.h"

#include "sgfd.hh"
#include "selector.hh"

#include "acoustic_control.hh"
#include "fixed_acoustic.hh"

#include "parallel_acousti.hh"

#include <xmmintrin.h>

#if USE_VAMPIR
int _VTUpdateClass_,_VTBCClass_;
int _VTVelUpdate_,_VTStressUpdate_;
int _VTBCSaveVel_,_VTBCDoVel_,_VTBCDoStress_;
int _VTBCTransStress_,_VTBCFreeStress_;

#endif

//
//Globals
int p_xargc;
char **p_xargv;

// standard
time_t StartTime;
char* CommandLine=nullptr;
int Verbose=TRUE;
char* CDFSourceName=nullptr;

#define COPYWRITE "* Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC (NTESS).\n \
* Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains\n \
* certain rights in this software.\n \
* \n\
* NOTICE:\n \
* For five (5) years from  the United States Government is granted for itself and others\n \
*  acting on its behalf a paid-up, nonexclusive, irrevocable worldwide license in this\n \
*  data to reproduce, prepare derivative works, and perform publicly and display\n \
*  publicly, by or on behalf of the Government. There is provision for the possible\n \
*  extension of the term of this license. Subsequent to that period or any extension\n \
*  granted, the United States Government is granted for itself and others acting on its\n \
*  behalf a paid-up, nonexclusive, irrevocable worldwide license in this data to reproduce,\n \
*  prepare derivative works, distribute copies to the public, perform publicly and display\n \
*  publicly, and to permit others to do so. The specific term of the license can be\n \
*  identified by inquiry made to National Technology and Engineering Solutions of Sandia,\n \
*  LLC or DOE.\n \
* NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT OF ENERGY, NOR\n \
*  NATIONAL TECHNOLOGY AND ENGINEERING SOLUTIONS OF SANDIA, LLC, NOR ANY OF THEIR EMPLOYEES,\n \
*  MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL RESPONSIBILITY FOR THE\n \
*  ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS\n \
*  DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.\n \
* Any licensee of this software has the obligation and responsibility to abide by the\n \
*  applicable export control laws, regulations, and general prohibitions relating to the\n \
*  export of technical data. Failure to obtain an export control license or other authority\n \
*  from the Government may result in criminal liability under U.S. laws.\n"

//
//MAIN PROGRAM
int main(int argc,char* argv[]){
  int oldMXCSR = _mm_getcsr(); //read the old MXCSR setting:  These turn denormals to zero.  Denormals were causing a big slowdown near startup
  int newMXCSR = oldMXCSR | 0x8040; // set DAZ and FZ bits
  _mm_setcsr( newMXCSR ); //write the new MXCSR setting to the MXCSR
  
#if NO_MMAP
#warning "Setting M_MMAP_MAX to 0; Better performance on Cplant"
  mallopt(M_MMAP_MAX,0); //Required to get decent performance on Cplant.
#endif
  
  //Always start a command parser, what to do this here so the -help option will
  // work outside of the parallel enviornment.
  parallelAcoustiCommandParser* theCommandParser=
  new parallelAcoustiCommandParser(argc,argv);
  
  StartTime=time(0);
  
  //Initialize Message Passing
  int tid=registerProcess(&Parent,&argc,&argv);
  if(tid>0){
    //This is a slave process.  Just loop taking action on messages.
    delete theCommandParser;
    acousticSlaveProcess* slaveLoop=new acousticSlaveProcess(tid,Parent);
    while(slaveLoop->processMessage());
    delete slaveLoop;
    return 1;
  }
  
  //Print out the copywrite assertion
  fprintf(stderr,"%s",COPYWRITE);
  
  //This is the master process; read the arguments. Note that the
  //modelDef is now initalized immediatly as the model is specified.
  theCommandParser->processArgs(argc,argv);
  if(! *(theCommandParser->_modelName)){
    assert(FALSE,
           "Single argument for the model not present in call.");
  }
  
  //Define some basic quantities
  DEF_MODEL_SIZE(theCommandParser->modelDef());
  
  //Define a structure for the parallelism then start the slave processes
  parallelDefStruct *parallelDef=
  newParallelDef(NX,NY,NZ,
                 theCommandParser->nxProc(),
                 theCommandParser->nyProc(),
                 theCommandParser->nzProc());
  
  NumProcs=
  theCommandParser->nxProc()*
  theCommandParser->nyProc()*
  theCommandParser->nzProc();
  startProcesses(argv[0],NumProcs,
                 !theCommandParser->_doCowork);
  tEprintf(Verbose,"MPI_Init complete at %i s\n",time(0)-StartTime);
  
  //allocate the buffer for sending and receiving messages from slaves
  //below is for general and source messaging.  Any selector messaging is
  //NOT accounted for.  Receiver re-allocation is done there.
  int bufferSize = MAX(10000,MIN(MAX_TIME_STEPS_SEND,NT)*sizeof(float)+100);
  setMessageBuffer(bufferSize);
  
  //Create the master model. This will be communicating with the
  //slaves to create the slave models
  masterSgfdModel* model;
  masterFixedMediaModel* mmModel=NULL;
  masterSgfdDependent* dependent;
  //Compact the holbergPcent into the fdCoeffs
  if(theCommandParser->_fdCoeffOrder==0) theCommandParser->_fdCoeffs[0]=theCommandParser->_holbergPcent;
  //We have two types of models: non-attenuative and attenuative, both fixed acoustic
  if(theCommandParser->_qSelectors) {
    model = new masterFixedAttenMediaModel(theCommandParser->_qSelectors,Tids,
                                           FIXED_MEDIA_ATTEN,
                                           theCommandParser->_surfaceBCMode,
                                           theCommandParser->_gridMultiplier,
                                           theCommandParser->modelDef(),parallelDef,
                                           theCommandParser->_modelName,
                                           theCommandParser->_fdCoeffOrder,
                                           theCommandParser->_fdCoeffs,
                                           theCommandParser->_printSlaveModelRead);
  } else {
    switch(theCommandParser->_mediaType){
      case FIXED_MEDIA:
        model=new masterFixedMediaModel(Tids,theCommandParser->_mediaType,
                                         theCommandParser->_surfaceBCMode,
                                         theCommandParser->_gridMultiplier,
                                         theCommandParser->modelDef(),parallelDef,
                                         theCommandParser->_modelName,
                                         theCommandParser->_fdCoeffOrder,
                                         theCommandParser->_fdCoeffs,
                                         theCommandParser->_printSlaveModelRead);
        break;
      default:
        assert(FALSE,
               "Unknown media type %i",
               theCommandParser->_mediaType);
    }
  }
  
  //The slave dependents get their info from the model.
  dependent=new masterSgfdDependent(model);

  //If required, modify the time vector that is read from the
  // model.
  if(theCommandParser->_newT){
    model->resetTimeVector(theCommandParser->_modelDef.minT,
                           theCommandParser->_modelDef.NT,
                           theCommandParser->_modelDef.dt);
  }else{
    char buffer[512]="NULL";
    initSend();
    sendMessage(AllProcesses,MESSAGE_GENERAL,"s",buffer);
  }
  
  
  //Create a loop to do the main work for the master process.
  acousticMasterProcess* mainLoop=
  new acousticMasterProcess(model,mmModel,dependent,theCommandParser);
  
  //Finish initialization.
  mainLoop->completeInit();
  
  //Do the run.
  mainLoop->doRun();
  mainLoop->runtime("Total Run Time is");
  
  //Exit.
  delete mainLoop;
  delete theCommandParser;
  _mm_setcsr( oldMXCSR );  //restore standard demnormal support
  return 0;
}

