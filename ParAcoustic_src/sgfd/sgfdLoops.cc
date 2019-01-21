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
 *  sgfdLoops.cc
 *
 *
 *  Defines class functions for process control of slaves and master
 *  Process classes.  These are the main control structures that
 *  call many other classes to do the work for each time step and
 *  for those things that happen every so many time steps.
 *
 *  Defines the following class functions:
 *  sgfdMasterProcess::sgfdMasterProcess
 *  sgfdMasterProcess::~sgfdMasterProcess
 *  sgfdMasterProcess::completeInit
 *  sgfdMasterProcess::completeInit
 *  sgfdMasterProcess::doRun
 *  sgfdMasterProcess::mainLoopWork
 *  sgfdMasterProcess::finalUpdate
 *  sgfdSlaveProcess::sgfdSlaveProcess
 *  sgfdSlaveProcess::~sgfdSlaveProcess
 *  sgfdSlaveProcess::processMessage
 *
 */

#include "sgfdLoops.hh"
#include "extra_output.hh"
#include "message_passing.h"
#include "sgfd.hh"
#include "sgfdReceiverNetwork.hh"

//Here are some tags for the various reasons that we might stop inside the main 
// loop.
#define NORMAL_SYNC       1
#define EXTRA_OUTPUT_STOP 2
#define TRACE_OUTPUT_STOP 4
#define CHECKPOINT_STOP   8

///Constructor, to be called after the model and dependent have been initilized.
sgfdMasterProcess::sgfdMasterProcess(int startTime,
      masterSgfdModel* model,masterSgfdDependent* dependent,
      sgfdCommandParser* runArguments){
  _model=model;
  _dependent=dependent;

  _traceOutputIterations=runArguments->_traceOutputIterations;
  _lastTraceOutputEndIndex=0;

  _iterationsPerStep=runArguments->_iterationsPerStep;

  _startTime=startTime;

  float dt=_model->_modelDef->dt;
  if(dt > 0.1){
    //Make the units for diagnostic output into seconds.
    strcpy(_tUnit,"s");
    _tMult=1.0;
  }else{
    //Make the units for diagnostic output into milli-seconds.
    strcpy(_tUnit,"ms");
    _tMult=1000.0;
  }      

  //Default is no checkpoints.
  _doCheckpoint=runArguments->_doCheckpoint;
  _nextCheckpoint=0;
  _checkpointIndicies=runArguments->_checkpointIndicies;
  strcpy(_checkpointDir,runArguments->_checkpointDir);
}

sgfdMasterProcess::~sgfdMasterProcess(){
  //Free memory and clean up
  //Shut down the slaves
  stopProcesses(AllProcesses,-1);
  if(_dependent){
    delete _dependent;
    _dependent=NULL;
  }
  if(_model){
    delete _model;
    _model=NULL;
  }
}

///Perform common initialization tasks that are done after the model
/// and the dependent variables have been initialized.
void sgfdMasterProcess::completeInit(char* modelName,
        int &iteration,sgfdCommandParser* runArguments){
  //Set initial conditions and boundary conditions.
  int readCP = 0;
  if(runArguments->_initialConditionName[0]!='\0')
    _dependent->setInitialConditions(runArguments->_initialConditionName,
             iteration,FALSE);
  else
    readCP = _dependent->setInitialConditions(runArguments->_checkpointDir,
             iteration,FALSE);
  if(readCP){
    //Acount for checkpoints that have already been written.
    for(int i=0;i<_doCheckpoint;i++){
      if(_checkpointIndicies[i] <= iteration){
        _nextCheckpoint++;
      }
    }
  }

  runArguments->sendBCMessage();

  //Set sources, receivers, and any extra output (slices, full wavefield).
  _dependent->setSources(runArguments->_extraSources);

  _dependent->setReceivers(runArguments->_extraReceivers,runArguments->_surfaceBCMode,!readCP);
  
  //send any extra info if needed
  runArguments->sendExtraInfo();
  
  if(runArguments->_extraOutput){
    _extraOutput=runArguments->_extraOutput;
  }else{
    _extraOutput=runArguments->allocateExtraOutput();
  }

  if(_extraOutput->_useModelExtraOutput)
    _extraOutput->readExtraOutput(modelName);
  if(!readCP){
    _extraOutput->initOutput(_model,_dependent->_sources,TRUE);
  }else{
    _extraOutput->readCheckpoint(_model,iteration,
         runArguments->_checkpointDir,0,0);
  }

  tEprintf(Verbose,"Initialization complete at %is\n",
     time(0)-_startTime);
  fflush(stdout); //flush startup messages from slaves
  fflush(stderr);

  DEF_MODEL_LIMITS(_dependent->modelDef());
  _extraOutput->doOutput(_dependent,minT,minT); //Check for zero time slice output.
}

///Another initialization method for programs that only want minimal initialization.
void sgfdMasterProcess::completeInit(char* modelName,int surfaceBCMode,
        extraOutputGroup* extraOutput){
  //Set initial conditions and boundary conditions.
  char initialConditionName[512]="";
  int iteration=0;
  _dependent->setInitialConditions(initialConditionName,iteration,FALSE);

  int nodes[6]={0,0,0,0,0,0};
  float vals[6]={1.0,1.0,1.0,1.0,1.0,1.0};
  _dependent->setBoundaryConditions(surfaceBCMode,nodes,vals);

  //Set sources, receivers, and any extra output (slices, full wavefield).
  _dependent->setSources(NULL);
  _dependent->setReceivers(NULL,surfaceBCMode,FALSE);

  if(extraOutput){
    _extraOutput=extraOutput;
  }else{
    _extraOutput=new extraOutputGroup(_model->modelDef());
  }
  //_extraOutput->readExtraOutput(modelName);
  _extraOutput->initOutput(_model,_dependent->_sources,TRUE,TRUE);

  tEprintf(Verbose,"Initialization complete at %is\n",
     time(0)-_startTime);
  fflush(stdout); //flush startup messages from slaves
  fflush(stderr);

  DEF_MODEL_LIMITS(_dependent->modelDef());
  _extraOutput->doOutput(_dependent,minT,minT); //Check for zero time slice output.
}

///Run the model and then do the finalUpdate.
int sgfdMasterProcess::doRun(int &iteration){
  while(mainLoopWork(iteration));
  tEprintf(Verbose,"Exited main loop\n");

  //Update velocity vector components to the final time t=tmax+(dt/2).
  finalUpdate(iteration);

  return iteration;
}

///What is done during the run, advance to the next important iteration.
/// That is where output needs to be done or the max number of iterations,
/// whichever is greater.
int sgfdMasterProcess::mainLoopWork(int &iteration){
  int loopStartTime=time(0);
  DEF_MODEL_SIZE(_dependent->modelDef());
  DEF_MODEL_LIMITS(_dependent->modelDef());
  float currTime=minT+iteration*dt;

  //Do the vel and stress updates
  int count=MIN(_iterationsPerStep,NT-iteration);
  int cause=NORMAL_SYNC;

  // When is the next output due.
  int eoCount=_extraOutput->iterationsToNext(currTime,dt);
  if(eoCount){
    if(eoCount==count){
      cause|=EXTRA_OUTPUT_STOP;
    }else if(eoCount<count){
      count=eoCount;
      cause=EXTRA_OUTPUT_STOP;
    }
  }

  // When is the next time to write trace output.
  if(iteration&&_traceOutputIterations>0 &&
     NT-iteration>0.5*_traceOutputIterations){
    int toCount=
    _traceOutputIterations-((iteration)%_traceOutputIterations);
    if(toCount==count){
      cause|=TRACE_OUTPUT_STOP;
    }else if(toCount<count){
      count=toCount;
      cause=TRACE_OUTPUT_STOP;
    }
  }

  // When is the next checkpoint due.
  if(_nextCheckpoint<_doCheckpoint){
    int cpCount=_checkpointIndicies[_nextCheckpoint]-iteration;
    if(cpCount==count){
      cause|=CHECKPOINT_STOP;
    }else if(cpCount<count){
      count=cpCount;
      cause|=CHECKPOINT_STOP;
    }
  }

  //Leave room for other things that derived classes might need to do.
  count=stepsToStop(count,iteration,cause);
  
  //Send the work out to the slaves
  tEprintf(Verbose,
     " Advancing to model time %7.2f%s (iteration %4i)\n",
     _tMult*(currTime+count*dt),_tUnit,iteration+count);
  _dependent->advance(iteration,count);
  int currWallTime=time(0);
  tEprintf(Verbose,
     "   Wall Time: %4is; %.2fs/iteration\n",
     currWallTime-_startTime,(float)(currWallTime-loopStartTime)/count);

  //Extra output.
  if(cause&EXTRA_OUTPUT_STOP){
    currTime=minT+iteration*dt;
    float tVel=currTime-0.5*dt;

    _extraOutput->doOutput(_dependent,tVel,currTime);
  }

  //Trace output.
  if(cause&TRACE_OUTPUT_STOP){
    _dependent->writeTraces(_lastTraceOutputEndIndex,iteration-1);
    _lastTraceOutputEndIndex=iteration-2;
  }

  //Checkpoint.
  if(cause&CHECKPOINT_STOP){
    _nextCheckpoint++;
    _dependent->doCheckpoint(iteration,_checkpointDir,0,0);
    _extraOutput->doCheckpoint(iteration,_checkpointDir,0,0);
  }

  //Leave room for other things that derived classes might need to do.
  mainLoopWorkExtraWork(iteration,cause);

  return iteration<NT;
}

void sgfdMasterProcess::finalUpdate(int& iteration){
  DEF_MODEL_LIMITS(_dependent->modelDef());

  //Update velocity vector components to the final time t=tmax+(dt/2).
  tEprintf(Verbose," Performing last velocity update\n");
  _dependent->advanceVel(iteration+1);

  //write slices for final time
  tEprintf(Verbose && _extraOutput->size(),"Checking final slices\n");
  _extraOutput->doOutput(_dependent,
       minT+iteration*dt+0.5*dt,minT+iteration*dt);

  //print output at reciever locations
  tEprintf(Verbose && _dependent->_receivers->size(),
     "Writing final trace output\n");
  _dependent->writeTraces(_lastTraceOutputEndIndex,iteration+1);
}

///Default constructor sets everything to NULL and initializes a buffer for message
/// passing.
sgfdSlaveProcess::sgfdSlaveProcess(int tid,int ptid){
  _tid=tid;
  _parent=ptid;
  
  _model=NULL;
  _dependent=NULL;
  _boundaries=NULL;
  
  //allocate an initial buffer to receive the startup message
  // from the parent process
  setMessageBuffer(1024+2048*sizeof(int));
}

sgfdSlaveProcess::~sgfdSlaveProcess(){
  if(_model){
    //       delete _model;
    _model=NULL; //So there is no double deletion.
  }
  if(_dependent){
    //        delete _dependent;
    _dependent=NULL; //So there is no double deletion.
  }
  if(_boundaries){
    //       delete _boundaries;
    _boundaries=NULL;
  }
  
  sendMessage(_parent,MESSAGE_EXIT,NULL);
  
  getMessage(_parent,MESSAGE_FINAL_EXIT,NULL);
  stopProcesses(_parent,0);
}

///Here is the method that receives a message from the parent and acts on it, this version
/// knows how to do the basic things. Derived class should pick off special messages and
/// call this version if they do not find a special match.
int sgfdSlaveProcess::processMessage(int msgtag){
  if(!msgtag){
    //get control message from the parent
    msgtag=getMessage(_parent,AnyMessageTag,NULL);  
  }

  //parse and process the message
  switch(msgtag){
  case MESSAGE_EXIT:
    //The work that was done here is now done in the deletion of the sgfdSlaveProcess.
    // Just return FALSE so we know to delete.
    return FALSE;

  case MESSAGE_SLICE:
    _dependent->slicer(_parent);
    break;
  case MESSAGE_FULL_FIELD:
    _dependent->fullFieldOutput(_parent);
    break;

  case MESSAGE_SEND_RECEIVER:
    {
      int index,startIndex,stopIndex;
      unpackMessage("iii",&index,&startIndex,&stopIndex);
      initSend();
      _dependent->packReceiverData(index,startIndex,stopIndex);
      sendMessage(_parent,MESSAGE_SEND_RECEIVER,NULL);
    }
    break;
  case MESSAGE_SEND_RECEIVER_GRID:
    {
      int index;
      unpackMessage("i",&index);
      _dependent->sendReceiverGridData(index);
    }
    break;

    //Update stress|velocity, advance 1/2 of a time step. These 
    // send and load the edges from neighbors before returning.
  case MESSAGE_ADVANCE_ONE_ITERATION:
    {
      int iteration,count;
      unpackMessage("ii",&iteration,&count);
      _dependent->advance(iteration,count);
    }
    break;
  case MESSAGE_UPDATE_VEL:
    {
      int iteration;
      unpackMessage("i",&iteration);
      _dependent->advanceVel(iteration,TRUE);
    }
    break;
  case MESSAGE_UPDATE_STRESS:
    {
      int iteration;
      unpackMessage("i",&iteration);
      _dependent->advanceStress(iteration,TRUE);
    }
    break;

    //
    //Default option
    // Unknown command
  default:
    tEprintf(Verbose,"PID[%i]: unknown control option %i\n",
       _tid,msgtag);
    initSend();
    sendMessage(_parent,MESSAGE_FAIL,"P",_model->procLim());
  }

  //Return TRUE, FALSE is returned only if the MESSAGE_EXIT message is processed.
  return TRUE;
}
