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
 *  sgfdLoops.hh
 *
 *
 *  Declares classes for process control of slaves and master
 *  Process classes.  These are the main control structures that
 *  call many other classes to do the work for each time step and
 *  for those things that happen every so many time steps.
 *
 *  Declares the following classes:
 *  sgfdMasterProcess
 *  sgfdSlaveProcess
 *
 */

#ifndef _sgfdLoops_hh_
#define _sgfdLoops_hh_

class masterSgfdDependent;
class masterSgfdModel;
class slaveSgfdDependent;
class slaveSgfdModel;
class extraOutputGroup;
class sgfdBoundaries;
class sgfdCommandParser;

#include <time.h>

#include "nstdutil.hh"

///Create a class for the master loop.
class sgfdMasterProcess{
public:
  masterSgfdModel* _model;
  masterSgfdDependent* _dependent;
  
  extraOutputGroup* _extraOutput;

  int _iterationsPerStep;
  int _startTime;
  int _traceOutputIterations,_lastTraceOutputEndIndex;
  //Checkpoint flags.
  int _doCheckpoint,_nextCheckpoint;
  int *_checkpointIndicies;
  char _checkpointDir[1024];
  
  char _tUnit[10];
  float _tMult;

  ///Constructor, to be called after the model and dependent have been initilized.
  sgfdMasterProcess(int startTime,
		    masterSgfdModel* model,masterSgfdDependent* dependent,
                    sgfdCommandParser* runArguments);

  virtual ~sgfdMasterProcess();

  ///Determine how long the current run has been going.
  time_t runtime(){return time(0)-_startTime;}

  ///Determine how long the current run has been going.
  time_t runtime(const char *printtag){
    int result=runtime();
    tEprintf(Verbose,"%s %is\n",printtag,result);
    return result;
  }

  ///Perform common initialization tasks that are done after the model
  /// and the dependent variables have been initialized.
  virtual void completeInit(char* modelName,
                            int &iteration,sgfdCommandParser* runArguments);

  ///Another initialization method for programs that only want minimal initialization.
  virtual void completeInit(char* modelName,int surfaceBCMode,
                            extraOutputGroup* extraOutput);

  ///Run the model and then do the finalUpdate.
  int doRun(int &iteration);

protected:
  ///Here is a hook for derived methods that need to stop in mainLoopWork for
  /// reasons to be determined.
  virtual int stepsToStop(int count,int iteration,int& cause){
    return count; //default version does nothing.
  }
  ///Here is a hook for derived methods that need to do something to be 
  /// determined in mainLoopWork.
  virtual void mainLoopWorkExtraWork(int iteration,int& cause){
    //default version does nothing.
  }

  ///What is done during the run, advance to the next important iteration. 
  /// That is where output needs to be done or the max number of iterations,
  /// whichever is greater.
  virtual int mainLoopWork(int &iteration);
  virtual void finalUpdate(int& iteration);
};

///Create another class for the slave loop. This hopefully encapsulates the common features
/// of the slaveLoop from my various codes. Derived classes look for their own specific
/// messages and call the parent class if there is no match.
class sgfdSlaveProcess{
public:
  int _tid,_parent; //Process id's for message passing.

  slaveSgfdModel *_model;
  slaveSgfdDependent *_dependent;
  sgfdBoundaries* _boundaries;

  ///Default constructor sets everything to NULL and initializes a buffer for message
  /// passing.
  sgfdSlaveProcess(int tid,int ptid);
  
  virtual ~sgfdSlaveProcess();

  ///Here is the initialization method, this should initialize _model and _dependent. NOTE:
  /// the subclasses should do the actual initialization on there own version of _model
  /// and _dependent (probably redefined to subclasses of slaveSgfdModel and slaveSgfdDependent)
  /// then call this version to set the variables within this class.
  virtual void doInit()=0;
  virtual void doInit(slaveSgfdModel *model,slaveSgfdDependent *dependent){
    _model=model;
    _dependent=dependent;
  }

  ///Here is the method that receives a message from the parent and acts on it, this version
  /// knows how to do the basic things. Derived class should pick off special messages and
  /// call this version if they do not find a special match.
  virtual int processMessage(int msgtag=FALSE);
};

#endif //#ifndef _sgfdLoops_hh_
