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
 *  parallel_acousti.hh
 *
 *
 *  This file contains declarations for classes parallelAcoustiCommandParser
 *  and acousticSlaveProcess.  The first of these is in charge of parsing
 *  the command line for paracousti-specific flags and options.  The second
 *  initializes slaves and handles some paracousti-specific messages.
 *
 *  Defines the following classes:
 *  parallelAcoustiCommandParser
 *  acousticMasterProcess
 *  acousticSlaveProcess
 *
 */
#ifndef _parallel_acousti_hh_
#define _parallel_acousti_hh_

class acousticSgfdBoundaries;
class masterFixedMediaModel;
class slaveAcousticModel;
class slaveAcousticDependent;
class selector;
class extraOutputGroup;

#include "sgfd.hh"
#include "sgfdLoops.hh"

//
//Here are references to the remaining global variables used for
//argument processing. Others are now state variables in the
//commandParser class.
extern char* CommandLine;
extern int Verbose;

/*Here is the real class that parses the arguments to the program.

  This class also holds most of the state variables used to run the program.
  I think this is an improvement over the way I used to do things with too
  many global variables.
*/
class parallelAcoustiCommandParser:public sgfdCommandParser{
public:
  //Options for initial conditions: -i[?]
  char _initialConditionName[512]; //-i
  int _iCNegateVelocities;      //-in
  int _usePML;                  //-bp
  float _xpmlFac; //-bpx
  float _alphaValues[6], _kValues[6]; //-bpc

  int _gridMultiplier;         //-D flag
  int _doCowork; //-pc (stationary media only).

  //And of course, flags having to do with media type.
  int _mediaType;
  
  //Attenuation based upon a series of standard linear fluid (spring and dashpot) input
  //Activate with the -Q family of flags
  selector* _qSelectors;

  //
  //Constructor
  parallelAcoustiCommandParser(int argc,char* argv[]):
    sgfdCommandParser(argc,argv){
    if(_flags)
      setDefaultValues();
  }

  ///Here is a specialized boundary condition setting method. The list of variables
  /// sent here needs to match what is read in elasticSlaveProcess doInit.
  virtual int sendBCMessage();
  
  //Here are the required virtual functions to process and argument or a flag.
  virtual int processNextArgument(int argc,char* argv[],int& i){
    if(!sgfdCommandParser::processNextArgument(argc,argv,i)){
      assert(FALSE,
	     "%s--requrires exactly 1 argument for model (already set to %s)",
	     argv[0],_modelName);
    }
    return i;
  }

  virtual int checkProcessFlag(int argc,char* argv[],int& i);

protected:
  virtual extraOutputGroup* allocateExtraOutput();

  virtual void setDefaultValues();
};

class acousticMasterProcess: public sgfdMasterProcess{
protected:
  int _iteration;
  parallelAcoustiCommandParser* _theCommandParser;
  masterFixedMediaModel* _mmModel;

public:
  acousticMasterProcess(masterSgfdModel* model,masterFixedMediaModel* mmModel,
			masterSgfdDependent* dependent,
			parallelAcoustiCommandParser* theCommandParser):
    sgfdMasterProcess(StartTime,model,dependent,theCommandParser){
    _theCommandParser=theCommandParser;
    _mmModel=mmModel;
    _iteration=0;

  }

  ///For convience redefine completeInit; call
  /// the sgfdMasterProcess version with most of the arguments derived from the 
  /// internally stored commandParser.
  void completeInit(){
    sgfdMasterProcess::completeInit(_theCommandParser->_modelName,
				    _iteration,_theCommandParser);
  }

  ///For convience redefine doRun, just calls the sgfdMasterProcess version with the
  /// internal _iteration variable for one of the arguments.
  int doRun(){
    return
      sgfdMasterProcess::doRun(_iteration);
  }

};

class acousticSlaveProcess: public sgfdSlaveProcess{
public:
  //Redefine model and dependent to the acoustic versions.
  acousticSgfdBoundaries* _boundaries;
  slaveAcousticModel *_model;
  slaveAcousticDependent *_dependent;

  /*Initialization info is read from messages so the constructor only
    needs the tid and the parent's tid.

    Note that derived classes might want to delay initialization.
  */
  acousticSlaveProcess(int tid,int ptid,int doInitNow=TRUE):
    sgfdSlaveProcess(tid,ptid){
    if(doInitNow)
      doInit();
  }

  virtual ~acousticSlaveProcess(){
     tEprintf(Verbose,
	     "[%i] Acoustic Slave Process Loop Exiting\n",
	     _tid);
  }

  ///Read details of how to initialize from messages from the master.
  virtual void doInit();

  /*Here is the method that receives a message from the parent and
    acts on it.

    This version is reimplemented from the the sgfdSlaveProcess. The parent knows how to process
    advance and most IO messages (the model or the dependent do the work).
  */
  virtual int processMessage(int msgtag=FALSE);

protected:
  /*Need to set slightly differenct BC's.

  4th order velocity material derivative message passing.
   
  -  x-5  x-5  x-4  x-4  x-3  x-3  x-2  x-2  x-1  x-1
  -   p    v    p    v    p    v    p    v    p    v
  -   C    C    C    C    C    C   P^   P^   P^   P^
  -            Pv   Pv   Pv   Pv    C    C    C    C    C    C
  -             p    v    p    v    p    v    p    v    p    v
  -             0    0    1    1    2    2    3    3    4    4
  */
  void setBoundaryConditions(int mediaType);
};

#endif //#ifndef _parallel_acousti_hh_
