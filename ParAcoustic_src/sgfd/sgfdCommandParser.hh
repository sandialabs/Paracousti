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
 *  sgfdCommandParser.hh
 *
 *
 *  Declares the sgfdCommandParser class.  This class
 *  forms the base class for all the finite difference codes and
 *  contains command line flags/options that are common to all of these
 *  codes.
 *
 *  Declares the following class:
 *  sgfdCommandParser
 *
 */

#ifndef _sgfdCommandParser_hh_
#define _sgfdCommandParser_hh_

class sourceNetwork;
class masterReceiverNetwork;

#include "commandParser.hh"

#include "sgfd.h"

///sgfdCommandParserBase localizes some of the commands that are
/// standardized across various sgfd codes. This should reduce some of the
/// redundancys.
///Because of the the final command parser needs to know about extraOutput, which
/// needs to know about the model and dependent classes this is a virtual class.
/// The final version is in sgfd.hh below the model and dependent definitions.
class sgfdCommandParserBase: public commandParser{
public:
  //Here is a variable required by the parser for this program but not by general
  // parsers.
  modelDefStruct _modelDef;
  
  //And here are a set of internal state variables to replace some of the global variables
  // used in the previous scheme.
  char _modelName[1024]; //argument
  float _holbergPcent; //Adjustments to finite-difference coefficients; -h flag
  float _fdCoeffs[MAX_NUM_SPACE_COEFFS]; //user input FD coeffs, allowing up to order 2*MAX_NUM_SPACE_COEFFS; -hc flag
  int _fdCoeffOrder;  //order of input FD coeffs; -hc flag
  float* _newT;  //Change the time vector, -T flag.
  
  //Flags to determine how the run looks.
  int _traceOutputIterations,_iterationsPerStep; //-t (-Rt), and -j.
  int _printSlaveModelRead; //-vr flag
  
  //Boundary conditions.
  // -bS option use spongy absorbing boundries
  int _spongeBCNodes[6];
  float _spongeBCValue[6];
  int _surfaceBCMode; //-bF option for free surface.
  
  //Initial condition.
  char _initialConditionName[1024];
  
  //Checkpoint flags.
  int _doCheckpoint;
  int *_checkpointIndicies;
  char _checkpointDir[1024];
  
  int _nxProc,_nyProc,_nzProc; //-p flag
  
  //May want to add extra receivers or sources that are not defined in the
  // earthmodel file
  masterReceiverNetwork* _extraReceivers; //add with -R[?] flags
  sourceNetwork* _extraSources;     //add with -S[?] flags
  
  //
  //Constructor
  sgfdCommandParserBase(int argc,char* argv[]):
  commandParser(argc,argv){
    if(argc<2) return;
    
    //Do an immediate check for the first argument of the form -help or --help.
    if(argc>1 &&
       (!strcmp(argv[1],"-help") || !strcmp(argv[1],"--help")))
      printUsage(TRUE);
  }
  
  //Here are accessor functions for the local state variables that replace the globals I
  // used to use to know about arguments.
  modelDefStruct* modelDef(){return &_modelDef;}
  int nxProc(){return _nxProc;}
  int nyProc(){return _nyProc;}
  int nzProc(){return _nzProc;}
  
  ///Here is a virtual method to pack the message for the boundary condition
  /// initialization.
  virtual int sendBCMessage();
  
  //New virtual function that can be tailored to send any extra info that is needed
  virtual int sendExtraInfo() {return 0;}
  
  ///Make a minor mod to processArgs; check that the model is initialized at the
  /// end of the call.
  virtual int processArgs(int argc,char* argv[],
                          int recursionDepth=0,int startArg=1){
    int returnVal=commandParser::processArgs(argc,argv,recursionDepth,startArg);
    
    if(recursionDepth==0)
      assert(*_modelName,
             "sgfdCommandParserBase::processArgs--model not set.");
    
    return returnVal;
  }
  
  ///Here is the required virtual function to process an argument.
  virtual int processNextArgument(int argc,char* argv[],int& i);
  
  ///Here is the required virtual function to process a flag.
  virtual int checkProcessFlag(int argc,char* argv[],int& i);
  
protected:
  ///Placeholder for a method to add the extra output flag description.
  virtual void addExtraOutputFlags()=0;
  
  ///Add flag descriptions pertaining to receivers.
  virtual void addExtraReceiverFlags();
  ///Add flag descriptions pertaining to sources.
  void addExtraSourceFlags();
  
  ///Add flag descriptions and set default values for internal variables.
  virtual void setDefaultValues();
};


#endif //#ifndef _sgfdCommandParser_hh_

