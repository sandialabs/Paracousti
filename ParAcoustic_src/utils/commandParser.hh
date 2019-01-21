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
 *  commandParser.hh
 *
 *
 *  Declares base class for parsing command line arguments and a helper class.
 *
 *  Declares the following classes:
 *  flagDescription
 *  commandParser
 *
 *  Defines the following macros:
 *  FLAG_REQUIRED
 *  FLAG_OPTIONAL
 *  FLAG_REPEATABLE
 *  FLAG_EXPERT
 *
 */
#ifndef _commandParser
#define _commandParser

#include "nstdutil.hh"
#include "array.hh"

extern int Verbose;
extern char *CommandLine;
extern char *CDFSourceName;

//Define a new class for parsing command line arguments. Subclass this to customize 
// acceptable args.

#define FLAG_REQUIRED 1
#define FLAG_OPTIONAL 2
#define FLAG_REPEATABLE 3
#define FLAG_EXPERT 4

///Define a class for the help output. This class holds a flag and some descriptive
/// info about the flag
class flagDescription{
public:
  int _flagType;
  char _flag[16],_fields[512];
  stringArray* _flagDescriptions;

  flagDescription(int flagType,const char* flag,const char* flagFields){
    _flagType=flagType;
    strcpy(_flag,flag);
    strcpy(_fields,flagFields);
    _flagDescriptions=new stringArray;
  }
  ~flagDescription(){
    for(int i=0;i<_flagDescriptions->size();free((*_flagDescriptions)[i++]));
    delete _flagDescriptions;
  }
  char* addDescription(const char* line){
    char* newString;
    assert((newString=(char*)malloc((strlen(line)+1)*sizeof(char)))!=NULL,
	   "addFlagDescription--unable to allocate %i chars for line %s: %s",
	   strlen(line)+1,_flag,line);
    strcpy(newString,line);
    _flagDescriptions->Add(newString);
    return newString;
  }

  void printUsage(char flagTypeDescription[5][256],int printExpert);
};
typedef flagDescription *flagDescriptionPtr;
typedef arrayI<flagDescriptionPtr> flagDescriptionArray;

///Define a new class for parsing command line arguments. Subclass
///this to customize acceptable args. Note this is a pure virtual
///class until processNextArgument, checkNextArgument and
///setDefaultValues are defined.
class commandParser{
protected:
  char _executable[1024];
  int _commandLineSize;

  //Trying to find a better way to build a useful help/diagnosis of bad call string. Try
  // maintaining 3 arrays with the flag, it's call arguments, and a description of what it
  // does.
  flagDescriptionArray *_flags;
  char _flagTypeDescription[5][256];

  char *_usageString;
public:
  ///Only 1 constructor, set the default values, this includes the usage string.
  commandParser(int argc,char* argv[]);

  virtual ~commandParser(){
    if(_usageString)
      free(_usageString);

    if(_flags){
      _flags->DeleteElems();
      delete _flags;
    }
  }
  ///Here are the two functions that must be redefined by subclasses to read the desired
  /// arguments and flags.
  virtual int processNextArgument(int argc,char* argv[],int& i)=0;
  virtual int checkProcessFlag(int argc,char* argv[],int& i)=0;

  ///Process the arguments, modifications to behaviour can easily be made 
  /// through methods processNextArgument or checkProcessFlag.
  virtual int processArgs(int argc,char* argv[],
                          int recursionDepth=0,int startArg=1);
protected:
  //Here is another pure virtual function. This one should be called by the initializer, it
  // needs to set default values for any local state variables and it should also build a
  // reasonable set of flag descriptions.
  virtual void setDefaultValues()=0;

  //Do a recursive call after reading arguments out of the given file.  
  int recursiveCallProcessArgs(char* filename,char* arg0,
                               int recursionDepth);

  //
  //Add a new flag description, assume the strings are all volatile so we need to make
  // local copies.
  int addFlagDescription(int flagType,
			 const char* flag,const char* flagField,
			 const char* description,
                         int extraDescriptionCount=0,...);

  //Here is a method to build a usage string from the flag descriptions. This should is called
  // at the beginning of process args.
  char* buildUsageString(char* argv0);

  ///Print the usage instructions on to stderr.
  void printUsage(int doExit,int printExpert=FALSE);
#ifdef STORE_RUN_HISTORY
  ///Write a history of calls.
  void printHistory(int nhistory);
  ///Repeat a command that is read from the history.
  void repeatCommand(int index);
#endif //#ifdef STORE_RUN_HISTORY
};

#endif //_commandParser
