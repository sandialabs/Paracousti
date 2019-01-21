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
 *  extra_output.hh
 *
 *
 *  Declares several classes that control slice and volume
 *  output for wavefield visualization.
 *
 *  Declares the following classes:
 *  extraOutput
 *  sliceOutput
 *  wavefieldOutput
 *  extraOutputGroup
 *
 *  Declares the following function:
 *  extraOutputCompFunc
 *
 */

#ifndef _extra_output_hh_
#define _extra_output_hh_

class masterSgfdDependent;
class masterSgfdModel;
class sourceNetwork;

#include "array.hh"
#include "sgfd.h"

///Start with a generic class for any type of extra output, there is always
/// a time when the extra output is to be produced and a name of the ouput
/// file. This is a pure virtual until type and doOutput methods are defined. This
/// contains methods that will generate an error to access the plane, comp, and coord
/// fields of a slice since that is by far the most common type of extra output.
class extraOutput{
public:
  float _t;
  float _dummy;

  extraOutput(){
    _t=-1;
  }
  extraOutput(float t){
    _t=t;
  }
  virtual ~extraOutput(){}

  float t(){return _t;}

  ///All derived classes must return a type.
  virtual int type()=0;
  ///All derived classes must return an associated time.
  float& outputTime(){return _t;}

  //Since most extra output is slices define some accessors for slice parameters
  // here in the base class.
  virtual int plane(){
    assert(FALSE,"extraOutput::plane--not defined for this derived class");
    return FALSE;
  }
  virtual int comp(){
    assert(FALSE,"extraOutput::comp--not defined for this derived class");
    return FALSE;
  }
  virtual float& coord(){
    assert(FALSE,"extraOutput::coord--not defined for this derived class");
    return _dummy=0.0;
  }

  ///Here is the routine that all extraOutput derived classes must implement to
  /// actually create the desired output file.
  virtual int doOutput(char *outfileName,masterSgfdDependent* forwardProblem,
		       const char* varName,float* work,
		       int **sliceIndices)=0;
#ifdef _serial_elasti_hh_
  virtual int doOutput(char *outfileName,serialDependent* forwardProblem,
		       const char* varName,float* work,
		       int **sliceIndices)=0;
#endif
};
typedef extraOutput* extraOutputPtr;
typedef arrayI<extraOutputPtr> extraOutputArray;
int extraOutputCompFunc(const void* t1,const void* t2);

//First major subclass of extra output is for slices. In addition to the time
// we need the component, the plane, and the coordinate.
#define SLICE_OUTPUT_TYPE 1

///First major subclass of extra output is for slices. In addition to the time
/// we need the component, the plane, and the coordinate.
class sliceOutput:public extraOutput{
public:
  int _comp,_plane;
  float _coord;

  sliceOutput(float t,int comp,int plane,float coord)
    :extraOutput(t){
    _comp=comp;
    _plane=plane;
    _coord=coord;
  }
  virtual int type(){return SLICE_OUTPUT_TYPE;}

  //Redefine plane, comp, and coord methods.
  virtual int plane(){return _plane;}
  virtual int comp(){return _comp;}
  virtual float& coord(){return _coord;}


#ifdef _serial_elasti_hh_
  virtual int doOutput(char *outfileName,serialDependent* forwardProblem,
		       const char* varName,float* work,
                       int **sliceIndices);
#endif
  virtual int doOutput(char *outfileName,masterSgfdDependent* forwardProblem,
		       const char* varName,float* work,
                       int **sliceIndices);
protected:
  virtual int writeCDFSlice(char* outfileName,int sliceNum,const char* varName,
			    int plane,int comp,float time,float position,
			    float* data,
                            int procLim[6]);
};

#define WAVEFIELD_OUTPUT_TYPE 2
///Second major subclass of extra output is for full wavefields. Since everything
/// is going to be written here we do not need much extra information.
class wavefieldOutput:public extraOutput{
public:
  char* _outfileName; //override outfileName in doOutput if set.

  ///Number of variables.
  int _numVars;
  ///Names of the variables need to be stored.
  char *_varNames[12];

  float _fieldScalar;

  ///Need a null constructor for derived subclasses.
  wavefieldOutput():extraOutput(0.0){
    _numVars=0;

    _fieldScalar=0.0;
  }

  ///Initializer takes a time a scalar value to apply to the variable before 
  /// it is ouput and an arbitrary number of variable names that are to
  /// written to the output file. The slave dependent that gets the message
  /// from this class must know how to get the data that is to be sent back.
  wavefieldOutput(float t,float fieldScalar,
		  int nVars,...);
  virtual ~wavefieldOutput(){
    for(int i=0;i<_numVars;i++){
      free(_varNames[i]);
    }
  }

  int addVar(char *varName){
    assert((_varNames[_numVars]=(char*)malloc(512*sizeof(char)))!=NULL,
	   "wavefieldOutput:::addVar--unable to allocate %i chars for variable name %i",
	   512,_numVars);
     strcpy(_varNames[_numVars],varName);
     return ++_numVars;
  }

  char* setFileName(char *outfileName){
    assert((_outfileName=(char*)malloc((strlen(outfileName)+1)*sizeof(char)))!=NULL,
	   "wavefieldOutput::setFileName--unable to allocate space for name %s",
	   outfileName);
    strcpy(_outfileName,outfileName);
    return _outfileName;
  }
  virtual int type(){return WAVEFIELD_OUTPUT_TYPE;}
#ifdef _serial_elasti_hh_
  virtual int doOutput(char *outfileName,serialDependent* forwardProblem,
		       const char* varName,float* work,
                       int **sliceIndices);
#endif
  virtual int doOutput(char *outfileName,masterSgfdDependent* forwardProblem,
		       const char* varName,float* work,
                       int **sliceIndices);
};

/*! \brief Here is container class that holds all of the extra output.

  This is the class that the program will allocate and use to check if it is time
  to create some extra output.
*/
class extraOutputGroup{
public:
  //Variables relating to slice output.
  long long _maxNumPerSliceFile;
  int _currFileIndex, _numSliceFiles;
  int* _sliceFileIndices;
  char _sliceOutputName[1024];
  int **_sliceIndices, **_sliceIndicesGlobal;
  float* _work;

  stringArray* _sliceCompNames;
  stringArray* _slicePlaneNames;

  //Variables relating to wavefield output.
  char _wavefieldOutputName[1024];

  //General variables.
  extraOutputArray* _theOutput;
  int _nextOutput;
  int _useModelExtraOutput; //-EM

  //Initializer.
  extraOutputGroup(modelDefStruct* modelDef);

  ///Method to write restart information required to correctly restart the extra
  /// output.
  virtual FILE* doCheckpoint(int iteration,char* cpDir,int cpID,
			     int callDepth,
                             FILE* cpFile=NULL);
  ///Here is a method to read a checkpoint file written by doCheckpoint.
  virtual FILE* readCheckpoint(masterSgfdModel* model,
			       int& iteration,char* cpDir,int cpID,
			       int callDepth,
                               FILE* cpFile=NULL);
 
  int size(){return _theOutput->size();}

  int type(int i){return (*_theOutput)[i]->type();}

  float t(int i){return (*_theOutput)[i]->t();}
  float& outputTime(int i){return (*_theOutput)[i]->outputTime();}
  float coord(int i){return (*_theOutput)[i]->coord();}
  int plane(int i){return (*_theOutput)[i]->plane();}
  int comp(int i){return (*_theOutput)[i]->comp();}

  virtual int numPlanes(){return _slicePlaneNames->size();}
  virtual char* planeText(int i){
    return (*_slicePlaneNames)[i];
  }
  virtual char* slicePlaneText(int i){
    return planeText((*_theOutput)[i]->plane());
  }

  virtual int numComps(){return _sliceCompNames->size();}
  virtual char* compText(int i){
    return (*_sliceCompNames)[i];
  }
  virtual char* sliceCompText(int i){
    return compText((*_theOutput)[i]->comp());
  }

  int iterationsToNext(float currT,float dt){
    if(_nextOutput<_theOutput->size()){
      int num=(int)rint(((*_theOutput)[_nextOutput]->t()-currT)/dt);
      return MAX(1,num);
    }
    return 0;
  }


  virtual sliceOutput* makeNewSlice(float t,int comp,int plane,float coord,
				    int verboseAdd=FALSE){
    sliceOutput* newSlice=new sliceOutput(t,comp,plane,coord);
    _theOutput->Add(newSlice);
    tEprintf(verboseAdd,
	     "Added output slice %s.%s; time %f; coord %.2f\n",
	     planeText(newSlice->plane()),
	     compText(newSlice->comp()),
	     newSlice->t(),newSlice->coord());
    return newSlice;
  }

  //Initialize output files (actully only the slice file needs to be initialized).
  int initOutput(masterSgfdModel* model,
		 sourceNetwork* sources,int allocateStorage,
                 int fileInit=TRUE);

  //Check if extra output is required at the current iteration.
#ifdef _serial_elasti_hh_
  int doOutput(serialDependent* forwardProblem,
               float t_vel,float t_press);
#endif
  int doOutput(masterSgfdDependent* forwardProblem,
               float t_vel,float t_press);

  //Write the extra output to a file for later reading, this should match the 
  // format that is read by readExtraOutput.
  int writeOutput(char* fileName);

  //Read extra output from a cdf file (model definition, this is the old way
  // of defining slices while generating the model).
  int readExtraOutput(char* infileName);

  ///Parse requests for additional output from the command line, note that some
  /// requests can not be honored unless this is a run call and not a generate
  /// model call.
  virtual int addExtraOutput(int& i,int argOffset,
			     int argc,char* argv[],
			     modelDefStruct* modelDef,
                             int runCall);

  void addWavefieldOutput(char* filename,
			  float t,float fieldScalar,
                          int nVars,...);
  void doWavefieldOutput(char* filename,
			 masterSgfdDependent* forwardProblem,
			  float t,float fieldScalar,
                         int nVars,...);
protected:
  int addNewSliceComp(const char* name){
    char *buffer;
    assert((buffer=(char*)malloc((strlen(name)+1)*sizeof(char)))!=NULL,
	   "addNewSliceComp--unable to allocate space for new slice type %s",
	   name);
    strcpy(buffer,name);
    _sliceCompNames->Add(buffer);
    return _sliceCompNames->size();
  }
  int addNewSlicePlane(const char* name){
    char *buffer;
    assert((buffer=(char*)malloc((strlen(name)+1)*sizeof(char)))!=NULL,
	   "addNewSlicePlane--unable to allocate space for new slice type %s",
	   name);
    strcpy(buffer,name);
    _slicePlaneNames->Add(buffer);
    return _slicePlaneNames->size();
  }
  void sortTimes(){
    _theOutput->sort(extraOutputCompFunc);
  }

  void initializeSliceVariables(masterSgfdModel* model,int allocateStorage);
  
  //Make these virtual so derived classes that need to define new slice types can 
  // do so.
  virtual int sliceComp2Num(char* compName);
  virtual float sliceCompTimeRasterOffset(int compNum){
    if(compNum!=4 && compNum!=5) 
      return 1./2;
    return 0.;
  }
  virtual int slicePlane2Num(char* planeName);
};
typedef extraOutputGroup* extraOutputGroupPtr;
#endif //#ifndef _extra_output_hh_
