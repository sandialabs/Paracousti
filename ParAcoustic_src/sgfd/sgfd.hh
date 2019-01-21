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
 *  sgfd.hh
 *
 *
 *  Declare several class functions that control the base class behavior
 *  for many aspects of the algorithm including models, slaves, masters,
 *  dependents and command line processing.
 *
 *  Declares the following classes:
 *  sgfdModel
 *  sgfdDependent
 *  masterSgfdModel
 *  slaveSgfdModel
 *  masterSgfdDependent
 *  slaveSgfdDependent
 *  sgfdCommandParser
 *
 *  Declares the following functions:
 *  newModelDef
 *  packModelDef
 *  unpackModelDef
 *  copyModelDef
 *  newParallelDef
 *  packParallelDef
 *  unpackParallelDef
 *  copyParallelDef
 *  readCDFRunParams
 *  readCDFRunParams
 *
 */
#ifndef _sgfd_hh_
#define _sgfd_hh_

class sourceNetwork;
class masterReceiverNetwork;
class slaveReceiverNetwork;
class extraOutputGroup;

#include "nstdutil.hh"

//The sgfd.h header is for the structures that are going to be
//required by both C and C++ files. In the original code this was
//elastic_structure.h
#include "sgfd.h"

//
//Procedure Declarations
//

//Subroutines for use with the modelDefStruct. 
// It would be better to define this as a class but it needs to be
// usable by c as well as c++ procedures.

///Make a new modelDef structure.
modelDefStruct* newModelDef(int xStart,int xStop,
			    int yStart,int yStop,
			    int zStart,int zStop,
			    int nt,
			    float minX,float dx,float minY,float dy,
			    float minZ,float dz,float minT,float dt,
			    float scalarSpeed,float scalarVel,
			    float scalarDen,float scalarStress,
			    int multiplier=1);
///Pack a message with an existing modelDef.
modelDefStruct* packModelDef(modelDefStruct* modelDef,int destroy=FALSE);
///Unpack a message into a new or existing modelDef struct.
modelDefStruct* unpackModelDef(modelDefStruct* modelDef=NULL);
///Copy a modelDef structure.
void copyModelDef(modelDefStruct* target,modelDefStruct *source);

///Make a new parallelDef structure.
parallelDefStruct* newParallelDef(int globalNX,int globalNY,int globalNZ,
				  int nxProc,int nyProc,int nzProc,
				  int procI=-1,int procJ=-1,int procK=-1,
				  int *neighbors=NULL);
///Pack a message with an existing parallelDef.
parallelDefStruct* packParallelDef(parallelDefStruct* parallelDef,int destroy=FALSE);
///Unpack a message into a new or existing parallelDef struct.
parallelDefStruct* unpackParallelDef(parallelDefStruct* parallelDef=NULL);
///Copy a parallelDef structure.
void copyParallelDef(parallelDefStruct* target,parallelDefStruct* source);

///Read the run parameters from a NetCDF file to modelDef struct.
int readCDFRunParams(const char* fileName,modelDefStruct* modelDef,
		     int print=FALSE);
///Read the run parameters from a NetCDF file to individual vars.
int readCDFRunParams(const char* fileName,
		     int& nx,int& ny,int& nz,int& nt,
		     int& nxy,int& nxyz,
		     float& minX,float& minY,float& minZ,float& minT,
		     float& dx,float& dy,float& dz,float& dt);

//The sgfdCommandParser reads the command line arguments that control the 
// program behaivior.
#include "sgfdCommandParser.hh"


//
//Classes

///Generic staggered-grid finite-difference model.

///This class encapsulates all features required for a model. There are subclasses for
/// master/slave and elastic/acoustic models.
class sgfdModel{
public:
  //Here are the variables needed by all types of model.
  int _gridMultiplier;

  modelDefStruct* _modelDef;
  char _modelName[1024];

  float _vMin,_vMax,_temp;

  //Sometimes we want to generate a unique filename for a given process. To pass
  // this back we need some permenant space, here it is.
  char _uniqueFilename[1024];

  /// Null initializer, does nothing but used by some subclasses.
  sgfdModel(){
    _gridMultiplier=1;

    _modelDef=NULL;
    _modelName[0]='\0';
  }    

  // Simple initilizer, slightly more than the null case.
  sgfdModel(modelDefStruct* modelDef,char* modelName=NULL){
    _gridMultiplier=1;

    _modelDef=modelDef;
    if(modelName){
      strcpy(_modelName,modelName);    
    }else{
      _modelName[0]='\0';
    } 
  }

  virtual ~sgfdModel() {}
  
  //! Name of the model.
  char* modelName(){return _modelName;}

  //! Definition of the model size, limits, etc.
  modelDefStruct* modelDef(){return _modelDef;}
  //! Definition of the parallel decomposition (program failure if this version called).
  virtual parallelDefStruct* parallelDef(){
    assert(FALSE,"sgfdModel::parallelDef--undefined in this class");
    return NULL;
  }
  //! Volume of 1 cell.
  float modelCellVolume(){return modelDef()->dx*modelDef()->dy*modelDef()->dz;}

  ///Min velocity (used for dispersion check).
  float& vMin(){return _vMin;}
  ///Max velocity (used for stability CFL check).
  float& vMax(){return _vMax;}
  ///Here is a function for "other" min and max values. In the acoustic problem
  /// these are the min and max values of vx, vy, and vz; indexed by i=1 to 6.
  virtual float& otherMinMax(int i){
    assert(FALSE,"sgfdModel::otherMinMax--has no possible calls");
    return _temp=0.0;
  }

  /// Easier access for various model parameters.
  inline int* procLim(){return _modelDef->procLim;}
  inline int nx(){return _modelDef->NX;}
  inline int ny(){return _modelDef->NY;}
  inline int nz(){return _modelDef->NZ;}
  inline int nxy(){return _modelDef->NXY;}
  inline int nxyz(){return _modelDef->NXYZ;}

  ///Include functions for the differentiator coefficients.
  // Centered (v o v x v o v).
  virtual float* px(){
    assert(FALSE,"sgfdModel::px--variable not defined");
    return NULL;}
  virtual float* py(){
    assert(FALSE,"sgfdModel::py--variable not defined");
    return NULL;}
  virtual float* pz(){
    assert(FALSE,"sgfdModel::pz--variable not defined");
    return NULL;}

  // Non-dimensionalizing diff coeffs.
  virtual float sx(){
    assert(FALSE,"sgfdModel::sx--variable not defined");
    return 0.0;}
  virtual float sy(){
    assert(FALSE,"sgfdModel::sy--variable not defined");
    return 0.0;}
  virtual float sz(){
    assert(FALSE,"sgfdModel::sz--variable not defined");
    return 0.0;}

  // Non-centered forms (v o v o x o v o v) for halved variables that are
  //  averaged.
  virtual float* rx(){
    assert(FALSE,"sgfdModel::rx--variable not defined");
    return NULL;}
  virtual float* ry(){
    assert(FALSE,"sgfdModel::ry--variable not defined");
    return NULL;}
  virtual float* rz(){
    assert(FALSE,"sgfdModel::rz--variable not defined");
    return NULL;}

  // Non-centered forms (v o v o x o v o v) for halved variables that are
  //  not averaged.
  virtual float* qx(){
    assert(FALSE,"sgfdModel::rx--variable not defined");
    return NULL;}
  virtual float* qy(){
    assert(FALSE,"sgfdModel::ry--variable not defined");
    return NULL;}
  virtual float* qz(){
    assert(FALSE,"sgfdModel::rz--variable not defined");
    return NULL;}

  // Another non-centered forms.
  virtual float* fx(){
    assert(FALSE,"sgfdModel::fx--variable not defined");
    return NULL;}
  virtual float* fy(){
    assert(FALSE,"sgfdModel::fy--variable not defined");
    return NULL;}

  // Gradient tensor differentiators.
  virtual float* g(){
    assert(FALSE,"sgfdModel::g--variable not defined");
    return NULL;}
  // 3 components of the 16 point 4th order interpolator.
  virtual float* a(){
    assert(FALSE,"sgfdModel::a--variable not defined");
    return NULL;}

  ///Include functions to access the model parameters.
  virtual float* rho(){
    assert(FALSE,"sgfdModel::rho--variable not defined");
    return NULL;
  }
  virtual unsigned char* vxfunc(){
    assert(FALSE,"sgfdModel::vxfunc--variable not defined");
    return NULL;
  }
  virtual unsigned char* vyfunc(){
    assert(FALSE,"sgfdModel::vyfunc--variable not defined");
    return NULL;
  }
  virtual unsigned char* vzfunc(){
    assert(FALSE,"sgfdModel::vzfunc--variable not defined");
    return NULL;
  }
  virtual unsigned char* ssfunc(){
    assert(FALSE,"sgfdModel::ssfunc--variable not defined");
    return NULL;
  }
  virtual unsigned char* xzfunc(){
    assert(FALSE,"sgfdModel::xzfunc--variable not defined");
    return NULL;
  }
  virtual unsigned char* yzfunc(){
    assert(FALSE,"sgfdModel::yzfunc--variable not defined");
    return NULL;
  }
  virtual unsigned char* xyfunc(){
    assert(FALSE,"sgfdModel::xyfunc--variable not defined");
    return NULL;
  }
  virtual unsigned char* muZero(){
    assert(FALSE,"sgfdModel::muZero--variable not defined");
    return NULL;
  }
  virtual bool useO24Only(){
    assert(FALSE,"sgfdModel::useO24Only--varibale not defined");
    return true;
  }
  virtual float* mu(){
    assert(FALSE,"sgfdModel::mu--variable not defined");
    return NULL;
  }
  virtual float* lambda(){
    assert(FALSE,"sgfdModel::lambda--variable not defined");
    return NULL;
  }
  virtual float* bulk(){
    assert(FALSE,"sgfdModel::bulk--variable not defined");
    return NULL;
  }
  virtual float* mmVx(){
    assert(FALSE,"sgfdModel::mmVx--variable not defined");
    return NULL;
  }
  virtual float* mmVy(){
    assert(FALSE,"sgfdModel::mmVy--variable not defined");
    return NULL;
  }
  virtual float* mmVz(){
    assert(FALSE,"sgfdModel::mmVz--variable not defined");
    return NULL;
  }


  ///Here is a virtual that can be redefined for extra model paramters. For
  /// instance, the poroelastic problems has about 8 extra. Instead of putting
  /// virtuals here just redefine this method to return the correct one based on
  /// the value of i.
  virtual float* extraVar(int i){
    assert(FALSE,"sgfdModel::extraVar--none defined");
    return NULL;
  }
    
  ///Here are some additional functions for derived model parameters at a given
  /// index.
  virtual float vp(int i,int j,int k){
    return vp(i+j*_modelDef->NX+k*_modelDef->NX*_modelDef->NY);
  }
  virtual float vs(int i,int j,int k){
    return vs(i+j*_modelDef->NX+k*_modelDef->NX*_modelDef->NY);
  }
  virtual float rho(int i,int j,int k){
    return rho(i+j*_modelDef->NX+k*_modelDef->NX*_modelDef->NY);
  }
  virtual float vp(int i){
    assert(FALSE,"sgfdModel::vp--variable not defined");
    return 0.0;
  }
  virtual float vs(int i){
    assert(FALSE,"sgfdModel::vs--variable not defined");
    return 0.0;
  }
  virtual float rho(int i){
    assert(FALSE,"sgfdModel::rho--variable not defined");
    return 0.0;
  }

  ///Here is a method to generate a unique filename based on the size of this
  /// process.
  char* uniqueFilename(const char* head,const char* dir=NULL);
protected:
  ///Here is an experimental method to increase the grid density by a specified
  /// integer multiplier. The user is responsible for changing the time vector to
  /// an appropriate value with the -T flag.
  virtual int multiplyGrid(int surfaceBCType,int gridMultiplier);

  //Subroutine HOLBERG_COEFFS calculates the two coefficients of the
  // 4th-order, staggered grid, finite-difference operator that approx-
  // imates a first derivative.  The formulae in Kindelan et al. 
  // (Geophysics, v. 55, p. 107-110) are used.  If the input quantity
  // rel_error=0, then these equal the conventional 4th-order staggered
  // grid numerical differentiator coefficients (c1=9/8 and c2=-1/24). 
  void holbergCoeffsC(float rel_error,
                      float* c1,float* c2);
  //calcHolbergCoef scales the coeficients for the x,y, and z dimensions
  // based on dx,dy,dz, and dt
  int calcHolbergCoef(float dx,float dy,float dz,float dt,
		      float& sx,float& sy,float& sz,
		      float cx[2],float cy[2],float cz[2],
                      float scalarSpeed,int fdCoeffOrder);

  ///Read material properties from a cdf file.
  int readCDFModel(char* fileName,
		   char* vsFileName,char* rhoFileName,
		   int* lim,
		   float* vp,float* vs,float* rho,
                   float& vMin,float& vMax);

  int minMaxVelocity(float* vp,float* vs,float *rho,float &vMin,float &vMax);

  
  ///And here is the even more complicated and longer method to read an interpolated
  /// grid.
  float* readCDFModelInterp(char* fileName,
			    int* lim,
			    const char* varName,float* varValues,
                            int noFindError=TRUE);
};

#define CHECKPOINT_START_FILE_FLAG 11234
///This class encapsulates all features required for dependent variables.

///This is pure virtual class until advance, advanceVel, and
///advanceStress methods are defined.
class sgfdDependent{
public:
  float* _dummy;
  sgfdModel* _model;

  int _checkpointID,_checkpointIteration;
  char _checkpointDir[1024];

  sgfdDependent(){
    _model=NULL;

    _checkpointID=-1;
    _checkpointIteration=0;
    *_checkpointDir='\0';
  }
  sgfdDependent(sgfdModel* model){
    setModel(model);

    _checkpointID=-1;
    _checkpointIteration=0;
    *_checkpointDir='\0';
  }

  ///Here is a method to let us know if the master has actual data or is just
  /// performing control.  
  virtual int doIHaveData(){return FALSE;}

  ///Advance the finite-difference algorithm count steps. Pure virtual here.
  virtual int advance(int& iteration,int count,int acknowledge=TRUE)=0;
  ///Advance the velocity one step. Pure virtual here.
  virtual int advanceVel(int iteration,int acknowledge=FALSE)=0;
  ///Advance the stress one step. Pure virtual here.
  virtual int advanceStress(int iteration,int acknowledge=FALSE)=0;

  ///Set the internal variable that holds the model to a new value.
  sgfdModel* setModel(sgfdModel* model){
    return _model=model;
  }

  ///Return the internal structure from the model that defines the space.
  modelDefStruct* modelDef(){return _model->modelDef();};
  ///Return the internal structure from the model that defines the parallel
  /// decomposition.
  parallelDefStruct* parallelDef(){return _model->parallelDef();}

  ///Return the internal variable that holds the model
  sgfdModel* model(){return _model;}
  ///Return the internal variable that is set when the model is read.
  char* modelName(){return _model->modelName();}
  virtual int isAnelastic(){return FALSE;}

  ///Here is a method that will look at a checkpoint directory and determine if
  /// it contains a valid checkpoint file.
  virtual int checkCheckpoint(char* cpDir,int cpID=0);
  
  ///Here is a method that will write a set ofcheckpoint (restart)
  /// files to the given directory.
  virtual FILE* doCheckpoint(int iteration,char* cpDir,int cpID,
			     int callDepth,
                             FILE* cpFile=NULL);
  ///Here is a method to read a checkpoint file written by doCheckpoint.
  virtual FILE* readCheckpoint(int& iteration,char* cpDir,int cpID,
			       int callDepth,
                               FILE* cpFile=NULL);

  //Sometimes need to use different start and stop indicies depending on the component.
  // Since these will have to be put into temps for sending to the Fortran subroutines
  // (pass by reference) I am going to store the C array values (0 start index).
  virtual int istart(){return 2;}
  virtual int istop(){return _model->modelDef()->NX-2;}

  virtual int jstart(){return 2;}
  virtual int jstop(){return _model->modelDef()->NY-2;}

  virtual int kstart(){return 2;}
  virtual int kstop(){return _model->modelDef()->NZ-2;}

  virtual int istop_wx(){return modelDef()->NX-3;}
  virtual int jstop_wy(){return modelDef()->NY-3;}
  virtual int kstop_wz(){return modelDef()->NZ-3;}

  //Put in functions to return the possible fields. All return NULL here.
  virtual floatPtr& vx(int i=0){return _dummy=NULL;}
  virtual floatPtr& vy(int i=0){return _dummy=NULL;}
  virtual floatPtr& vz(int i=0){return _dummy=NULL;}

  virtual floatPtr& P(int i=0){return _dummy=NULL;}

  //Here are some subroutines to access and modify the integer arrays used by the
  // Fortran code to determine which of the 4D arrays is the current time step
  // and which is the previous time step.
  virtual int* index_w(){return NULL;}
  virtual int* index_p(){return NULL;}

  //xx is sometimes used as a proxy for P
  virtual floatPtr& xx(int i=0){return _dummy=NULL;}
  virtual floatPtr& yy(int i=0){return _dummy=NULL;}
  virtual floatPtr& zz(int i=0){return _dummy=NULL;}

  virtual floatPtr& xy(int i=0){return _dummy=NULL;}
  virtual floatPtr& xz(int i=0){return _dummy=NULL;}
  virtual floatPtr& yz(int i=0){return _dummy=NULL;}

  ///Here is a virtual that can be redefined for extra model paramters. For
  /// instance, the poroelastic problems has about 8 extra. Instead of putting
  /// virtuals here just redefine this method to return the correct one based on
  /// the value of i.
  virtual float* extraVar(int i){
    assert(FALSE,"sgfdDependent::extraVar--none defined");
    return NULL;
  }

protected:
  ///Method to fill in a single variable slice on the YZ plane. This is 
  /// appropriate for values that can be developed by multiplication of 1 
  /// variable by a scalar.
  void doYZSlice(float* var,int yStart,int yEnd,int zStart,int zEnd,
		 float* slice,float scalar,
                 int i,float a1,float a2,float a3,float a4);

  ///Method to fill in a single variable slice on the XZ plane. This is 
  /// appropriate for values that can be developed by multiplication of 1 
  /// variable by a scalar.
  void doXZSlice(float* var,int xStart,int xEnd,int zStart,int zEnd,
		 float* slice,float scalar,
                 int j,float b1,float b2,float b3,float b4);
  ///Method to fill in a single variable slice on the YZ plane. This is 
  /// appropriate for variables that can only be derived by some mathematical
  /// combination of values. The function pointer valFunc should return the desired
  /// value at the index.
  void doXZSlice(float (*valFunc)(int,int,int),
		 int xStart,int xEnd,int zStart,int zEnd,
		 float* slice,
                 int j,float b1,float b2,float b3,float b4);

  ///Method to fill in a single variable slice on the XY plane. This is 
  /// appropriate for values that can be developed by multiplication of 1 
  /// variable by a scalar.
  void doXYSlice(float* var,int xStart,int xEnd,int yStart,int yEnd,
		 float* slice,float scalar,
                 int k,float c1,float c2,float c3,float c4);
};

#define BASIC_MODEL 0 //Numeric code for this type of model.
///Class for the master model.

///This class does not actually store anything but it sends messages
///to the slaves to do initialization.
class masterSgfdModel:public sgfdModel{
public:
  parallelDefStruct* _parallelDef;
  int _maxSlaveNX, _maxSlaveNY, _maxSlaveNZ;

  ///This initializer does not send the message, used by derived classes.
  masterSgfdModel(modelDefStruct* modelDef,parallelDefStruct* parallelDef,
		  char* modelName=NULL,
		  int surfaceBCType=0,
		  int gridMultiplier=1)
    :sgfdModel(modelDef,modelName){
    _parallelDef=parallelDef;
    _maxSlaveNX=_maxSlaveNY=_maxSlaveNZ=0;
    if(gridMultiplier!=1)
      multiplyGrid(surfaceBCType,gridMultiplier);
  }

  /// \brief This initializer actually sets up the subdomains.
  /// Determine breakdown and send out messages to get thing started
  /// on the slave side.
  masterSgfdModel(int tids[],
		  modelDefStruct* modelDef,parallelDefStruct* parallelDef,
		  char* modelName,int fdCoeffOrder,float* fdCoeffs,
		  int longRead=FALSE)
    :sgfdModel(modelDef,modelName){
    _parallelDef=parallelDef;

    _maxSlaveNX=_maxSlaveNY=_maxSlaveNZ=0;
    initializeSlaves(tids,
		     modelName,fdCoeffOrder,fdCoeffs,
		     longRead);
  }

  /// \brief This initializer actually sets up the subdomains.
  /// This is an experimental version that multiplies the grid by a sepecified
  ///  integer to allow propagation of higher frequency source without re-writing
  ///  the model.
  masterSgfdModel(int surfaceBCType,int gridMultiplier,
		  int tids[],
		  modelDefStruct* modelDef,parallelDefStruct* parallelDef,
		  char* modelName,int fdCoeffOrder,float* fdCoeffs,
		  int longRead=FALSE)
    :sgfdModel(modelDef,modelName){
    _parallelDef=parallelDef;
    _maxSlaveNX=_maxSlaveNY=_maxSlaveNZ=0;

    if(gridMultiplier!=1)
      multiplyGrid(surfaceBCType,gridMultiplier);
    initializeSlaves(tids,
		     modelName,fdCoeffOrder,fdCoeffs,
		     longRead);
  }

  virtual ~masterSgfdModel() {}
  virtual parallelDefStruct* parallelDef(){return _parallelDef;}

  virtual int resetTimeVector(float minT,int nt,float dt);
protected:
  virtual int multiplyGrid(int surfaceBCType,int gridMultiplier);

  ///This virtual function is called at the end of initializeSlaves to receive an
  /// acknoledgement that the slaves have completed initialization. This can be 
  /// redefined to give different behavior (maybe reporting more information on the
  /// material properties). The minimum information sent back is the min and max
  /// velocity and min and max density (since all codes using these modules save those).
  virtual void getInitAck(int tids[]);
    
  ///This procedure does almost all the work of the initialization. It calls the virtual
  /// sendInitMessage which should be re-defined in child classes that need to send additional
  /// information.
  virtual int initializeSlaves(int tids[],
			       char* modelName,int fdCoeffOrder,float* fdCoeffs,
                               int longRead);

  ///Here is a local function to build the model def portion of the domain
  /// decomposition.
  modelDefStruct** buildSlaveModelDef();

  ///Determine the neighbors a given subdomain.
  int setNeighbors(parallelDefStruct* parallelDef,
		   int tids[],modelDefStruct* slaveModelDefs[],
		   int nxProc,int nyProc,int nzProc,
                   int i,int j,int k);

  ///This small subroutine actually sends the initialization message. This is a virtual function
  /// so it can be changed, for instance a viscoElasticModel sends additional information about
  /// the definition of the anelastic parameters.
  virtual void sendInitMessage(char* modelName,int fdCoeffOrder,float* fdCoeffs,int longRead,
			       parallelDefStruct** slaveParallelDefs,
			       int tids[],int i,int j,int k,
			       int nxProc,int nxyProc,
                               int coworkParadyme=FALSE);
};

/// \brief Definition of a model for slaves.
/// This still does not actually define space for the model parameters
/// but includes some methods that are common to all slave models.
class slaveSgfdModel:public sgfdModel{
public:
  parallelDefStruct* _parallelDef;

  //Here are the differentiator coefficients.
  float _cx[2],_cy[2],_cz[2];
  float _sx,_sy,_sz;

  int _doLongRead;

  ///Certian derived classes require the null initializer.
  slaveSgfdModel():sgfdModel(){
    _doLongRead=FALSE;

    _sx=_cx[0]=_cx[1]=0.0;
    _sy=_cy[0]=_cy[1]=0.0;
    _sz=_cz[0]=_cz[1]=0.0;
  }

  ///The model is read from a cdf file. But the 
  /// model initialization routine receives the message
  /// from the parent containing the name (as well as this model limits)
  slaveSgfdModel(int convertParams);

  virtual ~slaveSgfdModel(){
    free(_modelDef);
    free(_parallelDef);
  }

  virtual int resetTimeVector();

  ///Include functions for the differentiator coefficients.
  virtual float* px(){return _cx;}
  virtual float* py(){return _cy;}
  virtual float* pz(){return _cz;}

  virtual float sx(){return _sx;}
  virtual float sy(){return _sy;}
  virtual float sz(){return _sz;}

  //Return the most important parameters
  inline parallelDefStruct* parallelDef(){return _parallelDef;}

  //Determine if a point is inside my domain.
  int isInside(float x,float y,float z,
               int &i,int &j,int &k);
  /*! \brief Determine if a point is inside the given domain. 

    If the domain is on edge of the model then points outside that 
    edge are considered inside. Use this with interp3.
  */ 
  int isInside3(float x,float y,float z);
  /*! \brief Determine an interpolated value of a given variable at an
    arbitrary location.

    If the domain is on edge of the model then points outside that 
    edge are considered inside. Note if the point is strictly inside
    the domain this just returns interp2(var,x,y,z); if the point is
    outside the domain it returns interp2 of a point on the domain edge.
    Use this with isInside3.
  */ 
  float interp3(float *var,
                float x,float y,float z);
protected:
  //Perform 2nd order 3D interpolation on some field.
  float interp2(float *var,
                float x,float y,float z);
};

/// A master dependent class, this is often sufficient.
class masterSgfdDependent:public sgfdDependent{
public:
  masterSgfdModel* _model;

  sourceNetwork* _sources;

#ifdef SGFD_INVERSION
  masterInversionReceiverNetwork* _receivers;
#else
  masterReceiverNetwork* _receivers;
#endif

  ///Constructor that sends the required messages to set up the slaveDependents.
  masterSgfdDependent(masterSgfdModel* model,int doSend=TRUE):sgfdDependent(model){
    _model=model;
    _receivers=NULL;
    _sources=NULL;

    if(doSend)
      sendInitialMessage();
    tEprintf(Verbose,"Initialized dependent variables\n");
  }

  virtual ~masterSgfdDependent();

  ///Write a checkpoint file.
  virtual FILE* doCheckpoint(int iteration,char* cpDir,int cpID,
			     int callDepth,
                             FILE* cpFile=NULL);
  ///Read the checkpoint file.
  virtual FILE* readCheckpoint(int& iteration,char* cpDir,int cpID,
			       int callDepth,
                               FILE* cpFile=NULL);
  
  /// \brief Set initial conditions.
  /// This has received only limited testing.
  virtual int setInitialConditions(char* initialConditionsName,int& iteration,
                                   int negateVelocities);

  /// \brief Set the boundary conditions.
  /// This is where a free/ridgid surface is set.
  virtual int setBoundaryConditions(int surfaceBCMode,
                                    int spongeBCNodes[6],float spongeBCValue[6]);
 
  /// \brief Set the reciever geometry.
  /// This method looks in the model file for receivers and then adds
  /// any receivers that where defined on the command line. Send this
  /// off to the slaves.
  virtual int setReceivers(masterReceiverNetwork* extraReceivers,
			   int surfaceMode,
                           int writeTraceFile=TRUE,int doCheckActive=TRUE);
  parallelDefStruct* parallelDef(){return _model->parallelDef();}

  /// \brief Set the source geometry.
  /// This method looks in the model file for sources and then adds
  /// any sources that where defined on the command line. Send this
  /// off to the slaves.
  virtual int setSources(sourceNetwork* extraSources);
  
  /// \brief Write trace data.
  /// Can be called with iteration < NT to write part way through a
  /// run.
  virtual int writeTraces(int startI,int stopI);

  ///\brief Advance the run by count iterations.
  /// Note that this just sends a message, the slaves do the actual work.
  virtual int advance(int& iteration,int count,int acknowledge=TRUE);

  ///Advance the velocity by 1/2 iteration.
  virtual int advanceVel(int iteration,int acknowledge=FALSE);

  ///Advance the stress by 1/2 iteration.
  virtual int advanceStress(int iteration,int acknowledge=FALSE);

  //IO Procs
  ///Write a time slice.
  virtual int slicer(int plane,int comp,float coord,
		     floatPtr& outData,int currProcLim[6]){
    return
      assert(FALSE,"masterSgfdDependent::slicer--no data in this class");
  }

protected:
  ///Here is an internal method to get the acknowledgment that the advance is
  /// complete. This default version just implements a barrier to make sure all
  /// the processes are synced. A more advanced version might get some additional
  /// information from the slaves.
  virtual void getAcknowledgement();

  int sendInitialMessage();

  ///Here is a method to allocate the receiver network. Since this is virtual it
  /// can be overridden in derived classes that want to use a different type of 
  /// network.
  virtual masterReceiverNetwork* allocateReceivers(modelDefStruct *modelDef,
						   char* modelname,int surfaceMode,
						   int receiverDecimate,
                                                   int useModelReceivers);
};

/// \brief A slave dependent class, this is still virtual.
/// Methods common to all slave processes. Still has not defined the
/// advance methods. Those depend on the specific type of modeling.
class slaveSgfdDependent:public sgfdDependent{
public:
  slaveSgfdModel* _model;

  sourceNetwork* _sources;

#ifdef SGFD_INVERSION
  slaveInversionReceiverNetwork* _receivers;
#else
  slaveReceiverNetwork* _receivers;
#endif

  //Keep local space for the slice around instead of allocating/deallocating
  // for every slice. I have found this makes for more efficient memory use
  // on many platforms.
  float *_slice;
  int _sliceLength;

  ///Need a constructor that does very little for certian subclasses
  slaveSgfdDependent(slaveSgfdModel* model)
    :sgfdDependent(model){
    _model=model;

    _sources=NULL;
    _receivers=NULL;

    _slice=NULL;
    _sliceLength=0;
  }

  ///Delete all allocated variables
  virtual ~slaveSgfdDependent();
  
  virtual int doIHaveData(){return TRUE;}

  ///Write a checkpoint file.
  virtual FILE* doCheckpoint(int iteration,char* cpDir,int cpID,
			     int callDepth,
                             FILE* cpFile=NULL);
  ///Read the checkpoint file.
  virtual FILE* readCheckpoint(int& iteration,char* cpDir,int cpID,
			       int callDepth,
                               FILE* cpFile=NULL);
  
  ///Set the sources from the message from the master.
  virtual int setSources(float* rho=NULL,int useCubicExtrap=FALSE);

  ///Set the receivers from the message from the master.
  virtual int setReceivers(int useCubicExtrap=TRUE);

  ///\brief Advance by count iterations.
  ///This is just a loop calling the virtual methods (still undefined)
  ///advanceVel and advanceStress. 
  virtual int advance(int& iteration,int count,int acknowledge=TRUE);

  //
  // IO subroutines
  //

  ///Pack the data from the receivers for sending to the master.
  int packReceiverData(int index,int startI,int stopI);

  ///Send the data from the receiver grids for sending to the master.
  int sendReceiverGridData(int index);

  /*! \brief Send a complete volume to the requestor.

  This is another method that must be defined in any real subclass.
   */
  virtual int fullFieldOutput(int recepient,
			      char* msgtag=NULL)=0;

  /*! \brief A Helper version of fullFieldOutput.

    Pass an argument to this version and it will send the data one layer at
    a time to the target.
  */
  virtual int fullFieldOutput(int recepient,float* vx);
    
  /// \brief Create a block of data containing a slice and send to the requestor.
  /// In serial mode this is written directly, in parallel mode this
  /// is sent as message to the parent. This will probably just call
  /// one of the other slicer methods with the correct arguments.
  virtual int slicer(int recepient)=0;

  /// Version 1 does everything including read the message).
  virtual int slicer(int recepient,
		     float *vx,float *vy,float *vz,
		     float *xx,float *yy,float *zz,
                     float *xy,float *xz,float *yz);

  ///Version 2, the message
  // has already been read (perhaps to pick off new types of slices.
  virtual int slicer(int recepient,
		     int plane,int comp,float coord,
		     float *vx,float *vy,float *vz,
		     float *xx,float *yy,float *zz,
                     float *xy,float *xz,float *yz);

protected:
  ///Here is a method to allocate the receiver network. Since this is virtual it
  /// can be overridden in derived classes that want to use a different type of 
  /// network.
  virtual slaveReceiverNetwork* allocateReceivers(modelDefStruct* modelDef,
						  int allocateData,
						  int mesgFromParent,
                                                  int useCubicExtrap);

  ///Simple internal method to map some method that returns vp (defined in the 
  ///model) to method in the dependent, this may make generating slices much easier.
  float modelVp(int i,int j,int k){return _model->vp(i,j,k);}

  ///Here is an internal method to send the acknowledgment that the advance is
  /// complete. This default version just implements a barrier to make sure all
  /// the processes are synced. A more advanced version might send some additional
  /// information to the master.
  virtual void sendAcknowledgement();

  void sliceXZCompValue(int comp,int j,
			int sliceXStart,int sliceXEnd,
			int sliceZStart,int sliceZEnd,
			float b1,float b2,float b3,float b4,
			float *vx,float *vy,float *vz,
			float *xx,float *yy,float *zz,
                        float *xy,float *xz,float *yz);

  void sliceYZCompValue(int comp,int i,
			int sliceYStart,int sliceYEnd,
			int sliceZStart,int sliceZEnd,
			float a1,float a2,float a3,float a4,
			float *vx,float *vy,float *vz,
			float *xx,float *yy,float *zz,
                        float *xy,float *xz,float *yz);

  void sliceXYCompValue(int comp,int k,
			int sliceXStart,int sliceXEnd,
			int sliceYStart,int sliceYEnd,
			float c1,float c2,float c3,float c4,
			float *vx,float *vy,float *vz,
			float *xx,float *yy,float *zz,
                        float *xy,float *xz,float *yz);
};

//And as promised in sgfdCommandParser here is the final version of the parser.
///The real sgfdCommandParser class with the extraOuput stuff included.
class sgfdCommandParser: public sgfdCommandParserBase{
public:
  //The final local variable.
  extraOutputGroup* _extraOutput; //add with -E[?] flags.

  ///The constructor does nothing extra from the parent. 
  sgfdCommandParser(int argc,char* argv[]):sgfdCommandParserBase(argc,argv){}

  ///Allocate an extraOutputGroup, this is a virtual so derived classes can
  /// create specialized versions if needed.
  virtual extraOutputGroup* allocateExtraOutput();

  ///Here is the required virtual function to process a flag.
  virtual int checkProcessFlag(int argc,char* argv[],int& i);

protected:
  ///Set any extra default values.
  virtual void setDefaultValues(){
    sgfdCommandParserBase::setDefaultValues();
    addExtraOutputFlags();
    _extraOutput=NULL; //add with -E[?] flags.
  }

  ///Add a description of flags pertaining to extra output.
  virtual void addExtraOutputFlags();
};


#endif //#ifndef _sgfd_hh_
