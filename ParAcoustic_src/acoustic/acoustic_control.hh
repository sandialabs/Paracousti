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
 *  acoustic_control.hh
 *
 *
 *  Declare slave model and dependent classes.  These classes live
 *  on specific domains and handle the actual model specific to that
 *  domain as well as the dependent variables for that domain.  The
 *  dependent class handles the time step updating for the dependent
 *  variables.
 *
 *  Declares classes: 
 *    slaveAcousticModel
 *    slaveAcousticReplacementModel
 *    slaveAcousticPseudoAttenuateModel
 *    slaveAcousticDependent
 *    slaveAcousticReplacementDependent
 *    slaveAcousticPseudoAttenuateDependent
 *
 */

#ifndef _acoustic_control_hh_
#define _acoustic_control_hh_

//Here are defines to determine the type of model. I now only need to know about still media.
#define FIXED_MEDIA         1
#define FIXED_MEDIA_ATTEN  2

class selector;
class acousticSgfdBoundaries;

#include "sgfd.hh"
#include "xtrautil.hh"

//Global variables and include file for Vampir parallel profiling API.
#if USE_VAMPIR
#if USE_WARNING
#warning "Using Vampir User API, Cplant only"
#endif
#include "/usr/local/share/vt/include/VT.h"
extern int _VTUpdateClass_,_VTBCClass_;
extern int _VTVelUpdate_,_VTStressUpdate_;
extern int _VTBCSaveVel_,_VTBCDoVel_,_VTBCDoStress_;
#endif

/// \brief Here is a real slave model for the stationary acoustic problem.
///The slave model needs to store bulk modulus instead of lambda and mu.
class slaveAcousticModel: public slaveSgfdModel{
public:
  float *_bulk,*_rho;
  float rhoRatioLimit, compLimit;
  unsigned char* _vxfunc, *_vyfunc, *_vzfunc, *_ssfunc;

  //Define some additional variables that correspond to dfa's formula notation
  float _rx[2],_ry[2],_rz[2],_fx[2],_fy[2];

  //The model is read from a cdf file. But the 
  // model initialization routine receives the message
  // from the parent containing the name (as well as this model limits)
  //The default rhoRatLim and cLimit control activation of order switching
  //for high contrast interfaces in the model
  slaveAcousticModel(float rhoRatLim=4.0, float cLimit=20.,bool finalizeInit=true);

  //I also need another initializer that does nothing at all for some derived classes. Use
  // a dummy int argument.
  slaveAcousticModel(int dummy):slaveSgfdModel(){
    _rho=_bulk=NULL;
  }

  //destructor
  ~slaveAcousticModel(){
    if(_bulk)
      free(_bulk);
    if(_rho)
      free(_rho);
    if(_vxfunc)
      free(_vxfunc);
    if(_vyfunc)
      free(_vyfunc);
    if(_vzfunc)
      free(_vzfunc);
    if(_ssfunc)
      free(_ssfunc);
    _vxfunc=_vyfunc=_vzfunc=_ssfunc=NULL;
  }
  //helper routines to retrieve model parameters
  virtual float vp(int i);
  virtual float rho(int i);

  virtual unsigned char* vxfunc() {
    return _vxfunc;
  }
  
  virtual unsigned char* vyfunc() {
    return _vyfunc;
  }
  
  virtual unsigned char* vzfunc() {
    return _vzfunc;
  }
  
  virtual unsigned char* ssfunc() {
    return _ssfunc;
  }
  
  virtual void convertVels(float &vMin, float &vMax, float &rhoMin, float &rhoMax,
                           float &bMMin, float &bMMax);
  virtual int resetTimeVector();
  virtual float* rx(){return _rx;}
  virtual float* ry(){return _ry;}
  virtual float* rz(){return _rz;}

  virtual float* fx(){return _fx;}
  virtual float* fy(){return _fy;}

  ///Include functions to access the model parameters.
  virtual float* rho(){return _rho;}
  virtual float* bulk(){return _bulk;}

  ///Here is a function to return the type of media being used (still here).
  virtual int mmType(){return FIXED_MEDIA;}

  //Some other minor changes because of the different material property definitions.
  //Return the min and maximum velocity in the model
  // this requires a conversion from lambda, mu, boy to vel
  // also calculates min and max den, vp and vs but only returns
  // min zero vel and max vel
  void minMax(float& min,float& max){
    float rhoMin,rhoMax;
    minMax(min,max,rhoMin,rhoMax);
  }
  void minMax(float& vpMin,float& vpMax,
	      float& rhoMin,float& rhoMax){
    DEF_MODEL_SIZE(modelDef());
    DEF_MODEL_SCALE(modelDef());

    for(int i=0;i<NXYZ;i++){
      float rho=_rho[i];

      float vp=sqrt(_bulk[i]/rho);

      rho*=scalarDen;
      vp*=scalarSpeed;

      SETMINMAX(vpMin,vpMax,vp,i);
      SETMINMAX(rhoMin,rhoMax,rho,i);
    }
  }
};
typedef slaveAcousticModel* slaveAcousticModelPtr;

/* This class adds standard linear fluid attenuation to the acoustic model.  It
  adds a few fields that must be specified in the input netcdf model file.
*/
class slaveAcousticAttenModel: public slaveAcousticModel {
public:
  //There is an array of all the selectors.
  selector* _selectors;
  
  //And an array with 1 pointer for each node aimed at the Q for that node.
  int* _Q;
  int* _nMechs;
  int _nM;
  float **_decayRates,**_ampP;
  float *_omegaAmpSum;
  
  //Keep one variable for the number of relaxation mechanisms, this is actually
  // the maximum number in any of the selectors.
  int _numRelaxation;
  
  //most of the machinery for this class can be handled by the superclass, but
  //it defines a few new variables that must be read from the input model file.
  slaveAcousticAttenModel(float rhoRatLim=4.0, float cLimit=20.);

  virtual ~slaveAcousticAttenModel();

  void completeFullSetup();
  
  int* readCDFQindex(char* fileName,
                     int* lim,
                     int* qindex) {
    DEF_PARALLEL(parallelDef());
    return readModelVariable(qindex,fileName,"Qindex",lim,
                             globalNX,globalNY,globalNZ,-1,FALSE);
  }
};

class slaveAcousticDependent:public slaveSgfdDependent{
public:
  slaveAcousticModel* _model; //Overrule the definition in slaveSgfdDependent.
  acousticSgfdBoundaries* _boundaries;

  //Need to change the starting point in the z-direction for special BC's.
  int _kstart;

  //The dependent variables.
  float *_vx,*_vy,*_vz;
  float* _pressure;
  
  ///Here is the constructor.
  slaveAcousticDependent(slaveAcousticModel* model,
                         int doAllocate=TRUE);
  //destructor
  virtual ~slaveAcousticDependent(){
    if(_vx){
      free(_vx);
      _vx=NULL;
    }
    if(_vy){
      free(_vy);
      _vy=NULL;
    }
    if(_vz){
      free(_vz);
      _vz=NULL;
    }
    if(_pressure){
      free(_pressure);
      _pressure=NULL;
    }
  }

  ///Here is a method that will write a set ofcheckpoint (restart)
  /// files to the given directory.
  virtual FILE* doCheckpoint(int iteration,char* cpDir,int cpID,
			     int callDepth,
                             FILE* cpFile=NULL);
  ///Here is a method to read a checkpoint file written by doCheckpoint.
  virtual FILE* readCheckpoint(int& iteration,char* cpDir,int cpID,
			       int callDepth,
                               FILE* cpFile=NULL);

  //Put in functions to return the possible fields. Only 1 time plane available.
  virtual floatPtr& vx(int i=0){
    assert(i==0,
	   "slaveAcousticDependent::vx--attempt access plane (%i) above 0",i);
    return _vx;
  }
  virtual floatPtr& vy(int i=0){
    assert(i==0,
	   "slaveAcousticDependent::vy--attempt access plane (%i) above 0",i);
    return _vy;
  }
  virtual floatPtr& vz(int i=0){
    assert(i==0,
	   "slaveAcousticDependent::vz--attempt access plane (%i) above 0",i);
    return _vz;
  }

  virtual floatPtr& P(int i=0){
    assert(i==0,
	   "slaveAcousticDependent::P--attempt access plane (%i) above 0",i);
    return _pressure;
  }
  virtual floatPtr& xx(int i=0){return P(i);}

  //Sometimes need to use different start and stop indicies depending on the component.
  // Since these will have to be put into temps for sending to the Fortran subroutines
  // (pass by reference) I am going to store the C array values (0 start index).
  virtual int kstart(){return _kstart;}

  void setInitialConditions();
  
  virtual void setBoundaryConditions(acousticSgfdBoundaries* boundaries);
  
  virtual int loadSeismogram(int iteration);

  //override super's advance to chech we are stable after every count iterations
  virtual int advance(int& iteration,int count,int acknowledge=TRUE){
    slaveSgfdDependent::advance(iteration,count,acknowledge);
    DEF_MODEL_SIZE(modelDef());
    for(int i=0;i<NXYZ;i++) {
      assert(isfinite(_vx[i]),"vx not finite at iteration %d, point %d\n",iteration,i);
      assert(isfinite(_vy[i]),"vy not finite at iteration %d, point %d\n",iteration,i);
      assert(isfinite(_vz[i]),"vz not finite at iteration %d, point %d\n",iteration,i);
    }
    
    return iteration;
  }
  // Advance the velocities one time step
  virtual int advanceVel(int iteration,int acknowledge);
  // Perform and update to the stress grid
  virtual int advanceStress(int iteration,int acknowledge);

  ///Define an easier to use version of slicer (required for a real class).
  virtual int slicer(int receipient){
    return 
      slaveSgfdDependent::slicer(receipient,
				 _vx,_vy,_vz,
				 _pressure,NULL,NULL,
				 NULL,NULL,NULL);
  }

  ///Define a real version of fullFieldOutput (required for a real class).
  virtual int fullFieldOutput(int recepient,
                              char* msgtag=NULL);
  
protected:
  //default initialization is zero for all dependent variables
  int initVariables(){
    DEF_MODEL_SIZE(modelDef());
    for(int i=0;i<NXYZ;i++){
      _vx[i]=_vy[i]=_vz[i]=0.0;
      _pressure[i]=0.0;
    }
    return NXYZ;
  }

};
typedef slaveAcousticDependent* slaveAcousticDependentPtr;

//This class overrides slaveAcousticDependent to set up and run attenuative models
//using standard linear fluid mechanisms
class slaveAcousticAttenDependent:public slaveAcousticDependent {
public:
  float **_rP;
  slaveAcousticAttenModel *_model;  //Override
  
  slaveAcousticAttenDependent(slaveAcousticAttenModel* model,int doAllocate=TRUE);
  
  ~slaveAcousticAttenDependent() {
    for(int i=0;i<_model->_numRelaxation;++i)
      free(_rP[i]);
    free(_rP);
  }
  
  //methods to advance one time step
  int advanceVel(int iteration,int acknowledge=FALSE);
  int advanceStress(int iteration,int acknowledge=FALSE);
};
#endif //#ifndef _acoustic_control_hh_
