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
 *  sgfdSources.hh
 *
 *
 *  Declares classes used to manage and represent sources.  These classes both handle
 *  individual sources and also the sourceNetwork class which acts a single point of contact 
 *  to handle all sources as a group.
 *
 *  Declares the following classes:
 *  source
 *  tdbcSource
 *  pointSource
 *  momentSource
 *  forceSource
 *  tractionSource
 *  sourceNetwork
 *
 *  Defines the following macros:
 *  DEFAULT_CHECK_BANDWIDTH_PCENT
 *  FORCE_SOURCE
 *  MOMENT_SOURCE
 *  TRACTION_SOURCE
 *  TDBC_SOURCE
 *
 */
//These utility functions are for sources in the SGFD problem

#ifndef _sgfdSources_hh_
#define _sgfdSources_hh_

#include "nstdutil.hh"
#include "array.hh"
#include "sgfd.h"

//Set some constants for default values of various parameters.
#define DEFAULT_CHECK_BANDWIDTH_PCENT 1.0

//
//A pure virtual parent class for all sources, actual sources will be derived from this.
// All sources will return a type.
#define FORCE_SOURCE    1
#define MOMENT_SOURCE   2
#define TRACTION_SOURCE 3
#define TDBC_SOURCE     4

extern const int AllProcesses;

class source{
public:
  int _active;
  float _dummyReturn,*_dummyPtrReturn;

  source(){
    _active=FALSE;
  }
  
  virtual ~source(){}

  //Access functions
  int active(){return _active;}

  //Here are the virtuals that must be defined by each real subclass.
  virtual int type()=0;
  virtual int setActive(modelDefStruct* modelDef)=0;
  virtual void packSourceMessage(modelDefStruct* modelDef,
				 int startIteration)=0;
  virtual void writeSource(int NT,float ampScalar,size_t index,
			   char* filename,int outfile,
			   int xVar,int yVar,int zVar,
			   int ampVar,int dataVar,
			   int* extraVars)=0;

  //These access functions are defined to throw an error.
  virtual float& x(){
    assert(FALSE,"source::x--function not valid");
    return _dummyReturn=0;
  }
  virtual float& y(){
    assert(FALSE,"source::y--function not valid");
    return _dummyReturn=0;
  }
  virtual float& z(){
    assert(FALSE,"source::z--function not valid");
    return _dummyReturn=0;
  }
  virtual float& amp(){
    assert(FALSE,"source::amp--function not valid");
    return _dummyReturn=0;
  }
  virtual float& iso(){
    assert(FALSE,"source::iso--function not valid");
    return _dummyReturn=0;
  }
  virtual floatPtr& data(){
    assert(FALSE,"source::data--function not valid");
    return _dummyPtrReturn=NULL;
  }

  virtual float& ax(){
    assert(FALSE,"source::ax--function not valid");
    return _dummyReturn=0;
  }
  virtual float& ay(){
    assert(FALSE,"source::ay--function not valid");
    return _dummyReturn=0;
  }
  virtual float& az(){
    assert(FALSE,"source::az--function not valid");
    return _dummyReturn=0;
  }

  virtual float& xxS(){
    assert(FALSE,"source::xxS--function not valid");
    return _dummyReturn=0;
  }
  virtual float& yyS(){
    assert(FALSE,"source::yyS--function not valid");
    return _dummyReturn=0;
  }
  virtual float& zzS(){
    assert(FALSE,"source::zzS--function not valid");
    return _dummyReturn=0;
  }

  virtual float& xyS(){
    assert(FALSE,"source::xyS--function not valid");
    return _dummyReturn=0;
  }
  virtual float& xzS(){
    assert(FALSE,"source::xzS--function not valid");
    return _dummyReturn=0;
  }
  virtual float& yzS(){
    assert(FALSE,"source::yzS--function not valid");
    return _dummyReturn=0;
  }

  virtual float& xyA(){
    assert(FALSE,"source::xyA--function not valid");
    return _dummyReturn=0;
  }
  virtual float& xzA(){
    assert(FALSE,"source::xzA--function not valid");
    return _dummyReturn=0;
  }
  virtual float& yzA(){
    assert(FALSE,"source::yzA--function not valid");
    return _dummyReturn=0;
  }

  //Need to apply the source to velocities and stresses.
  virtual int applyVel(modelDefStruct* modelDef,int iteration,
		       float *vx,float *vy,float *vz,
		       float* rho){
    //Default is to do nothing
    return FALSE;
  }
  virtual int applyVelMM(modelDefStruct* modelDef,int iteration,
		       float *vx,float *vy,float *vz,
		       float* rho){
    //Default is to do nothing
    return FALSE;
  }
  virtual int applyStress(modelDefStruct* modelDef,int iteration,
			  float *xx,float* yy,float *zz,
			  float *xy,float* xz,float *yz){
    //Default is to do nothing
    return FALSE;
  }

  //Need to apply the source to acoustic velocities and pressure.
  virtual int applyAcousticVel(modelDefStruct* modelDef,int iteration,
			       float *vx,float *vy,float *vz,
			       float* rho){
    //Default is to do nothing
    return FALSE;
  }
  virtual int applyAcousticVelMM(modelDefStruct* modelDef,int iteration,
			       float *vx,float *vy,float *vz,
			       float* rho){
    //Default is to do nothing
    return FALSE;
  }
  virtual int applyPressure(modelDefStruct* modelDef,int iteration,
			    float *pressure){
    //Default is to do nothing
    return FALSE;
  }
  virtual int applyPressureMM(modelDefStruct* modelDef,int iteration,
			    float *pressure){
    //Default is to do nothing
    return FALSE;
  }
};
typedef source* sourcePtr;
typedef arrayI<sourcePtr> sourceArray;
typedef sourceArray* sourceArrayPtr;

//
//Here are the distributed source types.
class tdbcSource: public source{
public:
  char _infileName[1024];

  //Variables for this source type. The number of points and their indicies.
  int _timeCutoffIndex;
  int _nPoints,*_pointIndicies;
  int *_gridIndicies,*_gridIndiciesMx,*_gridIndiciesMy,*_gridIndiciesMz;
  int _infile,_vxVar,_vyVar,_vzVar,_pressVar; //An open netcdf file

  //May implement as a real TDBC or a pseudo-TDBC were we add the difference
  // on each iteration.
  int _doDifferences;

  //So we don't have to read the entire
  // file at the start or do a read from the file on every iteration.
  int _bufferIndex;
  float** _vxBuffer,**_vyBuffer,**_vzBuffer,**_pressBuffer;

  //Some different variables if the file is already written on the correct
  // locations.
  int _nVxPoints,_nVyPoints,_nVzPoints,_nPressPoints;
  int *_vxPointIndicies,*_vyPointIndicies,*_vzPointIndicies,*_pressPointIndicies;
  int *_vxGridIndicies,*_vyGridIndicies,*_vzGridIndicies,*_pressGridIndicies;

  tdbcSource():source(){
    _doDifferences=FALSE;
    _infile=_vxVar=_vyVar=_vzVar=_pressVar=-1;
    _nPoints=0;
    _pointIndicies=NULL;
    _gridIndicies=_gridIndiciesMx=_gridIndiciesMy=_gridIndiciesMz=NULL;
    _vxBuffer=_vyBuffer=_vzBuffer=_pressBuffer=NULL;
  }
  
  tdbcSource(modelDefStruct* modelDef,char* filename);

  //Read a file corresponding to mapTDBC translation of John Holland's output
  // format.
  tdbcSource(modelDefStruct* modelDef);
    
  virtual ~tdbcSource();

  //Functions unique to a tdbcSource.
  int setDifference(){return _doDifferences=TRUE;}

  //Functions that all forces must define.
  virtual int type(){return TDBC_SOURCE;}
  virtual int setActive(modelDefStruct* modelDef){return TRUE;}
  virtual void packSourceMessage(modelDefStruct* modelDef,
                                 int startIteration);
  virtual void writeSource(int NT,float ampScalar,size_t index,
			   char* filename,int outfile,
			   int xVar,int yVar,int zVar,
			   int ampVar,int dataVar,
			   int* extraVars){}

  //This source type is only applied to velocities, but it can work for either
  // an elastic or an acoustic problem.
  virtual int applyAcousticVel(modelDefStruct* modelDef,int iteration,
		       float *vx,float *vy,float *vz,
		       float* rho){
    return applyVel(modelDef,iteration,
		    vx,vy,vz,
		    rho);
  }
  virtual int applyVel(modelDefStruct* modelDef,int iteration,
		       float *vx,float *vy,float *vz,
                       float* rho);
  
  virtual int applyPressure(modelDefStruct* modelDef,int iteration,
                            float *pressure);

protected:
  void initializeEdgePoints(modelDefStruct* modelDef,const char* flag,
			    int &nPoints,int &var,
			    intPtr &pointIndicies,intPtr &gridIndicies,
                            float*** buffer);

  void initializeUnitPoints(modelDefStruct* modelDef,
                            int numPointsDim,size_t numPoints);

  void readNextBuffer(modelDefStruct* modelDef,int index);
  void readNextEdgeBuffer(modelDefStruct* modelDef,int index);
};

class pointSource: public source{
public:
  float _x,_y,_z; //position variables
  float _amp;

  float* _data; //the source time function

  float _dummyReturn;

  //Constructors

  ///Null constructor.
  pointSource(){
    _x=_y=_z=_amp=0.0;
    _data=NULL;
  }

  ///Parameters passed as arguments.
  pointSource(modelDefStruct* modelDef,
              float x,float y,float z,float amp,float* data);

  ///Parameters read from message.
  pointSource(modelDefStruct* modelDef);

  ///Parameters read from CDF file.
  pointSource(modelDefStruct* modelDef,
	      char* filename,int infile,size_t index,
	      int xVar,int yVar,int zVar,
              int ampVar,int dataVar);

  virtual ~pointSource(){
    if(_data)
      free(_data);
    _data=NULL;
  }

  ///Additional message reading function. If there are a large number of 
  /// time-steps we use an additional function to read the rest of the data.
  virtual int readAuxMessage(int NT,int startIndex,
                             int startTag,int iGetMessg);

  //Access functions
  virtual float& x(){return _x;}
  virtual float& y(){return _y;}
  virtual float& z(){return _z;}
  virtual float& amp(){return _amp;}
  virtual floatPtr& data(){return _data;}

  virtual int setActive(modelDefStruct* modelDef);

  //Pack the basic parameters, this should be the same format as is read by the constructor.
  virtual void packSourceMessage(modelDefStruct* modelDef,
                                 int startIteration);

  //Write the basic parameters to a cdf file, this should be the same format as is read by 
  // the constructor.
  virtual void writeSource(int NT,float ampScalar,size_t index,
			   char* filename,int outfile,
			   int xVar,int yVar,int zVar,
			   int ampVar,int dataVar,
                           int* extraVars);
};

//
//Now define the derived classes that are actually used.
//
///The subclassed pointSource that implements a body moment.
class momentSource: public virtual pointSource{
public:
  //Want to use some shortcuts here is source is all anisotropic or isotropic.
  int _isAnisotropic,_isIsotropic;

  //Need some additional parameters.
  float _xxS,_yyS,_zzS;
  float _xyS,_xzS,_yzS;
  float _xyA,_xzA,_yzA;

  //Additional variable to save some time in the acoustic case.
  float _iso;

  //Needs to be able to extrapolate values to the 3 components of velocity. And the 6 components
  // of velocity, but the diagonal stress components can share an interpolator.
  interpStruct _xy_vInterp,_yz_vInterp,_xz_vInterp;
  interpStruct _xxInterp,_xyInterp,_xzInterp,_yzInterp;

  //Testing of cubic rather than linear extrapolators.
  int _useCubicExtrap, _useMomentRate;
  cubicInterpStruct _cubicXxInterp;

  //Constructors
  ///Parameters passed as arguments.
  momentSource(modelDefStruct* modelDef,
	       float x,float y,float z,
	       float xxS,float yyS,float zzS,
	       float xyS,float xzS,float yzS,
	       float xyA,float xzA,float yzA,
	       float amp,float* data,
               float* rho,int useMomentRate);
  ///Parameters read from message.
  momentSource(modelDefStruct* modelDef,int startTag,float* rho,
	       int useCubicExtrap=FALSE);

  // Parameters read from CDF file.
  momentSource(modelDefStruct* modelDef,
	       char* filename,int infile,size_t index,
	       int xVar,int yVar,int zVar,
	       int xxSVar,int yySVar,int zzSVar,
	       int xySVar,int xzSVar,int yzSVar,
	       int xyAVar,int xzAVar,int yzAVar,
	       int ampVar,int dataVar,
	       float* rho);

  //Now the rest of the functions.
  virtual int type(){return MOMENT_SOURCE;}

  virtual float& xxS(){return _xxS;}
  virtual float& yyS(){return _yyS;}
  virtual float& zzS(){return _zzS;}

  virtual float& xyS(){return _xyS;}
  virtual float& xzS(){return _xzS;}
  virtual float& yzS(){return _yzS;}

  virtual float& xyA(){return _xyA;}
  virtual float& xzA(){return _xzA;}
  virtual float& yzA(){return _yzA;}

  virtual float& iso(){return _iso;}

  virtual void packSourceMessage(modelDefStruct* modelDef,
                                 int startIteration);
  virtual void writeSource(int NT,float ampScalar,size_t index,
			   char* filename,int outfile,
			   int xVar,int yVar,int zVar,
			   int ampVar,int dataVar,
                           int* extraVars);

  //In the elastic/anelastic case this source is applied to velocities and stresses.
  virtual int applyVel(modelDefStruct* modelDef,int iteration,
		       float *vx,float *vy,float *vz,
                       float* rho);

  virtual int applyStress(modelDefStruct* modelDef,int iteration,
			  float *xx,float* yy,float *zz,
                          float *xy,float* xz,float *yz);

  //And in the acoustic case this source is applyed to both the velocities and
  // the pressure but in a slightly different way.
  virtual int applyAcousticVel(modelDefStruct* modelDef,int iteration,
			       float *vx,float *vy,float *vz,
			       float* rho){
    //Need to fix this later, for now just use elastic version.
    return applyVel(modelDef,iteration,
		    vx,vy,vz,
		    rho);
  }

  virtual int applyPressure(modelDefStruct* modelDef,int iteration,
                            float *pressure);

  virtual int applyPressureMM(modelDefStruct* modelDef,int iteration,
                              float *pressure);

protected:
  void calcInterpCoeffs(modelDefStruct* modelDef,float* rho);
};

class forceSource: public virtual pointSource{
public:
  //Need some additional parameters.
  float _ax,_ay,_az;
  float _rhox[8], _rhoy[8], _rhoz[8];

  //Needs to be able to extrapolate values to the 3 components of velocity. Also requires
  // the interpolated bouyancy for scaling.
  interpStruct _vxInterp,_vyInterp,_vzInterp;

  //Constructors

  ///Null constructor for derived classes.
  forceSource():pointSource(){
    _ax=_ay=_az=0.0;
    for(int i=0;i<8;++i) _rhox[i] = _rhoy[i] = _rhoz[i] = 1e13;
  }

  ///Parameters passed as arguments.
  forceSource(modelDefStruct* modelDef,
	      float x,float y,float z,
	      float ax,float ay,float az,
	      float amp,float* data,
              float rho);

  ///Values passed.
  forceSource(modelDefStruct* modelDef,
	      float x,float y,float z,
	      float ax,float ay,float az,
	      float amp,float* data,
              float* rho);

  ///Parameters read from message.
  forceSource(modelDefStruct* modelDef,int startTag,float* rho);

  ///Parameters read from CDF file.
  forceSource(modelDefStruct* modelDef,
	      char* filename,int infile,size_t index,
	      int xVar,int yVar,int zVar,
	      int axVar,int ayVar,int azVar,
	      int ampVar,int dataVar,
	      float* rho);

  //Now the rest of the functions.
  virtual int type(){return FORCE_SOURCE;}
  virtual float& ax(){return _ax;}
  virtual float& ay(){return _ay;}
  virtual float& az(){return _az;}

  virtual void packSourceMessage(modelDefStruct* modelDef,
                                 int startIteration);
  virtual void writeSource(int NT,float ampScalar,size_t index,
			   char* filename,int outfile,
			   int xVar,int yVar,int zVar,
			   int ampVar,int dataVar,
                           int* extraVars);

  //Need to apply this source to velocities only.
  virtual int applyVel(modelDefStruct* modelDef,int iteration,
		       float *vx,float *vy,float *vz,
                       float* rho);
  virtual int applyVelMM(modelDefStruct* modelDef,int iteration,
		       float *vx,float *vy,float *vz,
                         float* rho);
  //And in the acoustic case this source is applyed to just the velocities in exactly
  // the same way as the elastic case.
  virtual int applyAcousticVel(modelDefStruct* modelDef,int iteration,
			       float *vx,float *vy,float *vz,
			       float* rho){
    //Need to fix this later, for now just use elastic version.
    return applyVel(modelDef,iteration,
		    vx,vy,vz,
		    rho);
  }
  virtual int applyAcousticVelMM(modelDefStruct* modelDef,int iteration,
			       float *vx,float *vy,float *vz,
			       float* rho){
    //Need to fix this later, for now just use elastic version.
    return applyVelMM(modelDef,iteration,
		    vx,vy,vz,
		    rho);
  }


protected:
  void calcInterpCoeffs(modelDefStruct* modelDef,float* rho);
};

///Here is the subclassed pointSource that implements a surface traction. NOTE
/// the scaling here is not correct. We should be using some factors of lambda and
/// mu or the bulk modulus that I have yet figured out how to pass.
class tractionSource: public virtual pointSource{
public:
  //Need some additional parameters.
  float _ax,_ay,_az;
  float _rho;

  //Needs to be able to extrapolate values to the 3 components of velocity. Also requires
  // the interpolated bouyancy for scaling.
  interpStruct _vxInterp,_vyInterp,_vzInterp;

  //Constructors
  // Parameters passed as arguments.
  tractionSource(modelDefStruct* modelDef,
		 float x,float y,float z,
		 float ax,float ay,float az,
		 float amp,float* data,
                 float rho);
  tractionSource(modelDefStruct* modelDef,
		 float x,float y,float z,
		 float ax,float ay,float az,
		 float amp,float* data,
                 float* rho);

  // Parameters read from message.
  tractionSource(modelDefStruct* modelDef,int startTag,float* rho);

  // Parameters read from CDF file.
  tractionSource(modelDefStruct* modelDef,
		 char* filename,int infile,size_t index,
		 int xVar,int yVar,int zVar,
		 int axVar,int ayVar,int azVar,
		 int ampVar,int dataVar,
		 float* rho);

  //Now the rest of the functions.
  virtual int type(){return TRACTION_SOURCE;}
  virtual float& ax(){return _ax;}
  virtual float& ay(){return _ay;}
  virtual float& az(){return _az;}

  virtual void packSourceMessage(modelDefStruct* modelDef,
                                 int startIteration);
  virtual void writeSource(int NT,float ampScalar,size_t index,
			   char* filename,int outfile,
			   int xVar,int yVar,int zVar,
			   int ampVar,int dataVar,
                           int* extraVars);

  //Need to apply this source to velocities only.
  virtual int applyVel(modelDefStruct* modelDef,int iteration,
		       float *vx,float *vy,float *vz,
                       float* rho);
  //And in the acoustic case this source is applyed to just the velocities in exactly
  // the same way as the elastic case.
  virtual int applyAcousticVel(modelDefStruct* modelDef,int iteration,
			       float *vx,float *vy,float *vz,
			       float* rho){
    //Need to fix this later, for now just use elastic version.
    return applyVel(modelDef,iteration,
		    vx,vy,vz,
		    rho);
  }


protected:
  void calcInterpCoeffs(modelDefStruct* modelDef,float rho){
    calcInterpCoeffs(modelDef,(float*)NULL);
    _rho=rho;
  }
  void calcInterpCoeffs(modelDefStruct* modelDef,float* rho);
};

//
//Now a single class to hold all the sources in a given problem.
class sourceNetwork{
public:
  sourceArray* _sources;
  size_t _nfs,_nms,_nts,_ntdbc; 
  //These could be derived by it is handy to have them available.

  float* _wavelet; //The wavelet for new sources.
  float _bandCheckPCent; //Percentage level to check wavelet bandwidth.
  int _useModelSources; //-SM
  bool _useMomentRate;
  //Constructors
  // The null constructor, sources specified on the command line or added indivially.
  sourceNetwork(){
    _useModelSources=TRUE;
    _useMomentRate = false;
    _bandCheckPCent=DEFAULT_CHECK_BANDWIDTH_PCENT;
    _wavelet=NULL;
    _nfs=_nms=_nts=_ntdbc=0;
    _sources=new sourceArray;
  }

  // Sources passed in message.
  sourceNetwork(modelDefStruct* modelDef,float *rho=NULL,
                int useCubicExtrap=FALSE);

  // Sources read from cdf file
  sourceNetwork(modelDefStruct* modelDef,char* filename,
                int useModelSources=TRUE,float *rho=NULL);

  virtual ~sourceNetwork(){
    if(_wavelet)
      free(_wavelet);
    _wavelet=NULL;

    if(_sources){
      for(int i=0;i<size();i++){
	//source* currSource=(*_sources)[i];
	//delete currSource; //KLUDGE, this is yeilding an error for unknown reasons.
      }
      delete _sources;
    }
    _sources=NULL;
  }

  //A few helper routines.
  int size(){return _sources->size();}
  int nfs(){return _nfs;}
  int nms(){return _nms;}
  int nts(){return _nts;}
  int ntdbc(){return _ntdbc;}

  float& x(int i){return (*_sources)[i]->x();}
  float& y(int i){return (*_sources)[i]->y();}
  float& z(int i){return (*_sources)[i]->z();}
  float& amp(int i){return (*_sources)[i]->amp();}

  int type(int i){return (*_sources)[i]->type();}
  float* data(int i){return (*_sources)[i]->data();}

  float& ax(int i){return (*_sources)[i]->ax();}
  float& ay(int i){return (*_sources)[i]->ay();}
  float& az(int i){return (*_sources)[i]->az();}

  float& xxS(int i){return (*_sources)[i]->xxS();}
  float& yyS(int i){return (*_sources)[i]->yyS();}
  float& zzS(int i){return (*_sources)[i]->zzS();}

  float& xyS(int i){return (*_sources)[i]->xyS();}
  float& xzS(int i){return (*_sources)[i]->xzS();}
  float& yzS(int i){return (*_sources)[i]->yzS();}

  float& xyA(int i){return (*_sources)[i]->xyA();}
  float& xzA(int i){return (*_sources)[i]->xzA();}
  float& yzA(int i){return (*_sources)[i]->yzA();}

  source* s(int index){return (*_sources)[index];}

  float* wavelet(){return _wavelet;}

  int combine(sourceNetwork* extraSources);

  int sendSources(modelDefStruct* modelDef,
                  int target=AllProcesses);

  virtual int applyVel(modelDefStruct* modelDef,int iteration,
		       float *vx,float *vy,float *vz,
                       float* rho);
  virtual int applyStress(modelDefStruct* modelDef,int iteration,
			  float *xx,float* yy,float *zz,
                          float *xy,float* xz,float *yz);

  virtual int applyAcousticVel(modelDefStruct* modelDef,int iteration,
			       float *vx,float *vy,float *vz,
                               float* rho);
  virtual int applyAcousticVelMM(modelDefStruct* modelDef,int iteration,
			       float *vx,float *vy,float *vz,
                                 float* rho);
  virtual int applyPressure(modelDefStruct* modelDef,int iteration,
                            float *pressure);
  virtual int applyPressureMM(modelDefStruct* modelDef,int iteration,
                              float *pressure);

  //Bigger routines.
  // Write all the current sources to a cdf file.
  int writeSources(modelDefStruct* modelDef,char* filename,int createNew=FALSE);

  ///Set the source differentiation factor.
  virtual int setDispersionFactor(int i,int doPrint);

  ///Perform dispersion analysis and check CFL criteria.
  virtual int checkDispersion(modelDefStruct* modelDef,
			      float vMin,float vMax,
                              int doPrint=TRUE,float *wavelengths=NULL);

  // Add a source to an existing source network, this routine is called from inside 
  //  processArgs.
  forceSource* addForceSource(modelDefStruct* modelDef,
			      float x,float y,float z,
			      float ax,float ay,float az,
			      float amp,float* data,
			      float* rho){
    forceSource* newSource;
    _sources->Add(newSource=
		  new forceSource(modelDef,x,y,z,
				  ax,ay,az,amp,data,rho));
    _nfs++;
    return newSource;
  }
  forceSource* addForceSource(modelDefStruct* modelDef,
			      float x,float y,float z,
			      float ax,float ay,float az,
			      float amp,float* data,
			      float rho){
    forceSource* newSource;
    _sources->Add(newSource=
		  new forceSource(modelDef,x,y,z,
				  ax,ay,az,amp,data,rho));
    _nfs++;
    return newSource;
  }
  momentSource* addExplosiveSource(modelDefStruct* modelDef,
				   float x,float y,float z,float samp,
				   float* sourceTimeFunc){
    momentSource* newSource;
    _sources->Add(newSource=
		  new momentSource(modelDef,x,y,z,
				   1.0,1.0,1.0,
				   0.0,0.0,0.0,
				   0.0,0.0,0.0,
				   samp,
				   sourceTimeFunc,NULL,_useMomentRate));
    _nms++;
    return newSource;
  }
  virtual int addSources(sourceNetwork* extraSources);
  virtual int addSources(int& i,int argOffset,
			 int argc,char* argv[],
                         modelDefStruct* modelDef);

protected:
  void checkWaveletDispersion(float dx,float dy,float dz,
			      int nt,float dt,
                              float minVel,int derivative);

  float* readWavelet(int NT,char* filename);

  float* rickerWavelet(int nt,float dt,float fmode,
		       int normalizationDerivative=0,
                       int setPeakTime=0,float peakTime=0.0);

  //
  //Helper routines to write sources to CDF file.

  ///Write force sources to a NetCDF file.
  int writeForceSources(modelDefStruct* modelDef,char* filename,
                        int outfile,int ntDim);

  ///Write moment sources to a NetCDF file.
  int writeMomentSources(modelDefStruct* modelDef,char* filename,
                         int outfile,int ntDim);

  ///Write traction sources to a NetCDF file.
  int writeTractionSources(modelDefStruct* modelDef,char* filename,
                           int outfile,int ntDim);

  //Helper routines to read sources from CDF file.
  ///Read force sources.
  virtual int readForceSources(modelDefStruct* modelDef,
			       char* filename,int infile,
                               float* rho);

  ///Method to read any moment sources form the CDF file.
  virtual int readMomentSources(modelDefStruct* modelDef,
				char* filename,int infile,
                                float* rho);

  ///Read traction sources from a NetCDF file.
  virtual int readTractionSources(modelDefStruct* modelDef,
				  char* filename,int infile,
                                  float* rho);

  //
  //Direct wavelet generation subroutines. 
  void cTimesCExp(double* inOut,double r,double c);
  double* rickerSpec(int nf,double df,double fmode,double tmin,
                     double *spect);
  void buildGaussianWavelet(float freq,int nt,float dt);

  //
  //And some more internal functions to calculate the dispersion limits for a given source.
  int nFreqs(int nt);
  //calculate the frequency step (depends on nf and Dt)
  float dFreq(int nf,float dt){
    return 1.0/((float)2*(float)nf*dt);
  }
  //calculate the fourier aspec of the data
  // returns wrk array and allocates it as required
  float* sourceSpectrum(float* data,int scaleFactor,
			int nt,float dt,
			int& nf,float& df,
                        floatPtr& wrk);
  //Subroutine BANDIT determines the bandwidth of a frequency spectrum
  // between amplitude levels that are a specified percentage of the 
  // maximum amplitude.
  float bandit(float pcent,
	      float dx,float dy,float dz,
	      float df,int nf,float vmin,
               float *data,int doPrint=TRUE);
};
typedef sourceNetwork* sourceNetworkPtr;
#endif //#ifndef _sgfdSources_hh_
