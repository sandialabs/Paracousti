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
 *  movingSgfdSources.hh
 *
 *
 *  Declare several classes that control experimental moving
 *  sources.  Even stationary sources can call these classes because
 *  they simply override the standard stationary sources.  The
 *  moving sources override the specific standard stationary classes
 *  defined in sgfdSources.  The movingSource Network class controls
 *  groups of sources.
 *
 *  Declares the following classes:
 *  movingSource
 *  movingForceSource
 *  movingMomentSource
 *  movingSourceNetwork
 *
 */
#ifndef _movingSgfdSources_hh_
#define _movingSgfdSources_hh_

#include "sgfdSources.hh"

//Define some new types.
#define MOVING_FORCE_SOURCE    5
#define MOVING_MOMENT_SOURCE   6

///Here is a new class for moving sources. This will be used by derived classes
/// for multiple inheretance with a class derived from a point source (force,
/// moment, etc.).
class movingSource: public virtual pointSource{
public:
  float *_xMoving,*_yMoving,*_zMoving; 
  //Need to save a local copy of the bouyancy since it will be
  // required to calculate the interpolator coefficients at each
  // iteration.
  float* _rhoField;

  ///This constructor for this one is similar to a point source but it needs the
  /// location as a function of time.
  movingSource(modelDefStruct* modelDef,
	       float* x,float* y,float* z,float amp,
	       float* data,float *rho):
    pointSource(modelDef,x[0],y[0],z[0],amp,data){
    _xMoving=x;
    _yMoving=y;
    _zMoving=z;

    _rhoField=rho;
  }

  ///Here is constructor that reads the values from a message.
  movingSource(modelDefStruct* modelDef,int startTag,float *rho);

  ///Constructor that reads parameters read from CDF file.
  movingSource(modelDefStruct* modelDef,
	       char* filename,int infile,size_t index,
	       int xVar,int yVar,int zVar,
	       int ampVar,int dataVar);

  virtual ~movingSource(){
    if(_xMoving)
      free(_xMoving);
    if(_yMoving)
      free(_yMoving);
    if(_zMoving)
      free(_zMoving);
  }

  ///Additional message reading function. If there are a large number of 
  /// time-steps we use an additional function to read the rest of the data.
  virtual int readAuxMessage(int NT,int startIndex,
                             int startTag,int iGetMessg);

  ///Pack the correct source message. This is just type time variable location,
  /// we assume the other parameters are packed by other inhereited classes.
  virtual void packSourceMessage(modelDefStruct* modelDef,
                                 int startIteration);

  ///The applyVel and applyAcousticVel are meant to be called in
  /// conjunction with the same methods for the actual implemented
  /// source type. All this class method has to do is reset the 
  /// location. The derived class will then use the new location
  /// and call the calcInterpCoeffs method to reset the interpolators.
  virtual int applyVel(modelDefStruct* modelDef,int iteration,
		       float *vx,float *vy,float *vz,
		       float* rho){
    _x=_xMoving[iteration];
    _y=_yMoving[iteration];
    _z=_zMoving[iteration];
    return iteration;
  }

  ///Acoustic case is exactly the same for this class, just call
  /// applyVel.
  virtual int applyAcousticVel(modelDefStruct* modelDef,int iteration,
			       float *vx,float *vy,float *vz,
			       float* rho){
    //Need to fix this later, for now just use elastic version.
    return movingSource::applyVel(modelDef,iteration,
				  vx,vy,vz,
				  rho);
  }
  virtual int applyStress(modelDefStruct* modelDef,int iteration,
			  float *xx,float* yy,float *zz,
			  float *xy,float* xz,float *yz){
    _x=_xMoving[iteration];
    _y=_yMoving[iteration];
    _z=_zMoving[iteration];
    return iteration;
  }
  virtual int applyPressure(modelDefStruct* modelDef,int iteration,
			    float *pressure){
    return movingSource::applyStress(modelDef,iteration,
				     pressure,NULL,NULL,NULL,NULL,NULL);
  }
};

///Here is an actual moving-source class. This is the moving force
/// source so it inherits both the movingSource and the forceSource
/// classes.
class movingForceSource:public  forceSource, public movingSource{
public:
  //Do not need any new local variables, everything is either in
  // movingSource or forceSource.

  //
  //Constructors.

  ///Pass the values in as arguments.
  movingForceSource(modelDefStruct* modelDef,
		    float* x,float* y,float* z,
		    float ax,float ay,float az,
		    float amp,float* data,
		    float *rho):
    pointSource(modelDef,x[0],y[0],z[0],amp,data),
    forceSource(modelDef,x[0],y[0],z[0],ax,ay,az,amp,data,rho),
    movingSource(modelDef,x,y,z,amp,data,rho){}
  
  ///Parameters read from CDF file.
  movingForceSource(modelDefStruct* modelDef,
		    char* filename,int infile,size_t index,
		    int xVar,int yVar,int zVar,
		    int axVar,int ayVar,int azVar,
		    int ampVar,int dataVar,
		    float* rho):
    pointSource(modelDef,filename,infile,index,-1,-1,-1,ampVar,dataVar),
    forceSource(modelDef,filename,infile,index,-1,-1,-1,
              axVar,ayVar,azVar,ampVar,dataVar,rho),
    movingSource(modelDef,filename,infile,index,xVar,yVar,zVar,ampVar,dataVar){
    calcInterpCoeffs(modelDef,rho);
  }

  ///Parameters read from a message.
  movingForceSource(modelDefStruct* modelDef,int startTag,
                    float* rho);

  ///What about a destructor?
  virtual ~movingForceSource(){}

  ///Return the MOVING_FORCE_SOURCE type.
  virtual int type(){return MOVING_FORCE_SOURCE;}


  ///Pack the correct source message. This is just a combination of force source
  /// and moving source things in the correct order.
  virtual void packSourceMessage(modelDefStruct* modelDef,
				 int startIteration){
    forceSource::packSourceMessage(modelDef,startIteration);
    movingSource::packSourceMessage(modelDef,startIteration);
  }
  
  ///The applyVel method is easy, just call movingSource::applyVel
  /// and then forceSource::applyVel.
  virtual int applyVel(modelDefStruct* modelDef,int iteration,
		       float *vx,float *vy,float *vz,
                       float* rho);
  ///The applyAcousticVel method is easy, just call movingSource::applyAcousticVel
  /// and then forceSource::applyAcousticVel.
  virtual int applyAcousticVel(modelDefStruct* modelDef,int iteration,
			       float *vx,float *vy,float *vz,
                               float* rho);
};

///Here is an actual moving-source class. This is the moving moment
/// source so it inherits both the movingSource and the momentSource
/// classes.
class movingMomentSource:public  momentSource, public movingSource{
public:
  //Do not need any new local variables, everything is either in
  // movingSource or momentSource.

  //
  //Constructors.

  ///Pass the values in as arguments.
  movingMomentSource(modelDefStruct* modelDef,
		     float* x,float* y,float* z,
		     float xxS,float yyS,float zzS,
		     float xyS,float xzS,float yzS,
		     float xyA,float xzA,float yzA,
		     float amp,float* data,
		     float *rho, int useMomentRate):
    pointSource(modelDef,x[0],y[0],z[0],amp,data),
    momentSource(modelDef,x[0],y[0],z[0],
               xxS,yyS,zzS,xyS,xzS,yzS,xyA,xzA,yzA,amp,data,rho,useMomentRate),
    movingSource(modelDef,x,y,z,amp,data,rho){}
		 
  ///Parameters read from CDF file.
  movingMomentSource(modelDefStruct* modelDef,
		     char* filename,int infile,size_t index,
		     int xVar,int yVar,int zVar,
		     int xxSVar,int yySVar,int zzSVar,
		     int xySVar,int xzSVar,int yzSVar,
		     int xyAVar,int xzAVar,int yzAVar,
		     int ampVar,int dataVar,
		     float* rho):
    pointSource(modelDef,filename,infile,index,-1,-1,-1,ampVar,dataVar),
    momentSource(modelDef,filename,infile,index,-1,-1,-1,
               xxSVar,yySVar,zzSVar,xySVar,xzSVar,yzSVar,xyAVar,xzAVar,yzAVar,
               ampVar,dataVar,rho),
    movingSource(modelDef,filename,infile,index,xVar,yVar,zVar,ampVar,dataVar){
    calcInterpCoeffs(modelDef,rho);
  }
		 
  ///Parameters read from a message.
  movingMomentSource(modelDefStruct* modelDef,int startTag,float* rho);
		
  ///What about a destructor?
  virtual ~movingMomentSource(){}

  ///Return the MOVING_MOMENT_SOURCE type.
  virtual int type(){return MOVING_MOMENT_SOURCE;}


  ///Pack the correct source message. This is just a combination of moment source
  /// and moving source things in the correct order.
  virtual void packSourceMessage(modelDefStruct* modelDef,
				 int startIteration){
    momentSource::packSourceMessage(modelDef,startIteration);
    movingSource::packSourceMessage(modelDef,startIteration);
  }
  
  ///The applyVel method is easy, just call movingSource::applyVel
  /// and then momentSource::applyVel.
  virtual int applyVel(modelDefStruct* modelDef,int iteration,
		       float *vx,float *vy,float *vz,
                       float* rho);
  ///The applyAcousticVel method is easy, just call movingSource::applyAcousticVel
  /// and then momentSource::applyAcousticVel.
  virtual int applyAcousticVel(modelDefStruct* modelDef,int iteration,
			       float *vx,float *vy,float *vz,
                               float* rho);
  virtual int applyStress(modelDefStruct* modelDef,int iteration,
			  float *xx,float* yy,float *zz,
                          float *xy,float* xz,float *yz);
  virtual int applyPressure(modelDefStruct* modelDef,int iteration,
                            float *pressure);
};


///Now redefine the container class so it can add the moving sources.
class movingSourceNetwork: public sourceNetwork{
public:
  //Constructors
  ///The null constructor, sources specified on the command line or added 
  /// indivially.
  movingSourceNetwork():sourceNetwork(){}
  
  ///Sources passed in message.
  movingSourceNetwork(modelDefStruct* modelDef,float *rho=NULL,
                      int useCubicExtrap=FALSE);
  
  ///Sources read from cdf file
  movingSourceNetwork(modelDefStruct* modelDef,char* filename,
                      int useModelSources=TRUE,float *rho=NULL);

  ///What about a destructor?
  virtual ~movingSourceNetwork(){}

 protected:
  ///Add a force source to an existing source network, this routine is called from 
  /// inside processArgs.
  forceSource* addMovingForceSource(modelDefStruct* modelDef,
				    float* x,float* y,float* z,
				    float ax,float ay,float az,
				    float amp,float* data,
				    float* rho){
    movingForceSource* newSource;
    _sources->Add(newSource=
		  new movingForceSource(modelDef,x,y,z,
					ax,ay,az,amp,data,rho));
    _nfs++;
    return newSource;
  }
    
  ///Read force sources with special treatment for the moving ones.
  virtual int readForceSources(modelDefStruct* modelDef,
			       char* filename,int infile,
                               float* rho);

  ///Modification of the readMomentSources method. First strip off any 
  /// moving-moment sources, then call the sourceNetwork method.
  virtual int readMomentSources(modelDefStruct* modelDef,
				char* filename,int infile,
                                float* rho);
      
  ///Pick of the new types and set the factor. The dispersion check should really
  /// acount for the Doppler shift but I will ignore that for now.
  virtual int setDispersionFactor(int i,int doPrint);
};
#endif //#ifndef _movingSgfdSources_hh_
