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
 *  fixed_acoustic.hh
 *
 *
 *  This file contains class function declarations for the master
 *  fixed acoustic model.  It also defines a output group type.
 *
 *  Defines the following classes:
 *  acousticOutputGroup
 *  masterFixedMediaModel
 *
 */

#ifndef _fixed_acoustic_hh_
#define _fixed_acoustic_hh_

class selector;

#include "acoustic_control.hh"
#include "extra_output.hh"

/*! Output group for the acoustic problem

  This is were the messages and initialization methods are defined that 
  know about the new types of slices and acoustic volume output.
 Note:  This has not been tested recently
*/
class acousticOutputGroup:public extraOutputGroup{
public:
  ///Now need to change the contructor.
  acousticOutputGroup(modelDefStruct* modelDef):
    extraOutputGroup(modelDef){}

  ///Need to pick off some command line arguments.
  virtual int addExtraOutput(int& i,int argOffset,
			     int argc,char* argv[],
			     modelDefStruct* modelDef,
                             int runCall);
};

///class for master.  Need to see if this is really necessary
class masterFixedMediaModel:public masterSgfdModel{
public:

  //Use different initializers for the different types of
  // media.
  ///Here is the simplest constructor for derived classes.
  masterFixedMediaModel(int surfaceBCType,int gridMultiplier,
			 modelDefStruct* modelDef,parallelDefStruct* parallelDef,
			 char* modelName)
    :masterSgfdModel(modelDef,parallelDef,modelName,
		     surfaceBCType,gridMultiplier){
  }
  /// Constructor
  masterFixedMediaModel(int tids[],int mediaType,
			 int surfaceBCType,int gridMultiplier,
			 modelDefStruct* modelDef,parallelDefStruct* parallelDef,
			 char* modelName,int fdCoeffOrder,float* fdCoeffs,
			 int longRead=TRUE);

  ~masterFixedMediaModel(){
  }

protected:
  ///This virtual function is called at the end of initializeSlaves to receive an
  /// acknoledgement that the slaves have completed initialization. Redefine here to give some
  /// additional information on the acoustic model.
  virtual void getInitAck(int tids[]);

  ///Here is a revised (and non-virtual) acknoledgement method which is used
  /// by the masterSgfdModel constructors.
  void getAck(int tag,int* tids);
};

//class for master fixed attenuative master model
class masterFixedAttenMediaModel:public masterFixedMediaModel{
public:
  selector* _selectors;

  masterFixedAttenMediaModel(selector* selectors,int tids[],int mediaType,
                             int surfaceBCType,int gridMultiplier,
                             modelDefStruct* modelDef,parallelDefStruct* parallelDef,
                             char* modelName,int fdCoeffOrder,float* fdCoeffs,
                             int longRead=TRUE);
protected:
  virtual void sendInitMessage(char* modelName,int fdCoeffOrder,float* fdCoeffs,int longRead,
                               parallelDefStruct** slaveParallelDefs,
                               int tids[],int i,int j,int k,
                               int nxProc,int nxyProc,
                               int coworkParadyme=FALSE);
};

#endif //#ifndef _fixed_acoustic_hh_
