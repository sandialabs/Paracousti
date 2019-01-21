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
 *  acousticBoundary.hh
 *
 *
 *  Declare acoustic boundary related classes.  These are mainly
 *  classes that hold information about the boundaries instead of actually
 *  doing the boundary updates.  The neighbor process classes do the
 *  message passing among domains.
 *
 *  Declares classes:
 *    acousticBoundary
 *    acousticNeighborProcess
 *    acousticSgfdBoundaries
 *
 */

#ifndef _acousticBoundary_hh_
#define _acousticBoundary_hh_

#include "nstdutil.hh"

#include "sgfdBoundary.hh"

//
//Here are subclasses of the main classes, these are the versions that are actually going
// to be used for the acoustic code.
class acousticNeighborProcess: public neighborProcess{
public:
  ///Dumb constructor for derived subclasses.
  acousticNeighborProcess(sgfdModel* model,
			  int sideI,int sideJ,int sideK):
    neighborProcess(model,sideI,sideJ,sideK){
  }

  ///Constructor needs velocities and pressures, only 1 time step.
  acousticNeighborProcess(sgfdModel* model,sgfdDependent* dependent,
			  int sideI,int sideJ,int sideK);

protected:
  //Internal functions
  virtual void receiveNeighborVel(modelDefStruct* modelDef,
				  float* vx,float* vy,float* vz,
                                  int messageTagInc=0);
  virtual void sendNeighborVel(modelDefStruct* modelDef,
			       float* vx,float* vy,float* vz,
                               int messageTagInc=0);

  virtual void sendNeighborStress(modelDefStruct* modelDef,
				  float* xx,float* yy,float* zz,
				  float* xy,float* xz,float* yz,
                                  int messageTagInc=0);

  virtual void receiveNeighborStress(modelDefStruct* modelDef,
				     float* xx,float* yy,float* zz,
				     float* xy,float* xz,float* yz,
                                     int messageTagInc=0);

#if USE_MPI_SEND
  virtual void buildMPITypes(sgfdDependent* dependent);
#endif //#if USE_MPI_SEND
};

///Derived acoustic class for the complete static boundary condition.
class acousticSgfdBoundaries: public sgfdBoundaries{
public:
  int _usePML;
  int _freeSurfaceMode;
  bool minXb, maxXb, minYb, maxYb, minZb, maxZb;
  bool fminXb, fmaxXb, fminYb, fmaxYb, fminZb, fmaxZb;
  ///Dumb constructor for derived classes.
  acousticSgfdBoundaries():sgfdBoundaries(){}
  
  ///\brief Constructor for the acoustic parallel problem.
  acousticSgfdBoundaries(sgfdModel* model,sgfdDependent* dependent,
                         int freeSurfaceMode,int usePML);
};

#endif //#ifndef _acousticBoundary_hh_
