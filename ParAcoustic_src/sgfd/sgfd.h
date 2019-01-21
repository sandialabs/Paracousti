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
 *  sgfd.h
 *
 *
 *  Contains several definitions and declarations needed by many parts
 *  of the code.  Macros used to simplify writing code and structures
 *  that contain model and parallel information are declared.  Function
 *  prototypes for interpolation and extrapolation in the grid are
 *  also declared.
 *
 *  Declares the following functions:
 *  trilinCoeff
 *  trilinInterp
 *  trilinInterpOffset
 *  trilinExtrap
 *  triCubicCoeff
 *  triCubicInterp
 *  triCubicExtrap
 *
 *  Declares many macros and the modelDef and parallelDef structures as
 *  well.
 *
 */
/*Structures and MACROs that are used by C and C++ code

    In order to simplify subroutine calls there are some structures
    defined in this file which can be used by both C and C++
    codes. This is also where the numeric codes for the main sgfd
    messages are defined.
*/

#ifndef _sgfd_h_
#define _sgfd_h_

//There can be a problem with sending the messages for source/receiver info 
// when there are a very large number of time-steps.
#ifndef MAX_TIME_STEPS_SEND
#define MAX_TIME_STEPS_SEND    40002
#endif

/*Define the meaning of the controls used for interprocess communications*/
#define MESSAGE_UPDATE_VEL    10
#define MESSAGE_UPDATE_STRESS 11
#define MESSAGE_ADVANCE_ONE_ITERATION 12

#define MESSAGE_SLICE         100
#define MESSAGE_FULL_FIELD    101
#define MESSAGE_SEND_RECEIVER 102
#define MESSAGE_SEND_SOURCE   103
#define MESSAGE_SEND_RECEIVER_GRID 105

#define MESSAGE_SET_BC          1003
#define MESSAGE_SET_SOURCES     1004

#define MESSAGE_SLAVE_TO_SLAVE 1201

//Defines for output types
#define SLICE_YZ_PLANE 1
#define SLICE_XZ_PLANE 2
#define SLICE_XY_PLANE 3

#define SLICE_VX_COMP    1
#define SLICE_VY_COMP    2
#define SLICE_VZ_COMP    3
#define SLICE_PRESS_COMP 4
#define SLICE_RX_COMP    8
#define SLICE_RY_COMP    9
#define SLICE_RZ_COMP   10

#define SLICE_VP_COMP    5
#define SLICE_VS_COMP    6
#define SLICE_RHO_COMP   7

/*Macro's to access members of the main arrays (earth model, velocity, and stress)*/
#define RHO(i,j,k) (rho[(i)+(j)*NX+(k)*NXY])
#define LAMBDA(i,j,k) (lambda[(i)+(j)*NX+(k)*NXY])
#define OOMU(i,j,k) (mu[(i)+(j)*NX+(k)*NXY])
#define MU(i,j,k) (mu[(i)+(j)*NX+(k)*NXY]==-1?0.0:(1.0/mu[(i)+(j)*NX+(k)*NXY]))

#define VX(i,j,k) (vx[(i)+(j)*NX+(k)*NXY])
#define VY(i,j,k) (vy[(i)+(j)*NX+(k)*NXY])
#define VZ(i,j,k) (vz[(i)+(j)*NX+(k)*NXY])

#define XX(i,j,k) (xx[(i)+(j)*NX+(k)*NXY])
#define YY(i,j,k) (yy[(i)+(j)*NX+(k)*NXY])
#define ZZ(i,j,k) (zz[(i)+(j)*NX+(k)*NXY])
#define XY(i,j,k) (xy[(i)+(j)*NX+(k)*NXY])
#define XZ(i,j,k) (xz[(i)+(j)*NX+(k)*NXY])
#define YZ(i,j,k) (yz[(i)+(j)*NX+(k)*NXY])

#define PRESSURE(i,j,k) (pressure[(i)+(j)*NX+(k)*NXY])

///Model definition,

/// This structure contains the primary definition of the domain size
/// and limits. It is used by the master process to define the entire
/// region and by the slaves to define the subdomains. The
/// non-dimensionalizing scalars are also held in this structure. The
/// MACRO's DEF_MODEL_SIZE, DEF_MODEL_LIMITS, and DEF_MODEL_SCALE can
/// be called with a pointer to this structure to define variables
/// defining the relavent properties.
typedef struct _modelDef_{
  int NX,NY,NZ,NT,NXY,NXYZ;
  float minX,dx,minY,dy,minZ,dz,minT,dt;

  float scalarSpeed,scalarVel,scalarDen;
  float scalarStress;
  int procLim[6];
} modelDefStruct;
typedef modelDefStruct* modelDefStructPtr;

//define a macro to get the most commonly used modelDefStruct
// values into local vars for faster exection
#define DEF_MODEL_SIZE(md) int NX=(md)->NX,NY=(md)->NY,NZ=(md)->NZ,NT=(md)->NT, \
                               NXY=(md)->NXY,NXYZ=(md)->NXYZ, \
                               *procLim=(md)->procLim
#define DEF_MODEL_LIMITS(md) float minX=(md)->minX,dx=(md)->dx, \
                                   minY=(md)->minY,dy=(md)->dy, \
                                   minZ=(md)->minZ,dz=(md)->dz, \
                                   minT=(md)->minT,dt=(md)->dt
#define DEF_MODEL_SCALE(md) float scalarSpeed=(md)->scalarSpeed,scalarVel=(md)->scalarVel,\
                                  scalarDen=(md)->scalarDen,scalarStress=(md)->scalarStress


///Define the parallel decomposition,

///Like the modelDefStruct this is used by the master to define all
///the subdomains and by the slaves to define the local neighbors for
///message passing across boundarys. The MACRO DEF_PARALLEL defines
///variable extracted from this structure.
typedef struct _parallelDef_{
  int neighbors[27];
  modelDefStruct* neighborDefs[27];
  int nxProc,nyProc,nzProc,nxyProc;
  int procI,procJ,procK;
  int globalNX,globalNY,globalNZ;
} parallelDefStruct;
typedef parallelDefStruct* parallelDefStructPtr;

//define a macro to define the commonly used elements as local vars
#define DEF_PARALLEL(pd) int *neighbors=(pd)->neighbors, \
                             nxProc=(pd)->nxProc,nyProc=(pd)->nyProc,nzProc=(pd)->nzProc, \
                             nxyProc=(pd)->nxyProc, \
                             procI=(pd)->procI,procJ=(pd)->procJ,procK=(pd)->procK, \
                             globalNX=(pd)->globalNX,globalNY=(pd)->globalNY,globalNZ=(pd)->globalNZ

/*Macro's to access*/
#define SUBDOMAIN(elem,i,j,k) ((elem)[(i)+(j)*nxProc+(k)*nxyProc])
#define TID(i,j,k) SUBDOMAIN(Tids,i,j,k)

#define COWORKDOMAIN(pelem,elem,i,j,k) ((!i&&!j&&!k)?(pelem): \
				   ((elem)[(i)+(j)*nxProc+(k)*nxyProc-1]))
#define COWORKTID(i,j,k) COWORKDOMAIN(Parent,Tids,i,j,k)

#define BOUNDARY(bc,i,j,k) ((bc)[(i+1)+3*(j+1)+9*(k+1)])
#define NEIGHBOR(i,j,k) BOUNDARY(neighbors,i,j,k)

//
//Structures for interpolation.
//
///Cubic interplation/extraplation uses 64 nodes surrounding the point.
typedef struct _cubicInterpStruct_{
  int shortCircut; //If one coeff is 1, make iptr, jptr, kptr go to that index.
  float coeff[64]; //add this to the node
  int iptr,jptr,kptr;   //node indices of the lower left front corner.
} cubicInterpStruct;

///Linear interplation/extraplation uses 8 nodes surrounding the point.
typedef struct _interpStruct_{
  float coeff[8]; //add this to the node
  int iptr[8],jptr[8],kptr[8];   //node indices
} interpStruct;
  
//
//Also need to have prototypes for some utility functions that are defined in sgfd_util.c
#if defined(__cplusplus)
extern "C" {
#endif
  /*!Computes the eight coefficients of a 
    trilinear interpolator appropriate for a specified receiver
    location. */
  void trilinCoeff(float xr,float yr,float zr,
		   float xmin,float dx,float ymin,float dy,float zmin,float dz,
		   interpStruct* interpCoeff);

  /*! Computes an approximate value of an input
      3D array using trilinear interpolation. */
  float trilinInterp(modelDefStruct* modelDef,
		     interpStruct* interpCoeff,
		     float* vx);
  /*! Computes an approximate value of an input
      3D array using trilinear interpolation with offsets. */
  float trilinInterpOffset(modelDefStruct* modelDef,
			   interpStruct* interpCoeff,
			   int xoff,int yoff,int zoff,
			   float* vx);


  /*! Updates a 3D array by extrapolating an
      input scalar value to the eight surrounding gridpoints. */
  float trilinExtrap (modelDefStruct* modelDef,
		      float* vx,
		      interpStruct* interpCoeff,
		      float value);

  /*!Computes the 64 coefficients of a 
    tricubic interpolator appropriate for a specified receiver
    location. */
  void triCubicCoeff(float x,float y,float z,
		    float minX,float dx,float minY,float dy,float minZ,float dz,
		    cubicInterpStruct *cubicInterp);
    /*! Computes an approximate value of an input
      3D array using tricubic interpolation. */
  float triCubicInterp(modelDefStruct* modelDef,
		       cubicInterpStruct* interpCoeff,
		       float* vx);
  /*! Updates a 3D array by extrapolating an
      input scalar value to the 64 surrounding gridpoints. */
  float triCubicExtrap(modelDefStruct* modelDef,
		       cubicInterpStruct* interpCoeff,
		       float* vx,
		       float value);
#if defined(__cplusplus)
}
#endif

#endif //_sgfd_h_
