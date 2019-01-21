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
 *  constants.h
 *
 *
 *  Defines several basic macros that have general usability.
 *
 *  Declares the following macros:
 *  USE_DOUBLES_ONLY
 *  BOOL
 *  TRUE
 *  FALSE
 *  FLT_MAX
 *  SGN
 *  VAL
 *  ISEVEN
 *  ISODD
 *  MAX
 *  MAX3
 *  MIN
 *  MIN3
 *  ISMID
 *  ISMID_OC
 *  ISMID_CO
 *  ISOVERLAP
 *  SETMIN
 *  SETMAX
 *  SETMAXINDEX
 *  SETMINMAX
 *  ABS
 *  SQRNPS
 *  PI
 *  DEG_TO_RAD
 *  RAD_TO_DEG
 *
 *  Defines the following typedefs
 *  charPtr
 *  floatPtr
 *  doublePtr
 *  intPtr
 *
 */
/*
// MACROS AND CONSTANTS
*/
#ifndef _constants_h
#define _constants_h

#include <time.h>
#include <sys/time.h>

/*Globals used in all my programs*/
extern int Verbose;
extern char *CommandLine; /*This is the call string*/
extern time_t StartTime;  /*This is the starting time of the program*/

#ifdef USE_DOUBLES_ONLY
#define float double
#endif

#ifndef DELETE_WHEN_DONE
#define DELETE_WHEN_DONE 1
#endif

#ifndef BOOL
#define BOOL int
#endif

#ifndef TRUE
#define TRUE (1)
#endif /* TRUE */

#ifndef FALSE
#define FALSE (0)
#endif /* FALSE */

#ifndef FLT_MAX
#define FLT_MAX 1E+37
#endif

#ifndef SGN
#define SGN(val) (((val)<0)?(-1):(1))
#endif

#ifndef VAL
#define VAL(array,i,j,k) ((array)[(i)+(j)*NX+(k)*NXY])
#endif

#ifndef ISEVEN
#define ISEVEN(i)((2*(int)floor((double)(i)/2.))==i)
#endif

#ifndef ISODD
#define ISODD(i) (!ISEVEN(i))
#endif

#ifndef MAX
#define MAX(a,b) ( (a)>(b)?(a):(b))
#endif
#ifndef MAX3
#define MAX3(a,b,c) (MAX(MAX((a),(b)),(c)))
#endif

#ifndef MIN
#define MIN(a,b) ( (a)<(b)?(a):(b))
#endif
#ifndef MIN3
#define MIN3(a,b,c) (MIN(MIN((a),(b)),(c)))
#endif

#ifndef ISMID
#define ISMID(a,b,c) (((a)<=(b))&&((b)<=(c)))
#endif
#ifndef ISMID_OC
#define ISMID_OC(a,b,c) (((a)<=(b))&&((b)<(c)))
#endif
#ifndef ISMID_CO
#define ISMID_CO(a,b,c) (((a)<(b))&&((b)<=(c)))
#endif

#ifndef ISOVERLAP
#define ISOVERLAP(x1,x2,y1,y2) (ISMID_OC(x1,y1,x2) || ISMID_CO(x1,y2,x2) || ISMID_OC(y1,x1,y2) || ISMID_CO(y1,x2,y2))
#endif

#ifndef SETMIN
#define SETMIN(target,source,index) \
 ((target)=((index)?(((target)>(source))?(source):(target)):(source)))
#endif

#ifndef SETMAX
#define SETMAX(target,source,index) \
 ((target)=((index)?(((target)<(source))?(source):(target)):(source)))
#endif

#ifndef SETMAXINDEX
#define SETMAXINDEX(target,source,targetIndex,index) \
 ((!(index) || ((target)<(source)))?(((targetIndex)=(index)),((target)=(source))):(target))
#endif

#ifndef SETMINMAX
#define SETMINMAX(min,max,source,index) \
 (SETMIN(min,source,index),SETMAX(max,source,index))
#endif


#ifndef ABS
#define ABS(x) (((x)<0)?(-(x)):(x))
#endif

#ifndef SQRNPS
#define  SQRNPS(x)  ((x) * (x))
#endif

#ifndef PI
#define PI 3.1415926535897932385
#endif

#ifndef DEG_TO_RAD
#define DEG_TO_RAD 0.0174532925199431
#endif

#ifndef RAD_TO_DEG
#define RAD_TO_DEG 57.29577951308232
#endif

/*
// common typedefs
*/
typedef char* charPtr;
typedef float* floatPtr;
typedef double* doublePtr;
typedef int* intPtr;

#endif
