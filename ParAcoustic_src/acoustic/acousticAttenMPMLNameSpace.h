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
 *  acousticMPMLNameSpace.h
 *  
 *
 *  This file contains the updating kernels that actually form the core
 *  of the code for the MPML boundary.  Both O2 and O4 updating for the 
 *  MPML zones and O4 for the interior zone.  These are called from 
 *  updateAcousticMPMLFull.cc
 *
 *  The first section includes the varibles that are actually defined 
 *  in acousticMPMLNameSpace.cc.
 *
 *  Defines several inline functions used in updateAcousticMPMLFull.cc
 *  or internally:
 *    sqr
 *    reflToD
 *    updatePMLMemVar
 *    updatePMLOne
 *    updatePMLMemVarComp
 *    updateAcousticVx
 *    updateAcousticVxO2
 *    updateAcousticVy
 *    updateAcousticVyO2
 *    updateAcousticVz
 *    updateAcousticVzO2
 *    updateAcousticVxNoPML
 *    updateAcousticVyNoPML
 *    updateAcousticVzNoPML
 *    updateAcousticPressureNoPMLCore
 *    updateAcousticPressureNoPMLCoreO2
 *    updateAcousticPressureNoPML
 *    updateAcousticPressureNoPMLO2
 *    updateAcousticPressurePMLCore
 *    updateAcousticPressurePMLCoreO2
 *    updateAcousticPressurePML
 *    updateAcousticPressurePMLO2
 */
#ifndef _ACOUSTIC_ATTEN_MPML_NS
#define _ACOUSTIC_ATTEN_MPML_NS

#include<math.h>

#include<stdio.h>

namespace acousticAtten_MPML {
  extern float *__restrict__ kxTaper, *__restrict__ kyTaper, *__restrict__ kzTaper;
  extern float *__restrict__ alphaXTaper, *__restrict__ alphaYTaper, *__restrict__ alphaZTaper;
  extern float *__restrict__ xTaper, *__restrict__ yTaper, *__restrict__ zTaper;
  extern float *__restrict__ axTaper, *__restrict__ ayTaper, *__restrict__ azTaper;
  extern float *__restrict__ bxTaper, *__restrict__ byTaper, *__restrict__ bzTaper;
  extern float *__restrict__ kxTaperPH, *__restrict__ kyTaperPH, *__restrict__ kzTaperPH;
  extern float *__restrict__ alphaXTaperPH, *__restrict__ alphaYTaperPH, *__restrict__ alphaZTaperPH;
  extern float *__restrict__ xTaperPH, *__restrict__ yTaperPH, *__restrict__ zTaperPH;
  extern float *__restrict__ axTaperPH, *__restrict__ ayTaperPH, *__restrict__ azTaperPH;
  extern float *__restrict__ bxTaperPH, *__restrict__ byTaperPH, *__restrict__ bzTaperPH;
  extern float *__restrict__ xTaperRel, *__restrict__ yTaperRel;
  extern float *__restrict__ xTaperRelPH, *__restrict__ yTaperRelPH;
  
  extern float *__restrict__ vxxXmin, *__restrict__ vxxXmax, *__restrict__ vxxYmin, *__restrict__ vxxYmax, *__restrict__ vxxZmin, *__restrict__ vxxZmax;
  extern float *__restrict__ vyyXmin, *__restrict__ vyyXmax, *__restrict__ vyyYmin, *__restrict__ vyyYmax, *__restrict__ vyyZmin, *__restrict__ vyyZmax;
  extern float *__restrict__ vzzXmin, *__restrict__ vzzXmax, *__restrict__ vzzYmin, *__restrict__ vzzYmax, *__restrict__ vzzZmin, *__restrict__ vzzZmax;
  
  extern float *__restrict__ pxXmin, *__restrict__ pxXmax, *__restrict__ pxYmin, *__restrict__ pxYmax, *__restrict__ pxZmin, *__restrict__ pxZmax;
  extern float *__restrict__ pyXmin, *__restrict__ pyXmax, *__restrict__ pyYmin, *__restrict__ pyYmax, *__restrict__ pyZmin, *__restrict__ pyZmax;
  extern float *__restrict__ pzXmin, *__restrict__ pzXmax, *__restrict__ pzYmin, *__restrict__ pzYmax, *__restrict__ pzZmin, *__restrict__ pzZmax;

extern float *__restrict__*__restrict__ axE110TaperxPH;
extern float *__restrict__*__restrict__ bxE110TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ azC000Taper;
extern float *__restrict__*__restrict__*__restrict__ bzC000Taper;
extern float *__restrict__*__restrict__*__restrict__ ayC001TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC001TaperzPHxPH;
extern float *__restrict__*__restrict__ azE102TaperxPH;
extern float *__restrict__*__restrict__ bzE102TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ azC010TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC010TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ axC010TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC010TaperyPH;
extern float *__restrict__*__restrict__ azE120TaperxPH;
extern float *__restrict__*__restrict__ bzE120TaperxPH;
extern float *__restrict__*__restrict__ axE022Taper;
extern float *__restrict__*__restrict__ bxE022Taper;
extern float *__restrict__*__restrict__ azE202Taper;
extern float *__restrict__*__restrict__ bzE202Taper;
extern float *__restrict__*__restrict__ ayE022TaperzPH;
extern float *__restrict__*__restrict__ byE022TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ ayC110TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC110TaperzPHxPH;
extern float *__restrict__*__restrict__ axE220TaperxPH;
extern float *__restrict__*__restrict__ bxE220TaperxPH;
extern float *__restrict__*__restrict__ ayE022Taper;
extern float *__restrict__*__restrict__ byE022Taper;
extern float *__restrict__*__restrict__*__restrict__ ayC000TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ byC000TaperzPH;
extern float *__restrict__*__restrict__ ayE102TaperxPH;
extern float *__restrict__*__restrict__ byE102TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC111Taper;
extern float *__restrict__*__restrict__*__restrict__ byC111Taper;
extern float *__restrict__*__restrict__ axE011Taper;
extern float *__restrict__*__restrict__ bxE011Taper;
extern float *__restrict__*__restrict__ azE012TaperzPH;
extern float *__restrict__*__restrict__ bzE012TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ ayC011TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ byC011TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ ayC010Taper;
extern float *__restrict__*__restrict__*__restrict__ byC010Taper;
extern float *__restrict__*__restrict__*__restrict__ ayC011TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC011TaperyPHxPH;
extern float *__restrict__ ayX0TaperPH;
extern float *__restrict__ byX0TaperPH;
extern float *__restrict__*__restrict__*__restrict__ axC100Taper;
extern float *__restrict__*__restrict__*__restrict__ bxC100Taper;
extern float *__restrict__*__restrict__*__restrict__ azC000TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC000TaperyPHxPH;
extern float *__restrict__*__restrict__ ayE201TaperxPH;
extern float *__restrict__*__restrict__ byE201TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ axC101TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC101TaperyPH;
extern float *__restrict__*__restrict__ axE120TaperyPHxPH;
extern float *__restrict__*__restrict__ bxE120TaperyPHxPH;
extern float *__restrict__*__restrict__ ayE120TaperyPHxPH;
extern float *__restrict__*__restrict__ byE120TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC110TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ byC110TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ azC101TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC101TaperyPH;
extern float *__restrict__*__restrict__ axE220Taper;
extern float *__restrict__*__restrict__ bxE220Taper;
extern float *__restrict__*__restrict__*__restrict__ azC011TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC011TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC110TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC110TaperyPHxPH;
extern float *__restrict__*__restrict__ ayE202Taper;
extern float *__restrict__*__restrict__ byE202Taper;
extern float *__restrict__*__restrict__*__restrict__ axC101TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC101TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ azC000TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bzC000TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ azC010Taper;
extern float *__restrict__*__restrict__*__restrict__ bzC010Taper;
extern float *__restrict__*__restrict__*__restrict__ azC000TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC000TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ axC110TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC110TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ axC110TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC110TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ ayC110TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ byC110TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ ayC100TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ byC100TaperxPH;
extern float *__restrict__ ayZ0Taper;
extern float *__restrict__ byZ0Taper;
extern float *__restrict__*__restrict__*__restrict__ ayC010TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC010TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ azC010TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC010TaperzPHxPH;
extern float *__restrict__*__restrict__ ayE220TaperyPH;
extern float *__restrict__*__restrict__ byE220TaperyPH;
extern float *__restrict__*__restrict__ azE021Taper;
extern float *__restrict__*__restrict__ bzE021Taper;
extern float *__restrict__*__restrict__ ayE101TaperxPH;
extern float *__restrict__*__restrict__ byE101TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC111TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ byC111TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ axC011TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC011TaperxPH;
extern float *__restrict__*__restrict__ azE011Taper;
extern float *__restrict__*__restrict__ bzE011Taper;
extern float *__restrict__*__restrict__ axE210TaperyPH;
extern float *__restrict__*__restrict__ bxE210TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ axC001Taper;
extern float *__restrict__*__restrict__*__restrict__ bxC001Taper;
extern float *__restrict__*__restrict__*__restrict__ ayC101Taper;
extern float *__restrict__*__restrict__*__restrict__ byC101Taper;
extern float *__restrict__*__restrict__ ayE110TaperyPHxPH;
extern float *__restrict__*__restrict__ byE110TaperyPHxPH;
extern float *__restrict__ axY0TaperPH;
extern float *__restrict__ bxY0TaperPH;
extern float *__restrict__*__restrict__*__restrict__ azC010TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC010TaperyPH;
extern float *__restrict__*__restrict__ axE202TaperxPH;
extern float *__restrict__*__restrict__ bxE202TaperxPH;
extern float *__restrict__*__restrict__ azE022TaperzPH;
extern float *__restrict__*__restrict__ bzE022TaperzPH;
extern float *__restrict__*__restrict__ axE012TaperzPH;
extern float *__restrict__*__restrict__ bxE012TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ axC111TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC111TaperzPHyPH;
extern float *__restrict__*__restrict__ azE110TaperyPHxPH;
extern float *__restrict__*__restrict__ bzE110TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ azC110TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC110TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ axC000TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC000TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ ayC101TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC101TaperyPHxPH;
extern float *__restrict__*__restrict__ azE201TaperzPHxPH;
extern float *__restrict__*__restrict__ bzE201TaperzPHxPH;
extern float *__restrict__*__restrict__ axE202TaperzPH;
extern float *__restrict__*__restrict__ bxE202TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ azC110TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC110TaperyPHxPH;
extern float *__restrict__ azY1Taper;
extern float *__restrict__ bzY1Taper;
extern float *__restrict__*__restrict__*__restrict__ axC010TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC010TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ axC101TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bxC101TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ ayC101TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC101TaperzPHxPH;
extern float *__restrict__*__restrict__ axE012TaperyPH;
extern float *__restrict__*__restrict__ bxE012TaperyPH;
extern float *__restrict__*__restrict__ axE220TaperyPH;
extern float *__restrict__*__restrict__ bxE220TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ ayC001TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ byC001TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ axC001TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bxC001TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ axC100TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC100TaperzPHyPH;
extern float *__restrict__*__restrict__ axE012TaperzPHyPH;
extern float *__restrict__*__restrict__ bxE012TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ axC000TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC000TaperzPHxPH;
extern float *__restrict__*__restrict__ ayE012TaperzPH;
extern float *__restrict__*__restrict__ byE012TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ axC010TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bxC010TaperzPH;
extern float *__restrict__*__restrict__ ayE110TaperyPH;
extern float *__restrict__*__restrict__ byE110TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ axC110TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC110TaperxPH;
extern float *__restrict__*__restrict__ ayE202TaperzPHxPH;
extern float *__restrict__*__restrict__ byE202TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ axC000TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC000TaperxPH;
extern float *__restrict__*__restrict__ ayE022TaperyPH;
extern float *__restrict__*__restrict__ byE022TaperyPH;
extern float *__restrict__*__restrict__ ayE220TaperxPH;
extern float *__restrict__*__restrict__ byE220TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC011TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC011TaperzPHxPH;
extern float *__restrict__*__restrict__ ayE210TaperyPHxPH;
extern float *__restrict__*__restrict__ byE210TaperyPHxPH;
extern float *__restrict__*__restrict__ ayE110TaperxPH;
extern float *__restrict__*__restrict__ byE110TaperxPH;
extern float *__restrict__*__restrict__ ayE120TaperyPH;
extern float *__restrict__*__restrict__ byE120TaperyPH;
extern float *__restrict__*__restrict__ ayE021TaperyPH;
extern float *__restrict__*__restrict__ byE021TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ azC111TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC111TaperzPHyPH;
extern float *__restrict__ azX0Taper;
extern float *__restrict__ bzX0Taper;
extern float *__restrict__ ayX0Taper;
extern float *__restrict__ byX0Taper;
extern float *__restrict__*__restrict__*__restrict__ azC100Taper;
extern float *__restrict__*__restrict__*__restrict__ bzC100Taper;
extern float *__restrict__*__restrict__*__restrict__ azC001TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC001TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ axC001TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC001TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ azC111TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC111TaperzPHxPH;
extern float *__restrict__*__restrict__ ayE220Taper;
extern float *__restrict__*__restrict__ byE220Taper;
extern float *__restrict__*__restrict__ azE210Taper;
extern float *__restrict__*__restrict__ bzE210Taper;
extern float *__restrict__ axZ0Taper;
extern float *__restrict__ bxZ0Taper;
extern float *__restrict__*__restrict__ azE101TaperzPHxPH;
extern float *__restrict__*__restrict__ bzE101TaperzPHxPH;
extern float *__restrict__ axZ1TaperPH;
extern float *__restrict__ bxZ1TaperPH;
extern float *__restrict__*__restrict__ ayE102TaperzPHxPH;
extern float *__restrict__*__restrict__ byE102TaperzPHxPH;
extern float *__restrict__*__restrict__ azE120TaperyPHxPH;
extern float *__restrict__*__restrict__ bzE120TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ azC000TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC000TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ axC111TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC111TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ azC010TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC010TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ azC111TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bzC111TaperzPH;
extern float *__restrict__*__restrict__ azE012TaperyPH;
extern float *__restrict__*__restrict__ bzE012TaperyPH;
extern float *__restrict__*__restrict__ axE011TaperzPHyPH;
extern float *__restrict__*__restrict__ bxE011TaperzPHyPH;
extern float *__restrict__*__restrict__ ayE102TaperzPH;
extern float *__restrict__*__restrict__ byE102TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ axC010TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC010TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ azC000TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC000TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ axC000TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC000TaperyPHxPH;
extern float *__restrict__*__restrict__ axE101TaperzPH;
extern float *__restrict__*__restrict__ bxE101TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ azC111Taper;
extern float *__restrict__*__restrict__*__restrict__ bzC111Taper;
extern float *__restrict__*__restrict__ azE022TaperyPH;
extern float *__restrict__*__restrict__ bzE022TaperyPH;
extern float *__restrict__*__restrict__ azE102TaperzPH;
extern float *__restrict__*__restrict__ bzE102TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ axC100TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bxC100TaperzPH;
extern float *__restrict__ axY1Taper;
extern float *__restrict__ bxY1Taper;
extern float *__restrict__ ayZ1Taper;
extern float *__restrict__ byZ1Taper;
extern float *__restrict__*__restrict__*__restrict__ azC001TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC001TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ axC001TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC001TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ axC110Taper;
extern float *__restrict__*__restrict__*__restrict__ bxC110Taper;
extern float *__restrict__*__restrict__*__restrict__ azC101TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC101TaperyPHxPH;
extern float *__restrict__*__restrict__ azE101Taper;
extern float *__restrict__*__restrict__ bzE101Taper;
extern float *__restrict__*__restrict__*__restrict__ ayC000Taper;
extern float *__restrict__*__restrict__*__restrict__ byC000Taper;
extern float *__restrict__*__restrict__*__restrict__ azC101Taper;
extern float *__restrict__*__restrict__*__restrict__ bzC101Taper;
extern float *__restrict__*__restrict__*__restrict__ ayC011TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ byC011TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ axC111TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC111TaperxPH;
extern float *__restrict__*__restrict__ axE201TaperzPHxPH;
extern float *__restrict__*__restrict__ bxE201TaperzPHxPH;
extern float *__restrict__*__restrict__ azE220TaperyPH;
extern float *__restrict__*__restrict__ bzE220TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ ayC100Taper;
extern float *__restrict__*__restrict__*__restrict__ byC100Taper;
extern float *__restrict__*__restrict__*__restrict__ azC001TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC001TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ axC111TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bxC111TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ ayC110TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ byC110TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ azC110TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC110TaperxPH;
extern float *__restrict__ ayZ0TaperPH;
extern float *__restrict__ byZ0TaperPH;
extern float *__restrict__*__restrict__*__restrict__ azC101TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC101TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ axC110TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC110TaperyPHxPH;
extern float *__restrict__*__restrict__ axE110TaperyPH;
extern float *__restrict__*__restrict__ bxE110TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ axC000TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bxC000TaperzPH;
extern float *__restrict__*__restrict__ azE202TaperzPHxPH;
extern float *__restrict__*__restrict__ bzE202TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC010TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ byC010TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ axC100TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC100TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ azC001TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC001TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC111TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ byC111TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ ayC001TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ byC001TaperzPHyPH;
extern float *__restrict__*__restrict__ azE021TaperyPH;
extern float *__restrict__*__restrict__ bzE021TaperyPH;
extern float *__restrict__*__restrict__ ayE120Taper;
extern float *__restrict__*__restrict__ byE120Taper;
extern float *__restrict__*__restrict__*__restrict__ azC101TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC101TaperzPHyPH;
extern float *__restrict__*__restrict__ axE021Taper;
extern float *__restrict__*__restrict__ bxE021Taper;
extern float *__restrict__*__restrict__ axE102Taper;
extern float *__restrict__*__restrict__ bxE102Taper;
extern float *__restrict__*__restrict__*__restrict__ axC001TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC001TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ azC001TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC001TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ azC100TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC100TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ ayC010TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ byC010TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC100TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ byC100TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ ayC000TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ byC000TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ azC111TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC111TaperyPH;
extern float *__restrict__*__restrict__ axE101Taper;
extern float *__restrict__*__restrict__ bxE101Taper;
extern float *__restrict__ azX0TaperPH;
extern float *__restrict__ bzX0TaperPH;
extern float *__restrict__*__restrict__ ayE110Taper;
extern float *__restrict__*__restrict__ byE110Taper;
extern float *__restrict__*__restrict__ axE012Taper;
extern float *__restrict__*__restrict__ bxE012Taper;
extern float *__restrict__*__restrict__*__restrict__ azC101TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bzC101TaperzPH;
extern float *__restrict__*__restrict__ ayE101TaperzPHxPH;
extern float *__restrict__*__restrict__ byE101TaperzPHxPH;
extern float *__restrict__ ayX1TaperPH;
extern float *__restrict__ byX1TaperPH;
extern float *__restrict__*__restrict__ ayE022TaperzPHyPH;
extern float *__restrict__*__restrict__ byE022TaperzPHyPH;
extern float *__restrict__*__restrict__ ayE202TaperxPH;
extern float *__restrict__*__restrict__ byE202TaperxPH;
extern float *__restrict__*__restrict__ axE110TaperyPHxPH;
extern float *__restrict__*__restrict__ bxE110TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ azC010TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bzC010TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ ayC111TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC111TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC001TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC001TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ azC110TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bzC110TaperzPH;
extern float *__restrict__*__restrict__ axE120Taper;
extern float *__restrict__*__restrict__ bxE120Taper;
extern float *__restrict__*__restrict__ ayE012TaperyPH;
extern float *__restrict__*__restrict__ byE012TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ azC011Taper;
extern float *__restrict__*__restrict__*__restrict__ bzC011Taper;
extern float *__restrict__*__restrict__ azE012TaperzPHyPH;
extern float *__restrict__*__restrict__ bzE012TaperzPHyPH;
extern float *__restrict__*__restrict__ ayE102Taper;
extern float *__restrict__*__restrict__ byE102Taper;
extern float *__restrict__*__restrict__*__restrict__ ayC111TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC111TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ axC100TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC100TaperyPH;
extern float *__restrict__*__restrict__ azE220Taper;
extern float *__restrict__*__restrict__ bzE220Taper;
extern float *__restrict__*__restrict__ ayE021Taper;
extern float *__restrict__*__restrict__ byE021Taper;
extern float *__restrict__*__restrict__ axE210TaperxPH;
extern float *__restrict__*__restrict__ bxE210TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ axC011Taper;
extern float *__restrict__*__restrict__*__restrict__ bxC011Taper;
extern float *__restrict__ axY0Taper;
extern float *__restrict__ bxY0Taper;
extern float *__restrict__*__restrict__*__restrict__ ayC000TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ byC000TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC100TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC100TaperzPHxPH;
extern float *__restrict__*__restrict__ axE120TaperyPH;
extern float *__restrict__*__restrict__ bxE120TaperyPH;
extern float *__restrict__*__restrict__ azE012Taper;
extern float *__restrict__*__restrict__ bzE012Taper;
extern float *__restrict__*__restrict__ axE021TaperyPH;
extern float *__restrict__*__restrict__ bxE021TaperyPH;
extern float *__restrict__*__restrict__ azE022TaperzPHyPH;
extern float *__restrict__*__restrict__ bzE022TaperzPHyPH;
extern float *__restrict__ azY0Taper;
extern float *__restrict__ bzY0Taper;
extern float *__restrict__*__restrict__ azE220TaperxPH;
extern float *__restrict__*__restrict__ bzE220TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ axC000Taper;
extern float *__restrict__*__restrict__*__restrict__ bxC000Taper;
extern float *__restrict__*__restrict__ axE022TaperzPH;
extern float *__restrict__*__restrict__ bxE022TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ azC011TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC011TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ axC100TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC100TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ azC100TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC100TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC101TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ byC101TaperxPH;
extern float *__restrict__*__restrict__ ayE201Taper;
extern float *__restrict__*__restrict__ byE201Taper;
extern float *__restrict__*__restrict__ ayE210TaperyPH;
extern float *__restrict__*__restrict__ byE210TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ ayC100TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC100TaperyPHxPH;
extern float *__restrict__*__restrict__ azE201TaperzPH;
extern float *__restrict__*__restrict__ bzE201TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ azC101TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC101TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ axC001TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC001TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ azC110TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC110TaperzPHyPH;
extern float *__restrict__*__restrict__ axE102TaperxPH;
extern float *__restrict__*__restrict__ bxE102TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC000TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC000TaperyPHxPH;
extern float *__restrict__*__restrict__ azE210TaperyPH;
extern float *__restrict__*__restrict__ bzE210TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ axC011TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC011TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ ayC000TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ byC000TaperzPHyPH;
extern float *__restrict__*__restrict__ ayE012Taper;
extern float *__restrict__*__restrict__ byE012Taper;
extern float *__restrict__*__restrict__*__restrict__ azC000TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC000TaperzPHyPH;
extern float *__restrict__*__restrict__ azE210TaperxPH;
extern float *__restrict__*__restrict__ bzE210TaperxPH;
extern float *__restrict__*__restrict__ azE021TaperzPH;
extern float *__restrict__*__restrict__ bzE021TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ ayC100TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ byC100TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ axC100TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC100TaperxPH;
extern float *__restrict__*__restrict__ ayE210TaperxPH;
extern float *__restrict__*__restrict__ byE210TaperxPH;
extern float *__restrict__*__restrict__ ayE201TaperzPH;
extern float *__restrict__*__restrict__ byE201TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ ayC111TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ byC111TaperzPHyPH;
extern float *__restrict__*__restrict__ azE021TaperzPHyPH;
extern float *__restrict__*__restrict__ bzE021TaperzPHyPH;
extern float *__restrict__*__restrict__ axE110Taper;
extern float *__restrict__*__restrict__ bxE110Taper;
extern float *__restrict__*__restrict__*__restrict__ axC011TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC011TaperzPHxPH;
extern float *__restrict__*__restrict__ axE021TaperzPHyPH;
extern float *__restrict__*__restrict__ bxE021TaperzPHyPH;
extern float *__restrict__*__restrict__ azE011TaperzPH;
extern float *__restrict__*__restrict__ bzE011TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ ayC111TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ byC111TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ axC101TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC101TaperxPH;
extern float *__restrict__ ayX1Taper;
extern float *__restrict__ byX1Taper;
extern float *__restrict__*__restrict__ ayE011TaperzPH;
extern float *__restrict__*__restrict__ byE011TaperzPH;
extern float *__restrict__*__restrict__ azE220TaperyPHxPH;
extern float *__restrict__*__restrict__ bzE220TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC101TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ byC101TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ ayC101TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ byC101TaperyPH;
extern float *__restrict__*__restrict__ ayE220TaperyPHxPH;
extern float *__restrict__*__restrict__ byE220TaperyPHxPH;
extern float *__restrict__ azY1TaperPH;
extern float *__restrict__ bzY1TaperPH;
extern float *__restrict__*__restrict__ azE102Taper;
extern float *__restrict__*__restrict__ bzE102Taper;
extern float *__restrict__*__restrict__*__restrict__ ayC001TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ byC001TaperxPH;
extern float *__restrict__*__restrict__ ayE011Taper;
extern float *__restrict__*__restrict__ byE011Taper;
extern float *__restrict__*__restrict__ axE120TaperxPH;
extern float *__restrict__*__restrict__ bxE120TaperxPH;
extern float *__restrict__ azY0TaperPH;
extern float *__restrict__ bzY0TaperPH;
extern float *__restrict__*__restrict__*__restrict__ axC000TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC000TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ axC110TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bxC110TaperzPH;
extern float *__restrict__*__restrict__ axE201Taper;
extern float *__restrict__*__restrict__ bxE201Taper;
extern float *__restrict__*__restrict__*__restrict__ ayC110TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ byC110TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ azC001TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bzC001TaperzPH;
extern float *__restrict__ azX1Taper;
extern float *__restrict__ bzX1Taper;
extern float *__restrict__*__restrict__*__restrict__ ayC011Taper;
extern float *__restrict__*__restrict__*__restrict__ byC011Taper;
extern float *__restrict__*__restrict__ ayE201TaperzPHxPH;
extern float *__restrict__*__restrict__ byE201TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ axC101TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC101TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC110Taper;
extern float *__restrict__*__restrict__*__restrict__ byC110Taper;
extern float *__restrict__*__restrict__*__restrict__ azC011TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bzC011TaperzPH;
extern float *__restrict__*__restrict__ azE210TaperyPHxPH;
extern float *__restrict__*__restrict__ bzE210TaperyPHxPH;
extern float *__restrict__*__restrict__ azE110TaperyPH;
extern float *__restrict__*__restrict__ bzE110TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ ayC001Taper;
extern float *__restrict__*__restrict__*__restrict__ byC001Taper;
extern float *__restrict__*__restrict__ azE202TaperzPH;
extern float *__restrict__*__restrict__ bzE202TaperzPH;
extern float *__restrict__ ayZ1TaperPH;
extern float *__restrict__ byZ1TaperPH;
extern float *__restrict__*__restrict__ axE201TaperxPH;
extern float *__restrict__*__restrict__ bxE201TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ azC010TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC010TaperyPHxPH;
extern float *__restrict__ axY1TaperPH;
extern float *__restrict__ bxY1TaperPH;
extern float *__restrict__*__restrict__*__restrict__ axC010Taper;
extern float *__restrict__*__restrict__*__restrict__ bxC010Taper;
extern float *__restrict__*__restrict__ azE110TaperxPH;
extern float *__restrict__*__restrict__ bzE110TaperxPH;
extern float *__restrict__*__restrict__ axE011TaperyPH;
extern float *__restrict__*__restrict__ bxE011TaperyPH;
extern float *__restrict__*__restrict__ ayE101TaperzPH;
extern float *__restrict__*__restrict__ byE101TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ axC010TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC010TaperyPHxPH;
extern float *__restrict__*__restrict__ ayE202TaperzPH;
extern float *__restrict__*__restrict__ byE202TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ azC011TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC011TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ axC001TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC001TaperyPH;
extern float *__restrict__*__restrict__ azE011TaperyPH;
extern float *__restrict__*__restrict__ bzE011TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ axC011TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bxC011TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ azC111TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC111TaperxPH;
extern float *__restrict__*__restrict__ azE102TaperzPHxPH;
extern float *__restrict__*__restrict__ bzE102TaperzPHxPH;
extern float *__restrict__*__restrict__ azE120TaperyPH;
extern float *__restrict__*__restrict__ bzE120TaperyPH;
extern float *__restrict__*__restrict__ axE220TaperyPHxPH;
extern float *__restrict__*__restrict__ bxE220TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ axC011TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC011TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ ayC101TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ byC101TaperzPHyPH;
extern float *__restrict__*__restrict__ axE102TaperzPHxPH;
extern float *__restrict__*__restrict__ bxE102TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ azC011TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC011TaperzPHyPH;
extern float *__restrict__*__restrict__ azE201TaperxPH;
extern float *__restrict__*__restrict__ bzE201TaperxPH;
extern float *__restrict__*__restrict__ axE202TaperzPHxPH;
extern float *__restrict__*__restrict__ bxE202TaperzPHxPH;
extern float *__restrict__ azX1TaperPH;
extern float *__restrict__ bzX1TaperPH;
extern float *__restrict__*__restrict__ axE021TaperzPH;
extern float *__restrict__*__restrict__ bxE021TaperzPH;
extern float *__restrict__*__restrict__ axE101TaperxPH;
extern float *__restrict__*__restrict__ bxE101TaperxPH;
extern float *__restrict__*__restrict__ azE120Taper;
extern float *__restrict__*__restrict__ bzE120Taper;
extern float *__restrict__*__restrict__ azE101TaperzPH;
extern float *__restrict__*__restrict__ bzE101TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ azC100TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ bzC100TaperzPHyPH;
extern float *__restrict__*__restrict__ azE110Taper;
extern float *__restrict__*__restrict__ bzE110Taper;
extern float *__restrict__*__restrict__ ayE210Taper;
extern float *__restrict__*__restrict__ byE210Taper;
extern float *__restrict__*__restrict__*__restrict__ ayC011TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ byC011TaperzPH;
extern float *__restrict__*__restrict__ azE201Taper;
extern float *__restrict__*__restrict__ bzE201Taper;
extern float *__restrict__*__restrict__*__restrict__ axC110TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC110TaperyPH;
extern float *__restrict__*__restrict__ axE022TaperzPHyPH;
extern float *__restrict__*__restrict__ bxE022TaperzPHyPH;
extern float *__restrict__*__restrict__ ayE101Taper;
extern float *__restrict__*__restrict__ byE101Taper;
extern float *__restrict__*__restrict__ axE101TaperzPHxPH;
extern float *__restrict__*__restrict__ bxE101TaperzPHxPH;
extern float *__restrict__*__restrict__ azE011TaperzPHyPH;
extern float *__restrict__*__restrict__ bzE011TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ axC111Taper;
extern float *__restrict__*__restrict__*__restrict__ bxC111Taper;
extern float *__restrict__*__restrict__*__restrict__ azC110TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC110TaperzPHxPH;
extern float *__restrict__*__restrict__ ayE012TaperzPHyPH;
extern float *__restrict__*__restrict__ byE012TaperzPHyPH;
extern float *__restrict__*__restrict__ ayE011TaperyPH;
extern float *__restrict__*__restrict__ byE011TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ ayC100TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ byC100TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ axC101Taper;
extern float *__restrict__*__restrict__*__restrict__ bxC101Taper;
extern float *__restrict__ axZ1Taper;
extern float *__restrict__ bxZ1Taper;
extern float *__restrict__*__restrict__ azE022Taper;
extern float *__restrict__*__restrict__ bzE022Taper;
extern float *__restrict__*__restrict__*__restrict__ azC011TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC011TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC000TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC000TaperzPHxPH;
extern float *__restrict__*__restrict__ axE201TaperzPH;
extern float *__restrict__*__restrict__ bxE201TaperzPH;
extern float *__restrict__*__restrict__ ayE120TaperxPH;
extern float *__restrict__*__restrict__ byE120TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ azC001Taper;
extern float *__restrict__*__restrict__*__restrict__ bzC001Taper;
extern float *__restrict__*__restrict__*__restrict__ axC010TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC010TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ axC011TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC011TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC001TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ byC001TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ azC111TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC111TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC011TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ byC011TaperxPH;
extern float *__restrict__*__restrict__ axE210TaperyPHxPH;
extern float *__restrict__*__restrict__ bxE210TaperyPHxPH;
extern float *__restrict__*__restrict__ axE210Taper;
extern float *__restrict__*__restrict__ bxE210Taper;
extern float *__restrict__*__restrict__*__restrict__ axC111TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC111TaperzPHxPH;
extern float *__restrict__*__restrict__ ayE021TaperzPH;
extern float *__restrict__*__restrict__ byE021TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ azC100TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ bzC100TaperzPH;
extern float *__restrict__*__restrict__ ayE011TaperzPHyPH;
extern float *__restrict__*__restrict__ byE011TaperzPHyPH;
extern float *__restrict__*__restrict__*__restrict__ ayC010TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ byC010TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ axC111TaperyPH;
extern float *__restrict__*__restrict__*__restrict__ bxC111TaperyPH;
extern float *__restrict__*__restrict__ ayE021TaperzPHyPH;
extern float *__restrict__*__restrict__ byE021TaperzPHyPH;
extern float *__restrict__*__restrict__ axE202Taper;
extern float *__restrict__*__restrict__ bxE202Taper;
extern float *__restrict__*__restrict__ azE101TaperxPH;
extern float *__restrict__*__restrict__ bzE101TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC010TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ byC010TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ ayC010TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ byC010TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ azC100TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC100TaperxPH;
extern float *__restrict__*__restrict__*__restrict__ azC110Taper;
extern float *__restrict__*__restrict__*__restrict__ bzC110Taper;
extern float *__restrict__*__restrict__ axE022TaperyPH;
extern float *__restrict__*__restrict__ bxE022TaperyPH;
extern float *__restrict__*__restrict__ axE102TaperzPH;
extern float *__restrict__*__restrict__ bxE102TaperzPH;
extern float *__restrict__*__restrict__ axE011TaperzPH;
extern float *__restrict__*__restrict__ bxE011TaperzPH;
extern float *__restrict__*__restrict__*__restrict__ axC101TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bxC101TaperzPHxPH;
extern float *__restrict__*__restrict__*__restrict__ azC100TaperyPHxPH;
extern float *__restrict__*__restrict__*__restrict__ bzC100TaperyPHxPH;
extern float *__restrict__*__restrict__ azE202TaperxPH;
extern float *__restrict__*__restrict__ bzE202TaperxPH;
extern float *__restrict__ axZ0TaperPH;
extern float *__restrict__ bxZ0TaperPH;

  extern float *__restrict__ tcr, *__restrict__ tap;
  extern float *__restrict__ tdvx, *__restrict__ tdvy, *__restrict__ tdvz;
  
  extern int nyInt, nzInt;
  
  extern int xMinSize, xMaxSize, yMinSize, yMaxSize, zMinSize, zMaxSize;
  extern int _nXmin, _nXmax, _nYmin, _nYmax, _nZmin, _nZmax;
  extern int _xMaxStart, _yMaxStart, _zMaxStart;
  extern int xStart1, xStart2, yStart1, yStart2, zStart, zStart_z;
  extern int xEnd, yEnd, zEnd, xEnd1, yEnd1, zEnd1;
  extern int _surfaceBCMode;
  extern int _kstart;
  
  //calculates the square of a number, saving some typing
  inline float sqr(float x) {return x*x;}
  //this function converts the input desired reflection coefficient (r)
  //to the sigma at the flank of the model using the theoretical
  //reflection coefficient for MPMLs with the width (w, in meters) and
  //the maximum vp in the model anywhere (vMax)
  inline float reflToD(float r, float w, float vMax) {
    return -1.5f*logf(r)*vMax/w;
  }
  
  //This returns the new value of MPML memory variable at this time step
  //using the previous memory varibale value (prev), the MPML parameters
  //at this point (ataper and btaper) and the spatial derivative of the
  //corresponding dependent variable (derivative)
  inline float updatePMLMemVar(float prev,float ataper,float btaper,float derivative) {
    return btaper*prev+ataper*derivative;
  }
  
  //Acutally update the memory variable (xxx) and velocities for velocity
  //components (vx) using MPML parameters axtaper, bxtaper, and kxtaper
  //at index in the vx array (which corresponds to pindex in the memory
  //variable array and the density (rho) and spatial derivatives of
  //the appropriate stresses (derX)
  inline void updatePMLOne(float*__restrict__ vx, int index,float*__restrict__ xxx, int pindex,float axtaper,float bxtaper, 
                           float kxtaper,float derX, float rho) {
    xxx[pindex]=updatePMLMemVar(xxx[pindex],axtaper,bxtaper,derX);
    vx[index] -= (derX*kxtaper+xxx[pindex])/rho;  //note: kxtaper is one over the original
  }
  
  //Actually update the memory variables (xxx, xyy and xzz) at pindex
  //using MPML parameters axtaper, aytaper, aztaper, etc. and spatial
  //derivatives of velocties (derX, derY, derZ).  Used for pressure
  //updates.
  inline void updatePMLMemVarComp(float*__restrict__ xxx, float*__restrict__ xyy,float*__restrict__ xzz,int pindex,float axtaper,float aytaper, float aztaper, float bxtaper, float bytaper,
                                  float bztaper, float derX, float derY, float derZ) {
    xxx[pindex]=updatePMLMemVar(xxx[pindex],axtaper,bxtaper,derX);
    xyy[pindex]=updatePMLMemVar(xyy[pindex],aytaper,bytaper,derY);
    xzz[pindex]=updatePMLMemVar(xzz[pindex],aztaper,bztaper,derZ);
  }    
  
  //Update the vx for O4 MPML.  The c variables refer to O4 FD coefficients
  //while the d variables refer to O2 coefficients.  High contrast
  //interfaces trigger the change from O4 to O2.  This is determined
  //in acoustic_control.cc and a flag array vxfunc holds this info.  px
  //is the MPML memory variable and ii, jj and kk refer to the 3-D indices
  //within vx
  inline void updateAcousticVx(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float cx0, float cx1, float dcx0,
                               float cy0, float cy1, float dcy0, float cz0, float cz1, float dcz0, unsigned char*__restrict__ vxfunc, int pIndex, int index,
                               float*__restrict__ px, float*__restrict__ vx, float axtap, float aytap,
                               float aztap,float bxtap,float bytap,float bztap,int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+1]);
    float dfx;
    unsigned char dd = vxfunc[index];
    unsigned char dd1 = 1-dd;
    float cx0d = cx0*dd1+dcx0*dd;
    float cx1d = cx1*dd1;
    dfx = cx0d*(pressure[index+1]-pressure[index   ]) +
    cx1d*(pressure[index+2]-pressure[index-1]);
    updatePMLOne(vx,index,px,pIndex,axtap,bxtap,
                 kxTaperPH[ii],dfx,rhoave);
  }
  
  //Update the vx for O2 MPML.  This is only used near the very flanks of
  //the model
  inline void updateAcousticVxO2(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float dcx0,
                                 float dcy0,float dcz0, int pIndex, int index,
                                 float*__restrict__ px, float*__restrict__ vx, float axtap, float aytap,
                                 float aztap,float bxtap,float bytap,float bztap,int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+1]);
    float dfx;
    dfx = dcx0*(pressure[index+1]-pressure[index   ]);
    updatePMLOne(vx,index,px,pIndex,axtap,bxtap,
                 kxTaperPH[ii],dfx,rhoave);
  }
  
  //Update vy for O4 MPML.
  inline void updateAcousticVy(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float cx0, float cx1, float dcx0,
                               float cy0, float cy1, float dcy0, float cz0, float cz1, float dcz0, unsigned char*__restrict__ vyfunc,int pIndex, int index,
                               float*__restrict__ py, float*__restrict__ vy, float axtap, float aytap,
                               float aztap,float bxtap,float bytap,float bztap,int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+NX]);
    float dfy;
    unsigned char dd = vyfunc[index];
    unsigned char dd1 = 1-dd;
    float cy0d = cy0*dd1+dcy0*dd;
    float cy1d = cy1*dd1;
    dfy = cy0d*(pressure[index+NX]-pressure[index   ]) +
    cy1d*(pressure[index+2*NX]-pressure[index-NX]);
    updatePMLOne(vy,index,py,pIndex,aytap,bytap,
                 kyTaperPH[jj],dfy,rhoave);
  }
  
  //update vy for O2 MPML.  Only used on the very flanks of the model
  inline void updateAcousticVyO2(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float dcx0,
                                 float dcy0,float dcz0, int pIndex, int index,
                                 float*__restrict__ py, float*__restrict__ vy, float axtap, float aytap,
                                 float aztap,float bxtap,float bytap,float bztap,int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+NX]);
    float dfy;
    dfy = dcy0*(pressure[index+NX]-pressure[index   ]);
    updatePMLOne(vy,index,py,pIndex,aytap,bytap,
                 kyTaperPH[jj],dfy,rhoave);
  }
  
  //update vz for O4 MPML.
  inline void updateAcousticVz(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY,
                               float cx0, float cx1, float dcx0,
                               float cy0, float cy1, float dcy0, float cz0, float cz1, float dcz0, unsigned char*__restrict__ vzfunc,int pIndex, int index,
                               float*__restrict__ pz,float*__restrict__ vz, float axtap, float aytap,
                               float aztap,float bxtap,float bytap,float bztap,int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+NXY]);
    float dfz;
    unsigned char dd = vzfunc[index];
    unsigned char dd1 = 1-dd;
    float cz0d = cz0*dd1+dcz0*dd;
    float cz1d = cz1*dd1;
    dfz=cz0d*(pressure[index+NXY]-pressure[index   ]) +
    cz1d*(pressure[index+2*NXY]-pressure[index-NXY]);
    updatePMLOne(vz,index,pz,pIndex,aztap,bztap,
                 kzTaperPH[kk],dfz,rhoave);
  }
  
  //update vz for O2 MPML.  Only used on the very flanks of the model.
  inline void updateAcousticVzO2(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY,
                                 float dcx0,
                                 float dcy0,float dcz0,int pIndex, int index,
                                 float*__restrict__ pz,float*__restrict__ vz, float axtap, float aytap,
                                 float aztap,float bxtap,float bytap,float bztap,int ii, int jj, int kk) {
    float rhoave=(rho[index]+rho[index+NXY]);
    float dfz;
    dfz=dcz0*(pressure[index+NXY]-pressure[index   ]);
    updatePMLOne(vz,index,pz,pIndex,aztap,bztap,
                 kzTaperPH[kk],dfz,rhoave);
  }
  
  //update vx in the interior (no MPML)
  inline void updateAcousticVxNoPML(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float cx0, float cx1, float dcx0,
                                    float cy0, float cy1, float dcy0, float cz0, float cz1, float dcz0, unsigned char*__restrict__ vxfunc, int index,
                                    float*__restrict__ vx) {
    float rhoave=(rho[index]+rho[index+1]);
    unsigned char dd = vxfunc[index];
    unsigned char dd1 = 1-dd;
    float cx0d = cx0*dd1+dcx0*dd;
    float cx1d = cx1*dd1;
    vx[index]-=(cx0d*(pressure[index+1]-pressure[index   ]) +
                cx1d*(pressure[index+2]-pressure[index-1]))/rhoave;
  }
  
  //update vy in the interior (no MPML)
  inline void updateAcousticVyNoPML(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY, float cx0, float cx1, float dcx0,
                                    float cy0, float cy1, float dcy0, float cz0, float cz1, float dcz0, unsigned char*__restrict__ vyfunc, int index,
                                    float*__restrict__ vy) {
    float rhoave=(rho[index]+rho[index+NX]);
    unsigned char dd = vyfunc[index];
    unsigned char dd1 = 1-dd;
    float cy0d = cy0*dd1+dcy0*dd;
    float cy1d = cy1*dd1;
    vy[index]-=(cy0d*(pressure[index+NX]-pressure[index   ]) +
                cy1d*(pressure[index+2*NX]-pressure[index-NX]))/rhoave;
  }
  
  //update vz in the interior (no MPML)
  inline void updateAcousticVzNoPML(float*__restrict__ pressure, float*__restrict__ rho, int NX, int NXY,
                                    float cx0, float cx1, float dcx0,
                                    float cy0, float cy1, float dcy0, float cz0, float cz1, float dcz0, unsigned char*__restrict__ vzfunc,int index,
                                    float*__restrict__ vz) {
    float rhoave=(rho[index]+rho[index+NXY]);
    unsigned char dd = vzfunc[index];
    unsigned char dd1 = 1-dd;
    float cz0d = cz0*dd1+dcz0*dd;
    float cz1d = cz1*dd1;
    vz[index]-=(cz0d*(pressure[index+NXY]-pressure[index   ]) +
                cz1d*(pressure[index+2*NXY]-pressure[index-NXY]))/rhoave;
  }    

  inline void updateAcousticAttenNoPMLRest(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,
                                                 float *__restrict__  currRp,float *__restrict__  bulk,int NX,int NXY,int index, int i,
                                                 float *__restrict__  pressure) {
    float dvx = tdvx[i];
    float dvy = tdvy[i];
    float dvz = tdvz[i];
    float bb = bulk[index];
    
    float dvTot = dvx+dvy+dvz;
    float omega=tcr[i];
#if USE_ALT_MEM_VAR
    float ccAmpP=tap[i]*bb;
#else
    float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
    
    float tmoo2po = (2.0f-omega)/(2.0f+omega);
    float to2po = 2.0f/(2.0f+omega);
    float dVa = to2po*ccAmpP*dvTot;
    
    float prevRp=currRp[index];
    currRp[index]=tmoo2po*prevRp-dVa;
    
    //update the totals for the stress update
    float totalRp = omega*(currRp[index]+prevRp)*0.5f;
    
#if !USE_ALT_MEM_VAR
    totalRp*=bb;
#endif
    pressure[index]-=totalRp;
  }

  //Update the pressure in the interior (no MPML).  ssfunc fills the same
  //role as vxfunc, etc above.  This is just the kernel.
  inline void updateAcousticAttenPressureNoPMLCore(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,float *__restrict__  currRp,
                                              float *__restrict__  bulk,int NX,int NXY,float cx0,float cx1, 
                                              float cy0,float cy1,float cz0,float cz1,float dcx0,float dcy0,float dcz0,unsigned char *__restrict__  ssfunc, int index, int i,
                                              float *__restrict__  pressure) {
    float dvx, dvy, dvz;
    float bb = bulk[index];
    unsigned char dd = ssfunc[index];
    unsigned char dd1= 1-dd;
    float cx0d = cx0*dd1+dcx0*dd;
    float cx1d = cx1*dd1;
    float cy0d = cy0*dd1+dcy0*dd;
    float cy1d = cy1*dd1;
    float cz0d = cz0*dd1+dcz0*dd;
    float cz1d = cz1*dd1;
    tdvx[i] = dvx = cx0d*(vx[index]-vx[index-1])+
    cx1d*(vx[index+1]-vx[index-2]);
    tdvy[i] = dvy = cy0d*(vy[index]-vy[index-NX])+
    cy1d*(vy[index+NX]-vy[index-2*NX]);
    tdvz[i] = dvz = cz0d*(vz[index]-vz[index-NXY])+
    cz1d*(vz[index+NXY]-vz[index-2*NXY]);
    
    float dvTot = dvx+dvy+dvz;
    float omega=tcr[i];
#if USE_ALT_MEM_VAR
    float ccAmpP=tap[i]*bb;
#else
    float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
    
    float tmoo2po = (2.0f-omega)/(2.0f+omega);
    float to2po = 2.0f/(2.0f+omega);
    float dVa = to2po*ccAmpP*dvTot;
    
    float prevRp = currRp[index];
    currRp[index]=tmoo2po*prevRp-dVa;
    
    //update the totals for the pressure update
    float totalRp = omega*(currRp[index]+prevRp)*0.5f;
    
#if !USE_ALT_MEM_VAR
    totalRp*=bb;
#endif
    pressure[index]-=dvTot*bb+totalRp;
  }
  
  //Update the pressure in the interior (no MPML) for O2.  This is not actually ever called since O4 is used everywhere interior.  This is just the kernel.
  inline void updateAcousticAttenPressureNoPMLCoreO2(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,float *__restrict__  currRp,
                                                float *__restrict__  bulk,int NX,int NXY,float dcx0,float dcy0,float dcz0,
                                                int index, int i,
                                                float *__restrict__  pressure) {
    float dvx, dvy, dvz;
    float bb = bulk[index];
    tdvx[i] = dvx = dcx0*(vx[index]-vx[index-1]);
    tdvy[i] = dvy = dcy0*(vy[index]-vy[index-NX]);
    tdvz[i] = dvz = dcz0*(vz[index]-vz[index-NXY]);
    
    float dvTot = dvx+dvy+dvz;
    float omega=tcr[i];
#if USE_ALT_MEM_VAR
    float ccAmpP=tap[i]*bb;
#else
    float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
    
    float tmoo2po = (2.0f-omega)/(2.0f+omega);
    float to2po = 2.0f/(2.0f+omega);
    float dVa = to2po*ccAmpP*dvTot;
    
    float prevRp = currRp[index];
    currRp[index]=tmoo2po*prevRp-dVa;
    
    //update the totals for the pressure update
    float totalRp = omega*(currRp[index]+prevRp)*0.5f;
    
#if !USE_ALT_MEM_VAR
    totalRp*=bb;
#endif
    pressure[index]-=dvTot*bb+totalRp;
  }
  
  //This is what is actually called from updateAcousticMPMLFull.cc.  This performs the interior (no MPML) pressure update
  inline void updateAcousticAttenPressureNoPML(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz, float **__restrict__ rp, float **__restrict__ decayRates,float **__restrict__ ampP,
                                          int* __restrict__ Qindex, int nMechs,
                                          int i0, int ie,
                                          float *__restrict__  bulk,int NX,int NXY,float cx0,float cx1, 
                                          float cy0,float cy1,float cz0,float cz1,float dcx0,float dcy0,float dcz0,unsigned char *__restrict__  ssfunc, int jkind,
                                          float *__restrict__  pressure) {
    int ii=0;
    float* __restrict__ currRp=rp[ii];
    float* __restrict__ cRates=decayRates[ii];
    float* __restrict__ cAmpP=ampP[ii];
    for(int i=i0;i<ie;++i) {
      int qind = Qindex[i+jkind];
      tcr[i] = cRates[qind];
      tap[i] = cAmpP[qind];
    }
    for(int i=i0;i<ie;i++) {
      int index = i+jkind;
      updateAcousticAttenPressureNoPMLCore(vx,vy,vz,currRp,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,index,i,pressure);
    }
    for(ii=1;ii<nMechs;++ii) {
      currRp=rp[ii];
      cRates=decayRates[ii];
      cAmpP=ampP[ii];
      for(int i=i0;i<ie;++i) {
        int qind = Qindex[i+jkind];
        tcr[i] = cRates[qind];
        tap[i] = cAmpP[qind];
      }
      for(int i=i0;i<ie;i++) {
        int index = i+jkind;
        updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
      }
    }
  }
  
  //This is not ever called, but is analogous to updateAcousticPressureNoPML but for O2 updating
  inline void updateAcousticAttenPressureNoPMLO2(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz, float **__restrict__ rp, float **__restrict__ decayRates,float **__restrict__ ampP,
                                            int* __restrict__ Qindex, int nMechs,
                                            int i0, int ie,
                                            float *__restrict__  bulk,int NX,int NXY,float dcx0,float dcy0,float dcz0,int jkind,
                                            float *__restrict__  pressure) {
    int ii=0;
    float* __restrict__ currRp=rp[ii];
    float* __restrict__ cRates=decayRates[ii];
    float* __restrict__ cAmpP=ampP[ii];
    for(int i=i0;i<ie;++i) {
      int qind = Qindex[i+jkind];
      tcr[i] = cRates[qind];
      tap[i] = cAmpP[qind];
    }
    for(int i=i0;i<ie;i++) {
      int index = i+jkind;
      updateAcousticAttenPressureNoPMLCoreO2(vx,vy,vz,currRp,bulk,NX,NXY,dcx0,dcy0,dcz0,index,i,pressure);
    }
    for(ii=1;ii<nMechs;++ii) {
      currRp=rp[ii];
      cRates=decayRates[ii];
      cAmpP=ampP[ii];
      for(int i=i0;i<ie;++i) {
        int qind = Qindex[i+jkind];
        tcr[i] = cRates[qind];
        tap[i] = cAmpP[qind];
      }
      for(int i=i0;i<ie;i++) {
        int index = i+jkind;
        updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
      }
    }
  }
  
  //Update pressure with O4 updating in the MPML zone. MPML memory
  //variables are vxx, vyy and vzz, otherwise the variables have similar
  //meaning to above calls.  This is just the kernel.
  inline void updateAcousticAttenPressurePMLCore(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,float *__restrict__  currRp,
                                            float *__restrict__ vxx,float *__restrict__ vyy,float *__restrict__ vzz,
                                            float *__restrict__  bulk,int NX,int NXY,float cx0,float cx1, 
                                            float cy0,float cy1,float cz0,float cz1,float dcx0,float dcy0,float dcz0,unsigned char *__restrict__  ssfunc, int index, int pIndex, int i,
                                            float *__restrict__  pressure,
                                            float axtaper, float aytaper, float aztaper,
                                            float bxtaper, float bytaper, float bztaper,
                                            float kxtaper, float kytaper, float kztaper) {
    float dvx, dvy, dvz;
    float bb = bulk[index];
    unsigned char dd = ssfunc[index];
    unsigned char dd1= 1-dd;
    float cx0d = cx0*dd1+dcx0*dd;
    float cx1d = cx1*dd1;
    float cy0d = cy0*dd1+dcy0*dd;
    float cy1d = cy1*dd1;
    float cz0d = cz0*dd1+dcz0*dd;
    float cz1d = cz1*dd1;
    dvx = cx0d*(vx[index]-vx[index-1])+
    cx1d*(vx[index+1]-vx[index-2]);
    dvy = cy0d*(vy[index]-vy[index-NX])+
    cy1d*(vy[index+NX]-vy[index-2*NX]);
    dvz = cz0d*(vz[index]-vz[index-NXY])+
    cz1d*(vz[index+NXY]-vz[index-2*NXY]);

    //PML updating in group of lines below
    updatePMLMemVarComp(vxx,vyy,vzz,pIndex,axtaper,aytaper,aztaper,
                        bxtaper,bytaper,bztaper,
                        dvx,dvy,dvz);
    tdvx[i] = dvx = dvx*kxtaper+vxx[pIndex];
    tdvy[i] = dvy = dvy*kytaper+vyy[pIndex];
    tdvz[i] = dvz = dvz*kztaper+vzz[pIndex];
    
    float dvTot = dvx+dvy+dvz;
    float omega=tcr[i];
#if USE_ALT_MEM_VAR
    float ccAmpP=tap[i]*bb;
#else
    float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
    
    float tmoo2po = (2.0f-omega)/(2.0f+omega);
    float to2po = 2.0f/(2.0f+omega);
    float dVa = to2po*ccAmpP*dvTot;
    
    float prevRp = currRp[index];
    currRp[index]=tmoo2po*prevRp-dVa;
    
    //update the totals for the pressure update
    float totalRp = omega*(currRp[index]+prevRp)*0.5f;
    
#if !USE_ALT_MEM_VAR
    totalRp*=bb;
#endif
    pressure[index]-=dvTot*bb+totalRp;
  }
  
  //Update the pressure with O2 accuracy in the MPML.  This is only used on the very flank of the model.  This is just the kernel.
  inline void updateAcousticAttenPressurePMLCoreO2(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,float *__restrict__  currRp,
                                              float *__restrict__ vxx,float *__restrict__ vyy,float *__restrict__ vzz,
                                              float *__restrict__  bulk,int NX,int NXY, 
                                              float dcx0,float dcy0,float dcz0,int index, int pIndex, int i,
                                              float *__restrict__  pressure,
                                              float axtaper, float aytaper, float aztaper,
                                              float bxtaper, float bytaper, float bztaper,
                                              float kxtaper, float kytaper, float kztaper) {
    float dvx, dvy, dvz;
    float bb = bulk[index];
    dvx = dcx0*(vx[index]-vx[index-1]);
    dvy = dcy0*(vy[index]-vy[index-NX]);
    dvz = dcz0*(vz[index]-vz[index-NXY]);
    
    //PML updating in group of lines below
    updatePMLMemVarComp(vxx,vyy,vzz,pIndex,axtaper,aytaper,aztaper,
                        bxtaper,bytaper,bztaper,
                        dvx,dvy,dvz);
    tdvx[i] = dvx = dvx*kxtaper+vxx[pIndex];
    tdvy[i] = dvy = dvy*kytaper+vyy[pIndex];
    tdvz[i] = dvz = dvz*kztaper+vzz[pIndex];
    
    float dvTot = dvx+dvy+dvz;
    float omega=tcr[i];
#if USE_ALT_MEM_VAR
    float ccAmpP=tap[i]*bb;
#else
    float ccAmpP=tap[i];
#endif //USE_ALT_MEM_VAR
    
    float tmoo2po = (2.0f-omega)/(2.0f+omega);
    float to2po = 2.0f/(2.0f+omega);
    float dVa = to2po*ccAmpP*dvTot;
    
    float prevRp = currRp[index];
    currRp[index]=tmoo2po*prevRp-dVa;
    
    //update the totals for the pressure update
    float totalRp = omega*(currRp[index]+prevRp)*0.5f;
    
#if !USE_ALT_MEM_VAR
    totalRp*=bb;
#endif
    pressure[index]-=dvTot*bb+totalRp;
  }
  
  //This is what is actually called from updateAcousticMPMLFull.cc for updating pressure in the MPML zone with O4 accuracy
  inline void updateAcousticAttenPressurePML(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,
                                        float *__restrict__ vxx, float *__restrict__ vyy, float* __restrict__ vzz,float **__restrict__ rp, float **__restrict__ decayRates,float **__restrict__ ampP,
                                        int* __restrict__ Qindex, int nMechs,
                                        int i0, int ie, int di,
                                        float *__restrict__  bulk,int NX,int NXY,float cx0,float cx1, 
                                        float cy0,float cy1,float cz0,float cz1,float dcx0,float dcy0,float dcz0,unsigned char *__restrict__  ssfunc, int jkind, int jkpIndex,
                                        float *__restrict__  pressure, float axtaper,
                                        float aytaper, float aztaper, float bxtaper, float bytaper, float bztaper, float kytaper, float kztaper) {
    int ii=0;
    float* __restrict__ currRp=rp[ii];
    float* __restrict__ cRates=decayRates[ii];
    float* __restrict__ cAmpP=ampP[ii];
    for(int i=i0;i<ie;++i) {
      int qind = Qindex[i+jkind];
      tcr[i] = cRates[qind];
      tap[i] = cAmpP[qind];
    }
    for(int i=i0;i<ie;i++) {
      int index = i+jkind;
      int pIndex = i+jkpIndex;
      updateAcousticAttenPressurePMLCore(vx,vy,vz,currRp,vxx,vyy,vzz,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,index,pIndex,i,pressure,axtaper,aytaper,aztaper,bxtaper,bytaper,bztaper,kxTaper[i],kytaper,kztaper);
    }
    for(ii=1;ii<nMechs;++ii) {
      currRp=rp[ii];
      cRates=decayRates[ii];
      cAmpP=ampP[ii];
      for(int i=i0;i<ie;++i) {
        int qind = Qindex[i+jkind];
        tcr[i] = cRates[qind];
        tap[i] = cAmpP[qind];
      }
      for(int i=i0;i<ie;i++) {
        int index = i+jkind;
        updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
      }
    }
  }

  inline void updateAcousticAttenPressurePML(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,
                                        float *__restrict__ vxx, float *__restrict__ vyy, float* __restrict__ vzz,float **__restrict__ rp, float **__restrict__ decayRates,float **__restrict__ ampP,
                                        int* __restrict__ Qindex, int nMechs,
                                        int i0, int ie, int di,
                                        float *__restrict__  bulk,int NX,int NXY,float cx0,float cx1, 
                                        float cy0,float cy1,float cz0,float cz1,float dcx0,float dcy0,float dcz0,unsigned char *__restrict__  ssfunc, int jkind, int jkpIndex,
                                        float *__restrict__  pressure, float *__restrict__ axtaper,
                                        float *__restrict__ aytaper, float *__restrict__ aztaper, float *__restrict__ bxtaper, float *__restrict__ bytaper, float *__restrict__ bztaper, float kytaper, float kztaper) {
    int ii=0;
    float* __restrict__ currRp=rp[ii];
    float* __restrict__ cRates=decayRates[ii];
    float* __restrict__ cAmpP=ampP[ii];
    for(int i=i0;i<ie;++i) {
      int qind = Qindex[i+jkind];
      tcr[i] = cRates[qind];
      tap[i] = cAmpP[qind];
    }
    for(int i=i0;i<ie;i++) {
      int ti = i-di;
      int index = i+jkind;
      int pIndex = i+jkpIndex;
      updateAcousticAttenPressurePMLCore(vx,vy,vz,currRp,vxx,vyy,vzz,bulk,NX,NXY,cx0,cx1,cy0,cy1,cz0,cz1,dcx0,dcy0,dcz0,ssfunc,index,pIndex,i,pressure,axtaper[ti],aytaper[ti],aztaper[ti],bxtaper[ti],bytaper[ti],bztaper[ti],kxTaper[i],kytaper,kztaper);
    }
    for(ii=1;ii<nMechs;++ii) {
      currRp=rp[ii];
      cRates=decayRates[ii];
      cAmpP=ampP[ii];
      for(int i=i0;i<ie;++i) {
        int qind = Qindex[i+jkind];
        tcr[i] = cRates[qind];
        tap[i] = cAmpP[qind];
      }
      for(int i=i0;i<ie;i++) {
        int index = i+jkind;
        updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
      }
    }
  }


  
  //This is what is actually called from updateAcousticMPMLFull.cc for updating pressure in the MPML zone with O2 accuracy.  This is only used on the very flanks of the model.
  inline void updateAcousticAttenPressurePMLO2(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,
                                          float *__restrict__ vxx, float *__restrict__ vyy, float* __restrict__ vzz,float **__restrict__ rp, float **__restrict__ decayRates,float **__restrict__ ampP,
                                          int* __restrict__ Qindex, int nMechs,
                                          int i0, int ie, int di,
                                          float *__restrict__  bulk,int NX,int NXY, 
                                          float dcx0,float dcy0,float dcz0,int jkind, int jkpIndex,
                                          float *__restrict__  pressure, float axtaper,
                                          float aytaper, float aztaper, float bxtaper, float bytaper, float bztaper, float kytaper, float kztaper) {
    int ii=0;
    float* __restrict__ currRp=rp[ii];
    float* __restrict__ cRates=decayRates[ii];
    float* __restrict__ cAmpP=ampP[ii];
    for(int i=i0;i<ie;++i) {
      int qind = Qindex[i+jkind];
      tcr[i] = cRates[qind];
      tap[i] = cAmpP[qind];
    }
    for(int i=i0;i<ie;i++) {
      int index = i+jkind;
      int pIndex = i+jkpIndex;
      updateAcousticAttenPressurePMLCoreO2(vx,vy,vz,currRp,vxx,vyy,vzz,bulk,NX,NXY,dcx0,dcy0,dcz0,index,pIndex,i,pressure,axtaper,aytaper,aztaper,bxtaper,bytaper,bztaper,kxTaper[i],kytaper,kztaper);
    }
    for(ii=1;ii<nMechs;++ii) {
      currRp=rp[ii];
      cRates=decayRates[ii];
      cAmpP=ampP[ii];
      for(int i=i0;i<ie;++i) {
        int qind = Qindex[i+jkind];
        tcr[i] = cRates[qind];
        tap[i] = cAmpP[qind];
      }
      for(int i=i0;i<ie;i++) {
        int index = i+jkind;
        updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
      }
    }
  }

  inline void updateAcousticAttenPressurePMLO2(float *__restrict__  vx,float *__restrict__  vy,float *__restrict__  vz,
                                          float *__restrict__ vxx, float *__restrict__ vyy, float* __restrict__ vzz, float **__restrict__ rp, float **__restrict__ decayRates,float **__restrict__ ampP,
                                          int* __restrict__ Qindex, int nMechs,
                                          int i0, int ie, int di,
                                          float *__restrict__  bulk,int NX,int NXY, 
                                          float dcx0,float dcy0,float dcz0,int jkind, int jkpIndex,
                                          float *__restrict__  pressure, float *__restrict__ axtaper,
                                          float *__restrict__ aytaper, float *__restrict__ aztaper, float *__restrict__ bxtaper, float *__restrict__ bytaper, float *__restrict__ bztaper, float kytaper, float kztaper) {
    int ii=0;
    float* __restrict__ currRp=rp[ii];
    float* __restrict__ cRates=decayRates[ii];
    float* __restrict__ cAmpP=ampP[ii];
    for(int i=i0;i<ie;++i) {
      int qind = Qindex[i+jkind];
      tcr[i] = cRates[qind];
      tap[i] = cAmpP[qind];
    }
    for(int i=i0;i<ie;i++) {
      int ti = i-di;
      int index = i+jkind;
      int pIndex = i+jkpIndex;
      updateAcousticAttenPressurePMLCoreO2(vx,vy,vz,currRp,vxx,vyy,vzz,bulk,NX,NXY,dcx0,dcy0,dcz0,index,pIndex,i,pressure,axtaper[ti],aytaper[ti],aztaper[ti],bxtaper[ti],bytaper[ti],bztaper[ti],kxTaper[i],kytaper,kztaper);
    }
    for(ii=1;ii<nMechs;++ii) {
      currRp=rp[ii];
      cRates=decayRates[ii];
      cAmpP=ampP[ii];
      for(int i=i0;i<ie;++i) {
        int qind = Qindex[i+jkind];
        tcr[i] = cRates[qind];
        tap[i] = cAmpP[qind];
      }
      for(int i=i0;i<ie;i++) {
        int index = i+jkind;
        updateAcousticAttenNoPMLRest(vx,vy,vz,currRp,bulk,NX,NXY,index,i,pressure);
      }
    }
  }



}

#endif
