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
 *  updateAnelasticMPMLFull.h
 *  
 *
 *  Declarations of functions used for MPML updating
 *
 *  Declares the following functions:
 *  setupAcousticMPMLBoundsFull
 *  updateVxMPMLBounds
 *  updateVyMPMLBounds
 *  updateVzMPMLBounds
 *  updateVzPressFreeMPMLBounds
 *  updateAcousticPressureMPML
 *  updateAcousticPressurePressFreeMPML
 *
 */

class sgfdModel;

#include "sgfd.h"

void setupAcousticAttenMPMLBoundsFull(int nXmin, float xMinVal, float xMinAVal, float xMinKVal,int nXmax, float xMaxVal, float xMaxAval, float xMaxKVal, int nYmin, float yMinVal, float yMinAVal, float yMinKVal,
                                  int nYmax, float yMaxVal, float yMaxAVal, float yMaxKVal, int nZmin, float zMinVal, float zMinAVal, float zMinKVal,int nZmax, float zMaxVal, float zMaxAVal, float zMaxKVal,int kstart, int kstart_z, sgfdModel* model, int sbcmode, float xfac);
void updateVxAttenMPMLBounds(modelDefStruct* modelDef,
                        float* __restrict__ vx,
                        float* __restrict__ pressure,
                        float cx[2],float cy[2],float cz[2],
                        float* __restrict__ rho,unsigned char* __restrict__ vxfunc,
                        bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                        bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb);
void updateVyAttenMPMLBounds(modelDefStruct* modelDef,
                        float* __restrict__ vy,
                        float* __restrict__ yy,
                        float cx[2],float cy[2],float cz[2],
                        float* __restrict__ rho,
                        unsigned char* __restrict__ vyfunc,
                        bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                        bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb);
void updateVzAttenMPMLBounds(modelDefStruct* modelDef,
                        float* __restrict__ vz,
                        float* __restrict__ zz,
                        float cx[2],float cy[2],float cz[2],
                        float* __restrict__ rho,
                        unsigned char* __restrict__ vzfunc,
                        bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                        bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb);
void updateAcousticAttenPressureMPML(modelDefStruct* modelDef,
                                   float*__restrict__ pressure,
                                     float*__restrict__ vx,float*__restrict__ vy,float*__restrict__ vz,int*__restrict__ Qindex, int nMechs,
                                     float**__restrict__ decayRates, float**__restrict__ ampP,
                                     float**__restrict__ rp,
                                   float cx[2],float cy[2],float cz[2],
                                   float*__restrict__ bulk,
                                   unsigned char*__restrict__ ssfunc,
                                   bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                                   bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb);
void updateVzPressFreeAttenMPMLBounds(modelDefStruct* modelDef,
                        float* __restrict__ vz,
                        float* __restrict__ zz,
                        float cx[2],float cy[2],float cz[2],
                        float* __restrict__ rho,
                        unsigned char* __restrict__ vzfunc,
                        bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb);
void updateAcousticAttenPressurePressFreeMPML(modelDefStruct* modelDef,
                                float*__restrict__ pressure,
                                float*__restrict__ vx,float*__restrict__ vy,float*__restrict__ vz,
                                float cx[2],float cy[2],float cz[2],
                                float*__restrict__ bulk,
                                unsigned char*__restrict__ ssfunc,
                                bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb);
