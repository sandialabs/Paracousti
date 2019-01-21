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
 *  updateAnelasticSpongeFull.h
 *  
 *
 *  Declarations of function udes in sponge updating.
 *
 *  Declares the following functions:
 *  setupAcousticSpongeBoundsFull
 *  updateVxSpongeBounds
 *  updateVySpongeBounds
 *  updateVzSpongeBounds
 *  updateVzPressFreeSpongeBounds
 *  updateAcousticPressureSponge
 *  updateAcousticPressurePressFreeSponge
 *
 */

class sgfdModel;

#include "sgfd.h"

void setupAcousticSpongeBoundsFull(int nXmin, float xMinVal, int nXmax, float xMaxVal, int nYmin, float yMinVal,
                                  int nYmax, float yMaxVal, int nZmin, float zMinVal, int nZmax, float zMaxVal, int kstart, int kstart_z, sgfdModel* model, int sbcmode);
void updateVxSpongeBounds(modelDefStruct* modelDef,
                        float* __restrict__ vx,
                        float* __restrict__ pressure,
                        float cx[2],float cy[2],float cz[2],
                        float* __restrict__ rho,unsigned char* __restrict__ vxfunc,
                        bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                        bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb);
void updateVySpongeBounds(modelDefStruct* modelDef,
                        float* __restrict__ vy,
                        float* __restrict__ yy,
                        float cx[2],float cy[2],float cz[2],
                        float* __restrict__ rho,
                        unsigned char* __restrict__ vyfunc,
                        bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                        bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb);
void updateVzSpongeBounds(modelDefStruct* modelDef,
                        float* __restrict__ vz,
                        float* __restrict__ zz,
                        float cx[2],float cy[2],float cz[2],
                        float* __restrict__ rho,
                        unsigned char* __restrict__ vzfunc,
                        bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                        bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb);
void updateAcousticPressureSponge(modelDefStruct* modelDef,
                                   float*__restrict__ pressure,
                                   float*__restrict__ vx,float*__restrict__ vy,float*__restrict__ vz,
                                   float cx[2],float cy[2],float cz[2],
                                   float*__restrict__ bulk,
                                   unsigned char*__restrict__ ssfunc,
                                   bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb,
                                   bool fminXb, bool fmaxXb, bool fminYb, bool fmaxYb, bool fminZb, bool fmaxZb);
void updateVzPressFreeSpongeBounds(modelDefStruct* modelDef,
                        float* __restrict__ vz,
                        float* __restrict__ zz,
                        float cx[2],float cy[2],float cz[2],
                        float* __restrict__ rho,
                        unsigned char* __restrict__ vzfunc,
                        bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb);
void updateAcousticPressurePressFreeSponge(modelDefStruct* modelDef,
                                float*__restrict__ pressure,
                                float*__restrict__ vx,float*__restrict__ vy,float*__restrict__ vz,
                                float cx[2],float cy[2],float cz[2],
                                float*__restrict__ bulk,
                                unsigned char*__restrict__ ssfunc,
                                bool minXb, bool maxXb, bool minYb, bool maxYb, bool minZb, bool maxZb);
