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
 *  sgfd.cc
 *
 *
 *  Define several class functions that control the base class behavior
 *  for many aspects of the algorithm including models, slaves, masters,
 *  dependents and command line processing.
 *
 *  Defines the following class functions:
 *  sgfdModel::uniqueFilename
 *  sgfdModel::multiplyGrid
 *  sgfdModel::holbergCoeffsC
 *  sgfdModel::calcHolbergCoef
 *  sgfdModel::readCDFModel
 *  sgfdModel::minMaxVelocity
 *  sgfdModel::readCDFModelInterp
 *  sgfdDependent::checkCheckpoint
 *  sgfdDependent::doCheckpoint
 *  sgfdDependent::readCheckpoint
 *  sgfdDependent::doYZSlice
 *  sgfdDependent::doXZSlice
 *  sgfdDependent::doXZSlice
 *  sgfdDependent::doXYSlice
 *  masterSgfdModel::resetTimeVector
 *  masterSgfdModel::multiplyGrid
 *  masterSgfdModel::getInitAck
 *  masterSgfdModel::initializeSlaves
 *  masterSgfdModel::buildSlaveModelDef
 *  masterSgfdModel::setNeighbors
 *  masterSgfdModel::sendInitMessage
 *  slaveSgfdModel::slaveSgfdModel
 *  slaveSgfdModel::resetTimeVector
 *  slaveSgfdModel::isInside
 *  slaveSgfdModel::isInside3
 *  slaveSgfdModel::interp3
 *  slaveSgfdModel::interp2
 *  masterSgfdDependent::doCheckpoint
 *  masterSgfdDependent::readCheckpoint
 *  masterSgfdDependent::setInitialConditions
 *  masterSgfdDependent::setBoundaryConditions
 *  masterSgfdDependent::setReceivers
 *  masterSgfdDependent::setSources
 *  masterSgfdDependent::writeTraces
 *  masterSgfdDependent::advance
 *  masterSgfdDependent::advanceVel
 *  masterSgfdDependent::advanceStress
 *  masterSgfdDependent::getAcknowledgement
 *  masterSgfdDependent::sendInitialMessage
 *  masterSgfdDependent::allocateReceivers
 *  slaveSgfdDependent::doCheckpoint
 *  slaveSgfdDependent::readCheckpoint
 *  slaveSgfdDependent::setSources
 *  slaveSgfdDependent::setReceivers
 *  slaveSgfdDependent::advance
 *  slaveSgfdDependent::packReceiverData
 *  slaveSgfdDependent::sendReceiverGridData
 *  slaveSgfdDependent::fullFieldOutput
 *  slaveSgfdDependent::slicer
 *  slaveSgfdDependent::slicer
 *  slaveSgfdDependent::allocateReceivers
 *  slaveSgfdDependent::sendAcknowledgement
 *  slaveSgfdDependent::sliceXZCompValue
 *  slaveSgfdDependent::sliceYZCompValue
 *  slaveSgfdDependent::sliceXYCompValue
 *  sgfdCommandParser::allocateExtraOutput
 *  sgfdCommandParser::checkProcessFlag
 *   sgfdCommandParser::addExtraOutputFlags
 *
 *  Defines the following functions:
 *  extraOutputCompFunc
 *  newModelDef
 *  packModelDef
 *  unpackModelDef
 *  copyModelDef
 *  newParallelDef
 *  packParallelDef
 *  unpackParallelDef
 *  copyParallelDef
 *  readCDFRunParams
 *  readCDFRunParams
 *
 */

#include "sgfd.hh"
#include "xtrautil.hh"
#include "sgfdSources.hh"
#include "movingSgfdSources.hh"
#include "extra_output.hh"
#include "sgfdReceiverNetwork.hh"
#include "message_passing.h"

#ifdef SGFD_INVERSION
#include "inversionReceivers.hh"
#endif


//Actual definitions of extern variables from header files.
// sgfdReceivers.hh
int NumReceiverTypeNames=15;
const char* ReceiverTypeNames[]={"","Velocity","Pressure","3C","4C","Inversion",
                           "Vx","Vy","Vz",
                           "XX","YY","ZZ","XY","XZ","YZ","Complete"};

//Symbols from extra_output.hh
int extraOutputCompFunc(const void* t1,const void* t2){
  extraOutput *tt1=*((extraOutput**)t1),*tt2=*((extraOutput**)t2);
  if     (tt1->t()<tt2->t()) return -1;
  else if(tt1->t()>tt2->t()) return  1;
  return 0;
}

//
//Procedure Declarations
//

//Subroutines for use with the modelDefStruct. It would be
// better to define this as a class but it needs to be usable
// by c as well as c++ procedures.
modelDefStruct* newModelDef(int xStart,int xStop,
			    int yStart,int yStop,
			    int zStart,int zStop,
			    int nt,
			    float minX,float dx,float minY,float dy,
			    float minZ,float dz,float minT,float dt,
			    float scalarSpeed,float scalarVel,
			    float scalarDen,float scalarStress,
			    int multiplier){
  modelDefStruct* modelDef=(modelDefStruct*)malloc(sizeof(modelDefStruct));
  assert(modelDef!=NULL,
	 "newModelDef--unable to allocate %i bytes for modelDefStruct",
	 sizeof(modelDefStruct));
  
  modelDef->procLim[0]=xStart;
  modelDef->procLim[1]=xStop;
  modelDef->procLim[2]=yStart;
  modelDef->procLim[3]=yStop;
  modelDef->procLim[4]=zStart;
  modelDef->procLim[5]=zStop;

  modelDef->NX=xStop-xStart;
  modelDef->NY=yStop-yStart;
  modelDef->NZ=zStop-zStart;
  modelDef->NT=nt;

  modelDef->NXY=modelDef->NX*modelDef->NY;
  modelDef->NXYZ=modelDef->NXY*modelDef->NZ;

  modelDef->minX=minX+dx*xStart;
  modelDef->dx=dx;

  modelDef->minY=minY+dy*yStart;
  modelDef->dy=dy;

  modelDef->minZ=minZ+dz*zStart;
  modelDef->dz=dz;

  modelDef->minT=minT;
  modelDef->dt=dt;

  modelDef->scalarSpeed=scalarSpeed;
  modelDef->scalarVel=scalarVel;
  modelDef->scalarDen=scalarDen;
  modelDef->scalarStress=scalarStress;
  return modelDef;
}

modelDefStruct* packModelDef(modelDefStruct* modelDef,int destroy){
  if(!modelDef){
    packMessage("i",0);
    return modelDef;
  }else{
    packMessage("i",1);
  }

  packMessage("I",modelDef->procLim,6);
  //need to pack and unpack NX, NY, NZ and NT instead of deriving from
  // procLim as was done previously. Then changes with different multipliers
  // only need to be made in one place.
  packMessage("iiii",modelDef->NX,modelDef->NY,modelDef->NZ,modelDef->NT);

  packMessage("ffffffff",
	      modelDef->minX,modelDef->dx,modelDef->minY,modelDef->dy,
	      modelDef->minZ,modelDef->dz,modelDef->minT,modelDef->dt);

  packMessage("ffff",
	      modelDef->scalarSpeed,modelDef->scalarVel,
	      modelDef->scalarDen,modelDef->scalarStress);

  if(destroy){
    free(modelDef);
    return NULL;
  }
  return modelDef;
}

modelDefStruct* unpackModelDef(modelDefStruct* modelDef){
  int flag;
  unpackMessage("i",&flag);
  if(!flag) return NULL;

  if(!modelDef)
    assert((modelDef=(modelDefStruct*)malloc(sizeof(modelDefStruct)))!=NULL,
	   "unpackModelDef--unable to allocate %i bytes for new modelDefStruct",
	   sizeof(modelDefStruct));

  unpackMessage("I",modelDef->procLim,6);
  //need to pack and unpack NX, NY, NZ and NT instead of deriving from
  // procLim as was done previously. Then changes with different multipliers
  // only need to be made in one place.
  unpackMessage("iiii",&modelDef->NX,&modelDef->NY,&modelDef->NZ,&modelDef->NT);
  //still calculate NXY and NXYZ
  modelDef->NXY=modelDef->NX*modelDef->NY;
  modelDef->NXYZ=modelDef->NXY*modelDef->NZ;

  unpackMessage("ffffffff",
		&modelDef->minX,&modelDef->dx,&modelDef->minY,&modelDef->dy,
		&modelDef->minZ,&modelDef->dz,&modelDef->minT,&modelDef->dt);

  unpackMessage("ffff",
		&modelDef->scalarSpeed,&modelDef->scalarVel,
		&modelDef->scalarDen,&modelDef->scalarStress);

  return modelDef;
}
void copyModelDef(modelDefStruct* target,modelDefStruct *source){
  for(int i=0;i<6;i++)
    target->procLim[i]=source->procLim[i];

  target->NX=source->NX;
  target->NY=source->NY;
  target->NZ=source->NZ;
  target->NT=source->NT;

  target->NXY=target->NX*target->NY;
  target->NXYZ=target->NXY*target->NZ;

  target->minX=source->minX;
  target->dx=source->dx;

  target->minY=source->minY;
  target->dy=source->dy;

  target->minZ=source->minZ;
  target->dz=source->dz;

  target->minT=source->minT;
  target->dt=source->dt;

  target->scalarSpeed=source->scalarSpeed;
  target->scalarVel=source->scalarVel;
  target->scalarDen=source->scalarDen;
  target->scalarStress=source->scalarStress;
}

//Same for the parallelDefStruct
parallelDefStruct* newParallelDef(int globalNX,int globalNY,int globalNZ,
				  int nxProc,int nyProc,int nzProc,
				  int procI,int procJ,int procK,
				  int *neighbors){
  parallelDefStruct* parallelDef=(parallelDefStruct*)malloc(sizeof(parallelDefStruct));
  assert(parallelDef!=NULL,
	 "newParallelDef--unable to allocate %i bytes for new parallelDefStruct",
	 sizeof(parallelDefStruct));

  parallelDef->globalNX=globalNX;
  parallelDef->globalNY=globalNY;
  parallelDef->globalNZ=globalNZ;

  parallelDef->nxProc=nxProc;
  parallelDef->nyProc=nyProc;
  parallelDef->nzProc=nzProc;
  parallelDef->nxyProc=nxProc*nyProc;

  parallelDef->procI=procI;
  parallelDef->procJ=procJ;
  parallelDef->procK=procK;

  for(int i=0;neighbors && i<27;i++){
    parallelDef->neighbors[i]=neighbors[i];
  }

  return parallelDef;
}
parallelDefStruct* packParallelDef(parallelDefStruct* parallelDef,int destroy){
  if(!parallelDef){
    packMessage("i",0);
    return NULL;
  }else{
    packMessage("i",1);
  }

  packMessage("iii",
	      parallelDef->globalNX,parallelDef->globalNY,parallelDef->globalNZ);

  packMessage("iii",parallelDef->nxProc,parallelDef->nyProc,parallelDef->nzProc);

  packMessage("iii",parallelDef->procI,parallelDef->procJ,parallelDef->procK);

  packMessage("I",parallelDef->neighbors,27);

  for(int i=0;i<27;packModelDef(parallelDef->neighborDefs[i++]));

  if(destroy){
    free(parallelDef);
    return NULL;
  }
  return parallelDef;
}
parallelDefStruct* unpackParallelDef(parallelDefStruct* parallelDef){
  int flag;
  unpackMessage("i",&flag);
  if(!flag) return NULL;

  if(!parallelDef)
    assert((parallelDef=(parallelDefStruct*)malloc(sizeof(parallelDefStruct)))!=NULL,
	   "unpackParallelDef--unable to allocate %i bytes for new parallelDefStruct",
	   sizeof(parallelDefStruct));

  unpackMessage("iii",
		&(parallelDef->globalNX),
		&(parallelDef->globalNY),
		&(parallelDef->globalNZ));
  unpackMessage("iii",
		&(parallelDef->nxProc),&(parallelDef->nyProc),&(parallelDef->nzProc));
  parallelDef->nxyProc=parallelDef->nxProc*parallelDef->nyProc;

  unpackMessage("iii",&parallelDef->procI,&parallelDef->procJ,&parallelDef->procK);
  
  unpackMessage("I",parallelDef->neighbors,27);

  for(int i=0;i<27;i++){
    parallelDef->neighborDefs[i]=
      unpackModelDef();
  }

  return parallelDef;
}
void copyParallelDef(parallelDefStruct* target,parallelDefStruct* source){
  target->globalNX=source->globalNX;
  target->globalNY=source->globalNY;
  target->globalNZ=source->globalNZ;

  target->nxProc=source->nxProc;
  target->nyProc=source->nyProc;
  target->nzProc=source->nzProc;
  target->nxyProc=target->nxProc*target->nyProc;

  target->procI=source->procI;
  target->procJ=source->procJ;
  target->procK=source->procK;

  for(int i=0;i<27;i++){
    target->neighbors[i]=source->neighbors[i];
  }
}

//And subroutines to read the modelDef from a netCDF model file.
//read run parameters from a cdf file
// at this time the run parameters are the grid size and limits
// as well as the scale factors
int readCDFRunParams(const char* fileName,modelDefStruct* modelDef,
		     int print){
  int returnVal=readCDFRunParams(fileName,
				 modelDef->NX,modelDef->NY,modelDef->NZ,modelDef->NT,
				 modelDef->NXY,modelDef->NXYZ,
				 modelDef->minX,modelDef->minY,modelDef->minZ,modelDef->minT,
				 modelDef->dx,modelDef->dy,modelDef->dz,modelDef->dt);
  tEprintf(print,"Run Parameters from\n\t\"%s\"\n",fileName);
  tEprintf(print,
	   "\tX: start %.1f; dx %.1f; nx %i=>stop %.1f\n",
	   modelDef->minX,modelDef->dx,modelDef->NX,
	   modelDef->minX+modelDef->dx*(modelDef->NX-1));
  tEprintf(print,
	   "\tY: start %.1f; dy %.1f; ny %i=>stop %.1f\n",
	   modelDef->minY,modelDef->dy,modelDef->NY,
	   modelDef->minY+modelDef->dy*(modelDef->NY-1));
  tEprintf(print,
	   "\tZ: start %.1f; dz %.1f; nz %i=>stop %.1f\n",
	   modelDef->minZ,modelDef->dz,modelDef->NZ,
	   modelDef->minZ+modelDef->dz*(modelDef->NZ-1));
  tEprintf(print,
	   "\tT: start %.3fms; dt %.3fms; nt %i=>stop %.3fms\n",
	   1000*modelDef->minT,1000*modelDef->dt,modelDef->NT,
	   1000*(modelDef->minT+modelDef->dt*(modelDef->NT-1)));
  tEprintf(print,"\tModel is ~%.2f million nodes\n",
	   (float)modelDef->NXYZ/1e6);
  return returnVal;
}   
int readCDFRunParams(const char* fileName,
		     int& nx,int& ny,int& nz,int& nt,
		     int& nxy,int& nxyz,
		     float& minX,float& minY,float& minZ,float& minT,
		     float& dx,float& dy,float& dz,float& dt){
  //Make sure filename ends in .cdf
  char buffer[512];
  if(!strcmp(fileName+strlen(fileName)-4,".cdf")){
    sprintf(buffer,"%s",fileName);
  }else{
    sprintf(buffer,"%s.cdf",fileName);
  }

  //open cdf file for input
  int inFile;
  assert(nc_open(buffer,NC_NOWRITE,&inFile)==NC_NOERR,
	 "readCDFRunParams--unable to open file %s for reading\n",
	 buffer);

  //get the file dimensions
  size_t temp;
  int nxDim,nyDim,nzDim,ntDim;
  assert(nc_inq_dimid(inFile,"NX",&nxDim)==NC_NOERR,
	 "readCDFRunParams--unable to read NX from file %s",
	 buffer);
  nc_inq_dimlen(inFile,nxDim,&temp);
  nx=temp;
  assert(nc_inq_dimid(inFile,"NY",&nyDim)==NC_NOERR,
	 "readCDFRunParams--unable to read NY from file %s",
	 buffer);
  nc_inq_dimlen(inFile,nyDim,&temp);
  ny=temp;
  assert(nc_inq_dimid(inFile,"NZ",&nzDim)==NC_NOERR,
	 "readCDFRunParams--unable to read NZ from file %s",
	 buffer);
  nc_inq_dimlen(inFile,nzDim,&temp);
  nz=temp;
  assert(nc_inq_dimid(inFile,"NT",&ntDim)==NC_NOERR,
	 "readCDFRunParams--unable to read NT from file %s",
	 buffer);
  nc_inq_dimlen(inFile,ntDim,&temp);
  nt=temp;

  nxy=nx*ny;
  nxyz=nx*ny*nz;

  //get the file limits
  int minima;
  assert(nc_inq_varid(inFile,"minima",&minima)==NC_NOERR,
	 "readCDFRunParams--unable to open variable minima in file %s",
	 buffer);
  size_t index=0;
  nc_get_var1_float(inFile,minima,&index,&minX);
  index=1;
  nc_get_var1_float(inFile,minima,&index,&minY);
  index=2;
  nc_get_var1_float(inFile,minima,&index,&minZ);
  index=3;
  nc_get_var1_float(inFile,minima,&index,&minT);

  //get the increments
  int increments;
  assert(nc_inq_varid(inFile,"increments",&increments)==NC_NOERR,
	 "readCDFRunParams--unable to open variable increments in file %s",
	 buffer);
  index=0;
  nc_get_var1_float(inFile,increments,&index,&dx);
  index=1;
  nc_get_var1_float(inFile,increments,&index,&dy);
  index=2;
  nc_get_var1_float(inFile,increments,&index,&dz);
  index=3;
  nc_get_var1_float(inFile,increments,&index,&dt);
  
  nc_close(inFile);
  return nxyz;
}

char* sgfdModel::uniqueFilename(const char* head,const char* indir){
  char dir[512];
  DEF_MODEL_SIZE(modelDef());
  
  if(!indir){
    sprintf(dir,"");
  }else if(indir[strlen(indir)]!='/'){
    sprintf(dir,"%s/",indir);
  } else
    sprintf(dir,"%s",indir);
  
  if(head && *head){
    sprintf(_uniqueFilename,"%s%s.%i-%i.%i-%i.%i-%i",
            dir,head,
            procLim[0],procLim[1],
            procLim[2],procLim[3],
            procLim[4],procLim[5]);
  }else{
    sprintf(_uniqueFilename,"%s%i-%i.%i-%i.%i-%i",
            dir,
            procLim[0],procLim[1],
            procLim[2],procLim[3],
            procLim[4],procLim[5]);
  }
  return _uniqueFilename;
}

int sgfdModel::multiplyGrid(int surfaceBCType,int gridMultiplier){
  _gridMultiplier=gridMultiplier;
  
  _modelDef->NX= 1+gridMultiplier*(_modelDef->NX-1);
  _modelDef->NY= 1+gridMultiplier*(_modelDef->NY-1);
  if(!surfaceBCType){
    _modelDef->NZ= 1+gridMultiplier*(_modelDef->NZ-1);
  }else{
    _modelDef->NZ= 1+gridMultiplier*(_modelDef->NZ-2);
    _modelDef->minZ= _modelDef->minZ+2*_modelDef->dz-
    2*_modelDef->dz/(float)gridMultiplier;
  }
  
  _modelDef->NXY= _modelDef->NX*_modelDef->NY;
  _modelDef->NXYZ= _modelDef->NX*_modelDef->NY*_modelDef->NZ;
  
  _modelDef->dx/= (float)gridMultiplier;
  _modelDef->dy/= (float)gridMultiplier;
  _modelDef->dz/= (float)gridMultiplier;
  
  return _modelDef->NXYZ;
}

void sgfdModel::holbergCoeffsC(float rel_error,
                    float* c1,float* c2){
  double aa,bb,cc,amp,phase,part_real,part_imag;
  
  aa=1.0-rel_error;
  bb=1.0+rel_error;
  
  amp=pow(bb,(1.0/3.0));
  phase=atan2(2.0*sqrt(rel_error),aa)/3.0;
  part_real=amp*cos(phase);
  part_imag=amp*sin(phase);
  
  cc=pow(bb,(2.0/3.0));
  cc=cc*(part_real+sqrt(3.0)*part_imag);
  cc=cc-aa;
  
  *c1=(9.0*aa+3.0*cc)/8.0;
  *c2= -(aa+3.0*cc)/24.0;
}

int sgfdModel::calcHolbergCoef(float dx,float dy,float dz,float dt,
                    float& sx,float& sy,float& sz,
                    float cx[2],float cy[2],float cz[2],
                    float scalarSpeed,int fdCoeffOrder){
  //Define some constants.
  float c1,c2;
  if(fdCoeffOrder==0)
    holbergCoeffsC(cx[0],&c1,&c2);  //holbergPcent will be in cx[0] here
  else {
    //for now we still assume 4th order
    c1 = cx[0];
    c2 = cx[1];
  }
  sx=(dt/dx)*scalarSpeed;
  sy=(dt/dy)*scalarSpeed;
  sz=(dt/dz)*scalarSpeed;
  
  cx[0]=c1*sx;cx[1]=c2*sx;
  cy[0]=c1*sy;cy[1]=c2*sy;
  cz[0]=c1*sz;cz[1]=c2*sz;
  
  sx=0.5*sx;
  sy=0.5*sy;
  sz=0.5*sz;
  return 6;
}

int sgfdModel::readCDFModel(char* fileName,
                 char* vsFileName,char* rhoFileName,
                 int* lim,
                 float* vp,float* vs,float* rho,
                 float& vMin,float& vMax){
  //Just use readCDFModelInterp since it defaults to a simple read if the
  // _gridMultiplier is 1. Read vp, conditionally read vs, and read rho.
  if(!readCDFModelInterp(fileName,lim,"vp",vp,FALSE))
    assert(readCDFModelInterp(fileName,lim,"c",vp)!=NULL,
           "readCDFModel--unable to read vp from\n\t%s\n",
           fileName);
  if(vs){
    assert(readCDFModelInterp((vsFileName && *vsFileName)?vsFileName:fileName,
                              lim,"vs",vs)!=NULL,
           "readCDFModel--unable to read vs from\n\t%s\n",
           fileName);
  }
  assert(readCDFModelInterp((rhoFileName && *rhoFileName)?rhoFileName:fileName,
                            lim,"rho",rho)!=NULL,
         "readCDFModel--unable to read rho from\n\t%s\n",
         fileName);
  
  return minMaxVelocity(vp,vs,rho,vMin,vMax);
}

int sgfdModel::minMaxVelocity(float* vp,float* vs,float *rho,float &vMin,float &vMax){
  DEF_MODEL_SIZE(modelDef());
  //find the min and max of the parameters
  float vpMin=0.0,vpMax=0.0,vsMax=0.0,vsMin=0.0,rhoMin=0.0,rhoMax=0.0;
  for(int i=0;i<NXYZ;i++){
    if(vs){
      //This is the normal elastic model.
      float alfa=vp[i];
      float beta=vs[i];
      float rhoC=rho[i];
      SETMINMAX(vpMin,vpMax,alfa,i);
      SETMINMAX(vsMin,vsMax,beta,i);
      SETMINMAX(rhoMin,rhoMax,rhoC,i);
      SETMIN(vMin,alfa,i);
      if(beta>0.0)
        SETMIN(vMin,beta,i);
    }else{
      float alfa=vp[i];
      float rhoC=rho[i];
      SETMINMAX(vpMin,vpMax,alfa,i);
      SETMINMAX(rhoMin,rhoMax,rhoC,i);
      SETMIN(vMin,alfa,i);
    }      
  }
  
  vMax=vpMax;
  
  return NXYZ;
}

float* sgfdModel::readCDFModelInterp(char* fileName,
                          int* lim,
                          const char* varName,float* varValues,
                          int noFindError){
  DEF_MODEL_SIZE(modelDef());
  DEF_PARALLEL(parallelDef());
  
  //Check that space has been allocated for the variable.
  assert(varValues!=NULL,
         "readCDFModelInterp--%s not initialized",varName);
  
  if(_gridMultiplier==1){
    return
    readModelVariable(varValues,fileName,varName,
                      procLim,globalNX,globalNY,globalNZ,
                      -1,noFindError);
  }
  
  //Create a local version of the procLimits accounting for the multiplier.
  int rLim[6]={(int)floor((float)lim[0]/_gridMultiplier),(int)ceil((float)lim[1]/_gridMultiplier),
    (int)floor((float)lim[2]/_gridMultiplier),(int)ceil((float)lim[3]/_gridMultiplier),
    (int)floor((float)lim[4]/_gridMultiplier),(int)ceil((float)lim[5]/_gridMultiplier)};
  int rNX=rLim[1]-rLim[0],rNY=rLim[3]-rLim[2],rNZ=rLim[5]-rLim[4];
  int rNXY=rNX*rNY,rNXYZ=rNX*rNY*rNZ;
  
  //Allocate local space to read the values from the file, these values will
  // be interpolated to the fine grid later.
  float* readValues;
  assert((readValues=(float*)malloc(rNXYZ*sizeof(float)))!=NULL,
         "readCDFModelInterp--unable to allocate local array[%i]",
         rNXYZ);
  
  readModelVariable(readValues,fileName,varName,lim,
                    globalNX,globalNY,globalNZ);
  
  //Now do the interpolation.
  for(int k=0;k<NZ;k++){
    int rK=(int)floor((float)k/(float)_gridMultiplier);
    rK=MIN(rK,rNZ-2);
    float ddz=(float)(k-rK*_gridMultiplier)/(float)_gridMultiplier;
    float hdz=1.0-ddz;
    for(int j=0;j<NY;j++){
      int rJ=(int)floor((float)j/(float)_gridMultiplier);
      rJ=MIN(rJ,rNY-2);
      float ddy=(float)(j-rJ*_gridMultiplier)/(float)_gridMultiplier;
      float hdy=1.0-ddy;
      for(int i=0;i<NX;i++){
        int rI=(int)floor((float)i/(float)_gridMultiplier);
        rI=MIN(rI,rNX-2);
        float ddx=(float)(i-rI*_gridMultiplier)/(float)_gridMultiplier;
        float hdx=1.0-ddx;
        
        int index=i+j*NX+k*NXY;
        int rindex=rI+rJ*rNX+rK*rNXY;
        varValues[index]=
        (hdx*hdy*hdz*readValues[rindex           ] +
         ddx*hdy*hdz*readValues[rindex+1         ] +
         hdx*ddy*hdz*readValues[rindex  +rNX     ] +
         ddx*ddy*hdz*readValues[rindex+1+rNX     ] +
         
         hdx*hdy*ddz*readValues[rindex      +rNXY] +
         ddx*hdy*ddz*readValues[rindex+1    +rNXY] +
         hdx*ddy*ddz*readValues[rindex  +rNX+rNXY] +
         ddx*ddy*ddz*readValues[rindex+1+rNX+rNXY]);
      }
    }
  }
  //Free the local memory and return.
  free(readValues);
  
  return varValues;      
}

int sgfdDependent::checkCheckpoint(char* cpDir,int cpID){
  //Since this is the base the file has not been open yet,
  // open the binary file that will hold all the info need for
  // the restart.
  char filename[1024];
  sprintf(filename,"%s/checkpoint_%i.bin",cpDir,cpID);
  FILE* cpFile=fopen(filename,"rb");
  if(cpFile==NULL)
    return FALSE;
  
  //OK, this looks like a checkpoint file. Look for the magic number.
  int cpFlag,returnVal;
  fread(&cpFlag,sizeof(int),1,cpFile);
  if(cpFlag!=CHECKPOINT_START_FILE_FLAG){
    fprintf(stderr,"WARNING--file %s looks like checkpoint but does not start with magic number %i(%i)\n",
            filename,CHECKPOINT_START_FILE_FLAG,cpFlag);
    returnVal=FALSE;
  }else{
    //Now check for an exact command line match.
    int clLen;
    fread(&clLen,sizeof(int),1,cpFile);
    char *clCheck;
    assert((clCheck=(char*)malloc((1+clLen)*sizeof(char)))!=NULL,
           "checkCheckpoint--unable to allocate %i chars for checkpoint command line check",
           clLen);
    fread(clCheck,sizeof(char),clLen+1,cpFile);
    if(strcmp(clCheck,CommandLine)){
      tEprintf(Verbose,
               "Command line mismatch, not using checkpoint from\n\t%s\n",
               cpDir);
      returnVal=FALSE;
    }else{
      returnVal=TRUE;
      
      int iteration;
      fread(&iteration,sizeof(int),1,cpFile);
      fprintf(stderr,
              "Restarting from checkpoint file\n\t%s\n\tIteration: %i\n",
              filename,iteration);
    }
    free(clCheck);
  }
  fclose(cpFile);
  return returnVal;
}

FILE* sgfdDependent::doCheckpoint(int iteration,char* cpDir,int cpID,
                   int callDepth,
                   FILE* cpFile){
  //Save the checkpoint information.
  _checkpointID=cpID;
  strcpy(_checkpointDir,cpDir);
  
  //Since this is the base the file has not been open yet,
  // open the binary file that will hold all the info need for
  // the restart.
  char filename[1024];
  sprintf(filename,"%s/checkpoint_%i.bin",cpDir,cpID);
  cpFile=fopen(filename,"wb");
  assert(cpFile!=NULL,
         "sgfdDependent::doCheckpoint--unable to open file\n\t%s\n",
         filename);
  
  //Write the checkpoint flag.
  int cpFlag=CHECKPOINT_START_FILE_FLAG;
  fwrite(&cpFlag,sizeof(int),1,cpFile);
  
  //Write the command line length and the command line.
  int clLen=strlen(CommandLine);
  fwrite(&clLen,sizeof(int),1,cpFile);
  fwrite(CommandLine,sizeof(char),clLen+1,cpFile);
  
  //Write the iteration.
  fwrite(&iteration,sizeof(int),1,cpFile);
  
  if(!callDepth){
    fclose(cpFile);
    cpFile=NULL;
  }
  return cpFile;
}

FILE* sgfdDependent::readCheckpoint(int& iteration,char* cpDir,int cpID,
                     int callDepth,
                     FILE* cpFile){
  //Save the checkpoint information.
  _checkpointID=cpID;
  strcpy(_checkpointDir,cpDir);
  
  //Since this is the base the file has not been open yet,
  // open the binary file that will hold all the info need for
  // the restart.
  char filename[1024];
  sprintf(filename,"%s/checkpoint_%i.bin",cpDir,cpID);
  cpFile=fopen(filename,"rb");
  assert(cpFile!=NULL,
         "sgfdDependent::readCheckpoint--unable to open file\n\t%s\n",
         filename);
  
  //Read the checkpoint flag.
  int cpFlag;
  fread(&cpFlag,sizeof(int),1,cpFile);
  assert(cpFlag==CHECKPOINT_START_FILE_FLAG,
         "sgfdDependent::readCheckpoint--start flag mismatch %i should be %i in %s",
         cpFlag,CHECKPOINT_START_FILE_FLAG,filename);
  
  //Now check for an exact command line match.
  int clLen;
  fread(&clLen,sizeof(int),1,cpFile);
  char *clCheck;
  assert((clCheck=(char*)malloc((1+clLen)*sizeof(char)))!=NULL,
         "checkCheckpoint--unable to allocate %i chars for checkpoint command line check",
         clLen);
  fread(clCheck,sizeof(char),clLen+1,cpFile);
  assert(!strcmp(clCheck,CommandLine),
         "readCheckpoint--command line mismatch %s!=%s",
         clCheck,CommandLine);
  free(clCheck);
  
  //Read the iteration.
  fread(&iteration,sizeof(int),1,cpFile);
  _checkpointIteration=iteration;
  
  if(!callDepth){
    fclose(cpFile);
    cpFile=NULL;
  }
  return cpFile;
}

void sgfdDependent::doYZSlice(float* var,int yStart,int yEnd,int zStart,int zEnd,
               float* slice,float scalar,
               int i,float a1,float a2,float a3,float a4){
  DEF_MODEL_SIZE(modelDef());
  
  for(int k=zStart,sI=0;k<zEnd;k++){
    for(int j=yStart;j<yEnd;j++){
      int index=i+j*NX+k*NXY;
      float vel= scalar*
      (a1*var[index-1]+a2*var[index  ]+
       a3*var[index+1]+a4*var[index+2]);
      slice[sI++]=vel;
    }
  }
}

void sgfdDependent::doXZSlice(float* var,int xStart,int xEnd,int zStart,int zEnd,
               float* slice,float scalar,
               int j,float b1,float b2,float b3,float b4){
  DEF_MODEL_SIZE(modelDef());
  int NX2=2*NX;
  
  for(int k=zStart,sI=0;k<zEnd;k++){
    for(int i=xStart;i<xEnd;i++){
      int index=i+j*NX+k*NXY;
      float vel= scalar*
      (b1*var[index-NX]+b2*var[index    ]+
       b3*var[index+NX]+b4*var[index+NX2]);
      slice[sI++]=vel;
    }
  }
}

void sgfdDependent::doXZSlice(float (*valFunc)(int,int,int),
               int xStart,int xEnd,int zStart,int zEnd,
               float* slice,
               int j,float b1,float b2,float b3,float b4){
  for(int k=zStart,sI=0;k<zEnd;k++){
    for(int i=xStart;i<xEnd;i++){
      slice[sI++]=
      b1*valFunc(i-1,j,k)+b2*valFunc(i  ,j,k)+
      b3*valFunc(i+1,j,k)+b4*valFunc(i+2,j,k);
    }
  }
}

void sgfdDependent::doXYSlice(float* var,int xStart,int xEnd,int yStart,int yEnd,
               float* slice,float scalar,
               int k,float c1,float c2,float c3,float c4){
  DEF_MODEL_SIZE(modelDef());
  int NXY2=2*NXY;
  
  for(int j=yStart,sI=0;j<yEnd;j++){
    for(int i=xStart;i<xEnd;i++){
      int index=i+j*NX+k*NXY;
      float vel= scalar*
      (c1*var[index-NXY]+c2*var[index     ]+
       c3*var[index+NXY]+c4*var[index+NXY2]);
      slice[sI++]=vel;
    }
  }
}

int masterSgfdModel::resetTimeVector(float minT,int nt,float dt){
  _modelDef->NT=nt;
  _modelDef->dt=dt;
  _modelDef->minT=minT;
  
  char buffer[512];
  strcpy(buffer,"MESSAGE_RESET_TIME_VECTOR");
  initSend();
  sendMessage(AllProcesses,MESSAGE_GENERAL,"s fif",
              buffer,minT,nt,dt);
  return nt;
}

int masterSgfdModel::multiplyGrid(int surfaceBCType,int gridMultiplier){
  sgfdModel::multiplyGrid(surfaceBCType,gridMultiplier);
  
  _parallelDef->globalNX= 1+gridMultiplier*(_parallelDef->globalNX-1);
  _parallelDef->globalNY= 1+gridMultiplier*(_parallelDef->globalNY-1);
  if(!surfaceBCType){
    _parallelDef->globalNZ= 1+gridMultiplier*(_parallelDef->globalNZ-1);
  }else{
    _parallelDef->globalNZ= 1+gridMultiplier*(_parallelDef->globalNZ-2);
  }
  
  return _parallelDef->globalNX*_parallelDef->globalNY*_parallelDef->globalNZ;
}

void masterSgfdModel::getInitAck(int tids[]){
  float rhoMin=0.0,rhoMax=0.0;
  for(int i=0;i<NumProcs;i++){
    float localVMin,localVMax,localRMin,localRMax;
    getMessage(AllProcesses,MESSAGE_INITIALIZE,"ff ff",
               &localVMin,&localVMax,&localRMin,&localRMax);
    SETMIN(_vMin,localVMin,i);
    SETMAX(_vMax,localVMax,i);
    
    SETMIN(rhoMin,localRMin,i);
    SETMAX(rhoMax,localRMax,i);
  }
  //and return the global min and max to the slaves as well as printing a diagnostic
  sendMessage(AllProcesses,MESSAGE_INITIALIZE,"ff",_vMin,_vMax);
  tEprintf(Verbose,
           "Model read by slaves; global min vel %.1f; max vel %.1f\n",
           _vMin,_vMax);
  tEprintf(Verbose,
           "  Rho %.2f-%.2f\n",
           rhoMin,rhoMax);
}

int masterSgfdModel::initializeSlaves(int tids[],
                     char* modelName,int fdCoeffOrder,float* fdCoeffs,
                     int longRead){
  //Define a whole lot of locals
  DEF_MODEL_SIZE(_modelDef);
  
  DEF_PARALLEL(_parallelDef);
  int nxyzProc=nxyProc*nzProc;
  
  //need to fill out modelDefs and parallelDefs for each domain
  modelDefStruct** slaveModelDefs=buildSlaveModelDef();
  
  parallelDefStruct** slaveParallelDefs;
  assert((slaveParallelDefs=
          (parallelDefStruct**)malloc(nxyzProc*sizeof(parallelDefStruct*)))!=NULL,
         "initializeSlaves--unable to allocate %i parallelDefStruct*'s",
         nxyzProc);
  
  //Now a second trip to fill in the slaveParallelDefs; this requires that
  // the modelDefs are completed so it must be done in two passes
  for(int i=0;i<nxProc;i++){
    for(int j=0;j<nyProc;j++){
      for(int k=0;k<nzProc;k++){
        //initialize my parallelDef
        SUBDOMAIN(slaveParallelDefs,i,j,k)=
        newParallelDef(NX,NY,NZ,nxProc,nyProc,nzProc,i,j,k);
        
        //fill in the neighbor and neighborDef fields
        // this is the part that requires a second pass through
        // the nested loops
        setNeighbors(SUBDOMAIN(slaveParallelDefs,i,j,k),
                     tids,slaveModelDefs,
                     nxProc,nyProc,nzProc,i,j,k);
        
        //now send the actual data required for initialization
        sendInitMessage(modelName,fdCoeffOrder,fdCoeffs,longRead,
                        slaveParallelDefs,
                        tids,i,j,k,
                        nxProc,nxyProc);
      }
    }
  }
  
  //don't need the slave modelDefs and parallelDefs anymore
  for(int i=0;i<NumProcs;i++){
    free(slaveModelDefs[i]);
    free(slaveParallelDefs[i]);
  }
  free(slaveModelDefs);
  free(slaveParallelDefs);
  
  //And get an acknoledgement from the slaves that (this portion?) of the
  // initialization is complete.
  getInitAck(tids);
  return NumProcs;
}

modelDefStruct** masterSgfdModel::buildSlaveModelDef(){
  DEF_MODEL_SIZE(_modelDef);
  DEF_MODEL_LIMITS(_modelDef);
  DEF_MODEL_SCALE(_modelDef);
  
  DEF_PARALLEL(_parallelDef);
  int nxyzProc=nxyProc*nzProc;
  modelDefStruct** slaveModelDefs;
  assert((slaveModelDefs=(modelDefStruct**)malloc(nxyzProc*sizeof(modelDefStruct*)))!=NULL,
         "initializeSlaves--unable to allocate %i modelDefStruct*'s",
         nxyzProc);
  //First trip through the nested loop fill in the slaveModelDefs
  for(int i=0;i<nxProc;i++){
    for(int j=0;j<nyProc;j++){
      for(int k=0;k<nzProc;k++){
        //calculate the breaks between the processors
        int xStart=MAX(0,i*NX/nxProc-2);
        int xStop =MIN(NX,(i+1)*NX/nxProc+2);
        int yStart=MAX(0,j*NY/nyProc-2);
        int yStop =MIN(NY,(j+1)*NY/nyProc+2);
        int zStart=MAX(0,k*NZ/nzProc-2);
        int zStop =MIN(NZ,(k+1)*NZ/nzProc+2);
        
        //Allocate the modelDefStruct for this process.
        SUBDOMAIN(slaveModelDefs,i,j,k)=
        newModelDef(xStart,xStop,yStart,yStop,zStart,zStop,NT,
                    minX,dx,minY,dy,minZ,dz,minT,dt,
                    scalarSpeed,scalarVel,scalarDen,scalarStress);
        _maxSlaveNX = SUBDOMAIN(slaveModelDefs,i,j,k)->NX>_maxSlaveNX?SUBDOMAIN(slaveModelDefs,i,j,k)->NX:_maxSlaveNX;
        _maxSlaveNY = SUBDOMAIN(slaveModelDefs,i,j,k)->NY>_maxSlaveNY?SUBDOMAIN(slaveModelDefs,i,j,k)->NY:_maxSlaveNY;
        _maxSlaveNZ = SUBDOMAIN(slaveModelDefs,i,j,k)->NZ>_maxSlaveNZ?SUBDOMAIN(slaveModelDefs,i,j,k)->NZ:_maxSlaveNZ;
      }
    }
  }
  return slaveModelDefs;
}

int masterSgfdModel::setNeighbors(parallelDefStruct* parallelDef,
                 int tids[],modelDefStruct* slaveModelDefs[],
                 int nxProc,int nyProc,int nzProc,
                 int i,int j,int k){
  int nxyProc=nxProc*nyProc; //can't use DEF_PARALLEL here (shadows neighbors)
  int count=0;
  for(int n=k-1;n<=k+1;n++){
    for(int m=j-1;m<=j+1;m++){
      for(int l=i-1;l<=i+1;l++){
        int index=(l-i+1)+3*(m-j+1)+9*(n-k+1);
        //neighbors should only be directly adjacent for elastic/anelastic/FM acoustic
        //index==13 will be self.  Commented out and put this logic in acousticBoundary
        //if((abs(l-i)+abs(m-j)+abs(n-k)==1 || index==13) &&
        if(ISMID(0,l,nxProc-1) &&
           ISMID(0,m,nyProc-1) &&
           ISMID(0,n,nzProc-1)){
          parallelDef->neighbors[index]=
          SUBDOMAIN(tids,l,m,n);
          parallelDef->neighborDefs[index]=
          SUBDOMAIN(slaveModelDefs,l,m,n);
          
          count++;
        }else{
          parallelDef->neighbors[index]=-1;
          parallelDef->neighborDefs[index]=NULL;
        }
      }
    }
  }
  return count;
}

void masterSgfdModel::sendInitMessage(char* modelName,int fdCoeffOrder,float* fdCoeffs,int longRead,
                     parallelDefStruct** slaveParallelDefs,
                     int tids[],int i,int j,int k,
                     int nxProc,int nxyProc,
                     int coworkParadyme){
  initSend();
  packMessage("s iiF ii",
              modelName,Verbose,fdCoeffOrder,fdCoeffs,fdCoeffOrder>0?fdCoeffOrder:1,_gridMultiplier,longRead);
  
  //pack the parallelDef, since I am the center of the
  // neighbors and neighborDefs arrays this is all the info
  // required by the slave process
  packParallelDef(SUBDOMAIN(slaveParallelDefs,i,j,k));
  //send the message
  sendMessage(coworkParadyme?COWORKTID(i,j,k):TID(i,j,k),
              MESSAGE_INITIALIZE,NULL);
}

slaveSgfdModel::slaveSgfdModel(int convertParams):sgfdModel(){
  //get the model and initial condition filenames, these will provide
  // much of the information I need about myself.
  int fdCoeffOrder;
  unpackMessage("s ii",
                _modelName,&Verbose,&fdCoeffOrder);
  unpackMessage("F ii",
                _cx,fdCoeffOrder>0?fdCoeffOrder:1,&_gridMultiplier,&_doLongRead);
  
  //unpack parallelDef from the message
  _parallelDef=
  unpackParallelDef();
  
  //my modeldef is the center of the _parallelDef->neighborDefs
  _modelDef=BOUNDARY(_parallelDef->neighborDefs,0,0,0);
  deassert(modelDef()->NX<8 || modelDef()->NY<8 || modelDef()->NZ<8,
           "Grid size %ix%ix%i too small must be at least 8x8x8\n",
           modelDef()->NX,modelDef()->NY,modelDef()->NZ);
  
  //Now define variables for the rest of the initialization
  DEF_MODEL_LIMITS(modelDef());
  
  //Calculate the numerical differentiator coefficients.
  calcHolbergCoef(dx,dy,dz,dt,_sx,_sy,_sz,_cx,_cy,_cz,
                  modelDef()->scalarSpeed,fdCoeffOrder);
}

int slaveSgfdModel::resetTimeVector(){
  float minT,dt;
  int nt;
  unpackMessage("fif",&minT,&nt,&dt);
  
  _cx[0]*=dt/_modelDef->dt;
  _cx[1]*=dt/_modelDef->dt;
  _cy[0]*=dt/_modelDef->dt;
  _cy[1]*=dt/_modelDef->dt;
  _cz[0]*=dt/_modelDef->dt;
  _cz[1]*=dt/_modelDef->dt;
  
  _modelDef->NT=nt;
  _modelDef->dt=dt;
  _modelDef->minT=minT;
  
  return nt;
}

int slaveSgfdModel::isInside(float x,float y,float z,
             int &i,int &j,int &k){
  DEF_MODEL_SIZE(modelDef());
  DEF_MODEL_LIMITS(modelDef());
  
  float maxX=minX+dx*(NX-1),maxY=minY+dy*(NY-1),maxZ=minZ+dz*(NZ-1);
  if(!ISMID(minX,x,maxX) ||
     !ISMID(minY,y,maxY) ||
     !ISMID(minZ,z,maxZ))
    return FALSE;
  
  i=(int)floor((float)(x-minX)/dx);
  j=(int)floor((float)(y-minY)/dy);
  k=(int)floor((float)(z-minZ)/dz);
  return TRUE;
}

int slaveSgfdModel::isInside3(float x,float y,float z){
  DEF_MODEL_SIZE(modelDef());
  DEF_MODEL_LIMITS(modelDef());
  DEF_PARALLEL(_parallelDef);
  float maxX=minX+(NX-1)*dx,maxY=minY+(NY-1)*dy,maxZ=minZ+(NZ-1)*dz;
  
  if(NEIGHBOR(-1,0,0)!=-1 && x<minX)
    return FALSE;
  if(NEIGHBOR( 1,0,0)!=-1 && x>maxX)
    return FALSE;
  
  if(NEIGHBOR(0,-1,0)!=-1 && y<minY)
    return FALSE;
  if(NEIGHBOR(0, 1,0)!=-1 && y>maxY)
    return FALSE;
  
  if(NEIGHBOR(0,0,-1)!=-1 && z<minZ)
    return FALSE;
  if(NEIGHBOR(0,0, 1)!=-1 && z>maxZ)
    return FALSE;
  
  return TRUE;
}

float slaveSgfdModel::interp3(float *var,
              float x,float y,float z){
  DEF_MODEL_SIZE(modelDef());
  DEF_MODEL_LIMITS(modelDef());
  float maxX=minX+(NX-1)*dx,maxY=minY+(NY-1)*dy,maxZ=minZ+(NZ-1)*dz;
  if(x<=minX)
    return interp3(var,minX+dx,y,z);
  if(x>=maxX)
    return interp3(var,maxX-dx,y,z);
  
  if(y<minY)
    return interp3(var,x,minY+dy,z);
  if(y>=maxY)
    return interp3(var,x,maxY-dy,z);
  
  if(z<minZ)
    return interp3(var,x,y,minZ+dz);
  if(z>=maxZ)
    return interp3(var,x,y,maxZ-dz);
  
  return interp2(var,x,y,z);
}

float slaveSgfdModel::interp2(float *var,
              float x,float y,float z){
  DEF_MODEL_SIZE(modelDef());
  DEF_MODEL_LIMITS(modelDef());
  
  int i,j,k;
  assert(isInside(x,y,z,i,j,k),
         "slaveSgfdModel::interp2--point %.2f,%.2f,%.2f; not in region");
  int index=i+j*NX+k*NXY;
  
  float ddx=x-(minX+dx*i);
  float hdx=dx-ddx;
  float ddy=y-(minY+dy*j);
  float hdy=dy-ddy;
  float ddz=z-(minZ+dz*k);
  float hdz=dz-ddz;
  
  float sum=0.0;
  sum+=hdx*hdy*hdz*var[index         ];
  sum+=ddx*hdy*hdz*var[index+1       ];
  sum+=hdx*ddy*hdz*var[index  +NX    ];
  sum+=ddx*ddy*hdz*var[index+1+NX    ];
  
  sum+=hdx*hdy*ddz*var[index     +NXY];
  sum+=ddx*hdy*ddz*var[index+1   +NXY];
  sum+=hdx*ddy*ddz*var[index  +NX+NXY];
  sum+=ddx*ddy*ddz*var[index+1+NX+NXY];
  
  return sum/dx/dy/dz;
}

masterSgfdDependent::~masterSgfdDependent(){
  if(_receivers){
    delete _receivers;
    _receivers=NULL;
  }
  if(_sources){
    delete _sources;
    _sources=NULL;
  }
}

FILE* masterSgfdDependent::doCheckpoint(int iteration,char* cpDir,int cpID,
                   int callDepth,
                   FILE* cpFile){
  //Use the parent to open the file and required whatever information
  // it can write.
  if(!cpFile){
    for(int i=0;i<NumProcs;i++){
      initSend();
      sendMessage(Tids[i],MESSAGE_GENERAL,
                  "s isi","DO_CHECKPOINT",
                  iteration,cpDir,cpID+i+1);
    }
    
    cpFile=sgfdDependent::doCheckpoint(iteration,cpDir,cpID,
                                       callDepth+1);
  }
  if(_receivers)
    _receivers->doCheckpoint(iteration,cpDir,cpID,0);
  
  if(!callDepth){
    fclose(cpFile);
    cpFile=NULL;
  }
  return cpFile;
}

FILE* masterSgfdDependent::readCheckpoint(int& iteration,char* cpDir,int cpID,
                     int callDepth,
                     FILE* cpFile){
  //Use the parent to open the file and required whatever information
  // it can write.
  if(!cpFile){
    for(int i=0;i<NumProcs;i++){
      initSend();
      sendMessage(Tids[i],MESSAGE_GENERAL,
                  "s isi","READ_CHECKPOINT",
                  iteration,cpDir,cpID+i+1);
    }
    
    cpFile=sgfdDependent::readCheckpoint(iteration,cpDir,cpID,
                                         callDepth+1);
  }
  
  if(_receivers)
    _receivers->readCheckpoint(iteration,cpDir,cpID,0);
  
  if(!callDepth){
    fclose(cpFile);
    cpFile=NULL;
  }
  return cpFile;
}

int masterSgfdDependent::setInitialConditions(char* initialConditionsName,int& iteration,
                         int negateVelocities){
  //Send initial condition message.
  initSend();
  
  if(*initialConditionsName &&
     !strstr(initialConditionsName,".cdf\0")){
    char buffer[512]="CHECKPOINT_RESTART";
    sendMessage(AllProcesses,MESSAGE_INITIALIZE,"s",buffer);
    readCheckpoint(iteration,initialConditionsName,
                   0,0);
    return TRUE;
  }
  
  //Name of ic, there will be more here later.
  packMessage("s i",initialConditionsName,
              negateVelocities);
  //Other optional things to do to the IC's.
  packMessage("i",negateVelocities);
  
  sendMessage(AllProcesses,MESSAGE_INITIALIZE,NULL);
  
  //Look for an iteration in the IC's, if present then these were saved as
  // a checkpoint.
  if(initialConditionsName[0]){
    int infile=openCDFFile(initialConditionsName,FALSE,NC_NOWRITE);
    int varID;
    if(nc_inq_varid(infile,"iteration",&varID)==NC_NOERR){
      //For now just read the iteration, in the future also want to save/read the
      // receiver data so we can do a restart and get the same results.
      assert(nc_get_var_int(infile,varID,&iteration)==NC_NOERR,
             "setInitialConditions--unable to get var iteration from %s(%i,%i)",
             initialConditionsName,infile,varID);
      tEprintf(Verbose,"Starting at iteration %i\n",
               iteration);
    }
    assert(nc_close(infile)==NC_NOERR,
           "setInitialConditions--unable to close %s(%i)",
           initialConditionsName,infile);
  }
  
  //Wait for acknowlegement from slaves.
  for(int i=0;i<NumProcs;i++)
    getMessage(AllProcesses,MESSAGE_INITIALIZE,NULL);
  
  return initialConditionsName[0]!='\0';
}

int masterSgfdDependent::setBoundaryConditions(int surfaceBCMode,
                          int spongeBCNodes[6],float spongeBCValue[6]){
  tEprintf(Verbose,"Setting boundary conditions\n");
  initSend();
  
  //set any special boundry conditions
  // free surface and TDBC
  packMessage("i",
              surfaceBCMode);
  //sponge
  packMessage("IF",
              spongeBCNodes,6,spongeBCValue,6);
  
  sendMessage(AllProcesses,MESSAGE_SET_BC,NULL);
  for(int i=0;i<NumProcs;i++)
    getMessage(AllProcesses,MESSAGE_SET_BC,NULL);
  
  return NumProcs;
}

int masterSgfdDependent::setReceivers(masterReceiverNetwork* extraReceivers,
                 int surfaceMode,
                 int writeTraceFile,int doCheckActive){
  //Read receivers
  allocateReceivers(modelDef(),modelName(),surfaceMode,
                    extraReceivers?(extraReceivers->_timeSubsample):1,
                    extraReceivers?(extraReceivers->_useModelReceivers):TRUE);
  
  //add any receivers specified on the command line.
  _receivers->addReceivers(extraReceivers);
  _receivers->initializeNetwork(modelDef(),TRUE);
  if(doCheckActive)
    _receivers->checkActive(TRUE);
  
  if(writeTraceFile && _receivers->size()){
    _receivers->writeHeader(modelDef(),TRUE);
    if(_sources &&  _sources->size() && _sources->_bandCheckPCent>0 &&
       !_receivers->_splitFiles){
      _sources->writeSources(modelDef(),_receivers->_traceOutputName);
    }
  }
  
  return _receivers->size();
}

int masterSgfdDependent::setSources(sourceNetwork* extraSources){
  tEprintf(Verbose,"Setting sources\n");
  //read the sources from the model and send to out. After being sent the sources
  // are no longer required and can be deleted.
  _sources=
  new movingSourceNetwork(modelDef(),modelName(),
                          extraSources?(extraSources->_useModelSources):TRUE);
  if(extraSources)
    _sources->addSources(extraSources);
  
  _sources->checkDispersion(modelDef(),_model->vMin(),_model->vMax());
  
  int num=_sources->size();
  //Prepare to send info to the slave processes
  initSend();
  _sources->sendSources(modelDef());
  tEprintf(Verbose,"Set %i sources (%i force, %i moment)\n",
           _sources->size(),_sources->nfs(),_sources->nms());
  
  return num;
}

int masterSgfdDependent::writeTraces(int startI,int stopI){
  return _receivers->writeTraces(modelDef(),startI,stopI);
}

int masterSgfdDependent::advance(int& iteration,int count,int acknowledge){
  //send out the work
  initSend();
  sendMessage(AllProcesses,MESSAGE_ADVANCE_ONE_ITERATION,"ii",
              iteration,count);
  
  if(acknowledge)
    getAcknowledgement();
  
  return iteration+=count;
}

int masterSgfdDependent::advanceVel(int iteration,int acknowledge){
  //send out the work
  initSend();
  sendMessage(AllProcesses,MESSAGE_UPDATE_VEL,"i",iteration);
  
  // (2) Recieve the results
  for(int i=0;i<NumProcs;i++){
    checkForError();
    getMessage(AllProcesses,MESSAGE_UPDATE_VEL,NULL);
  }
  return iteration;
}

int masterSgfdDependent::advanceStress(int iteration,int acknowledge){
  //send out the work
  initSend();
  sendMessage(AllProcesses,MESSAGE_UPDATE_STRESS,"i",iteration);
  
  // (2) Recieve the results
  for(int i=0;i<NumProcs;i++){
    checkForError();
    getMessage(AllProcesses,MESSAGE_UPDATE_STRESS,NULL);
  }
  return iteration;
}

void masterSgfdDependent::getAcknowledgement(){
  doBarrier();
}

int masterSgfdDependent::sendInitialMessage(){
  //Send init message and wait for acknowledgment
  initSend();
  sendMessage(AllProcesses,MESSAGE_INITIALIZE,NULL);
  //     for(int i=0;i<NumProcs;i++)
  //       getMessage(AllProcesses,MESSAGE_INITIALIZE,NULL);
  doBarrier();
  return NumProcs;
}

masterReceiverNetwork* masterSgfdDependent::allocateReceivers(modelDefStruct *modelDef,
                                         char* modelname,int surfaceMode,
                                         int receiverDecimate,
                                         int useModelReceivers){
#ifdef SGFD_INVERSION
  _receivers=new
  masterInversionReceiverNetwork(modelDef,modelname,surfaceMode,
                                 useModelReceivers,receiverDecimate);
  masterSgfdDependent::_receivers=_receivers;
#else
  _receivers=new
  masterReceiverNetwork(modelDef,modelname,surfaceMode,
                        useModelReceivers,receiverDecimate);
#endif //#ifdef SGFD_INVERSION
  return _receivers;
}

slaveSgfdDependent::~slaveSgfdDependent(){
  if(_slice) free(_slice);
  
  if(_sources)
    delete _sources;
  if(_receivers)
    delete _receivers;
}

FILE* slaveSgfdDependent::doCheckpoint(int iteration,char* cpDir,int cpID,
                   int callDepth,
                   FILE* cpFile){
  //Use the parent class to open the file and required whatever information
  // it can write.
  if(!cpFile){
    cpFile=sgfdDependent::doCheckpoint(iteration,cpDir,cpID,
                                       callDepth+1);
  }
  if(_receivers)
    _receivers->doCheckpoint(iteration,cpDir,cpID,0);
  
  if(!callDepth){
    fclose(cpFile);
    cpFile=NULL;
  }
  return cpFile;
}

FILE* slaveSgfdDependent::readCheckpoint(int& iteration,char* cpDir,int cpID,
                     int callDepth,
                     FILE* cpFile){
  //Use the parent class to open the file and required whatever information
  // it can write.
  if(!cpFile){
    cpFile=sgfdDependent::readCheckpoint(iteration,cpDir,cpID,
                                         callDepth+1);
  }
  
  if(_receivers)
    _receivers->readCheckpoint(iteration,cpDir,cpID,0);
  
  if(!callDepth){
    fclose(cpFile);
    cpFile=NULL;
  }
  return cpFile;
}

int slaveSgfdDependent::setSources(float* rho,int useCubicExtrap){
  getMessage(Parent,MESSAGE_SET_SOURCES,NULL);
  if(_sources) delete _sources;
  _sources=new movingSourceNetwork(_model->modelDef(),rho,
                                   useCubicExtrap);
  
  //reply complete
  sendMessage(Parent,MESSAGE_SET_SOURCES,NULL);
  return _sources->size();
}

int slaveSgfdDependent::setReceivers(int useCubicExtrap){
  if(_receivers) delete _receivers;
  allocateReceivers(modelDef(),TRUE,TRUE,useCubicExtrap);
  
  return _receivers->size();
}

int slaveSgfdDependent::advance(int& iteration,int count,int acknowledge){
  for(int i=0;i<count;i++){
    advanceVel(iteration);
    advanceStress(iteration,FALSE);
    iteration++;
  }
  if(acknowledge) sendAcknowledgement();
  return iteration;
}

int slaveSgfdDependent::packReceiverData(int index,int startI,int stopI){
  assert(_receivers->active(index),
         "slaveDependent::packReceiver--receiver %i is not active in this process",
         index);
  (*_receivers)[index]->packData(modelDef(),_receivers->work(modelDef()),
                                 startI,stopI,_receivers->timeSubsample());
  return stopI-startI;
}

int slaveSgfdDependent::sendReceiverGridData(int index){
  int startI=0,stopI=0;
  (*_receivers->_receiverGrids)[index]->
  packData(modelDef(),_receivers->work(modelDef()),
           startI,stopI,_receivers->timeSubsample());
  return stopI-startI;
}

int slaveSgfdDependent::fullFieldOutput(int recepient,float* vx){
  DEF_MODEL_SIZE(modelDef());
  
  for(int k=0;k<NZ;k++){
    initSend();
    sendMessage(recepient,MESSAGE_FULL_FIELD,"F",&vx[k*NXY],NXY);
  }
  return NZ;
}

int slaveSgfdDependent::slicer(int recepient,
           float *vx,float *vy,float *vz,
           float *xx,float *yy,float *zz,
           float *xy,float *xz,float *yz){
  //read the parameters for the slice and prepare to send the return
  // message
  int plane,comp;
  float coord;
  unpackMessage("iif",&plane,&comp,&coord);
  return
  slicer(recepient,plane,comp,coord,
         vx,vy,vz,
         xx,yy,zz,xy,xz,yz);
}

int slaveSgfdDependent::slicer(int recepient,
           int plane,int comp,float coord,
           float *vx,float *vy,float *vz,
           float *xx,float *yy,float *zz,
           float *xy,float *xz,float *yz){
  DEF_MODEL_SIZE(_model->modelDef());
  DEF_MODEL_LIMITS(_model->modelDef());
  DEF_PARALLEL(_model->parallelDef());
  
  //Calculate the slice size.
  int sliceXStart=procI?2:0,sliceXEnd=NX-((procI==nxProc-1)?0:2);
  int sliceYStart=procJ?2:0,sliceYEnd=NY-((procJ==nyProc-1)?0:2);
  int sliceZStart=procK?2:0,sliceZEnd=NZ-((procK==nzProc-1)?0:2);
  
  int
  sliceNX=sliceXEnd-sliceXStart,
  sliceNY=sliceYEnd-sliceYStart,
  sliceNZ=sliceZEnd-sliceZStart;
  
  //create variables for the local grid limits
  float xmin_vx=minX+0.5*dx;
  float ymin_vy=minY+0.5*dy;
  float zmin_vz=minZ+0.5*dz;
  
  int currSliceLength;
  //Slice wavefield along an yz-plane (fixed x-coordinate).
  if (plane==SLICE_YZ_PLANE){
    //set coords and check limits
    int i;
    float xi;
    if(comp!=SLICE_VX_COMP){
      i=(int)floor((coord-minX)/dx);
      xi=minX+i*dx;
    }else{
      //reset coords for 1/2 node offset on vx
      i=(int)floor((coord-xmin_vx)/dx);
      xi=xmin_vx+i*dx;
    }
    
    //check if in bounds for this process
    if(i-1 < 0 ||
       i+2 >= NX-1){
      //send the failure message to the slice recipent
      getMessage(Parent,MESSAGE_SLICE,NULL);
      initSend();
      sendMessage(recepient,MESSAGE_SLICE,"i",FALSE);
      return FALSE;
    }
    
    float p=(coord-xi)/dx;
    
    float a1=-p*(p-1.0)*(p-2.0)/6.0;
    float a2=(p-1.0)*(p+1.0)*(p-2.0)/2.0;
    float a3=-p*(p+1.0)*(p-2.0)/2.0;
    float a4= p*(p-1.0)*(p+1.0)/6.0;
    
    //Allocate (if required) and zero memory for the slice.
    currSliceLength=sliceNY*sliceNZ;
    _slice=(float*)allocateEnough(_slice,_sliceLength,
                                  currSliceLength,sizeof(float));
    for(int ii=0;ii<currSliceLength;_slice[ii++]=0.0);
    
    sliceYZCompValue(comp,i,
                     sliceYStart,sliceYEnd,sliceZStart,sliceZEnd,
                     a1,a2,a3,a4,
                     vx,vy,vz,xx,yy,zz,xy,xz,yz);
  }else if(plane==SLICE_XZ_PLANE){
    //Slice wavefield along an xz-plane (fixed y-coordinate).
    //calculate coordinates
    int j;
    float yj;
    if(comp!=SLICE_VY_COMP){
      j=(int)floor((coord-minY)/dy);
      yj=minY+j*dy;
    }else{
      //reset the coords to account for the 1/2 node offset
      j=(int)floor((coord-ymin_vy)/dy);
      yj=ymin_vy+j*dy;
    }
    
    if(j-1 < 0 ||
       j+2 >= NY-1){
      //send the failure message to the slice recipent
      getMessage(Parent,MESSAGE_SLICE,NULL);
      initSend();
      sendMessage(recepient,MESSAGE_SLICE,"i",FALSE);
      return FALSE;
    }
    
    float q=(coord-yj)/dy;
    
    float b1=-q*(q-1.0)*(q-2.0)/6.0;
    float b2=(q-1.0)*(q+1.0)*(q-2.0)/2.0;
    float b3=-q*(q+1.0)*(q-2.0)/2.0;
    float b4= q*(q-1.0)*(q+1.0)/6.0;
    
    //allocate and zero memory for the slice
    currSliceLength=sliceNX*sliceNZ;
    _slice=(float*)allocateEnough(_slice,_sliceLength,
                                  currSliceLength,sizeof(float));
    for(int ii=0;ii<currSliceLength;_slice[ii++]=0.0);
    
    //And fill the slice with the correct values for this component.
    sliceXZCompValue(comp,j,
                     sliceXStart,sliceXEnd,sliceZStart,sliceZEnd,
                     b1,b2,b3,b4,
                     vx,vy,vz,xx,yy,zz,xy,xz,yz);
  }else if(plane==SLICE_XY_PLANE){
    //Slice the wavefield along an xy-plane (fixed z-coordinate).
    //Set the coords for all elems expect vz
    int k;
    float zk;
    if(comp!=SLICE_VZ_COMP){
      k=(int)floor((coord-minZ)/dz);
      zk=minZ+k*dz;
    }else{
      //reset coords to account for 1/2 node offset
      k=(int)floor((coord-zmin_vz)/dz);
      zk=zmin_vz+k*dz;
    }
    
    if(k-1 < 0 ||
       k+2 >= NZ-1){
      //send the failure message to the slice recipent
      getMessage(Parent,MESSAGE_SLICE,NULL);
      initSend();
      sendMessage(recepient,MESSAGE_SLICE,"i",FALSE);
      return FALSE;
    }
    
    float r=(coord-zk)/dz;
    
    float c1=-r*(r-1.0)*(r-2.0)/6.0;
    float c2=(r-1.0)*(r+1.0)*(r-2.0)/2.0;
    float c3=-r*(r+1.0)*(r-2.0)/2.0;
    float c4= r*(r-1.0)*(r+1.0)/6.0;
    
    //allocate and zero memory for the slice
    currSliceLength=sliceNX*sliceNY;
    _slice=(float*)allocateEnough(_slice,_sliceLength,
                                  currSliceLength,sizeof(float));
    for(int ii=0;ii<currSliceLength;_slice[ii++]=0.0);
    
    sliceXYCompValue(comp,k,
                     sliceXStart,sliceXEnd,
                     sliceYStart,sliceYEnd,
                     c1,c2,c3,c4,
                     vx,vy,vz,xx,yy,zz,xy,xz,yz);
  }
  
  //The _slice should now be filled in with the correct info; send it
  // to the recepient
  getMessage(Parent,MESSAGE_SLICE,NULL);
  initSend();
  sendMessage(recepient,MESSAGE_SLICE,"i iiiiii F",TRUE,
              procLim[0]+sliceXStart,procLim[0]+sliceXEnd,
              procLim[2]+sliceYStart,procLim[2]+sliceYEnd,
              procLim[4]+sliceZStart,procLim[4]+sliceZEnd,
              _slice,currSliceLength);
  
  return currSliceLength;
}

slaveReceiverNetwork* slaveSgfdDependent::allocateReceivers(modelDefStruct* modelDef,
                                        int allocateData,
                                        int mesgFromParent,
                                        int useCubicExtrap){
#ifdef SGFD_INVERSION
  _receivers=
  new slaveInversionReceiverNetwork(modelDef,allocateData,mesgFromParent,
                                    useCubicExtrap);
  slaveSgfdDependent::_receivers=_receivers;
#else
  _receivers=
  new slaveReceiverNetwork(modelDef,allocateData,
                           mesgFromParent,useCubicExtrap);
#endif
  return _receivers;
}

void slaveSgfdDependent::sendAcknowledgement(){
  doBarrier();
}

void slaveSgfdDependent::sliceXZCompValue(int comp,int j,
                      int sliceXStart,int sliceXEnd,
                      int sliceZStart,int sliceZEnd,
                      float b1,float b2,float b3,float b4,
                      float *vx,float *vy,float *vz,
                      float *xx,float *yy,float *zz,
                      float *xy,float *xz,float *yz){
  DEF_MODEL_SIZE(modelDef());
  DEF_MODEL_SCALE(modelDef());
  DEF_MODEL_LIMITS(modelDef());
  
  switch(comp){
    case SLICE_VX_COMP:
      //x-component of velocity.
      doXZSlice(vx,sliceXStart,sliceXEnd,sliceZStart,sliceZEnd,
                _slice,scalarVel,j,b1,b2,b3,b4);
      break;
      
    case SLICE_VY_COMP:
      //y-component of velocity.
      doXZSlice(vy,sliceXStart,sliceXEnd,sliceZStart,sliceZEnd,
                _slice,scalarVel,j,b1,b2,b3,b4);
      break;
      
    case SLICE_VZ_COMP:
      //z-component of velocity.
      doXZSlice(vz,sliceXStart,sliceXEnd,sliceZStart,sliceZEnd,
                _slice,scalarVel,j,b1,b2,b3,b4);
      break;
      
    case SLICE_PRESS_COMP:
      //Acoustic pressure.
      for(int k=sliceZStart,sI=0;k<sliceZEnd;k++){
        for(int i=sliceXStart;i<sliceXEnd;i++){
          if(yy){
            //Elastic problem.
            float aa=XX(i,j-1,k)+YY(i,j-1,k)+ZZ(i,j-1,k);
            float bb=XX(i,j  ,k)+YY(i,j  ,k)+ZZ(i,j  ,k);
            float cc=XX(i,j+1,k)+YY(i,j+1,k)+ZZ(i,j+1,k);
            float dd=XX(i,j+2,k)+YY(i,j+2,k)+ZZ(i,j+2,k);
            
            float press=b1*aa+b2*bb+b3*cc+b4*dd;
            _slice[sI++]=-press*scalarStress/3.0;
          }else{
            //Acoustic problem xx is already the pressure.
            _slice[sI++]= scalarStress*
            (b1*XX(i,j-1,k)+b2*XX(i,j  ,k)+
             b3*XX(i,j+1,k)+b4*XX(i,j+2,k));
          }
        }
      }
      break;
      
    case SLICE_RX_COMP:
      if(j<NY-3) {
        if(sliceZEnd>NZ-1) sliceZEnd = NZ-1;
        for(int k=sliceZStart,sI=0;k<sliceZEnd;k++) {
          for(int i=sliceXStart;i<sliceXEnd;i++) {
            float aa=((VZ(i,j  ,k)-VZ(i,j-1,k))/dy-(VY(i,j-1,k+1)-VY(i,j-1,k))/dz);
            float bb=((VZ(i,j+1,k)-VZ(i,j  ,k))/dy-(VY(i,j  ,k+1)-VY(i,j  ,k))/dz);
            float cc=((VZ(i,j+2,k)-VZ(i,j+1,k))/dy-(VY(i,j+1,k+1)-VY(i,j+1,k))/dz);
            float dd=((VZ(i,j+3,k)-VZ(i,j+2,k))/dy-(VY(i,j+2,k+1)-VY(i,j+2,k))/dz);
            _slice[sI++]=scalarVel*(b1*aa+b2*bb+b3*cc+b4*dd);
          }
        }
      }
      break;
      
    case SLICE_RY_COMP:
      if(sliceZEnd>NZ-1) sliceZEnd = NZ-1;
      if(sliceXEnd>NX-1) sliceXEnd = NX-1;
      for(int k=sliceZStart,sI=0;k<sliceZEnd;k++) {
        for(int i=sliceXStart;i<sliceXEnd;i++) {
          float aa = ((VX(i,j-1,k+1)-VX(i,j-1,k))/dz-(VZ(i+1,j-1,k)-VZ(i,j-1,k))/dx);
          float bb = ((VX(i,j  ,k+1)-VX(i,j  ,k))/dz-(VZ(i+1,j  ,k)-VZ(i,j  ,k))/dx);
          float cc = ((VX(i,j+1,k+1)-VX(i,j+1,k))/dz-(VZ(i+1,j+1,k)-VZ(i,j+1,k))/dx);
          float dd = ((VX(i,j+2,k+1)-VX(i,j+2,k))/dz-(VZ(i+1,j+2,k)-VZ(i,j+2,k))/dx);
          _slice[sI++]=scalarVel*(b1*aa+b2*bb+b3*cc+b4*dd);
        }
        if(sliceXEnd==NX-1) sI++;
      }
      break;
      
    case SLICE_RZ_COMP:
      if(j<NY-3) {
        if(sliceXEnd>NX-1) sliceXEnd = NX-1;
        for(int k=sliceZStart,sI=0;k<sliceZEnd;k++) {
          for(int i=sliceXStart;i<sliceXEnd;i++) {
            float aa=((VY(i+1,j-1,k)-VY(i,j-1,k))/dx-(VX(i,j  ,k)-VX(i,j-1,k))/dy);
            float bb=((VY(i+1,j  ,k)-VY(i,j  ,k))/dx-(VX(i,j+1,k)-VX(i,j  ,k))/dy);
            float cc=((VY(i+1,j+1,k)-VY(i,j+1,k))/dx-(VX(i,j+2,k)-VX(i,j+1,k))/dy);
            float dd=((VY(i+1,j+2,k)-VY(i,j+2,k))/dx-(VX(i,j+3,k)-VX(i,j+2,k))/dy);
            _slice[sI++]=scalarVel*(b1*aa+b2*bb+b3*cc+b4*dd);
          }
          if(sliceXEnd==NX-1) sI++;
        }
      }
      break;
      
    case SLICE_VP_COMP:
      for(int k=sliceZStart,sI=0;k<sliceZEnd;k++){
        for(int i=sliceXStart;i<sliceXEnd;i++){
          float val=
          b1*model()->vp(i,j-1,k)+b2*model()->vp(i,j  ,k)+
          b3*model()->vp(i,j+1,k)+b4*model()->vp(i,j+2,k);
          _slice[sI++]=val;
        }
      }
      break;
      
    case SLICE_VS_COMP:
      for(int k=sliceZStart,sI=0;k<sliceZEnd;k++){
        for(int i=sliceXStart;i<sliceXEnd;i++){
          float val=
          b1*model()->vs(i,j-1,k)+b2*model()->vs(i,j  ,k)+
          b3*model()->vs(i,j+1,k)+b4*model()->vs(i,j+2,k);
          _slice[sI++]=val;
        }
      }
      break;
      
    case SLICE_RHO_COMP:
      for(int k=sliceZStart,sI=0;k<sliceZEnd;k++){
        for(int i=sliceXStart;i<sliceXEnd;i++){
          float val=
          b1*model()->rho(i,j-1,k)+b2*model()->rho(i,j  ,k)+
          b3*model()->rho(i,j+1,k)+b4*model()->rho(i,j+2,k);
          _slice[sI++]=val;
        }
      }
      break;
      
    default:
      assert(FALSE,
             "sliceXZCompValue--do not know how to make slice of type %i",comp);
  }
}

void slaveSgfdDependent::sliceYZCompValue(int comp,int i,
                      int sliceYStart,int sliceYEnd,
                      int sliceZStart,int sliceZEnd,
                      float a1,float a2,float a3,float a4,
                      float *vx,float *vy,float *vz,
                      float *xx,float *yy,float *zz,
                      float *xy,float *xz,float *yz){
  DEF_MODEL_SIZE(modelDef());
  DEF_MODEL_SCALE(modelDef());
  DEF_MODEL_LIMITS(modelDef());
  
  switch(comp){
    case SLICE_VX_COMP:
      //x-component of velocity.
      doYZSlice(vx,sliceYStart,sliceYEnd,sliceZStart,sliceZEnd,
                _slice,scalarVel,i,a1,a2,a3,a4);
      break;
      
    case SLICE_VY_COMP:
      //y-component of velocity.
      doYZSlice(vy,sliceYStart,sliceYEnd,sliceZStart,sliceZEnd,
                _slice,scalarVel,i,a1,a2,a3,a4);
      break;
      
    case SLICE_VZ_COMP:
      //z-component of velocity.
      doYZSlice(vz,sliceYStart,sliceYEnd,sliceZStart,sliceZEnd,
                _slice,scalarVel,i,a1,a2,a3,a4);
      break;
      
    case SLICE_PRESS_COMP:
      //Acoustic pressure.
      for(int k=sliceZStart,sI=0;k<sliceZEnd;k++){
        for(int j=sliceYStart;j<sliceYEnd;j++){
          if(yy){
            //Elastic problem.
            float aa=XX(i-1,j,k)+YY(i-1,j,k)+ZZ(i-1,j,k);
            float dd=XX(i+2,j,k)+YY(i+2,j,k)+ZZ(i+2,j,k);
            float bb=XX(i  ,j,k)+YY(i  ,j,k)+ZZ(i  ,j,k);
            float cc=XX(i+1,j,k)+YY(i+1,j,k)+ZZ(i+1,j,k);
            
            float press=a1*aa+a2*bb+a3*cc+a4*dd;
            _slice[sI++]=-press*scalarStress/3.0;
          }else{
            //Acoustic problem xx is already the pressure.
            _slice[sI++]= scalarStress*
            (a1*XX(i-1,j,k)+a2*XX(i  ,j,k)+
             a3*XX(i+1,j,k)+a4*XX(i+2,j,k));
          }
        }
      }
      break;
      
    case SLICE_RX_COMP:
      if(sliceZEnd>NZ-1) sliceZEnd = NZ-1;
      if(sliceYEnd>NY-1) sliceYEnd = NY-1;
      for(int k=sliceZStart,sI=0;k<sliceZEnd;k++) {
        for(int j=sliceYStart;j<sliceYEnd;j++) {
          float aa=((VZ(i-1,j+1,k)-VZ(i-1,j,k))/dy-(VY(i-1,j,k+1)-VY(i-1,j,k))/dz);
          float bb=((VZ(i  ,j+1,k)-VZ(i  ,j,k))/dy-(VY(i  ,j,k+1)-VY(i  ,j,k))/dz);
          float cc=((VZ(i+1,j+1,k)-VZ(i+1,j,k))/dy-(VY(i+1,j,k+1)-VY(i+1,j,k))/dz);
          float dd=((VZ(i+2,j+1,k)-VZ(i+2,j,k))/dy-(VY(i+2,j,k+1)-VY(i+2,j,k))/dz);
          _slice[sI++]=scalarVel*(a1*aa+a2*bb+a3*cc+a4*dd);
        }
        if(sliceYEnd==NY-1) sI++;
      }
      break;
      
    case SLICE_RY_COMP:
      if(i<NX-3) {
        if(sliceZEnd>NZ-1) sliceZEnd = NZ-1;
        for(int k=sliceZStart,sI=0;k<sliceZEnd;k++) {
          for(int j=sliceYStart;j<sliceYEnd;j++) {
            float aa = ((VX(i-1,j,k+1)-VX(i-1,j,k))/dz-(VZ(i  ,j,k)-VZ(i-1,j,k))/dx);
            float bb = ((VX(i  ,j,k+1)-VX(i  ,j,k))/dz-(VZ(i+1,j,k)-VZ(i  ,j,k))/dx);
            float cc = ((VX(i+1,j,k+1)-VX(i+1,j,k))/dz-(VZ(i+2,j,k)-VZ(i+1,j,k))/dx);
            float dd = ((VX(i+2,j,k+1)-VX(i+2,j,k))/dz-(VZ(i+3,j,k)-VZ(i+2,j,k))/dx);
            _slice[sI++]=scalarVel*(a1*aa+a2*bb+a3*cc+a4*dd);
          }
        }
      }
      break;
      
    case SLICE_RZ_COMP:
      if(i<NX-3) {
        if(sliceYEnd>NY-1) sliceYEnd = NY-1;
        for(int k=sliceZStart,sI=0;k<sliceZEnd;k++) {
          for(int j=sliceYStart;j<sliceYEnd;j++) {
            float aa=((VY(i  ,j,k)-VY(i-1,j,k))/dx-(VX(i-1,j+1,k)-VX(i-1,j,k))/dy);
            float bb=((VY(i+1,j,k)-VY(i  ,j,k))/dx-(VX(i  ,j+1,k)-VX(i  ,j,k))/dy);
            float cc=((VY(i+2,j,k)-VY(i+1,j,k))/dx-(VX(i+1,j+1,k)-VX(i+1,j,k))/dy);
            float dd=((VY(i+3,j,k)-VY(i+2,j,k))/dx-(VX(i+2,j+1,k)-VX(i+2,j,k))/dy);
            _slice[sI++]=scalarVel*(a1*aa+a2*bb+a3*cc+a4*dd);
          }
          if(sliceYEnd==NY-1) sI++;
        }
      }
      break;
      
    case SLICE_VP_COMP:
      for(int k=sliceZStart,sI=0;k<sliceZEnd;k++){
        for(int j=sliceYStart;j<sliceYEnd;j++){
          float val=
          a1*model()->vp(i-1,j,k)+a2*model()->vp(i  ,j,k)+
          a3*model()->vp(i+1,j,k)+a4*model()->vp(i+2,j,k);
          _slice[sI++]=val;
        }
      }
      break;
      
    case SLICE_VS_COMP:
      for(int k=sliceZStart,sI=0;k<sliceZEnd;k++){
        for(int j=sliceYStart;j<sliceYEnd;j++){
          float val=
          a1*model()->vs(i-1,j,k)+a2*model()->vs(i  ,j,k)+
          a3*model()->vs(i+1,j,k)+a4*model()->vs(i+2,j,k);
          _slice[sI++]=val;
        }
      }
      break;
      
    case SLICE_RHO_COMP:
      for(int k=sliceZStart,sI=0;k<sliceZEnd;k++){
        for(int j=sliceYStart;j<sliceYEnd;j++){
          float val=
          a1*model()->rho(i-1,j,k)+a2*model()->rho(i  ,j,k)+
          a3*model()->rho(i+1,j,k)+a4*model()->rho(i+2,j,k);
          _slice[sI++]=val;
        }
      }
      break;
      
    default:
      assert(FALSE,
             "sliceYZCompValue--do not know how to make slice of type %i",comp);
  }
}

void slaveSgfdDependent::sliceXYCompValue(int comp,int k,
                      int sliceXStart,int sliceXEnd,
                      int sliceYStart,int sliceYEnd,
                      float c1,float c2,float c3,float c4,
                      float *vx,float *vy,float *vz,
                      float *xx,float *yy,float *zz,
                      float *xy,float *xz,float *yz){
  DEF_MODEL_SIZE(modelDef());
  DEF_MODEL_SCALE(modelDef());
  DEF_MODEL_LIMITS(modelDef());
  
  switch(comp){
    case SLICE_VX_COMP:
      //x-component of velocity.
      doXYSlice(vx,sliceXStart,sliceXEnd,sliceYStart,sliceYEnd,
                _slice,scalarVel,k,c1,c2,c3,c4);
      break;
      
    case SLICE_VY_COMP:
      //y-component of velocity.
      doXYSlice(vy,sliceXStart,sliceXEnd,sliceYStart,sliceYEnd,
                _slice,scalarVel,k,c1,c2,c3,c4);
      break;
      
    case SLICE_VZ_COMP:
      //z-component of velocity.
      doXYSlice(vz,sliceXStart,sliceXEnd,sliceYStart,sliceYEnd,
                _slice,scalarVel,k,c1,c2,c3,c4);
      break;
      
    case SLICE_PRESS_COMP:
      //Acoustic pressure.
      for(int j=sliceYStart,sI=0;j<sliceYEnd;j++){
        for(int i=sliceXStart;i<sliceXEnd;i++){
          if(yy){
            //Elastic problem.
            float aa=XX(i,j,k-1)+YY(i,j,k-1)+ZZ(i,j,k-1);
            float bb=XX(i,j,k)  +YY(i,j,k)  +ZZ(i,j,k);
            float cc=XX(i,j,k+1)+YY(i,j,k+1)+ZZ(i,j,k+1);
            float dd=XX(i,j,k+2)+YY(i,j,k+2)+ZZ(i,j,k+2);
            
            float press=c1*aa+c2*bb+c3*cc+c4*dd;
            _slice[sI++]=-press*scalarStress/3.0;
          }else{
            //Acoustic problem xx is already the pressure.
            _slice[sI++]= scalarStress*
            (c1*XX(i,j,k-1)+c2*XX(i,j,k  )+
             c3*XX(i,j,k+1)+c4*XX(i,j,k+2));
          }
        }
      }
      break;
      
    case SLICE_RX_COMP:
      if(k<NZ-3) {
        if(sliceYEnd>NY-1) sliceYEnd = NY-1;
        for(int j=sliceYStart,sI=0;j<sliceYEnd;j++) {
          for(int i=sliceXStart;i<sliceXEnd;i++) {
            float aa=((VZ(i,j+1,k-1)-VZ(i,j,k-1))/dy-(VY(i,j,k  )-VY(i,j,k-1))/dz);
            float bb=((VZ(i,j+1,k  )-VZ(i,j,k  ))/dy-(VY(i,j,k+1)-VY(i,j,k  ))/dz);
            float cc=((VZ(i,j+1,k+1)-VZ(i,j,k+1))/dy-(VY(i,j,k+2)-VY(i,j,k+1))/dz);
            float dd=((VZ(i,j+1,k+2)-VZ(i,j,k+2))/dy-(VY(i,j,k+3)-VY(i,j,k+2))/dz);
            _slice[sI++]=scalarVel*(c1*aa+c2*bb+c3*cc+c4*dd);
          }
        }
      }
      break;
      
    case SLICE_RY_COMP:
      if(k<NZ-3) {
        if(sliceXEnd>NX-1) sliceXEnd = NX-1;
        for(int j=sliceYStart,sI=0;j<sliceYEnd;j++) {
          for(int i=sliceXStart;i<sliceXEnd;i++) {
            float aa = ((VX(i,j,k  )-VX(i,j,k-1))/dz-(VZ(i+1,j,k-1)-VZ(i,j,k-1))/dx);
            float bb = ((VX(i,j,k+1)-VX(i,j,k  ))/dz-(VZ(i+1,j,k  )-VZ(i,j,k  ))/dx);
            float cc = ((VX(i,j,k+2)-VX(i,j,k+1))/dz-(VZ(i+1,j,k+1)-VZ(i,j,k+1))/dx);
            float dd = ((VX(i,j,k+3)-VX(i,j,k+2))/dz-(VZ(i+1,j,k+2)-VZ(i,j,k+2))/dx);
            _slice[sI++]=scalarVel*(c1*aa+c2*bb+c3*cc+c4*dd);
          }
          if(sliceXEnd==NX-1) sI++;
        }
      }
      break;
      
    case SLICE_RZ_COMP:
      if(sliceXEnd>NX-1) sliceXEnd = NX-1;
      if(sliceYEnd>NY-1) sliceYEnd = NY-1;
      for(int j=sliceYStart,sI=0;j<sliceYEnd;j++) {
        for(int i=sliceXStart;i<sliceXEnd;i++) {
          float aa=((VY(i+1,j,k-1)-VY(i,j,k-1))/dx-(VX(i,j+1,k-1)-VX(i,j,k-1))/dy);
          float bb=((VY(i+1,j,k  )-VY(i,j,k  ))/dx-(VX(i,j+1,k  )-VX(i,j,k  ))/dy);
          float cc=((VY(i+1,j,k+1)-VY(i,j,k+1))/dx-(VX(i,j+1,k+1)-VX(i,j,k+1))/dy);
          float dd=((VY(i+1,j,k+2)-VY(i,j,k+2))/dx-(VX(i,j+1,k+2)-VX(i,j,k+2))/dy);
          _slice[sI++]=scalarVel*(c1*aa+c2*bb+c3*cc+c4*dd);
        }
        if(sliceXEnd==NX-1) sI++;
      }
      break;
      
    case SLICE_VP_COMP:
      for(int j=sliceYStart,sI=0;j<sliceYEnd;j++){
        for(int i=sliceXStart;i<sliceXEnd;i++){
          float val=
          c1*model()->vp(i,j,k-1)+c2*model()->vp(i,j,k  )+
          c3*model()->vp(i,j,k+1)+c4*model()->vp(i,j,k+2);
          _slice[sI++]=val;
        }
      }
      break;
      
    case SLICE_VS_COMP:
      for(int j=sliceYStart,sI=0;j<sliceYEnd;j++){
        for(int i=sliceXStart;i<sliceXEnd;i++){
          float val=
          c1*model()->vs(i,j,k-1)+c2*model()->vs(i,j,k  )+
          c3*model()->vs(i,j,k+1)+c4*model()->vs(i,j,k+2);
          _slice[sI++]=val;
        }
      }
      break;
      
    case SLICE_RHO_COMP:
      for(int j=sliceYStart,sI=0;j<sliceYEnd;j++){
        for(int i=sliceXStart;i<sliceXEnd;i++){
          float val=
          c1*model()->rho(i,j,k-1)+c2*model()->rho(i,j,k  )+
          c3*model()->rho(i,j,k+1)+c4*model()->rho(i,j,k+2);
          _slice[sI++]=val;
        }
      }
      break;
      
    default:
      assert(FALSE,
             "sliceXYCompValue--do not know how to make slice of type %i",comp);
  }
}

extraOutputGroup* sgfdCommandParser::allocateExtraOutput(){
  return new extraOutputGroup(modelDef());
}

int sgfdCommandParser::checkProcessFlag(int argc,char* argv[],int& i){
  if(sgfdCommandParserBase::checkProcessFlag(argc,argv,i))
    return TRUE;
  
  switch (argv[i][1]){
    case 'E':
      if(!_extraOutput)
        _extraOutput=allocateExtraOutput();
      _extraOutput->addExtraOutput(i,2,argc,argv,modelDef(),TRUE);
      break;
    default:
      return FALSE;
  }
  return TRUE;
}

void sgfdCommandParser::addExtraOutputFlags(){
  addFlagDescription(FLAG_OPTIONAL,"-E[o]","s",
                     "Set slice output to file.");
  addFlagDescription(FLAG_OPTIONAL,"-Es","f s s f",
                     "Add 1 slice, time, component (Vx, Vy, Vz), plane (XY,XZ,YZ) at coordinate.");
  addFlagDescription(FLAG_OPTIONAL,"-En","i s s f",
                     "Add n slices, component (Vx, Vy, Vz), plane (XY,XZ,YZ) at coordinate.");
  addFlagDescription(FLAG_EXPERT,"-EA","s s f",
                     "Add slice at every time step, component (Vx, Vy, Vz), plane (XY,XZ,YZ) at",
                     1,"coordinate.");
  addFlagDescription(FLAG_OPTIONAL,"","","");
  
  addFlagDescription(FLAG_OPTIONAL,"-Ew","f",
                     "Add full waveform output at time.");
  addFlagDescription(FLAG_OPTIONAL,"-Ec","i",
                     "Add full waveform output checkpoints every n iterations");
  addFlagDescription(FLAG_OPTIONAL,"","","");
}
