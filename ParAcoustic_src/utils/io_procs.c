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
 *  io_procs.c
 *
 *
 *  Defines some functions used in reading and writing to NetCDF files.
 *
 *  Defines the following functions:
 *  openCDFFile
 *  writeCDFHeader
 *  writeCDFVarFloat
 *  readCDFVarBlock
 *
 */
#include "io_procs.h"
#include "nstdutil.h"
#include "netcdf.h"

/*open a file, this first checks if the file name ends in cdf
// then opens the file and makes sure it is really open*/
int openCDFFile(const char* fileName,int create,int mode){
  int cdfFile;

  /*Make sure filename ends in .cdf*/
  char buffer[1024];
  if(!strcmp(fileName+strlen(fileName)-4,".cdf")){
    sprintf(buffer,"%s",fileName);
  }else{
    sprintf(buffer,"%s.cdf",fileName);
  }

  if(create){
    assert(nc_create(buffer,NC_CLOBBER|NC_64BIT_OFFSET,&cdfFile)==NC_NOERR,
	   "openCDFFile--unable to create file %s\n",buffer);
  }else{
    if(nc_open(buffer,mode,&cdfFile)!=NC_NOERR)
      return -1;
  }
  return cdfFile;
}

/*Write a generic header. This includes all four dimensions, limits and increments
// Variables for x,y,z and t and the scalar multipliers
// if numVars is non-zero then x-y-z variables are defined, names are expected
// in varNames to correspond with numVars*/
int writeCDFHeader(const char* fileName,
		   int nx,int ny,int nz,int nt,
		   float* origin,
		   float dx,float minX,
		   float dy,float minY,
		   float dz,float minZ,
		   float dt,float minT,
		   int numVars,char** varNames){
  int numCoord,numSpatialCoord;
  int nxDim,nyDim,nzDim,ntDim;
  int minima,increments,originVar;
  int xVar,yVar,zVar,tVar;
  size_t index=0;
  int maxElems=MAX(MAX(nx,ny),MAX(nz,nt));
  float* nodeCenters;
  int i;
  int outFile=openCDFFile(fileName,TRUE,NC_CLOBBER);
  char buffer[512];

  sprintf(buffer,"parallel_elasti generic NETcdf file");
  nc_put_att_text(outFile,NC_GLOBAL,"title",strlen(buffer)+1,buffer);
  nc_put_att_text(outFile,NC_GLOBAL,"history",strlen(CommandLine)+1,CommandLine);

  nc_def_dim(outFile,"numCoord",4L,&numCoord);
  nc_def_dim(outFile,"numSpatialCoord",3L,&numSpatialCoord);
  nc_def_dim(outFile,"NX",(long)nx,&nxDim);
  nc_def_dim(outFile,"NY",(long)ny,&nyDim);
  nc_def_dim(outFile,"NZ",(long)nz,&nzDim);
  if(nt)
    nc_def_dim(outFile,"NT",(long)nt,&ntDim);

  /*Define variables
  // model limits*/
  nc_def_var(outFile,"minima",NC_FLOAT,1,&numCoord,&minima);
  nc_def_var(outFile,"increments",NC_FLOAT,1,&numCoord,&increments);
  if(origin)
    nc_def_var(outFile,"origin",NC_FLOAT,1,&numCoord,&originVar);

  nc_def_var(outFile,"x",NC_FLOAT,1,&nxDim,&xVar);
  nc_def_var(outFile,"y",NC_FLOAT,1,&nyDim,&yVar);
  nc_def_var(outFile,"z",NC_FLOAT,1,&nzDim,&zVar);
  if(nt)
    nc_def_var(outFile,"time",NC_FLOAT,1,&ntDim,&tVar);

  /*
  //Fill the data
  */
  nc_enddef(outFile);

  if(origin)
    nc_put_var_float(outFile,originVar,origin);

  index=0;
  nc_put_var1_float(outFile,minima,&index,&minX);
  nc_put_var1_float(outFile,increments,&index,&dx);

  index=1;
  nc_put_var1_float(outFile,minima,&index,&minY);
  nc_put_var1_float(outFile,increments,&index,&dy);

  index=2;
  nc_put_var1_float(outFile,minima,&index,&minZ);
  nc_put_var1_float(outFile,increments,&index,&dz);
  
  if(nt){
    index=3;
    nc_put_var1_float(outFile,minima,&index,&minT);
    nc_put_var1_float(outFile,increments,&index,&dt);
  }

  nodeCenters=(float*)malloc(maxElems*sizeof(float));
  assert(nodeCenters!=NULL,
	 "writeCDFHeader--unable to allocate %i floats for nodeCenters",
	 maxElems);

  for(i=0;i<nx;i++)
    nodeCenters[i]=minX+i*dx;
  nc_put_var_float(outFile,xVar,nodeCenters);

  for(i=0;i<ny;i++)
    nodeCenters[i]=minY+i*dy;
  nc_put_var_float(outFile,yVar,nodeCenters);

  for(i=0;i<nz;i++)
    nodeCenters[i]=minZ+i*dz;
  nc_put_var_float(outFile,zVar,nodeCenters);

  if(nt){
    for(i=0;i<nt;i++)
      nodeCenters[i]=minT+i*dt;
    nc_put_var_float(outFile,tVar,nodeCenters);
  }
  free(nodeCenters);

  /*Now define any variables with dimensions NX*NY*NZ
  // many calls will do this themselves so it is optional*/
  if(numVars){
    int newVar;
    int spatialDims[3]; /*last dimension varies fastest*/
    spatialDims[0]=nzDim;spatialDims[1]=nyDim;spatialDims[2]=nxDim;
    nc_redef(outFile); 
    for(i=0;i<numVars;i++){
      char* varName=varNames[i];
      assert(varName !=NULL && varName[0],
	     "writeCDFHeader--bad name for variable %i",
	     i);

      nc_def_var(outFile,varName,NC_FLOAT,3,spatialDims,&newVar);
    }
  }

  assert(nc_close(outFile)==NC_NOERR,
	 "writeCDFHeader--nc_close of %s unsuccessful",fileName);
  return 1;
}

/*Read and write basic elements, each of these will open and close
 the file, so it is not the best way to make alot of reads or writes
 but it is simple if only 1 or 2 values are needed*/
float writeCDFVarFloat(const char* fileName,
		       const char* varName,float data){
  size_t index=0;
  int outFile=openCDFFile(fileName,FALSE,NC_WRITE);
  int varID;
  if(nc_inq_varid(outFile,varName,&varID)!=NC_NOERR){
    nc_redef(outFile);
    assert(nc_def_var(outFile,varName,NC_FLOAT,
		      0,NULL,&varID)==NC_NOERR,
	   "writeCDFVarFloat--unable to create float variable %s in %s",
	   varName,fileName);
  }

  nc_enddef(outFile);
  nc_put_var1_float(outFile,varID,&index,&data);

  assert(nc_close(outFile)!=-1,
	 "writeCDFVarFloat--nc_close of %s unsuccessful",fileName);
  return data;
}

int readCDFVarBlock(const char* fileName,
		    const char* varName,float** varData,
		    int* procLim){
  size_t start[3];
  size_t count[3];

  int inFile=openCDFFile(fileName,FALSE,NC_NOWRITE);
  int varId;
  if(nc_inq_varid(inFile,varName,&varId)!=NC_NOERR)
    return FALSE;

  if(procLim){
    start[0]=procLim[4];
    start[1]=procLim[2];
    start[2]=procLim[0];

    count[0]=procLim[5]-procLim[4];
    count[1]=procLim[3]-procLim[2];
    count[2]=procLim[1]-procLim[0];

    if(!*varData)
      assert((*varData=(float*)malloc(count[0]*count[1]*count[2]*sizeof(float)))!=NULL,
	     "readCDFVarBlock--unable to allocate %i floats for varData",
	     count[0]*count[1]*count[2]);

    assert(nc_get_vara_float(inFile,varId,start,count,*varData)==NC_NOERR,
	   "readCDFVarBlock--unable to read (%i,%i,%i) block from %s at (%i,%i,%i)",
	   count[0],count[1],count[2],
	   fileName,
	   start[0],start[1],start[2]);
  }else{
    int dimID;
    assert(nc_inq_dimid(inFile,"NX",&dimID)==NC_NOERR,
	   "readCDFVarBlock--unable to get dimension NX from %s",fileName);
    assert(nc_inq_dimlen(inFile,dimID,&count[0])==NC_NOERR,
	   "readCDFVarBlock--unable to read dimension NX from %s",fileName);
    assert(nc_inq_dimid(inFile,"NY",&dimID)==NC_NOERR,
	   "readCDFVarBlock--unable to get dimension NY from %s",fileName);
    assert(nc_inq_dimlen(inFile,dimID,&count[1])==NC_NOERR,
	   "readCDFVarBlock--unable to read dimension NY from %s",fileName);
    assert(nc_inq_dimid(inFile,"NZ",&dimID)==NC_NOERR,
	   "readCDFVarBlock--unable to get dimension NZ from %s",fileName);
    assert(nc_inq_dimlen(inFile,dimID,&count[2])==NC_NOERR,
	   "readCDFVarBlock--unable to read dimension NZ from %s",fileName);
    
    if(!*varData)
      assert((*varData=(float*)malloc(count[0]*count[1]*count[2]*sizeof(float)))!=NULL,
	     "readCDFVarBlock--unable to allocate %i floats for varData",
	     count[0]*count[1]*count[2]);

    assert(nc_get_var_float(inFile,varId,*varData)==NC_NOERR,
	   "readCDFVarBlock--unable to read var %s from %s at (%i,%i,%i)",
	   varName,fileName);
  }

  assert(nc_close(inFile)!=-1,
	 "readCDFVarBlock--nc_close of %s unsuccessful",fileName);
  return count[0]*count[1]*count[2];
}
