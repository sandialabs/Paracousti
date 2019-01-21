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
 *  xtrautil.cc
 *
 *
 *  Defines functions used to read NetCDF files.
 *
 *  Defines the following functions:
 *  readModelVariable3D
 *  readModelVariable2DXZ
 *  readModelVariableSubdomain
 *  readModelVariableIndexSet
 *  readModelVariable
 *
 */
#include "xtrautil.hh"

#include <stdlib.h>


//
//coordinate transformations
//

/*
    Utility function returns the arc distance in degrees (delta),
    the azimuth from event to station (az; degress clockwise from north),
    and the backazimuth from station to event (baz; degrees clockwise from
    north) between two points: the "evt" and "sta". The event and station
    positions are specified by lat/lon pairs which are given as the
    decimal (double) values evtlat, evtlon, and stlat, stlon.
    The azimuth returned is from the event to the station (position 1
    to position 2).
*/
//
//Here are some new subroutines to read a variable from a model. The variable
// name is checked for upper and lower case versions. A possible 1D version of the
// variable is checked. The file is checked for versions that have been written
// with the variable name appended. And then we look for files that have been
// broken into subdomains.
static int readModelVariable3D(float* var,
			       int modelFile,char* modelName,const char* varName,
			       size_t start[3], size_t size[3],
			       int iOpenFile){
  int varId;
  char varNameL[64],varNameU[64];
  strcpy(varNameL,varName);strcpy(varNameU,varName);
  varNameL[0]=tolower(varName[0]);
  varNameU[0]=toupper(varName[0]);

  if(nc_inq_varid(modelFile,varNameL,&varId)==NC_NOERR ||
     nc_inq_varid(modelFile,varNameU,&varId)==NC_NOERR){
    //Success, read the variable, close the file and exit.
    assert(nc_get_vara_float(modelFile,varId,start,size,var)==NC_NOERR,
	   "readModelVariable3D--unable to read %s[%i,%i,%i=>%i,%i,%i] from %s (%i,%i)",
	   varName,start[2],start[1],start[0],size[2],size[1],size[0],
	   modelName,modelFile,varId);

    //Close the file and exit.
    if(iOpenFile)
      assert(nc_close(modelFile)==NC_NOERR,
	     "readModelVariable3D--unable to close %s (%i)",modelName,modelFile);
    return TRUE;
  }
  return FALSE;
}
static int readModelVariable3D(int* var,
			       int modelFile,char* modelName,const char* varName,
			       size_t start[3], size_t size[3],
			       int iOpenFile){
  int varId;
  char varNameL[64],varNameU[64];
  strcpy(varNameL,varName);strcpy(varNameU,varName);
  varNameL[0]=tolower(varName[0]);
  varNameU[0]=toupper(varName[0]);

  if(nc_inq_varid(modelFile,varNameL,&varId)==NC_NOERR ||
     nc_inq_varid(modelFile,varNameU,&varId)==NC_NOERR){
    //Success, read the variable, close the file and exit.
    assert(nc_get_vara_int(modelFile,varId,start,size,var)==NC_NOERR,
	   "readModelVariable3D--unable to read %s[%i,%i,%i=>%i,%i,%i] from %s (%i,%i)",
	   varName,start[2],start[1],start[0],size[2],size[1],size[0],
	   modelName,modelFile,varId);

    //Close the file and exit.
    if(iOpenFile)
      assert(nc_close(modelFile)==NC_NOERR,
	     "readModelVariable3D--unable to close %s (%i)",modelName,modelFile);
    return TRUE;
  }
  if(iOpenFile)
    assert(nc_close(modelFile)==NC_NOERR,
     "readModelVariable3D--unable to close %s (%i)",modelName,modelFile);
  return FALSE;
}
static int readModelVariable2DXZ(float* var,
				 int modelFile,char* modelName,const char* varName,
				 size_t start[3], size_t size[3],
				 int iOpenFile){
  int varId;
  char varNameL[64],varNameU[64];
  sprintf(varNameL,"%sXZ",varName);sprintf(varNameU,"%sXZ",varName);
  varNameL[0]=tolower(varName[0]);
  varNameU[0]=toupper(varName[0]);

  if(nc_inq_varid(modelFile,varNameL,&varId)==NC_NOERR ||
     nc_inq_varid(modelFile,varNameU,&varId)==NC_NOERR){
    int nx=size[2],ny=size[1],nz=size[0];

    //Success, read the variable.
    size_t start2[2],size2[2];
    start2[0]=start[0];start2[1]=start[2];
    size2[0]=size[0];size2[1]=size[2];
    float *xzPlane;
    assert((xzPlane=(float*)malloc(nx*nz*sizeof(float)))!=NULL,
	   "readModelVariable2DXZ--unable to allocate 2D buffer %ix%i",
	   nx,nz);
    assert(nc_get_vara_float(modelFile,varId,start2,size2,xzPlane)==NC_NOERR,
	   "readModelVariable2DXZ--unable to read %s[%i,%i=>%i,%i] from %s (%i,%i)",
	   varName,start2[1],start2[0],size2[1],size2[0],
	   modelName,modelFile,varId);

    //Propagate along the y direction.
    for(int k=0;k<nz;k++){
      for(int j=0;j<ny;j++){
	for(int i=0;i<nx;i++){
	  int index2=i+k*nx,index3=i+j*nx+k*nx*ny;
	  var[index3]=xzPlane[index2];
	}
      }
    }

    //Free the local space.
    free(xzPlane);

    //Close the file and exit.
    if(iOpenFile)
      assert(nc_close(modelFile)==NC_NOERR,
	     "readModelVariable3D--unable to close %s (%i)",modelName,modelFile);
    return TRUE;
  }
  return FALSE;
}

static int readModelVariableSubdomain(int* var,
				      char* modelName,const char* varName,
				      int *lim,
				      size_t start[3],size_t size[3],
				      int nx,int ny,int nz,
				      int modelFile,int iOpenFile){
  int subdomainFilenameVar;
  if(nc_inq_varid(modelFile,"subdomainFilename",&subdomainFilenameVar)!=NC_NOERR)
    return FALSE;

  //Read some important dimensions for reading these files.
  int dimID;
  size_t nSubdomains,subdomainFileNameLength;
  assert(nc_inq_dimid(modelFile,"nSubdomains",&dimID)==NC_NOERR,
	 "readModelVariableSubdomain--unable to read dim nSubdomains from %s(%i)",
	 modelName,modelFile);
  assert(nc_inq_dimlen(modelFile,dimID,&nSubdomains)==NC_NOERR,
	 "readModelVariableSubdomain--unable to get dim nSubdomains from %s(%i)",
	 modelName,modelFile);

  assert(nc_inq_dimid(modelFile,"subdomainFileNameLength",&dimID)==NC_NOERR,
	 "readModelVariableSubdomain--unable to read dim subdomainFileNameLength from %s(%i)",
	 modelName,modelFile);
  assert(nc_inq_dimlen(modelFile,dimID,&subdomainFileNameLength)==NC_NOERR,
	 "readModelVariableSubdomain--unable to get dim subdomainFileNameLength from %s(%i)",
	 modelName,modelFile);

  //Get and additional variable that defines the proclimits within each
  // subdomain file.
  int subdomainLimitVar;
  assert(nc_inq_varid(modelFile,"subdomainLimit",&subdomainLimitVar)==NC_NOERR,
	 "readModelVariableSubdomain--unable to get var subdomainLimit from %s(%i)",
	 modelName,modelFile);

  //Define start and size variables for reading the subdomain components.
  size_t subdomainStart[2]={0,0};
  size_t subdomainLimitSize[2]={1,6};
  size_t subdomainNameSize[2]={1,subdomainFileNameLength};
    
  //Fill the buffer with the path to the current file, prepend this path to the
  // subdomain filenames that are to be read from the file.
  char buffer[1024];
  strcpy(buffer,modelName);
  char* bufferRootPos=strrchr(buffer,'/');
	if(bufferRootPos==NULL)
		bufferRootPos=buffer;
	else
		bufferRootPos+=1;

  //Loop through the subdomains and read the required parts to fill in the 
  // variable. Note that to check if every element of the variable is filled in
  // we are pre-filling with -1.
  int nxy=nx*ny,nxyz=nx*ny*nz;
  for(int i=0;i<nxyz;var[i++]=-1);
  for(int i=0;i<nSubdomains;i++){
    subdomainStart[0]=i;

    //Read the current subdomain limits and check for overlap.
    int subdomainLimit[6];
    assert(nc_get_vara_int(modelFile,subdomainLimitVar,
			   subdomainStart,subdomainLimitSize,
			   subdomainLimit)==NC_NOERR,
	   "readModelVariableSubdomain--unable to read subdomainLimit[%i,%i,%i,%i] from %s(%i,%i)",
	   subdomainStart[0],subdomainStart[1],
	   subdomainLimitSize[0],subdomainLimitSize[1],
	   modelName,modelFile,subdomainLimitVar);
      
    if(ISOVERLAP(lim[0],lim[1],subdomainLimit[0],subdomainLimit[1]) &&
       ISOVERLAP(lim[2],lim[3],subdomainLimit[2],subdomainLimit[3]) &&
       ISOVERLAP(lim[4],lim[5],subdomainLimit[4],subdomainLimit[5])){
      //Read the portion that has overlap only.
      int readLimit[6]={MAX(lim[0],subdomainLimit[0]),
			MIN(lim[1],subdomainLimit[1]),
			MAX(lim[2],subdomainLimit[2]),
			MIN(lim[3],subdomainLimit[3]),
			MAX(lim[4],subdomainLimit[4]),
			MIN(lim[5],subdomainLimit[5])};
      size_t subdomainReadStart[3]={readLimit[4]-subdomainLimit[4],
				    readLimit[2]-subdomainLimit[2],
				    readLimit[0]-subdomainLimit[0]};
      size_t subdomainReadSize[3]={readLimit[5]-readLimit[4],
				   readLimit[3]-readLimit[2],
				   readLimit[1]-readLimit[0]};

      int subdomainNX=subdomainReadSize[2],subdomainNY=subdomainReadSize[1],subdomainNZ=subdomainReadSize[0];
      int subdomainNXY=subdomainNX*subdomainNY;
      int subdomainNXYZ=subdomainNX*subdomainNY*subdomainNZ;

      //Get the filename and open for reading. 
      assert(nc_get_vara_text(modelFile,subdomainFilenameVar,
			      subdomainStart,subdomainNameSize,
			      bufferRootPos)==NC_NOERR,
	     "readModelVariableSubdomain--unable to read subdomainName[%i,%i,%i,%i] from %s(%i,%i)",
	     subdomainStart[0],subdomainStart[1],
	     subdomainNameSize[0],subdomainNameSize[1],
	     modelName,modelFile,subdomainFilenameVar);
      bufferRootPos[subdomainNameSize[1]]='\0';
      int subdomainFile;
      assert((subdomainFile=openCDFFile(buffer,FALSE,NC_NOWRITE))>=0,
	     "readModelVariableSubdomain--unable to open %s for reading",buffer);

      //Allocate a variable to read the required portion of the subdomain into.
      // Should be able to do this in place with one of the "funky" nc_get
      // subroutines but that will have to wait for a later time.
      int* subdomainVar;
      assert((subdomainVar=(int*)malloc(subdomainNXYZ*sizeof(int)))!=NULL,
	     "readModelVariableSubdomain--unable to allocate %i (%i,%i,%i) ints for local space",
	     subdomainReadSize[2],subdomainReadSize[1],subdomainReadSize[0],
	     subdomainNXYZ);

      //Read the data into the local space and move to the correct part of the
      // var.
      if(!readModelVariable3D(subdomainVar,subdomainFile,buffer,varName,
			  subdomainReadStart,subdomainReadSize,TRUE)) {
          free(subdomainVar);
          if(iOpenFile)
            assert(nc_close(modelFile)==NC_NOERR,
             "readModelVariableSubdomain--unable to close %s(%i)",
             modelName,modelFile);
          return FALSE;
      }
      for(int kk=0;kk<subdomainNZ;kk++){
        for(int jj=0;jj<subdomainNY;jj++){
          for(int ii=0;ii<subdomainNX;ii++){
            int subdomainIndex=ii+jj*subdomainNX+kk*subdomainNXY;
            int varI=ii+readLimit[0]-lim[0],varJ=jj+readLimit[2]-lim[2],varK=kk+readLimit[4]-lim[4];
            int varIndex=varI+varJ*nx+varK*nxy;
            assert(varIndex<nxyz,
             "readModelVariableSubdomain--%s[%i,%i,%i] (%i) out of range (%i,%i,%i)",
             varName,varI,varJ,varK,varIndex,nx,ny,nz);
            var[varIndex]=subdomainVar[subdomainIndex];
          }
        }
      }
      free(subdomainVar);
    }
  }
  //Check that the var has been completely filled (no -1's remain).
  for(int k=0;k<nz;k++)
    for(int j=0;j<ny;j++)
      for(int i=0;i<nx;i++)
	assert(var[i+j*nx+k*nxy]!=-1,
	       "readModelVariableSubdomain--%s[%i,%i,%i(%i)] not filled",
	       varName,i,j,k,i+j*nx+k*nxy);

  if(iOpenFile)
    assert(nc_close(modelFile)==NC_NOERR,
	   "readModelVariableSubdomain--unable to close %s(%i)",
	   modelName,modelFile);
  return TRUE;
}
static int readModelVariableSubdomain(float* var,
				      char* modelName,const char* varName,
				      int *lim,
				      size_t start[3],size_t size[3],
				      int nx,int ny,int nz,
				      int modelFile,int iOpenFile){
  int subdomainFilenameVar;
  if(nc_inq_varid(modelFile,"subdomainFilename",&subdomainFilenameVar)!=NC_NOERR)
    return FALSE;

  //Read some important dimensions for reading these files.
  int dimID;
  size_t nSubdomains,subdomainFileNameLength;
  assert(nc_inq_dimid(modelFile,"nSubdomains",&dimID)==NC_NOERR,
	 "readModelVariableSubdomain--unable to read dim nSubdomains from %s(%i)",
	 modelName,modelFile);
  assert(nc_inq_dimlen(modelFile,dimID,&nSubdomains)==NC_NOERR,
	 "readModelVariableSubdomain--unable to get dim nSubdomains from %s(%i)",
	 modelName,modelFile);

  assert(nc_inq_dimid(modelFile,"subdomainFileNameLength",&dimID)==NC_NOERR,
	 "readModelVariableSubdomain--unable to read dim subdomainFileNameLength from %s(%i)",
	 modelName,modelFile);
  assert(nc_inq_dimlen(modelFile,dimID,&subdomainFileNameLength)==NC_NOERR,
	 "readModelVariableSubdomain--unable to get dim subdomainFileNameLength from %s(%i)",
	 modelName,modelFile);

  //Get and additional variable that defines the proclimits within each
  // subdomain file.
  int subdomainLimitVar;
  assert(nc_inq_varid(modelFile,"subdomainLimit",&subdomainLimitVar)==NC_NOERR,
	 "readModelVariableSubdomain--unable to get var subdomainLimit from %s(%i)",
	 modelName,modelFile);

  //Define start and size variables for reading the subdomain components.
  size_t subdomainStart[2]={0,0};
  size_t subdomainLimitSize[2]={1,6};
  size_t subdomainNameSize[2]={1,subdomainFileNameLength};
    
  //Fill the buffer with the path to the current file, prepend this path to the
  // subdomain filenames that are to be read from the file.
  char buffer[1024];
  strcpy(buffer,modelName);
  char* bufferRootPos=strrchr(buffer,'/');
	if(bufferRootPos==NULL)
		bufferRootPos=buffer;
	else
		bufferRootPos+=1;

  //Loop through the subdomains and read the required parts to fill in the 
  // variable. Note that to check if every element of the variable is filled in
  // we are pre-filling with -1e10.
  int nxy=nx*ny,nxyz=nx*ny*nz;
  for(int i=0;i<nxyz;var[i++]=-1e10);
  for(int i=0;i<nSubdomains;i++){
    subdomainStart[0]=i;

    //Read the current subdomain limits and check for overlap.
    int subdomainLimit[6];
    assert(nc_get_vara_int(modelFile,subdomainLimitVar,
			   subdomainStart,subdomainLimitSize,
			   subdomainLimit)==NC_NOERR,
	   "readModelVariableSubdomain--unable to read subdomainLimit[%i,%i,%i,%i] from %s(%i,%i)",
	   subdomainStart[0],subdomainStart[1],
	   subdomainLimitSize[0],subdomainLimitSize[1],
	   modelName,modelFile,subdomainLimitVar);
      
    if(ISOVERLAP(lim[0],lim[1],subdomainLimit[0],subdomainLimit[1]) &&
       ISOVERLAP(lim[2],lim[3],subdomainLimit[2],subdomainLimit[3]) &&
       ISOVERLAP(lim[4],lim[5],subdomainLimit[4],subdomainLimit[5])){
      //Read the portion that has overlap only.
      int readLimit[6]={MAX(lim[0],subdomainLimit[0]),
			MIN(lim[1],subdomainLimit[1]),
			MAX(lim[2],subdomainLimit[2]),
			MIN(lim[3],subdomainLimit[3]),
			MAX(lim[4],subdomainLimit[4]),
			MIN(lim[5],subdomainLimit[5])};
      size_t subdomainReadStart[3]={readLimit[4]-subdomainLimit[4],
				    readLimit[2]-subdomainLimit[2],
				    readLimit[0]-subdomainLimit[0]};
      size_t subdomainReadSize[3]={readLimit[5]-readLimit[4],
				   readLimit[3]-readLimit[2],
				   readLimit[1]-readLimit[0]};

      int subdomainNX=subdomainReadSize[2],subdomainNY=subdomainReadSize[1],subdomainNZ=subdomainReadSize[0];
      int subdomainNXY=subdomainNX*subdomainNY;
      int subdomainNXYZ=subdomainNX*subdomainNY*subdomainNZ;

      //Get the filename and open for reading. 
      assert(nc_get_vara_text(modelFile,subdomainFilenameVar,
			      subdomainStart,subdomainNameSize,
			      bufferRootPos)==NC_NOERR,
	     "readModelVariableSubdomain--unable to read subdomainName[%i,%i,%i,%i] from %s(%i,%i)",
	     subdomainStart[0],subdomainStart[1],
	     subdomainNameSize[0],subdomainNameSize[1],
	     modelName,modelFile,subdomainFilenameVar);
      bufferRootPos[subdomainNameSize[1]]='\0';
      int subdomainFile;
      assert((subdomainFile=openCDFFile(buffer,FALSE,NC_NOWRITE))>=0,
	     "readModelVariableSubdomain--unable to open %s for reading",buffer);

      //Allocate a variable to read the required portion of the subdomain into.
      // Should be able to do this in place with one of the "funky" nc_get
      // subroutines but that will have to wait for a later time.
      float* subdomainVar;
      assert((subdomainVar=(float*)malloc(subdomainNXYZ*sizeof(float)))!=NULL,
	     "readModelVariableSubdomain--unable to allocate %i (%i,%i,%i) floats for local space",
	     subdomainReadSize[2],subdomainReadSize[1],subdomainReadSize[0],
	     subdomainNXYZ);

      //Read the data into the local space and move to the correct part of the
      // var.
      if(!readModelVariable3D(subdomainVar,subdomainFile,buffer,varName,
                              subdomainReadStart,subdomainReadSize,TRUE)) {
        free(subdomainVar);
        if(iOpenFile)
          assert(nc_close(modelFile)==NC_NOERR,
                 "readModelVariableSubdomain--unable to close %s(%i)",
                 modelName,modelFile);
        return FALSE;
      }
      for(int kk=0;kk<subdomainNZ;kk++){
        for(int jj=0;jj<subdomainNY;jj++){
          for(int ii=0;ii<subdomainNX;ii++){
            int subdomainIndex=ii+jj*subdomainNX+kk*subdomainNXY;
            int varI=ii+readLimit[0]-lim[0],varJ=jj+readLimit[2]-lim[2],varK=kk+readLimit[4]-lim[4];
            int varIndex=varI+varJ*nx+varK*nxy;
            assert(varIndex<nxyz,
             "readModelVariableSubdomain--%s[%i,%i,%i] (%i) out of range (%i,%i,%i)",
             varName,varI,varJ,varK,varIndex,nx,ny,nz);
            var[varIndex]=subdomainVar[subdomainIndex];
          }
        }
      }
      free(subdomainVar);
    }
  }
  //Check that the var has been completely filled (no -1's remain).
  for(int k=0;k<nz;k++)
    for(int j=0;j<ny;j++)
      for(int i=0;i<nx;i++)
	assert(var[i+j*nx+k*nxy]!=-1e10,
	       "readModelVariableSubdomain--%s[%i,%i,%i(%i)] not filled",
	       varName,i,j,k,i+j*nx+k*nxy);

  if(iOpenFile)
    assert(nc_close(modelFile)==NC_NOERR,
	   "readModelVariableSubdomain--unable to close %s(%i)",
	   modelName,modelFile);
  return TRUE;
}
static int readModelVariableIndexSet(float* var,
				     char* modelName,const char* varName,
				     int *lim,
				     int nx,int ny,int nz,
				     int globalNX,int globalNY,int globalNZ,
				     int modelFile,int iOpenFile){
  //New 8/28/03 allow vp, vs, and rho to be set with indexed value sets.
  // If a given value contains all indicies set all nodes to that value
  // set.
  int indexModelDim;
  if(nc_inq_dimid(modelFile,"indexModelDim",&indexModelDim)!=NC_NOERR)
    return FALSE;

  int globalNXY=globalNX*globalNY;
  int nxy=nx*ny,nxyz=nx*ny*nz;

  size_t nIndexModel;
  int indexModelVar;
  assert(nc_inq_dimlen(modelFile,indexModelDim,&nIndexModel)==NC_NOERR,
	 "readModelVariableIndexSet--unable to read length of indexModelDim from %s(%i,%i)",
	 modelName,modelFile,indexModelDim);
  char varNameL[64],varNameU[64];
  sprintf(varNameL,"indexModel%s",varName);
  sprintf(varNameU,"indexModel%s",varName);
  varNameL[10]=tolower(varName[0]);
  varNameU[10]=toupper(varName[0]);
  
  assert(nc_inq_varid(modelFile,varNameL,&indexModelVar)==NC_NOERR ||
	 nc_inq_varid(modelFile,varNameU,&indexModelVar)==NC_NOERR,
	 "readModelVariableIndexSet--unable to get var %s from %s(%i)",
	 varNameU,modelName,modelFile);

  //Step through each model value set.
  for(size_t i=0;i<nIndexModel;i++){
    float currV;
    assert(nc_get_var1_float(modelFile,indexModelVar,&i,&currV)==NC_NOERR,
	   "readModelVariableIndexSet--unable to read var[%i] %s from %s(%i,%i)",
	   i,varNameU,modelName,modelFile,indexModelVar);

    //Check for a variable corresponding to indicies for this model index.
    char dimName[512],varName[512];
    sprintf(dimName,"indexModel%ldDim",i);
    int indexVarDim;
    if(nc_inq_dimid(modelFile,dimName,&indexVarDim)==NC_NOERR){
      size_t nIndexVar;
      int indexVar;
      sprintf(varName,"indexModel%ldIndicies",i);
      assert(nc_inq_dimlen(modelFile,indexVarDim,&nIndexVar)==NC_NOERR,
	     "readModelVariableIndexSet--unable to read length of %s from %s(%i,%i)",
	     dimName,modelName,modelFile,indexModelDim);
      assert(nc_inq_varid(modelFile,varName,&indexVar)==NC_NOERR,
	     "readModelVariableIndexSet--unable to get var %s from %s(%i)",
	     varName,modelName,modelFile);

      int *indicies;
      assert((indicies=(int*)malloc(nIndexVar*sizeof(int)))!=NULL,
	     "readModelVariableIndexSet--unable to allocate %i ints for indices",
	     nIndexVar);
      assert(nc_get_var_int(modelFile,indexVar,indicies)==NC_NOERR,
	     "readModelVariableIndexSet--unable to read var %s from %s(%i,%i)",
	     varName,modelName,modelFile,indexVar);
      for(int i=0;i<nIndexVar;i++){
	int kk=(int)floor((float)indicies[i]/(float)globalNXY);
	int jj=(int)floor((float)(indicies[i]-kk*globalNXY)/globalNX);
	int ii=indicies[i]-kk*globalNXY-jj*globalNX;
	if(ISMID(lim[0],ii,lim[1]-1) &&
	   ISMID(lim[2],jj,lim[3]-1) &&
	   ISMID(lim[4],kk,lim[5]-1)){
	  int index=(ii-lim[0])+(jj-lim[2])*nx+(kk-lim[4])*nxy;
	  assert(ISMID(0,index,nxyz-1),
		 "readModelVariableIndexSet--index %i=>%i out of bounds (%i,%i,%i)=>(%i,%i,%i)",
		 indicies[i],index,ii,jj,kk,ii-lim[0],jj-lim[2],kk-lim[4]);
	  var[index]=currV;
	}
      }
      free(indicies);
    }else{
      for(int i=0;i<nxyz;i++){
	var[i]=currV;
      }
    }
  }

  if(iOpenFile)
    assert(nc_close(modelFile)==NC_NOERR,
	   "readModelVariableSubdomain--unable to close %s(%i)",
	   modelName,modelFile);
  return TRUE;
}

float* readModelVariable(floatPtr &var,
			 char* modelName,const char* varName,
			 int *lim,int globalNX,int globalNY,int globalNZ,
			 int modelFile,int noFindError){
  char buffer[512];

  //Check if the file is provided, if not then open it.
  int iOpenFile=FALSE;
  if(modelFile<0){
    iOpenFile=TRUE;
    assert((modelFile=openCDFFile(modelName,FALSE,NC_NOWRITE))>=0,
	   "readModelVariable--unable to open %s for reading",modelName);
  }
  
  //Set up the limits for reading.
  size_t start[3]={lim[4],lim[2],lim[0]}; //last dimension varies fastest
  size_t size[3]={lim[5]-lim[4],lim[3]-lim[2],lim[1]-lim[0]};
  int nx=size[2],ny=size[1],nz=size[0];
  int nxy=nx*ny,nxyz=nx*ny*nz;

  //Check if the variable needs to be allocated.
  if(!var)
    assert((var=(float*)malloc(nxyz*sizeof(float)))!=NULL,
	   "readModelVariable--unable to allocate %i (%ix%ix%i) for %s",
	   nxyz,nx,ny,nz,varName);

  //Check for the existance of the variable, check for the variable name with 
  // the first character in upper case and lower case at the same time.
  if(readModelVariable3D(var,modelFile,modelName,varName,
			 start,size,iOpenFile)){
    return var;
  }

  //Try for a 2D of the variable.
  if(readModelVariable2DXZ(var,modelFile,modelName,varName,
			   start,size,iOpenFile)){
    return var;
  }

  //Try for a 1D version of the variable name.
  int varId;
  char varNameL[64],varNameU[64];
  sprintf(varNameL,"oneDModel%s",varName);sprintf(varNameU,"oneDModel%s",varName);
  varNameL[9]=tolower(varName[0]);
  varNameU[9]=toupper(varName[0]);
  if(nc_inq_varid(modelFile,varNameL,&varId)==NC_NOERR ||
     nc_inq_varid(modelFile,varNameU,&varId)==NC_NOERR){
    //Success, read the variable into locally allocated 1D space.
    float* oneDVar;
    assert((oneDVar=(float*)malloc(globalNZ*sizeof(float)))!=NULL,
	   "readModelVariable--unable to allocate %i for %s",nz,varNameU);
    assert(nc_get_var_float(modelFile,varId,oneDVar)==NC_NOERR,
	   "readModelVariable--unable to read %s in file %s(%i,%i)",
	   varNameU,modelName,modelFile,varId);

    //Use for loops to fill in the 3D variable with these 1D values.
    for(int k=0;k<nz;k++){
      for(int i=0;i<nxy;i++){
	int index=i+k*nxy;
	var[index]=oneDVar[k+lim[4]];
      }
    }

    //Close the file and exit.
    free(oneDVar);
    if(iOpenFile)
      assert(nc_close(modelFile)==NC_NOERR,
	     "readModelVariable--unable to close %s (%i)",modelName,modelFile);
    return var;
  }

  //Try and index set approach.
  if(readModelVariableIndexSet(var,modelName,varName,
			       lim,nx,ny,nz,
			       globalNX,globalNY,globalNZ,
			       modelFile,iOpenFile))
    return var;

  //Try for a different file with the variable name appended.
  int modelFile2;
  sprintf(buffer,modelName);
  char* bPtr=strstr(buffer,".cdf");
  if(!bPtr) bPtr=buffer+strlen(buffer);
  sprintf(bPtr,"_%s.cdf",varName);
  if((modelFile2=openCDFFile(buffer,FALSE,NC_NOWRITE))>0){
    assert(readModelVariable3D(var,modelFile2,buffer,varName,
			       start,size,TRUE),
	   "readModelVariable--unable to read %s from %s(%i)",
	   varName,buffer,modelFile2);
    if(iOpenFile)
      assert(nc_close(modelFile)==NC_NOERR,
	     "readModelVariable--unable to close %s (%i)",modelName,modelFile);
    return var;    
  }
   
  //Look for a variable in the model that defines files with subdomains.
  if(readModelVariableSubdomain(var,modelName,varName,lim,
				start,size,nx,ny,nz,
				modelFile,iOpenFile))
    return var;
  
  assert(!noFindError,"readModelVariable--unable to read any form of %s from %s",
	 varName,modelName);
  return NULL; //Only executed if noFindError is false.
}
 
int* readModelVariable(intPtr &var,
                       char* modelName,const char* varName,
                       int *lim,int globalNX,int globalNY,int globalNZ,
                       int modelFile,int noFindError){
  char buffer[512];
  bool iAlloc = false;
  
  //Check if the file is provided, if not then open it.
  int iOpenFile=FALSE;
  if(modelFile<0){
    iOpenFile=TRUE;
    assert((modelFile=openCDFFile(modelName,FALSE,NC_NOWRITE))>=0,
           "readModelVariable--unable to open %s for reading",modelName);
  }
  
  //Set up the limits for reading.
  size_t start[3]={lim[4],lim[2],lim[0]}; //last dimension varies fastest
  size_t size[3]={lim[5]-lim[4],lim[3]-lim[2],lim[1]-lim[0]};
  int nx=size[2],ny=size[1],nz=size[0];
  int nxyz=nx*ny*nz;
  
  //Check if the variable needs to be allocated.
  if(!var) {
    assert((var=(int*)malloc(nxyz*sizeof(int)))!=NULL,
           "readModelVariable--unable to allocate %i (%ix%ix%i) for %s",
           nxyz,nx,ny,nz,varName);
    iAlloc = true;
  }
  
  //Check for the existance of the variable, check for the variable name with
  // the first character in upper case and lower case at the same time.
  if(readModelVariable3D(var,modelFile,modelName,varName,
                         start,size,iOpenFile)){
    return var;
  }
  
  //Try for a different file with the variable name appended.
  int modelFile2;
  sprintf(buffer,modelName);
  char* bPtr=strstr(buffer,".cdf");
  if(!bPtr) bPtr=buffer+strlen(buffer);
  sprintf(bPtr,"_%s.cdf",varName);
  if((modelFile2=openCDFFile(buffer,FALSE,NC_NOWRITE))>0){
    assert(readModelVariable3D(var,modelFile2,buffer,varName,
                               start,size,TRUE),
           "readModelVariable--unable to read %s from %s(%i)",
           varName,buffer,modelFile2);
    if(iOpenFile)
      assert(nc_close(modelFile)==NC_NOERR,
             "readModelVariable--unable to close %s (%i)",modelName,modelFile);
    return var;
  }
  
  //Look for a variable in the model that defines files with subdomains.
  if(readModelVariableSubdomain(var,modelName,varName,lim,
                                start,size,nx,ny,nz,
                                modelFile,iOpenFile))
    return var;
  
  assert(!noFindError,"readModelVariable--unable to read any form of %s from %s",
         varName,modelName);
  if(iAlloc) {
    free(var);
    var = NULL;
  }
  return NULL; //Only executed if noFindError is false.
}

