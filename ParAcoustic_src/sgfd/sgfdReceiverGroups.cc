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
 *  sgfdReceiverGroups.cc
 *
 *
 *  Defines class functions used in groups of receivers, mainly
 *  of class receiverGrid.  These handle common types of receivers
 *  that vary in location, but otherwise are the same.
 *
 *  Defines the following class functions:
 *  receiverGrid::receiverGrid
 *  receiverGrid::receiverGrid
 *  receiverGrid::packInit
 *  receiverGrid::readData
 *  receiverGrid::packData
 *  receiverGrid::fill
 *  pressureReceiverGrid::fillvalue
 *
 */

#include "sgfdReceiverGroups.hh"
#include "message_passing.h"
///Here is the slave constructor, get information from a message that has
/// already been received.
receiverGrid::receiverGrid(modelDefStruct* modelDef,int gridIndex):receiver(modelDef){
  _gridIndex=gridIndex;

  //The parent has already read x, y, z, and amp. Now read dx, nx, dy, ny, dz,
  // and nz.
  unpackMessage("ii ii ii ii",&_gridTimeSubsample,&_rawOutput,
    &_dx,&_nx,&_dy,&_ny,&_dz,&_nz);
  if(_rawOutput)
    unpackMessage("s",_rawOutputDir);

  _nt=1+(int)floor((modelDef->NT-1)/_gridTimeSubsample);
  
  //Now modify _x, _y, _z, _nx, _ny, and _nz as required to fit the current 
  // subdomain.
  DEF_MODEL_SIZE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  float maxX=minX+dx*(NX-1),maxY=minY+dy*(NY-1),maxZ=minZ+dz*(NZ-1);

  _ix=_iy=_iz=0;
  if(!ISMID_CO(minX+dx,_x,maxX-2*dx) || 
     !ISMID_CO(minX+dx,_x+dx*_dx*(_nx-1),maxX-2*dx)){
    if(_dx>0 && _x<=minX+dx){
_ix=1+(int)floor((minX+dx-_x)/(dx*_dx));
_x+=dx*_dx*_ix;
_nx-=_ix;
    }else if(_dx<0 && _x>maxX-2*dx){
_ix=(int)ceil((_x-(maxX-2.*dx))/(dx*_dx));
_x-=dx*_dx*_ix;
_nx-=_ix;
    }

    if(_dx>0 && (_x+dx*_dx*(_nx-1))>maxX-2*dx){
_nx=1+(int)floor((maxX-2*dx-_x)/(dx*_dx));
    }else if(_dx<0 && (_x+dx*_dx*(_nx-1))<(minX+dx)){
_nx=(int)ceil((_x-(minX+dx))/(dx*_dx));
    }

    _nx=MAX(_nx,0);
  }
  if(!ISMID_CO(minY+dy,_y,maxY-2*dy) || 
     !ISMID_CO(minY+dy,_y+dy*_dy*(_ny-1),maxY-2*dy)){
    if(_dy>0 && _y<=minY+dy){
_iy=1+(int)floor((minY+dy-_y)/(dy*_dy));
_y+=dy*_dy*_iy;
_ny-=_iy;
    }else if(_dy<0 && _y>maxY-2*dy){
_iy=(int)ceil((_y-(maxY-2.*dy))/(dy*_dy));
_y-=dy*_dy*_iy;
_ny-=_iy;
    }

    if(_dy>0 && (_y+dy*_dy*(_ny-1))>maxY-2*dy){
_ny=1+(int)floor((maxY-2*dy-_y)/(dy*_dy));
    }else if(_dy<0 && (_y+dy*_dy*(_ny-1))<(minY+dy)){
_ny=(int)ceil((_y-(minY+dy))/(dy*_dy));
    }

    _ny=MAX(_ny,0);
  }
  if(!ISMID_CO(minZ+dz,_z,maxZ-2*dz) || 
     !ISMID_CO(minZ+dz,_z+dz*_dz*(_nz-1),maxZ-2*dz)){
    if(_dz>0 && _z<=minZ+dz){
_iz=1+(int)floor((minZ+dz-_z)/(dz*_dz));
_z+=dz*_dz*_iz;
_nz-=_iz;
    }else if(_dz<0 && _z>maxZ-2*dz){
_iz=(int)ceil((_z-(maxZ-2.*dz))/(dz*_dz));
_z-=dz*_dz*_iz;
_nz-=_iz;
    }

    if(_dz>0 && (_z+dz*_dz*(_nz-1))>maxZ-2*dz){
_nz=1+(int)floor((maxZ-2*dz-_z)/(dz*_dz));
    }else if(_dz<0 && (_z+dz*_dz*(_nz-1))<=(minZ+dz)){
_nz=(int)ceil((_z-(minZ+dz))/(dz*_dz));
    }

    _nz=MAX(_nz,0);
  }

  if(_rawOutput==2){
    assert((_xcorr=(float*)malloc(_nx*_ny*_nz*sizeof(float)))!=NULL,
     "receiverGrid--unable to allocate %ix%ix%i grid for xcorrelation",
     _nx,_ny,_nz);
    for(int i=0;i<_nx*_ny*_nz;_xcorr[i++]=0);
  }
}

///Here is the master version, get the info from arguments that are passed in.
receiverGrid::receiverGrid(modelDefStruct* modelDef,int allocate,float amp,
       float x,float rdx,int nx,
       float y,float rdy,int ny,
       float z,float rdz,int nz,
       int timeSubsample,
       int rawOutput,char* rawOutputDir):
  receiver(modelDef,x,y,z,amp){
  _rawOutput=rawOutput;
  if(!rawOutputDir){
    _rawOutputDir[0]='\0';
  }else{
    strcpy(_rawOutputDir,rawOutputDir);
  }

  //Set counts.
  _nx=nx;
  _ny=ny;
  _nz=nz;

  //Convert dx, dy, dz from a floating point to and integer #of grid points.
  _dx=(int)round(rdx/modelDef->dx);
  _dy=(int)round(rdy/modelDef->dy);
  _dz=(int)round(rdz/modelDef->dz);

  //Do some range checking.
  DEF_MODEL_SIZE(modelDef);
  DEF_MODEL_LIMITS(modelDef);
  float maxX=minX+dx*(NX-1),maxY=minY+dy*(NY-1),maxZ=minZ+dz*(NZ-1);

  assert(ISMID(minX+2.*dx,x,maxX-2*dx),
   "receiverGrid--start X (%.2f) out of range [%.2f-%.2f]\n",
   x,minX+2.*dx,maxX-2*dx);
  assert(ISMID(minX+2.*dx,x+dx*_dx*(_nx-1),maxX-2*dx),
   "receiverGrid--end X (%.2f) out of range [%.2f-%.2f]\n",
   x+dx*_dx*(_nx-1),minX+2.*dx,maxX-2*dx);

  assert(ISMID(minY+2.*dy,y,maxY-2*dy),
   "receiverGrid--start Y (%.2f) out of range [%.2f-%.2f]\n",
   y,minY+2.*dy,maxY-2*dy);
  assert(ISMID(minY+2.*dy,y+dy*_dy*(_ny-1),maxY-2*dy),
   "receiverGrid--end Y (%.2f) out of range [%.2f-%.2f]\n",
   y+dy*_dy*(_ny-1),minY+2.*dy,maxY-2*dy);

  assert(ISMID(minZ+2.*dz,z,maxZ-2*dz),
   "receiverGrid--start Z (%.2f) out of range [%.2f-%.2f]\n",
   z,minZ+2.*dz,maxZ-2*dz);
  assert(ISMID(minZ+2.*dz,z+dz*_dz*(_nz-1),maxZ-2*dz),
   "receiverGrid--end Z (%.2f) out of range [%.2f-%.2f]\n",
   z+dz*_dz*(_nz-1),minZ+2.*dz,maxZ-2*dz);

  //Setup the time stuff.
  _gridTimeSubsample=timeSubsample;
  _nt=1+(int)floor((modelDef->NT-1)/_gridTimeSubsample);
}

///Pack a message with the parameters for this receiver. Assumes the buffer is
/// ready to be filled.
void receiverGrid::packInit(){
  receiver::packInit();

  packMessage("ii ii ii ii",_gridTimeSubsample,_rawOutput,
  _dx,_nx,_dy,_ny,_dz,_nz);
  if(_rawOutput)
    packMessage("s",_rawOutputDir);
}

///Read the data in preperation for performing RTM cross-correlation.
int receiverGrid::readData(modelDefStruct* modelDef){
  char buffer[1024];
  sprintf(buffer,"%s/Grid%i_%i-%i_%i-%i_%i-%i_Raw.bin",
    _rawOutputDir,_gridIndex,
    modelDef->procLim[0],modelDef->procLim[1],
    modelDef->procLim[2],modelDef->procLim[3],
    modelDef->procLim[4],modelDef->procLim[5]);
  FILE* in=fopen(buffer,"rb");
  assert(buffer!=NULL,"receiverGrid::readData--unable to open %s.",buffer);

  //Check all the sizes match.
  int temp;
  fread(&temp,sizeof(int),1,in);
  assert(temp==_nt,"receiverGrid::readData--NT mismatch %i!=%i.",temp,_nt);

  fread(&temp,sizeof(int),1,in);
  assert(temp==_ix,"receiverGrid::readData--ix mismatch %i!=%i.",temp,_ix);
  fread(&temp,sizeof(int),1,in);
  assert(temp==_nx,"receiverGrid::readData--nx mismatch %i!=%i.",temp,_nx);

  fread(&temp,sizeof(int),1,in);
  assert(temp==_iy,"receiverGrid::readData--iy mismatch %i!=%i.",temp,_iy);
  fread(&temp,sizeof(int),1,in);
  assert(temp==_ny,"receiverGrid::readData--ny mismatch %i!=%i.",temp,_ny);

  fread(&temp,sizeof(int),1,in);
  assert(temp==_iz,"receiverGrid::readData--iz mismatch %i!=%i.",temp,_iz);
  fread(&temp,sizeof(int),1,in);
  assert(temp==_nz,"receiverGrid::readData--nz mismatch %i!=%i.",temp,_nz);

  int nxyz=_nx*_ny*_nz;
  float *wrk=(float*)malloc(_nt*sizeof(float));
  assert(wrk!=NULL,"receiverGrid::readData--unable to allocate wrk (%i).",_nt);
  for(int i=0;i<nxyz;i++){
    fread(&temp,sizeof(int),1,in);
    assert(temp==i,"receiverGrid::readData--count mismatch %i!=%i.",temp,i);
    fread(wrk,sizeof(float),(size_t)_nt,in);
    //Flip in time.
    for(int j=0;j<_nt;j++)
_data[j+i*_nt]=wrk[_nt-1-j];
  }
  fclose(in);
  free(wrk);

  return TRUE;
}

///Pack the data, since this one does multiple sends do our own buffer init.
int receiverGrid::packData(modelDefStruct* modelDef,float* wrk,
         int startI,int stopI,int decimation){
  if(_rawOutput==2){
    initSend();
    packMessage("i ii ii ii",_nt,_ix,_nx,_iy,_ny,_iz,_nz);
    int nxyz=_nx*_ny*_nz;
    if(nxyz>0)
packMessage("Fi",_xcorr,nxyz,_gridIndex);
    sendMessage(Parent,MESSAGE_SEND_RECEIVER_GRID,NULL);
  }else if(_rawOutput){
    char buffer[1024];
    sprintf(buffer,"%s/Grid%i_%i-%i_%i-%i_%i-%i_Raw.bin",
      _rawOutputDir,_gridIndex,
      modelDef->procLim[0],modelDef->procLim[1],
      modelDef->procLim[2],modelDef->procLim[3],
      modelDef->procLim[4],modelDef->procLim[5]);
    FILE* out=fopen(buffer,"wb");
    fwrite(&_nt,sizeof(int),1,out);
    fwrite(&_ix,sizeof(int),1,out);fwrite(&_nx,sizeof(int),1,out);
    fwrite(&_iy,sizeof(int),1,out);fwrite(&_ny,sizeof(int),1,out);
    fwrite(&_iz,sizeof(int),1,out);fwrite(&_nz,sizeof(int),1,out);

    int nxyz=_nx*_ny*_nz;
    for(int i=0;i<nxyz;i++){
fwrite(&i,sizeof(int),1,out);
fwrite(_data+i*_nt,sizeof(float),(size_t)_nt,out);
    }
    fclose(out);
    sendMessage(Parent,MESSAGE_SEND_RECEIVER_GRID,"i",_gridIndex);
  }else{
    initSend();
    sendMessage(Parent,MESSAGE_SEND_RECEIVER_GRID,
    "i ii ii ii",_nt,_ix,_nx,_iy,_ny,_iz,_nz);

    int nxyz=_nx*_ny*_nz;
    if(nxyz>0){
scaleData(modelDef,startI,stopI,NULL,NULL);
int currPacked=0,currSend=1;
initSend(); 
for(int i=0;i<nxyz;i++){
  packMessage("iF",i,_data+i*_nt,_nt);
  if(++currPacked==MAX_R_PER_SEND){
    sendMessage(Parent,MESSAGE_SEND_RECEIVER_GRID+currSend++,NULL);
    initSend();
    currPacked=0;
  }
}
if(currPacked)
  sendMessage(Parent,MESSAGE_SEND_RECEIVER_GRID+currSend++,NULL);
    }
  }
  return _nx*_ny*_nz*_nt;
}

///Define the fill method as an iterator, then the real subclasses just need
/// to define a method to actually fill the data.
float receiverGrid::fill(modelDefStruct* modelDef,int iteration,
       float* vx,float* vy,float* vz,
       float* xx,float* yy,float* zz,
       float* xy,float* xz,float* yz){
  if(!active()) 
    return -10;

  if(!(iteration%_gridTimeSubsample)){
    int tindex=(int)round(iteration/_gridTimeSubsample);
    if(tindex<0||tindex>=_nt)
return -10;

    for(int k=0;k<_nz;k++){
int koffset=k*_nx*_ny;
for(int j=0;j<_ny;j++){
  int joffset=j*_nx;
  for(int i=0;i<_nx;i++){
    int gindex=i+joffset+koffset;
    int index=tindex+_nt*gindex;
    assert(index>=0 && index<_n,
     "receverGrid::fill--index %i out of bounds 0 to %i.",index,_n);

    if(_rawOutput==2){
      _xcorr[gindex]+=_data[index]*
  fillvalue(modelDef,i,j,k,vx,vy,vz,xx,yy,zz,xy,xz,yz);
    }else{
      _data[index]=
  fillvalue(modelDef,i,j,k,vx,vy,vz,xx,yy,zz,xy,xz,yz);
    }
  }
}
    }
  }
  return 1.0;
}

float pressureReceiverGrid::fillvalue(modelDefStruct* modelDef,int i,int j,int k,
      float* vx,float* vy,float* vz,
      float* xx,float* yy,float* zz,
      float* xy,float* xz,float* yz){
  if(!yy){
    return _amp*
trilinInterpOffset(modelDef,&_interp,i*_dx,j*_dy,k*_dz,xx);
  }
  //Interpolate diagonal components of stress tensor onto receiver location.
  float ampx=trilinInterpOffset(modelDef,&_interp,i*_dx,j*_dy,k*_dz,xx);
  float ampy=trilinInterpOffset(modelDef,&_interp,i*_dx,j*_dy,k*_dz,yy);
  float ampz=trilinInterpOffset(modelDef,&_interp,i*_dx,j*_dy,k*_dz,zz);
    
  return -1.0/3.0*_amp*(ampx+ampy+ampz);
}  
