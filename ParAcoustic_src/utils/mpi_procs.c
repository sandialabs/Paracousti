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
 *  mpi_procs.c
 *
 *
 *  Defines functions used in interprocess communication.  Most of these are defined
 *  in message_passing.h
 *
 *  Defines the following functions:
 *  messageRank
 *  setMessageBuffer
 *  initSend
 *  currMessageBufferSize
 *  doBarrier
 *  packMessage
 *  sendMessage
 *  iSendMessage
 *  unpackMessage
 *  getMessage
 *  setBlocksType
 *  sendBlocks
 *  getBlocks
 *  immediateSendBlocks
 *  setBlockType
 *  packBlock
 *  unpackBlock
 *  checkForError
 *  sendError
 *  killProcesses
 *  stopProcesses
 *  registerProcess
 *  startProcesses
 *
 *  Defines the following macros:
 *  MESSAGE_SEND_BLOCKS
 *
 */

#include <stdio.h>
#include <signal.h>

#include "nstdutil.h"

#include "message_passing.h"
#include "mpi.h"

/*Global variable used in other processes*/
const int AnyMessageTag=MPI_ANY_TAG;
const int AllProcesses=MPI_ANY_SOURCE;
int NumProcs=-1;
int* Tids=NULL;
int Parent;
int _Rank_=-1;

/*Global variables for packing and unpacking data*/
int BufferSize=0,BufferPos=0;
char* Buffer=NULL;

//Extra variables required for general immediate sends.
MPI_Request ISendRequests[1000];
int ISendRequestIndex=0;

//Function to return the rank
int messageRank(){return _Rank_;}

/* Basic send and receive routines*/
int setMessageBuffer(int size){
  if(!Buffer){
    assert((Buffer=(char*)malloc(size))!=NULL,
           "setMessageBuffer--unable to allocate Buffer of size %i",
           size);
    BufferSize=size;
  }else if(size>BufferSize){
    assert((Buffer=(char*)realloc(Buffer,size))!=NULL,
           "setMessageBuffer--unable to re-allocate Buffer to size %i",
           size);
    BufferSize=size;
  }
  BufferPos=0;
  return BufferSize;
}

int initSend(void){
  BufferPos=0;
  return TRUE;
}
int currMessageBufferSize(){
  return BufferPos;
}

void doBarrier(){
  MPI_Barrier(MPI_COMM_WORLD);
}

/*send a message, initSend MUST BE CALLED FIRST
 format string for additional arguments
 i for int, f for float, s for string, P sends ProcLim (6 ints)
 F sends an array of floats requires two arguments a float* and int (length)
 I sends an array of ints requires two arguments a int* and int (length)
 B[n] sends  n blocks of floats from a 3d array using NX,NY, and NZ
 arguments are processed in order from left to right*/
int packMessage(const char* format,...){
  int count = 0;
  int ival,*iPtrVal;
  float fval,*fPtrVal;
  double dval,*dPtrVal;
  char* sval;
  int* procLim;
  
  va_list ap;
  va_start(ap,format);
  
  for(const char* p=format;*p;p++,count++){
    switch(*p){
      case ' ':
        //allow space to make messages more readable to humans.
        break;
      case 'B':
      {
        float *dummySpace[3];
        float **data=dummySpace;
        int nx,nxy,xMin,xMax,yMin,yMax,zMin,zMax;
        char* endptr;
        int i;
        
        //read the number of blocks to write
        int n=1;
        if(isdigit(*(p+1))){
          n=strtol(p+1,&endptr,10);
          p=endptr-1; //need to step back 1 so increment is next char
          if(n>3)
            data=(float**)malloc(n*sizeof(float*));
        }
        //read the grid size
        nx=va_arg(ap,int);
        nxy=va_arg(ap,int);
        
        //read the limits
        xMin=va_arg(ap,int);
        xMax=va_arg(ap,int);
        yMin=va_arg(ap,int);
        yMax=va_arg(ap,int);
        zMin=va_arg(ap,int);
        zMax=va_arg(ap,int);
        
        //fill the data
        for(i=0;i<n;i++)
          data[i]=va_arg(ap,float*);
        
        //pack the block
        packBlock(xMin,xMax,yMin,yMax,zMin,zMax,
                  n,nx,nxy,data);
        
        //free data if required
        if(n>3) free(data);
      }
        break;
      case 'i':
        ival=va_arg(ap,int);
        MPI_Pack(&ival,1,MPI_INT,
                 Buffer,BufferSize,&BufferPos,
                 MPI_COMM_WORLD);
        break;
      case 'I':
        iPtrVal=va_arg(ap,int*);
        ival=va_arg(ap,int);
        MPI_Pack(iPtrVal,ival,MPI_INT,
                 Buffer,BufferSize,&BufferPos,
                 MPI_COMM_WORLD);
        break;
      case 'f':
        fval=va_arg(ap,double);
        MPI_Pack(&fval,1,MPI_FLOAT,
                 Buffer,BufferSize,&BufferPos,
                 MPI_COMM_WORLD);
        break;
      case 'F':
        fPtrVal=va_arg(ap,float*);
        ival=va_arg(ap,int);
        MPI_Pack(fPtrVal,ival,MPI_FLOAT,
                 Buffer,BufferSize,&BufferPos,
                 MPI_COMM_WORLD);
        break;
      case 'd':
        dval=va_arg(ap,double);
        MPI_Pack(&dval,1,MPI_DOUBLE,
                 Buffer,BufferSize,&BufferPos,
                 MPI_COMM_WORLD);
        break;
      case 'D':
        dPtrVal=va_arg(ap,double*);
        ival=va_arg(ap,int);
        MPI_Pack(dPtrVal,ival,MPI_DOUBLE,
                 Buffer,BufferSize,&BufferPos,
                 MPI_COMM_WORLD);
        break;
      case 's':
        sval=va_arg(ap,char*);
        MPI_Pack(sval,512,MPI_CHAR,
                 Buffer,BufferSize,&BufferPos,
                 MPI_COMM_WORLD);
        break;
      case 'P':
        procLim=va_arg(ap,int*);
        MPI_Pack(procLim,6,MPI_INT,
                 Buffer,BufferSize,&BufferPos,
                 MPI_COMM_WORLD);
        break;
      default:
        assert(FALSE,"packMessage--unknown flag %s in format string %s",
               *p,format);
    }
  }
  va_end(ap);
  return count;
}

int sendMessage(int target,int tag,const char* format,...){
  if(format && format[0]){
    int ival,*iPtrVal;
    float fval,*fPtrVal;
    double dval,*dPtrVal;
    char* sval;
    int* procLim;
    
    va_list ap;
    va_start(ap,format);
    
    for(const char* p=format;*p;p++){
      switch(*p){
        case ' ':
          //allow space to make messages more readable to humans.
          break;
        case 'B':
        {
          float *dummySpace[3];
          float **data=dummySpace;
          int nx,nxy,xMin,xMax,yMin,yMax,zMin,zMax;
          char* endptr;
          int i;
          
          //read the number of blocks to write
          int n=1;
          if(isdigit(*(p+1))){
            n=strtol(p+1,&endptr,10);
            p=endptr-1; //need to step back 1 so increment is next char
            if(n>3)
              data=(float**)malloc(n*sizeof(float*));
          }
          //read the grid size
          nx=va_arg(ap,int);
          nxy=va_arg(ap,int);
          
          //read the limits
          xMin=va_arg(ap,int);
          xMax=va_arg(ap,int);
          yMin=va_arg(ap,int);
          yMax=va_arg(ap,int);
          zMin=va_arg(ap,int);
          zMax=va_arg(ap,int);
          
          //fill the data
          for(i=0;i<n;i++)
            data[i]=va_arg(ap,float*);
          
          //pack the block
          packBlock(xMin,xMax,yMin,yMax,zMin,zMax,
                    n,nx,nxy,data);
          
          //free data if required
          if(n>3) free(data);
        }
          break;
        case 'i':
          ival=va_arg(ap,int);
          MPI_Pack(&ival,1,MPI_INT,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        case 'I':
          iPtrVal=va_arg(ap,int*);
          ival=va_arg(ap,int);
          MPI_Pack(iPtrVal,ival,MPI_INT,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        case 'f':
          fval=va_arg(ap,double);
          MPI_Pack(&fval,1,MPI_FLOAT,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        case 'F':
          fPtrVal=va_arg(ap,float*);
          ival=va_arg(ap,int);
          MPI_Pack(fPtrVal,ival,MPI_FLOAT,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        case 'd':
          dval=va_arg(ap,double);
          MPI_Pack(&dval,1,MPI_DOUBLE,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        case 'D':
          dPtrVal=va_arg(ap,double*);
          ival=va_arg(ap,int);
          MPI_Pack(dPtrVal,ival,MPI_DOUBLE,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        case 's':
          sval=va_arg(ap,char*);
          MPI_Pack(sval,512,MPI_CHAR,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        case 'P':
          procLim=va_arg(ap,int*);
          MPI_Pack(procLim,6,MPI_INT,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        default:
          assert(FALSE,"sendMessage--unknown flag %s in format string %s",
                 *p,format);
      }
    }
    va_end(ap);
  }
  
  if(target==AllProcesses){
    int i;
    for(i=0;i<NumProcs;i++){
      MPI_Send(Buffer,BufferPos,
               MPI_PACKED,Tids[i],tag,
               MPI_COMM_WORLD);
    }
  }else{
    MPI_Send(Buffer,BufferPos,
             MPI_PACKED,target,tag,
             MPI_COMM_WORLD);
  }
  return tag;
}
int iSendMessage(int target,int tag,const char* format,...){
  if(format && format[0]){
    int ival,*iPtrVal;
    float fval,*fPtrVal;
    double dval,*dPtrVal;
    char* sval;
    int* procLim;
    
    va_list ap;
    va_start(ap,format);
    
    for(const char* p=format;*p;p++){
      switch(*p){
        case ' ':
          //allow space to make messages more readable to humans.
          break;
        case 'B':
        {
          float *dummySpace[3];
          float **data=dummySpace;
          int nx,nxy,xMin,xMax,yMin,yMax,zMin,zMax;
          char* endptr;
          int i;
          
          //read the number of blocks to write
          int n=1;
          if(isdigit(*(p+1))){
            n=strtol(p+1,&endptr,10);
            p=endptr-1; //need to step back 1 so increment is next char
            if(n>3)
              data=(float**)malloc(n*sizeof(float*));
          }
          //read the grid size
          nx=va_arg(ap,int);
          nxy=va_arg(ap,int);
          
          //read the limits
          xMin=va_arg(ap,int);
          xMax=va_arg(ap,int);
          yMin=va_arg(ap,int);
          yMax=va_arg(ap,int);
          zMin=va_arg(ap,int);
          zMax=va_arg(ap,int);
          
          //fill the data
          for(i=0;i<n;i++)
            data[i]=va_arg(ap,float*);
          
          //pack the block
          packBlock(xMin,xMax,yMin,yMax,zMin,zMax,
                    n,nx,nxy,data);
          
          //free data if required
          if(n>3) free(data);
        }
          break;
        case 'i':
          ival=va_arg(ap,int);
          MPI_Pack(&ival,1,MPI_INT,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        case 'I':
          iPtrVal=va_arg(ap,int*);
          ival=va_arg(ap,int);
          MPI_Pack(iPtrVal,ival,MPI_INT,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        case 'f':
          fval=va_arg(ap,double);
          MPI_Pack(&fval,1,MPI_FLOAT,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        case 'F':
          fPtrVal=va_arg(ap,float*);
          ival=va_arg(ap,int);
          MPI_Pack(fPtrVal,ival,MPI_FLOAT,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        case 'd':
          dval=va_arg(ap,double);
          MPI_Pack(&dval,1,MPI_DOUBLE,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        case 'D':
          dPtrVal=va_arg(ap,double*);
          ival=va_arg(ap,int);
          MPI_Pack(dPtrVal,ival,MPI_DOUBLE,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        case 's':
          sval=va_arg(ap,char*);
          MPI_Pack(sval,512,MPI_CHAR,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        case 'P':
          procLim=va_arg(ap,int*);
          MPI_Pack(procLim,6,MPI_INT,
                   Buffer,BufferSize,&BufferPos,
                   MPI_COMM_WORLD);
          break;
        default:
          assert(FALSE,"sendMessage--unknown flag %s in format string %s",
                 *p,format);
      }
    }
    va_end(ap);
  }
  
  if(target==AllProcesses){
    int i;
    for(i=0;i<NumProcs;i++){
      MPI_Isend(Buffer,BufferPos,
                MPI_PACKED,Tids[i],tag,
                MPI_COMM_WORLD,
                &ISendRequests[++ISendRequestIndex]);
    }
  }else{
    MPI_Isend(Buffer,BufferPos,
              MPI_PACKED,target,tag,
              MPI_COMM_WORLD,
              &ISendRequests[++ISendRequestIndex]);
  }
  
  return tag;
}

/*receive a message same format as sendMessage*/
int unpackMessage(const char* format,...){
  int count=0;
  int* ival,iDirectVal;
  float* fval;
  double* dval;
  char* sval;
  int* procLim;
  
  va_list ap;
  va_start(ap,format);
  
  for(const char* p=format;*p;p++){
    switch(*p){
      case ' ':
        //allow space to make messages more readable to humans.
        break;
      case 'B':
      {
        float *dummySpace[3];
        float **data=dummySpace;
        int nx,nxy,xMin,xMax,yMin,yMax,zMin,zMax;
        char* endptr;
        int i;
        
        //read the number of blocks to write
        int n=1;
        if(isdigit(*(p+1))){
          n=strtol(p+1,&endptr,10);
          p=endptr-1; //need to step back 1 so increment is next char
          if(n>3)
            data=(float**)malloc(n*sizeof(float*));
        }
        //read the grid size
        nx=va_arg(ap,int);
        nxy=va_arg(ap,int);
        
        //read the limits
        xMin=va_arg(ap,int);
        xMax=va_arg(ap,int);
        yMin=va_arg(ap,int);
        yMax=va_arg(ap,int);
        zMin=va_arg(ap,int);
        zMax=va_arg(ap,int);
        
        //fill the data
        for(i=0;i<n;i++)
          data[i]=va_arg(ap,float*);
        
        //pack the block
        unpackBlock(xMin,xMax,yMin,yMax,zMin,zMax,
                    n,nx,nxy,data);
        
        //free data if required
        if(n>3) free(data);
      }
        break;
      case 'i':
        ival=va_arg(ap,int*);
        MPI_Unpack(Buffer,BufferSize,&BufferPos,
                   ival,1,MPI_INT,MPI_COMM_WORLD);
        break;
      case 'I':
        ival=va_arg(ap,int*);
        iDirectVal=va_arg(ap,int);
        MPI_Unpack(Buffer,BufferSize,&BufferPos,
                   ival,iDirectVal,MPI_INT,MPI_COMM_WORLD);
        break;
      case 'f':
        fval=va_arg(ap,float*);
        MPI_Unpack(Buffer,BufferSize,&BufferPos,
                   fval,1,MPI_FLOAT,MPI_COMM_WORLD);
        break;
      case 'F':
        fval=va_arg(ap,float*);
        iDirectVal=va_arg(ap,int);
        MPI_Unpack(Buffer,BufferSize,&BufferPos,
                   fval,iDirectVal,MPI_FLOAT,MPI_COMM_WORLD);
        break;
      case 'd':
        dval=va_arg(ap,double*);
        MPI_Unpack(Buffer,BufferSize,&BufferPos,
                   dval,1,MPI_DOUBLE,MPI_COMM_WORLD);
        break;
      case 'D':
        dval=va_arg(ap,double*);
        iDirectVal=va_arg(ap,int);
        MPI_Unpack(Buffer,BufferSize,&BufferPos,
                   dval,iDirectVal,MPI_DOUBLE,MPI_COMM_WORLD);
        break;
      case 's':
        sval=va_arg(ap,char*);
        MPI_Unpack(Buffer,BufferSize,&BufferPos,
                   sval,512,MPI_CHAR,MPI_COMM_WORLD);
        break;
      case 'P':
        procLim=va_arg(ap,int*);
        MPI_Unpack(Buffer,BufferSize,&BufferPos,
                   procLim,6,MPI_INT,MPI_COMM_WORLD);
        break;
      default:
        assert(FALSE,"unpackMessage--unknown flag %s in format string %s",
               *p,format);
    }
  }
  va_end(ap);
  return count;
}

int getMessage(int source,int tag,const char* format,...){
  MPI_Status status;
  
  MPI_Recv(Buffer,BufferSize,MPI_PACKED,
           source,tag,
           MPI_COMM_WORLD,&status);
  BufferPos=0;
  //parse the format string
  if(format && format[0]){
    int* ival,iDirectVal;
    float* fval;
    double* dval;
    char* sval;
    int* procLim;
    
    va_list ap;
    va_start(ap,format);
    
    for(const char* p=format;*p;p++){
      switch(*p){
        case ' ':
          //allow space to make messages more readable to humans.
          break;
        case 'B':
        {
          float *dummySpace[3];
          float **data=dummySpace;
          int nx,nxy,xMin,xMax,yMin,yMax,zMin,zMax;
          char* endptr;
          int i;
          
          //read the number of blocks to write
          int n=1;
          if(isdigit(*(p+1))){
            n=strtol(p+1,&endptr,10);
            p=endptr-1; //need to step back 1 so increment is next char
            if(n>3)
              data=(float**)malloc(n*sizeof(float*));
          }
          //read the grid size
          nx=va_arg(ap,int);
          nxy=va_arg(ap,int);
          
          //read the limits
          xMin=va_arg(ap,int);
          xMax=va_arg(ap,int);
          yMin=va_arg(ap,int);
          yMax=va_arg(ap,int);
          zMin=va_arg(ap,int);
          zMax=va_arg(ap,int);
          
          //fill the data
          for(i=0;i<n;i++)
            data[i]=va_arg(ap,float*);
          
          //pack the block
          unpackBlock(xMin,xMax,yMin,yMax,zMin,zMax,
                      n,nx,nxy,data);
          
          //free data if required
          if(n>3) free(data);
        }
          break;
        case 'i':
          ival=va_arg(ap,int*);
          MPI_Unpack(Buffer,BufferSize,&BufferPos,
                     ival,1,MPI_INT,MPI_COMM_WORLD);
          break;
        case 'I':
          ival=va_arg(ap,int*);
          iDirectVal=va_arg(ap,int);
          MPI_Unpack(Buffer,BufferSize,&BufferPos,
                     ival,iDirectVal,MPI_INT,MPI_COMM_WORLD);
          break;
        case 'f':
          fval=va_arg(ap,float*);
          MPI_Unpack(Buffer,BufferSize,&BufferPos,
                     fval,1,MPI_FLOAT,MPI_COMM_WORLD);
          break;
        case 'F':
          fval=va_arg(ap,float*);
          iDirectVal=va_arg(ap,int);
          MPI_Unpack(Buffer,BufferSize,&BufferPos,
                     fval,iDirectVal,MPI_FLOAT,MPI_COMM_WORLD);
          break;
        case 'd':
          dval=va_arg(ap,double*);
          MPI_Unpack(Buffer,BufferSize,&BufferPos,
                     dval,1,MPI_DOUBLE,MPI_COMM_WORLD);
          break;
        case 'D':
          dval=va_arg(ap,double*);
          iDirectVal=va_arg(ap,int);
          MPI_Unpack(Buffer,BufferSize,&BufferPos,
                     dval,iDirectVal,MPI_DOUBLE,MPI_COMM_WORLD);
          break;
        case 's':
          sval=va_arg(ap,char*);
          MPI_Unpack(Buffer,BufferSize,&BufferPos,
                     sval,512,MPI_CHAR,MPI_COMM_WORLD);
          break;
        case 'P':
          procLim=va_arg(ap,int*);
          MPI_Unpack(Buffer,BufferSize,&BufferPos,
                     procLim,6,MPI_INT,MPI_COMM_WORLD);
          break;
        default:
          assert(FALSE,"getMessage--unknown flag %s in format string %s",
                 *p,format);
      }
    }
    va_end(ap);
  }
  
  return status.MPI_TAG;
}


//
//MPI Only; Advanced Block packing and unpacking routines, use the derived datatype
// constructors.
//
#define MESSAGE_SEND_BLOCKS 12347
void setBlocksType(MPI_Datatype* newtype,
                   int n,float** data,MPI_Datatype* types){
  int i;
  
  int *block_lengths;
  MPI_Aint *displacements;
  MPI_Aint startAddress,currAddress;
  
  //Allocate variables for the lengths (all 1) and displacements.
  assert((block_lengths=(int*)malloc(n*sizeof(int)))!=NULL,
         "setBlocksType--unable to allocate %i ints for block_lengths",
         n);
  assert((displacements=(MPI_Aint*)malloc(n*sizeof(MPI_Aint)))!=NULL,
         "setBlocksType--unable to allocate %i MPI_Aint for displacments",
         n);
  
  //Fill in the values (I could probably guess the displacents are 4*i).
  MPI_Address(data,&startAddress);
  for(i=0;i<n;i++){
    block_lengths[i]=1;
    MPI_Address(data[i],&currAddress);
    displacements[i]=currAddress-startAddress;
  }
  
  //Create the type.
  MPI_Type_struct(n,block_lengths,displacements,types,newtype);
  MPI_Type_commit(newtype);
}
void sendBlocks(int target,MPI_Datatype* group,float** data){
  MPI_Send(data,1,*group,target,
           MESSAGE_SEND_BLOCKS,MPI_COMM_WORLD);
}
void getBlocks(int source,MPI_Datatype* group,float** data){
  MPI_Status status;
  MPI_Recv(data,1,*group,source,
           MESSAGE_SEND_BLOCKS,MPI_COMM_WORLD,&status);
}

void immediateSendBlocks(MPI_Request* request,
                         int target,MPI_Datatype* group,float** data){
  MPI_Isend(data,1,*group,target,
            MESSAGE_SEND_BLOCKS,MPI_COMM_WORLD,request);
}

void setBlockType(MPI_Datatype* newtype,int *offset,
                  int nx,int ny,int nz,
                  int xStart,int xStop,
                  int yStart,int yStop,
                  int zStart,int zStop){
  int i=xStop-xStart,j=yStop-yStart,k=zStop-zStart;
  int nxy=nx*ny;
  *offset=xStart+yStart*nx+zStart*nxy;
  if(k!=nz){
    //This is the best case, since we are using all the xs and ys we can
    // send a contingous block collapse the
    // x and y. Then the block is defined by a single MPI_Type_vector.
    MPI_Type_contiguous(k*nx*ny,MPI_FLOAT,newtype);
  }else if(j!=ny){
    //This is the next best case, NZ entries each i*j long.
    MPI_Type_vector(k,i*j,nxy,MPI_FLOAT,newtype);
  }else if(i!=nx){
    //Worst case, each entry is i long and there are ny*nz of them.
    MPI_Type_vector(ny*nz,i,nx,MPI_FLOAT,newtype);
  }else{
    assert(FALSE,
           "setBlockType failure: new version requires at least 2 indicies to fill rnage");
  }
  
  //Commit the new type
  MPI_Type_commit(newtype);
}
int packBlock(int xStart,int xStop,
              int yStart,int yStop,
              int zStart,int zStop,
              int n,int nx,int nxy,
              float** data){
  int i=xStop-xStart,j,k;
  int ii;
  
  for(ii=0;ii<n;ii++){
    float* xx=data[ii];
    for(k=zStart;k<zStop;k++){
      for(j=yStart;j<yStop;j++){
        MPI_Pack(&xx[xStart+j*nx+k*nxy],i,MPI_FLOAT,
                 Buffer,BufferSize,&BufferPos,MPI_COMM_WORLD);
      }
    }
  }
  return i*j*k;
}
int unpackBlock(int xStart,int xStop,
                int yStart,int yStop,
                int zStart,int zStop,
                int n,int nx,int nxy,
                float** data){
  int i=xStop-xStart,j,k;
  int ii;
  
  for(ii=0;ii<n;ii++){
    float* xx=data[ii];
    for(k=zStart;k<zStop;k++){
      for(j=yStart;j<yStop;j++){
        MPI_Unpack(Buffer,BufferSize,&BufferPos,
                   &xx[xStart+j*nx+k*nxy],i,MPI_FLOAT,MPI_COMM_WORLD);
      }
    }
  }
  
  return i*j*k;
}

//error handling routines
// inform the master there is an error
// or for the master kill the rest of the slaves that do not have and error
void checkForError(void){
  //not yet implemented
}
static void sendError(int message){
  if(message==SIGSEGV)
    fprintf(stderr,
            "Seg Fault--");
  fprintf(stderr,
          "Process Terminating\n");
  sendMessage(Parent,MESSAGE_FAIL,NULL);
  MPI_Finalize();
}

//Subroutines to start and stop he slave processes
// if val is non-zero this exits the program
void killProcesses(int val){
  stopProcesses(-1,val);
}
void stopProcesses(int parent,int val){
  if(parent<0){
    int i;
    //kill slave processes
    initSend();
    sendMessage(AllProcesses,MESSAGE_EXIT,NULL);
    
    for(i=0;i<NumProcs;i++){
      getMessage(Tids[i],MESSAGE_EXIT,NULL);
    }
    
    sendMessage(AllProcesses,MESSAGE_FINAL_EXIT,NULL);
    fprintf(stderr,"Killed %i MPI processes\n",
            NumProcs);
  }
  
  MPI_Finalize();
  free(Tids);
  if(val>=0)
    exit(val);
}

int registerProcess(int* parent,
                    int* argc,char ***argv){
  int rank;
  MPI_Init(argc,argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  
  *parent=0;
  _Rank_=rank;
  Parent=*parent;
  
  if(rank==0){ //this is the master process
    //set the behavior for an error
    //signal(SIGINT,killProcesses);
    //signal(SIGSEGV,killProcesses);
    setErrorFunc(killProcesses);
  }else{
    //set the behavior for an error
    //signal(SIGINT,sendError);
    //signal(SIGSEGV,sendError);
    setErrorFunc(sendError);
    
    //Set the Tids array to have one item which is mytid.
    Tids=(int*)malloc(sizeof(int));
    *Tids=rank;
  }
  
  return rank;
}

int startProcesses(char* executableName,int numProcs,
                   int masterSlave){
  int i;
  int size;
  NumProcs=numProcs-(masterSlave?0:1);
  
  /*Initialize processes that do the work of updating grids*/
  Tids=(int*)malloc(NumProcs*sizeof(int));
  
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if((size-1)!=NumProcs){
    printf("Unable to start %i tasks, running %i procs\n",
           NumProcs,size);
    killProcesses(TRUE);
  }
  /*Set the Tids*/
  for(i=1;i<size;i++){
    Tids[i-1]=i;
  }
  
  tEprintf(Verbose,"Initialized %i processes into MPI\n",
           NumProcs);
  return NumProcs;
}

