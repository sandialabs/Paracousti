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
 *  message_passing.h
 *
 *
 *  Declares functions and macros having to do with general message passing.  Many
 *  MPI-specific functions are declared in mpi_procs.c.
 *
 *  Declares the following functions:
 *  initSend
 *  setMessageBuffer
 *  currMessageBufferSize
 *  sendMessage
 *  iSendMessage
 *  getMessage
 *  packMessage
 *  unpackMessage
 *  doBarrier
 *  setBlockType
 *  setBlocksType
 *  sendBlocks
 *  getBlocks
 *  immediateSendBlocks
 *  packBlock
 *  unpackBlock
 *  startProcesses
 *  registerProcess
 *  stopProcesses
 *  checkForError
 *  messageRank
 *
 *  Defines the following macros:
 *  MESSAGE_EXIT
 *  MESSAGE_FINAL_EXIT
 *  MESSAGE_INITIALIZE
 *  MESSAGE_DEBUG
 *  MESSAGE_FAIL
 *  MESSAGE_GENERAL
 *  MESSAGE_BARRIER
 *  MESSAGE_TEST
 *  MPIPP_H
 *
 */
/*These are subroutines used for message passing.

    Note that the definition of these routines is passing
    architecture independent. There are c implementations of these
    subroutines for MPI in the file mpi_procs.c and for PVM in the
    file pvm_procs.c. There are limited number of subroutines that can
    only be implemented in MPI.
*/

#ifndef _message_passing_h_
#define _message_passing_h_

/*There are some global variables*/
extern const int AnyMessageTag;
extern const int AllProcesses;
extern int NumProcs;
extern int* Tids;
extern int Parent;

/*Standard Message tags*/
#define MESSAGE_EXIT          1
#define MESSAGE_FINAL_EXIT    2
#define MESSAGE_INITIALIZE    3
#define MESSAGE_DEBUG         4
#define MESSAGE_FAIL          5
#define MESSAGE_GENERAL       6
#define MESSAGE_BARRIER       9

//Instead of making new messages for things that are not done very often this
// can be used by the application for a wide variety of possible messages. 
// For instance on the send end:
// char buffer[512]; //MPI implementation assumes 512 chars in a string
// sprintf(buffer,"APP_MESSAGE_1");
// packMessage("s",buffer);

// Matched by on the receive end:
// case MESSAGE_GENERAL
//   unpackMessage("s",string);
//   if(!strcmp(string,"APP_MESSAGE_1")){
//     .
//     .
//     .

/*Procedure declarations*/
#if defined(__cplusplus)
extern "C" {
#endif
  ///Call before starting a new message to clear the buffer.
  int initSend(void); 

  /// Call before first message to initialize the buffer.

  /// Also call this if the message buffer size has
  /// increased. Required for MPI only (null function for PVM).
  int setMessageBuffer(int size); 
  int currMessageBufferSize();

  /*! \brief Send a message with user defined parameters.

    Note initSend MUST BE CALLED FIRST.
    Format string for additional arguments
      -# i for int, f for float, s for string, P sends ProcLim (6 ints)
      -# F sends an array of floats requires two arguments a float* and int (length)
      -# d sends a double; D an array of Doubles
      -# I sends an array of ints requires two arguments a int* and int (length)
      -# B[n] sends  n blocks of floats from a 3d array requires 8+n
         arguments nx nxy xMin xMax yMin yMax zMin zMax data[0] ... 

    arguments are processed in order from left to right*/
  int sendMessage(int target,int tag,const char* format,...);

  ///Immediate send (MPI only PVM is always immediate).
  int iSendMessage(int target,int tag,const char* format,...);

  /// Receive a message same format as sendMessage.
  int getMessage(int source,int tag,const char* format,...);

  /// Pack a message but dont send.
  int packMessage(const char* format,...);

  // Unpack a message but dont receive.
  int unpackMessage(const char* format,...);

  // Implement a barrier (just call MPI_Barrier for mpi, more complex for pvm).
  void doBarrier();

#if USE_MPI_SEND
#define MPIPP_H 
  //This is a KLUDGE, under RH7 the mpi++ routines fail to compile, this keeps them
  // from being included.
#include "mpi.h"
  //
  //MPI Only; Advanced Block packing and unpacking routines, use the derived 
  //datatype constructors.
  //
  void setBlockType(MPI_Datatype* newtype,int *offset,
		    int nx,int ny,int nz,
		    int xStart,int xStop,
		    int yStart,int yStop,
		    int zStart,int zStop);

  void setBlocksType(MPI_Datatype* newtype,
		     int n,float** data,MPI_Datatype* types);

  //These are the blocking send and receives that I have used pre 2/14/03
  void sendBlocks(int target,MPI_Datatype* group,float** data);
  void getBlocks(int source,MPI_Datatype* group,float** data);

  //These are immediate versions, need to test for completion with 
  // MPI_Wait or MPI_test before proceding. The return values are the request
  // ids.
  //Checks are:
  // MPI_Wait(request,status);
  // MPI_Test(request,flag,status)
  void immediateSendBlocks(MPI_Request* request,
			   int target,MPI_Datatype* group,float** data);
#endif

  /* Pack and unpack n blocks of data*/  
  int packBlock(int xStart,int xStop, 
		int yStart,int yStop,
		int zStart,int zStop,
		int n,int nx,int nxy,
		float** data);

  //and the analogous unpack routines
  int unpackBlock(int xStart,int xStop,
		  int yStart,int yStop,
		  int zStart,int zStop,
		  int n,int nx,int nxy,
		  float** data);

  /*Start the slaves; parameters come from the global variables 
    // NXProc,NX,Dx,MinX, etc*/
  int startProcesses(char* executableName,int numProcs,
		     int masterSlave);
  /*if you are a slave; check in with the message passing protocol*/
  int registerProcess(int* parent,
		      int* argc,char ***argv);

  /*exit message passing, if called by the master, kill the slaves too*/
  void stopProcesses(int parent,int val);

  /*check for an error in one of the slaves; if it occurs
    // shutdown the rest of the slaves and exit*/
  void checkForError(void);

  /*Return a unique process id*/
  int messageRank();
#if defined(__cplusplus)
}
#endif

#endif /*#ifndef _message_passing_h_*/
