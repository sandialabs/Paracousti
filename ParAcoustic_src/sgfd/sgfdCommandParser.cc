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
 *  sgfdCommandParser.cc
 *
 *
 *  Defines class functions for the sgfdCommandParser class.  This class
 *  forms the base class for all the finite difference codes and
 *  contains command line flags/options that are common to all of these
 *  codes.
 *
 *  Defines the following class functions:
 *  sgfdCommandParserBase::sendBCMessage
 *  sgfdCommandParserBase::processNextArgument
 *  sgfdCommandParserBase::checkProcessFlag
 *  sgfdCommandParserBase::addExtraReceiverFlags
 *  sgfdCommandParserBase::addExtraSourceFlags
 *  sgfdCommandParserBase::setDefaultValues
 *
 */

#include "sgfdCommandParser.hh"
#include "sgfd.hh"
#include "message_passing.h"

//These are the generic versions of receivers, sources, and extra output.
#include "sgfdSources.hh"
#include "movingSgfdSources.hh"

#include "sgfdReceiverNetwork.hh"

#ifdef SGFD_INVERSION
#include "inversionReceivers.hh"
#endif


///Here is a virtual method to pack the message for the boundary condition
/// initialization.
int sgfdCommandParserBase::sendBCMessage(){
  tEprintf(Verbose,"Setting boundary conditions\n");
  initSend();

  //set any special boundry conditions
  // free surface and TDBC
  packMessage("i",
  _surfaceBCMode);  
  //sponge
  packMessage("IF",
  _spongeBCNodes,6,_spongeBCValue,6);

  sendMessage(AllProcesses,MESSAGE_SET_BC,NULL);
  for(int i=0;i<NumProcs;i++)
    getMessage(AllProcesses,MESSAGE_SET_BC,NULL);

  return NumProcs;
}

///Here is the required virtual function to process an argument.
int sgfdCommandParserBase::processNextArgument(int argc,char* argv[],int& i){
  if(!*_modelName){
    strcpy(_modelName,argv[i]);
    //run parameters are now read immediatly
    readCDFRunParams(_modelName,modelDef(),TRUE);
    return TRUE;
  }
  return FALSE;
}

///Here is the required virtual function to process a flag.
int sgfdCommandParserBase::checkProcessFlag(int argc,char* argv[],int& i){
  switch (argv[i][1]){
  case 'a':
    //Ignore this flag, it is used by totalview.
    tEprintf(Verbose,"Ignoring totalview -a flag\n");
    break;

  case 'v':
    if(strlen(argv[i])>2 && argv[i][2]=='r'){
_printSlaveModelRead= !_printSlaveModelRead;
tEprintf(Verbose,"Slaves %sprinting model read information\n",
   _printSlaveModelRead?"":"not ");
    }else{
Verbose=!Verbose;
if(strlen(argv[i])>2 && isdigit(argv[i][2]))
  Verbose=atoi(argv[i]+2);
tEprintf(Verbose,"Verbose flag set to %i\n",Verbose);
    }
    break;
  
  case 'j':
    assert(argc>i+1,_usageString,i,argv[i],argv[0]);
    _iterationsPerStep=atoi(argv[++i]);
    tEprintf(Verbose,"Running maximum of %i iterations between output steps\n",
       _iterationsPerStep);
    break;

  case 'p':
    assert(argc>i+3,_usageString,i,argv[i],argv[0]);
    _nxProc=atoi(argv[++i]);
    _nyProc=atoi(argv[++i]);
    _nzProc=atoi(argv[++i]);
    tEprintf(Verbose,"Using %ix%ix%i processes\n",_nxProc,_nyProc,_nzProc);
    break;
  case 'h':
      if(strlen(argv[i])>2) {
        if(argv[i][2]=='c') {
          assert(argc>i+1,_usageString,i,argv[i],argv[0]);
          _fdCoeffOrder = atoi(argv[++i]);
          assert(_fdCoeffOrder<=2*MAX_NUM_SPACE_COEFFS && _fdCoeffOrder%2==0,"Error: FD coefficient order must be <= %d and even",2*MAX_NUM_SPACE_COEFFS);
          _fdCoeffOrder/=2;
          assert(argc>i+_fdCoeffOrder,_usageString,i,argv[i],argv[0]);
          for(int ii=0;ii<_fdCoeffOrder;++ii)
            _fdCoeffs[ii] = atof(argv[++i]);
          tEprintf(Verbose,"Set %d FD coefficients to user specified values\n",_fdCoeffOrder);
        } else {
          assert(FALSE,_usageString,i,argv[i],argv[0]);
        }
      } else {
        assert(argc>i+1,_usageString,i,argv[i],argv[0]);
        _holbergPcent=atof(argv[++i]);
        if(_holbergPcent>10.0) _holbergPcent/= 100.0;
        tEprintf(Verbose,
           "Differentiators set using Holberg formula, %%%.2f level\n",
           100.0*_holbergPcent);
      }
    break;
  case 's':
    assert(argc>i+3,_usageString,i,argv[i],argv[0]);
    _modelDef.scalarSpeed=atof(argv[++i]);
    _modelDef.scalarVel=atof(argv[++i]);
    _modelDef.scalarDen=atof(argv[++i]);
    _modelDef.scalarStress=
_modelDef.scalarDen*_modelDef.scalarSpeed*_modelDef.scalarVel;
    tEprintf(Verbose,
       "Set scalar values: speed %.3g; vel %.3g; den %.3g => stress %.3g\n",
       _modelDef.scalarSpeed,_modelDef.scalarVel,
       _modelDef.scalarDen,_modelDef.scalarStress);
    break;
  case 't':
    assert(argc>i+1,_usageString,i,argv[i],argv[0]);
    _traceOutputIterations=atoi(argv[++i]);
    tEprintf(Verbose,"Trace output every %i iterations\n",
       _traceOutputIterations);
    assert(_traceOutputIterations>0,"Illegal value must be >0");
    break;

  case 'i':
    //Initial condition -i
    if(strlen(argv[i])<3){
assert(argc>i+1,_usageString,i,argv[i],argv[0]);
strcpy(_initialConditionName,argv[++i]);
tEprintf(Verbose,"Using initial conditions \"%s\"\n",
   _initialConditionName);
    }else{
return FALSE;
    }
    break;

  case 'C':
    assert(argc>2+i,_usageString,i,argv[i],argv[0]);
    setVectorValues(argc,argv,i,_doCheckpoint,_checkpointIndicies);
    strcpy(_checkpointDir,argv[++i]);

    tEprintf(Verbose,
       "Checkpoints at %i points to\n\t%s\n\t",
       _doCheckpoint,_checkpointDir);
    for(int ii=0;ii<_doCheckpoint;ii++)
tEprintf(Verbose,"%i%s",_checkpointIndicies[ii],((ii+1)%10)?" ":"\n\t");
    tEprintf(Verbose&&(_doCheckpoint%10),"\n");
    break;

  case 'T':
    assert(argc>i+1,_usageString,i,argv[i],argv[0]);
    assert(*_modelName,
     "-T option must only be called after model is set");
    {
float tmin,dt;
int nt;
setVectorValues(argc,argv,i,tmin,dt,nt);
assert((_newT=(float*)malloc(nt*sizeof(float)))!=NULL,
       "Unable to allocate %i floats for new t vector",
       nt);
for(int iii=0;iii<nt;++iii) _newT[iii]=tmin+dt*iii;
_modelDef.NT=nt;
_modelDef.minT=tmin;
_modelDef.dt=dt;
tEprintf(Verbose,
   "Set new tvector[%i]: [%.3f:%.5g:%.3f]\n",
   nt,tmin,dt,tmin+dt*nt);
    }
    break;

  case 'R':
    if(strlen(argv[i])>2 && (argv[i][2]=='t' || argv[i][2]=='w')){
assert(argc>i+1,_usageString,i,argv[i],argv[0]);
_traceOutputIterations=atoi(argv[++i]);
tEprintf(Verbose,"Trace output every %i iterations\n",
   _traceOutputIterations);
assert(_traceOutputIterations>0,"Illegal value must be >0");
    }else{
if(!_extraReceivers){
#ifdef SGFD_INVERSION
  _extraReceivers=new masterInversionReceiverNetwork();
#else
  _extraReceivers=new masterReceiverNetwork();
#endif
}
_extraReceivers->addReceivers(i,2,argc,argv,modelDef());
    }
    break;
  case 'S':
    if(!_extraSources)
_extraSources=new movingSourceNetwork();
    _extraSources->addSources(i,2,argc,argv,modelDef());
    break;

    //
    //boundry condition options -b
    case 'b':
if(strlen(argv[i])>2 && argv[i][2]!='S'){
  return FALSE;
}else{
  if(strlen(argv[i])<4){
    //This is my most common usage, set all flanks to use the same values.
    assert(argc>i+2,_usageString,i,argv[i],argv[0]);
    int nodes=atoi(argv[++i]);
    float value=atof(argv[++i]);
    if(value>=1.0) value/=100.0; //value must be in %
    assert(ISMID(0.0,value,1.0),
     "Sponge Boundry Conditions (%i,%f) value out of range",
     nodes,value);
    tEprintf(!ISMID(0.50,value,0.9995),
       "Sponge Boundry Conditions--recomended values in range 0.50-0.999\n");
    
    _spongeBCNodes[0]=_spongeBCNodes[1]=_spongeBCNodes[2]=nodes;
    _spongeBCNodes[3]=_spongeBCNodes[4]=_spongeBCNodes[5]=nodes;
    _spongeBCValue[0]=_spongeBCValue[1]=_spongeBCValue[2]=value;
    _spongeBCValue[3]=_spongeBCValue[4]=_spongeBCValue[5]=value;

    tEprintf(Verbose,
       "%s spongy absorbing boundry (%i nodes->%.2f%%)\n",
       nodes?"Using":"Not using",
       nodes,100.0*value);
  }else if(argv[i][3]=='3'){
    //Here is the intermediate version to set individual values for the 3
    // axes. NOTE that surface BC's will override any z-min value set
    // here and force it to 0 nodes=>value 1.
    assert(argc>i+6,_usageString,i,argv[i],argv[0]);

    //X-axis.
    int nodes=atoi(argv[++i]);
    float value=atof(argv[++i]);
    if(value>=1.0) value/=100.0; //value must be in %
    assert(ISMID(0.0,value,1.0),
     "Sponge Boundry Conditions (%i,%f) value out of range",
     nodes,value);
    tEprintf(!ISMID(0.50,value,0.9995),
       "Sponge Boundry Conditions--recomended values in range 0.50-0.999\n");
    _spongeBCNodes[0]=_spongeBCNodes[1]=nodes;
    _spongeBCValue[0]=_spongeBCValue[1]=value;
    tEprintf(Verbose,
       "%s x-axis spongy absorbing boundry (%i nodes->%.2f%%)\n",
       nodes?"Using":"Not using",
       nodes,100.0*value);

    //Y-axis.
    nodes=atoi(argv[++i]);
    value=atof(argv[++i]);
    if(value>=1.0) value/=100.0; //value must be in %
    assert(ISMID(0.0,value,1.0),
     "Sponge Boundry Conditions (%i,%f) value out of range",
     nodes,value);
    tEprintf(!ISMID(0.50,value,0.9995),
       "Sponge Boundry Conditions--recomended values in range 0.50-0.999\n");
    _spongeBCNodes[2]=_spongeBCNodes[3]=nodes;
    _spongeBCValue[2]=_spongeBCValue[3]=value;
    tEprintf(Verbose,
       "\ty-axis spongy absorbing boundry (%i nodes->%.2f%%)\n",
       nodes,100.0*value);

    //Z-axis.
    nodes=atoi(argv[++i]);
    value=atof(argv[++i]);
    if(value>=1.0) value/=100.0; //value must be in %
    assert(ISMID(0.0,value,1.0),
     "Sponge Boundry Conditions (%i,%f) value out of range",
     nodes,value);
    tEprintf(!ISMID(0.50,value,0.9995),
       "Sponge Boundry Conditions--recomended values in range 0.50-0.999\n");
    _spongeBCNodes[4]=_spongeBCNodes[5]=nodes;
    _spongeBCValue[4]=_spongeBCValue[5]=value;
    tEprintf(Verbose,
       "\tz-axis spongy absorbing boundry (%i nodes->%.2f%%)\n",
       nodes,100.0*value);
  }else if(argv[i][3]=='6'){
    //Here is the longest version to set individual values for all six 
    // flanks. NOTE that surface BC's will override any z-min value set
    // here and force it to 0 nodes=>value 1.
    assert(argc>i+12,_usageString,i,argv[i],argv[0]);

    //X min
    _spongeBCNodes[0]=atoi(argv[++i]);
    _spongeBCValue[0]=atof(argv[++i]);
    if(_spongeBCValue[0]>=1.0) _spongeBCValue[0]/=100.0; //value must be in %
    assert(ISMID(0.0,_spongeBCValue[0],1.0),
     "Sponge Boundry Conditions (%i,%f) value out of range",
     _spongeBCNodes[0],_spongeBCValue[0]);
    tEprintf(!ISMID(0.50,_spongeBCValue[0],0.9995),
       "Sponge Boundry Conditions--recomended values in range 0.50-0.999\n");
    tEprintf(Verbose,
       "%s x-min spongy absorbing boundry (%i nodes->%.2f%%)\n",
       _spongeBCNodes[0]?"Using":"Not using",
       _spongeBCNodes[0],100.0*_spongeBCValue[0]);

    //X max
    _spongeBCNodes[1]=atoi(argv[++i]);
    _spongeBCValue[1]=atof(argv[++i]);
    if(_spongeBCValue[1]>=1.0) _spongeBCValue[1]/=100.0; //value must be in %
    assert(ISMID(0.0,_spongeBCValue[1],1.0),
     "Sponge Boundry Conditions (%i,%f) value out of range",
     _spongeBCNodes[1],_spongeBCValue[1]);
    tEprintf(!ISMID(0.50,_spongeBCValue[1],0.9995),
       "Sponge Boundry Conditions--recomended values in range 0.50-0.999\n");
    tEprintf(Verbose,
       "\tx-max spongy absorbing boundry (%i nodes->%.2f%%)\n",
       _spongeBCNodes[1],100.0*_spongeBCValue[1]);

    //Y min
    _spongeBCNodes[2]=atoi(argv[++i]);
    _spongeBCValue[2]=atof(argv[++i]);
    if(_spongeBCValue[2]>=1.0) _spongeBCValue[2]/=100.0; //value must be in %
    assert(ISMID(0.0,_spongeBCValue[2],1.0),
     "Sponge Boundry Conditions (%i,%f) value out of range",
     _spongeBCNodes[2],_spongeBCValue[2]);
    tEprintf(!ISMID(0.50,_spongeBCValue[2],0.9995),
       "Sponge Boundry Conditions--recomended values in range 0.50-0.999\n");
    tEprintf(Verbose,
       "\ty-min spongy absorbing boundry (%i nodes->%.2f%%)\n",
       _spongeBCNodes[2],100.0*_spongeBCValue[2]);

    //Y max
    _spongeBCNodes[3]=atoi(argv[++i]);
    _spongeBCValue[3]=atof(argv[++i]);
    if(_spongeBCValue[3]>=1.0) _spongeBCValue[3]/=100.0; //value must be in %
    assert(ISMID(0.0,_spongeBCValue[3],1.0),
     "Sponge Boundry Conditions (%i,%f) value out of range",
     _spongeBCNodes[3],_spongeBCValue[3]);
    tEprintf(!ISMID(0.50,_spongeBCValue[3],0.9995),
       "Sponge Boundry Conditions--recomended values in range 0.50-0.999\n");
    tEprintf(Verbose,
       "\ty-max spongy absorbing boundry (%i nodes->%.2f%%)\n",
       _spongeBCNodes[3],100.0*_spongeBCValue[3]);

    //Z min
    _spongeBCNodes[4]=atoi(argv[++i]);
    _spongeBCValue[4]=atof(argv[++i]);
    if(_spongeBCValue[4]>=1.0) _spongeBCValue[4]/=100.0; //value must be in %
    assert(ISMID(0.0,_spongeBCValue[4],1.0),
     "Sponge Boundry Conditions (%i,%f) value out of range",
     _spongeBCNodes[4],_spongeBCValue[4]);
    tEprintf(!ISMID(0.50,_spongeBCValue[4],0.9995),
       "Sponge Boundry Conditions--recomended values in range 0.50-0.999\n");
    tEprintf(Verbose,
       "\tz-min spongy absorbing boundry (%i nodes->%.2f%%)\n",
       _spongeBCNodes[4],100.0*_spongeBCValue[4]);

    //Z max
    _spongeBCNodes[5]=atoi(argv[++i]);
    _spongeBCValue[5]=atof(argv[++i]);
    if(_spongeBCValue[5]>=1.0) _spongeBCValue[5]/=100.0; //value must be in %
    assert(ISMID(0.0,_spongeBCValue[5],1.0),
     "Sponge Boundry Conditions (%i,%f) value out of range",
     _spongeBCNodes[5],_spongeBCValue[5]);
    tEprintf(!ISMID(0.50,_spongeBCValue[5],0.9995),
       "Sponge Boundry Conditions--recomended values in range 0.50-0.999\n");
    tEprintf(Verbose,
       "\tz-max spongy absorbing boundry (%i nodes->%.2f%%)\n",
       _spongeBCNodes[5],100.0*_spongeBCValue[5]);
  }else{
    return FALSE;
  }
}
break;
  default:
    return FALSE;
  }
  return TRUE;
}

///Add flag descriptions pertaining to receivers.
void sgfdCommandParserBase::addExtraReceiverFlags(){
  addFlagDescription(FLAG_OPTIONAL,"-R[1]","s f f f f",
         "Add a single receiver, type (Vx, Vy, Vz, Pressure, 3C, 4C), x y z location",
         1," and amplitude scalar.");
  addFlagDescription(FLAG_OPTIONAL,"-Ro","s",
         "Set receiver output filename.");
  addFlagDescription(FLAG_OPTIONAL,"-R(d|a)","",
         "Integrate/differentiate velocity receivers to displacement/acceleration.");
  addFlagDescription(FLAG_EXPERT,"-Rs","i",
         "Subsample receiver output by n times (NO ANTI-ALIAS FILTER).");
  addFlagDescription(FLAG_OPTIONAL,"-Rg","s f f i f f i f f i",
         "Add a grid of receivers, type, ix dx nx, iy dy ny, iz dz nz.");
  addFlagDescription(FLAG_OPTIONAL,"-Rf[3]","s s",
         "Add a receivers rrom a file, type, filename. File is 4 column x y z amp",
         1," unless 3 is present then just x y z and amps are all 1.");
  addFlagDescription(FLAG_OPTIONAL,"","","");
}
///Add flag descriptions pertaining to sources.
void sgfdCommandParserBase::addExtraSourceFlags(){
  addFlagDescription(FLAG_OPTIONAL,"-S","s",
         "Read source wavelet from file.");
  addFlagDescription(FLAG_OPTIONAL,"-S(R|H|D)[t]","i|f",
         "Create ramp, heavyside, or delta function wavelet with offset in samples",
         1," or time.");
  addFlagDescription(FLAG_OPTIONAL,"-Sb","f",
         "Set bandwith for dispersion check (1%).");
  addFlagDescription(FLAG_EXPERT,"-SN","f",
         "Time reverse the source waveform with offset.");
  addFlagDescription(FLAG_OPTIONAL,"-Sr","f [0|1|2]",
         "Build a Ricker wavelet with given peak frequency. Optionally normalize",
         1," amplitude by 0th, 1st, or 2nd derivative.");
  addFlagDescription(FLAG_EXPERT,"-ST","s",
         "Use time-dependent boundary condition source from given file.",
         1," Use mapTDBC to build TDBC format.");
  addFlagDescription(FLAG_OPTIONAL,"-Sf","f f f f f f",
         "Add a force source, x y z amp, theta phi (0 0 is z, 90 0 is x, 90 90 is y).");
  addFlagDescription(FLAG_OPTIONAL,"-Se","f f f f",
         "Add an explosive source, x y z amp.");
  addFlagDescription(FLAG_OPTIONAL,"-Sd","f f f f f f f",
         "Add a double-couple moment source, x y z amp, strike dip rake of fault",
         1," plane and slip vector.");
  addFlagDescription(FLAG_EXPERT,"-Sm","f f f f f f f f f f f f f",
         "Add an arbitary momentsource, x y z amp, 9 components of the moment tensor.");
  addFlagDescription(FLAG_OPTIONAL,"","","");
}

///Add flag descriptions and set default values for internal variables.
void sgfdCommandParserBase::setDefaultValues(){
  char tmpSpace[512];
  addFlagDescription(FLAG_REQUIRED,"ModelName","",
         "Name of the model, must be fully qualified since this is read in parallel.");
  addFlagDescription(FLAG_OPTIONAL,"-i","s",
         "Load initial conditions from file.");

  addFlagDescription(FLAG_OPTIONAL,"-a","",
         "Ignored (argument flag for totalview debugger).");

  addFlagDescription(FLAG_OPTIONAL,"-vr","",
         "Print model information as subdomains are read by slaves.");  
  addFlagDescription(FLAG_OPTIONAL,"-p","i i i",
         "Domain decomposition, n procs is nx*ny*nz+1 (1 1 1).");
  addFlagDescription(FLAG_OPTIONAL,"-t","i",
         "Write partial receiver output at this many iterations.");

  addFlagDescription(FLAG_OPTIONAL,"-T","vector specification",
         "Reset the time vector.");

  addFlagDescription(FLAG_OPTIONAL,"-b[S]","i f",
         "Set spongy BC to this many nodes with this final taper (20, 95%).");

  sprintf(tmpSpace,"Manually set FD coefficients with order i (<=%d and even) and coefficients given by order/2 floats starting at central value moving outwards.",2*MAX_NUM_SPACE_COEFFS);
  addFlagDescription(FLAG_OPTIONAL,"-hc","i f ... f",
                     tmpSpace);

  addFlagDescription(FLAG_OPTIONAL,"","","");
  addExtraReceiverFlags();
  addExtraSourceFlags();
 //Model not read yet
  strcpy(_modelName,"");

  //Set the default values for the non-dimensionalizing scalars here.
  _modelDef.scalarSpeed=1000.0;
  _modelDef.scalarVel=1.0e-6;
  _modelDef.scalarDen=1000.0;
  _modelDef.scalarStress=
    _modelDef.scalarDen*_modelDef.scalarSpeed*_modelDef.scalarVel;

  //New advanced options to redefine the time vector.
  _newT=NULL;

  //And properties of internal state variables that are replacing the globals.
  // parallel processes
  _holbergPcent=0.0; //-h flag
  _fdCoeffOrder=0; //-hc flag
  for(int ii=0;ii<MAX_NUM_SPACE_COEFFS;++ii) _fdCoeffs[ii] = 0.0; //-hc flag
  _nxProc=1;_nyProc=1;_nzProc=1; //-p flag

  //Flags to determine how the run looks.
  _traceOutputIterations=MAX_TIME_STEPS_SEND;//-t or -R(w|t)
  _iterationsPerStep=25; //-j.
  _printSlaveModelRead=FALSE; //-vr flag

  //Boundary conditions.
  // -bS option use spongy absorbing boundries
  _spongeBCNodes[0]=_spongeBCNodes[1]=_spongeBCNodes[2]=20;
  _spongeBCNodes[3]=_spongeBCNodes[4]=_spongeBCNodes[5]=20;
  _spongeBCValue[0]=_spongeBCValue[1]=_spongeBCValue[2]=0.95;
  _spongeBCValue[3]=_spongeBCValue[4]=_spongeBCValue[5]=0.95;

  _surfaceBCMode=FALSE; //-bF option for free surface; -bV for v_z free.

  //Options for initial conditions: -i
  _initialConditionName[0]='\0'; //-i

  //Optional checkpoint every %i iterations.
  _doCheckpoint=FALSE;
  strcpy(_checkpointDir,"");

  //May want to add extra receivers or sources that are not defined in the 
  // earthmodel file
  _extraReceivers=NULL; //add with -R[?] flags
  _extraSources=NULL;     //add with -S[?] flags
}
