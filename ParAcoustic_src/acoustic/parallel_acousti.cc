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
 *  parallel_acousti.cc
 *
 *
 *  This file contains definitions for classes parallelAcoustiCommandParser
 *  and acousticSlaveProcess.  The first of these is in charge of parsing
 *  the command line for paracousti-specific flags and options.  The second
 *  initializes slaves and handles some paracousti-specific messages.
 *
 *  Defines the following class functions:
 *  parallelAcoustiCommandParser::sendBCMessage
 *  parallelAcoustiCommandParser::checkProcessFlag
 *  parallelAcoustiCommandParser::allocateExtraOutput
 *  parallelAcoustiCommandParser::setDefaultValues
 *  acousticSlaveProcess::doInit
 *  acousticSlaveProcess::processMessage
 *  acousticSlaveProcess::setBoundaryConditions
 *
 */


#define USAGE "IMPROPER CALL\nUSAGE: %s [-v[int_digit]][-p n n n]"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nstdutil.hh"

#include "message_passing.h"
#include "io_procs.h"

#include "parallel_acousti.hh"

#include "acoustic_control.hh"

#include "acousticBoundary.hh"
#include "selector.hh"
#include "updateAcousticCPMLFull.h"
#include "updateAcousticMPMLFull.h"
#include "updateAcousticSpongeFull.h"
#include "updateAcousticAttenMPMLFull.h"
#include "sgfdSources.hh"
#include "sgfdReceiverNetwork.hh"
#include "fixed_acoustic.hh"

//sends boundary specific (CPML, sponge) and pressure-free surface information to the slaves
int parallelAcoustiCommandParser::sendBCMessage(){
  tEprintf(Verbose,"Setting elastic boundary conditions\n");
  initSend();
  
  //set any special boundry conditions
  // free surface and TDBC
  packMessage("iif",
              _surfaceBCMode,_usePML,_xpmlFac);
  
  //sponge
  packMessage("IFFF",
              _spongeBCNodes,6,_spongeBCValue,6,_alphaValues,6,_kValues,6);
  
  sendMessage(AllProcesses,MESSAGE_SET_BC,NULL);
  for(int i=0;i<NumProcs;i++)
    getMessage(AllProcesses,MESSAGE_SET_BC,NULL);
  
  return NumProcs;
}

//parse the command line for acoustic-specific flags
int parallelAcoustiCommandParser::checkProcessFlag(int argc,char* argv[],int& i){
  if(sgfdCommandParser::checkProcessFlag(argc,argv,i))
    return TRUE;
  
  switch (argv[i][1]){
    case 'i':
      //
      //Initial condition options -i[?]
      if(strlen(argv[i])<3){
        assert(argc>i+1,_usageString,i,argv[i],argv[0]);
        strcpy(_initialConditionName,argv[++i]);
        tEprintf(Verbose,"Using initial conditions \"%s\"\n",
                 _initialConditionName);
      }else{
        assert(FALSE,_usageString,i,argv[i],argv[0]);
      }
      break;
    case 'D':
      //Multiply the grid by the specified amount by interpolating.
      assert(argc>i+1,USAGE,argv[0]);
      _gridMultiplier=atoi(argv[++i]);
      tEprintf(Verbose,"Multipling grid by %i times\n",_gridMultiplier);
      break;
      
      //
      //boundry condition options -b
    case 'b':
      //Note that -b or -bS are stripped off by sgfdComandParser.
      switch(argv[i][2]){
          //acoustic PML
        case 'p':
        {_usePML = 1;
          int argLen = strlen(argv[i]);
          if(argLen>3 && argv[i][3]=='c') {
            _usePML = 3;
          } else if(argLen>3 && argv[i][3]=='m') {
            _usePML = 4;
          }
          if(argLen<4 || (_usePML>1 && argLen<5)){
            //This is my most common usage, set all flanks to use the same values.
            assert(argc>i+2+(_usePML==2 || _usePML==4)+2*(_usePML==3),_usageString,i,argv[i],argv[0]);
            int nodes=atoi(argv[++i]);
            assert(nodes>1,"PML requires at least 2 nodes");
            float value=atof(argv[++i]);
            value/=100.0; //value must be in %
            assert(ISMID(0.0,value,1.0),
                   "PML Boundry Conditions (%i,%f) value out of range",
                   nodes,value);
            
            _spongeBCNodes[0]=_spongeBCNodes[1]=_spongeBCNodes[2]=nodes;
            _spongeBCNodes[3]=_spongeBCNodes[4]=_spongeBCNodes[5]=nodes;
            _spongeBCValue[0]=_spongeBCValue[1]=_spongeBCValue[2]=value;
            _spongeBCValue[3]=_spongeBCValue[4]=_spongeBCValue[5]=value;
            
            if(_usePML==2) {
              _xpmlFac = atof(argv[++i]);
              assert(ISMID(0.0,_xpmlFac,1.0),
                     "PML cross-factor Boundry Conditions (%f) value out of range",
                     _xpmlFac);
              tEprintf(Verbose,
                       "%s xPML absorbing boundry (%i nodes->%.2f%%, cross-factor=%f)\n",
                       nodes?"Using":"Not using",
                       nodes,100.0*value,_xpmlFac);
            } else if(_usePML==3 || _usePML==4) {
              float avalue = 0., kvalue = 1.;
              avalue = atof(argv[++i]);
              assert(avalue>=0.,"cPML alpha value %f out of range (>=0)",avalue);
              kvalue = atof(argv[++i]);
              assert(kvalue>=1.,"cPML kappa value %f out of range (>=1)",kvalue);
              _alphaValues[0]=_alphaValues[1]=_alphaValues[2]=_alphaValues[3]=_alphaValues[4]=_alphaValues[5]=avalue;
              _kValues[0]=_kValues[1]=_kValues[2]=_kValues[3]=_kValues[4]=_kValues[5]=kvalue;
              tEprintf(Verbose,
                       "%s cPML absorbing boundry (%i nodes->%.2f%%, alpha=%f, kappa=%f)\n",
                       nodes?"Using":"Not using",
                       nodes,100.0*value,avalue,kvalue);
              if(_usePML==4) {
                _xpmlFac = atof(argv[++i]);
                assert(ISMID(0.0,_xpmlFac,1.0),
                       "MPML cross-factor Boundry Conditions (%f) value out of range",
                       _xpmlFac);
                tEprintf(Verbose,
                         "MPML absorbing boundry cross-factor=%f)\n",_xpmlFac);
              }
            } else {
              tEprintf(Verbose,
                       "%s PML absorbing boundry (%i nodes->%.2f%%)\n",
                       nodes?"Using":"Not using",
                       nodes,100.0*value);
              float avalue = 0., kvalue = 1.;
              _alphaValues[0]=_alphaValues[1]=_alphaValues[2]=
              _alphaValues[3]=_alphaValues[4]=_alphaValues[5]=avalue;
              _kValues[0]=_kValues[1]=_kValues[2]=_kValues[3]=_kValues[4]=_kValues[5]=kvalue;
            }
          }else if(argv[i][argLen-1]=='3'){
            //Here is the intermediate version to set individual values for the 3
            // axes. NOTE that surface BC's will override any z-min value set
            // here and force it to 0 nodes=>value 1.
            assert(argc>i+6+(_usePML==2 || _usePML==4)+6*(_usePML==3),_usageString,i,argv[i],argv[0]);
            
            //X-axis.
            int nodes=atoi(argv[++i]);
            assert(nodes>1,"PML requires at least 2 nodes");
            float value=atof(argv[++i]);
            value/=100.0; //value must be in %
            assert(ISMID(0.0,value,1.0),
                   "PML Boundry Conditions (%i,%f) value out of range",
                   nodes,value);
            _spongeBCNodes[0]=_spongeBCNodes[1]=nodes;
            _spongeBCValue[0]=_spongeBCValue[1]=value;
            if(_usePML==3 ||_usePML==4) {
              float avalue = 0., kvalue = 1.;
              avalue = atof(argv[++i]);
              assert(avalue>=0.,"cPML alpha value %f out of range (>=0)",avalue);
              kvalue = atof(argv[++i]);
              assert(kvalue>=1.,"cPML kappa value %f out of range (>=1)",kvalue);
              _alphaValues[0]=_alphaValues[1]=avalue;
              _kValues[0]=_kValues[1]=kvalue;
              tEprintf(Verbose,
                       "%s x-axis cPML absorbing boundry (%i nodes->%.2f%%, alpha=%f, kappa=%f)\n",
                       nodes?"Using":"Not using",
                       nodes,100.0*value,avalue,kvalue);
            } else {
              tEprintf(Verbose,
                       "%s x-axis PML absorbing boundry (%i nodes->%.2f%%)\n",
                       nodes?"Using":"Not using",
                       nodes,100.0*value);
              float avalue = 0., kvalue = 1.;
              _alphaValues[0]=_alphaValues[1]=avalue;
              _kValues[0]=_kValues[1]=kvalue;
            }
            
            //Y-axis.
            nodes=atoi(argv[++i]);
            assert(nodes>1,"PML requires at least 2 nodes");
            value=atof(argv[++i]);
            value/=100.0; //value must be in %
            assert(ISMID(0.0,value,1.0),
                   "PML Boundry Conditions (%i,%f) value out of range",
                   nodes,value);
            _spongeBCNodes[2]=_spongeBCNodes[3]=nodes;
            _spongeBCValue[2]=_spongeBCValue[3]=value;
            if(_usePML==3 || _usePML==4) {
              float avalue = 0., kvalue = 1.;
              avalue = atof(argv[++i]);
              assert(avalue>=0.,"cPML alpha value %f out of range (>=0)",avalue);
              kvalue = atof(argv[++i]);
              assert(kvalue>=1.,"cPML kappa value %f out of range (>=1)",kvalue);
              _alphaValues[2]=_alphaValues[3]=avalue;
              _kValues[2]=_kValues[3]=kvalue;
              tEprintf(Verbose,
                       "\ty-axis cPML absorbing boundry (%i nodes->%.2f%%, alpha=%f, kappa=%f)\n",
                       nodes,100.0*value,avalue,kvalue);
            } else {
              tEprintf(Verbose,
                       "\ty-axis PML absorbing boundry (%i nodes->%.2f%%)\n",
                       nodes,100.0*value);
              float avalue = 0., kvalue = 1.;
              _alphaValues[2]=_alphaValues[3]=avalue;
              _kValues[2]=_kValues[3]=kvalue;
            }
            
            //Z-axis.
            nodes=atoi(argv[++i]);
            assert(nodes>1,"PML requires at least 2 nodes");
            value=atof(argv[++i]);
            value/=100.0; //value must be in %
            assert(ISMID(0.0,value,1.0),
                   "PML Boundry Conditions (%i,%f) value out of range",
                   nodes,value);
            _spongeBCNodes[4]=_spongeBCNodes[5]=nodes;
            _spongeBCValue[4]=_spongeBCValue[5]=value;
            if(_usePML==3 || _usePML==4) {
              float avalue = 0., kvalue = 1.;
              avalue = atof(argv[++i]);
              assert(avalue>=0.,"cPML alpha value %f out of range (>=0)",avalue);
              kvalue = atof(argv[++i]);
              assert(kvalue>=1.,"cPML kappa value %f out of range (>=1)",kvalue);
              _alphaValues[4]=_alphaValues[5]=avalue;
              _kValues[4]=_kValues[5]=kvalue;
              tEprintf(Verbose,
                       "\tz-axis cPML absorbing boundry (%i nodes->%.2f%%, alpha=%f, kappa=%f)\n",
                       nodes,100.0*value,avalue,kvalue);
            } else {
              tEprintf(Verbose,
                       "\tz-axis PML absorbing boundry (%i nodes->%.2f%%)\n",
                       nodes,100.0*value);
              float avalue = 0., kvalue = 1.;
              _alphaValues[4]=_alphaValues[5]=avalue;
              _kValues[4]=_kValues[5]=kvalue;
            }
            if(_usePML==2 || _usePML==4) {
              _xpmlFac = atof(argv[++i]);
              assert(ISMID(0.0,_xpmlFac,1.0),
                     "PML cross-factor Boundry Conditions (%f) value out of range",
                     _xpmlFac);
              tEprintf(Verbose,"\t cross-factor = %f\n",_xpmlFac);
            }
          }else if(argv[i][argLen-1]=='6'){
            //Here is the longest version to set individual values for all six
            // flanks. NOTE that surface BC's will override any z-min value set
            // here and force it to 0 nodes=>value 1.
            assert(argc>i+12+(_usePML==2 || _usePML==4)+12*(_usePML==3),_usageString,i,argv[i],argv[0]);
            
            //X min
            _spongeBCNodes[0]=atoi(argv[++i]);
            assert(_spongeBCNodes[0]>1,"PML requires at least 2 nodes");
            _spongeBCValue[0]=atof(argv[++i]);
            _spongeBCValue[0]/=100.0; //value must be in %
            assert(ISMID(0.0,_spongeBCValue[0],1.0),
                   "PML Boundry Conditions (%i,%f) value out of range",
                   _spongeBCNodes[0],_spongeBCValue[0]);
            if(_usePML==3 || _usePML==4) {
              float avalue = 0., kvalue = 1.;
              avalue = atof(argv[++i]);
              assert(avalue>=0.,"cPML alpha value %f out of range (>=0)",avalue);
              kvalue = atof(argv[++i]);
              assert(kvalue>=1.,"cPML kappa value %f out of range (>=1)",kvalue);
              _alphaValues[0]=avalue;
              _kValues[0]=kvalue;
              tEprintf(Verbose,
                       "%s x-min cPML absorbing boundry (%i nodes->%.2f%%, alpha=%f, kappa=%f)\n",
                       _spongeBCNodes[0]?"Using":"Not using",
                       _spongeBCNodes[0],100.0*_spongeBCValue[0],avalue,kvalue);
            } else {
              tEprintf(Verbose,
                       "%s x-min PML absorbing boundry (%i nodes->%.2f%%)\n",
                       _spongeBCNodes[0]?"Using":"Not using",
                       _spongeBCNodes[0],100.0*_spongeBCValue[0]);
              _alphaValues[0] = 0.;
              _kValues[0] = 1.;
           }
            
            //X max
            _spongeBCNodes[1]=atoi(argv[++i]);
            assert(_spongeBCNodes[1]>1,"PML requires at least 2 nodes");
            _spongeBCValue[1]=atof(argv[++i]);
            _spongeBCValue[1]/=100.0; //value must be in %
            assert(ISMID(0.0,_spongeBCValue[1],1.0),
                   "PML Boundry Conditions (%i,%f) value out of range",
                   _spongeBCNodes[1],_spongeBCValue[1]);
            if(_usePML==3 || _usePML==4) {
              float avalue = 0., kvalue = 1.;
              avalue = atof(argv[++i]);
              assert(avalue>=0.,"cPML alpha value %f out of range (>=0)",avalue);
              kvalue = atof(argv[++i]);
              assert(kvalue>=1.,"cPML kappa value %f out of range (>=1)",kvalue);
              _alphaValues[1]=avalue;
              _kValues[1]=kvalue;
              tEprintf(Verbose,
                       "\tx-max cPML absorbing boundry (%i nodes->%.2f%%, alpha=%f, kappa=%f)\n",
                       _spongeBCNodes[1],100.0*_spongeBCValue[1],avalue,kvalue);
            } else {
              tEprintf(Verbose,
                       "\tx-max PML absorbing boundry (%i nodes->%.2f%%)\n",
                       _spongeBCNodes[1],100.0*_spongeBCValue[1]);
              _alphaValues[1] = 0.;
              _kValues[1] = 1.;
            }
            
            //Y min
            _spongeBCNodes[2]=atoi(argv[++i]);
            assert(_spongeBCNodes[2]>1,"PML requires at least 2 nodes");
            _spongeBCValue[2]=atof(argv[++i]);
            _spongeBCValue[2]/=100.0; //value must be in %
            assert(ISMID(0.0,_spongeBCValue[2],1.0),
                   "PML Boundry Conditions (%i,%f) value out of range",
                   _spongeBCNodes[2],_spongeBCValue[2]);
            if(_usePML==3 || _usePML==4) {
              float avalue = 0., kvalue = 1.;
              avalue = atof(argv[++i]);
              assert(avalue>=0.,"cPML alpha value %f out of range (>=0)",avalue);
              kvalue = atof(argv[++i]);
              assert(kvalue>=1.,"cPML kappa value %f out of range (>=1)",kvalue);
              _alphaValues[2]=avalue;
              _kValues[2]=kvalue;
              tEprintf(Verbose,
                       "\ty-min cPML absorbing boundry (%i nodes->%.2f%%, alpha=%f, kappa=%f)\n",
                       _spongeBCNodes[2],100.0*_spongeBCValue[2],avalue,kvalue);
            } else {
              tEprintf(Verbose,
                       "\ty-min PML absorbing boundry (%i nodes->%.2f%%)\n",
                       _spongeBCNodes[2],100.0*_spongeBCValue[2]);
              _alphaValues[2] = 0.;
              _kValues[2] = 1.;
           }
            
            //Y max
            _spongeBCNodes[3]=atoi(argv[++i]);
            assert(_spongeBCNodes[3]>1,"PML requires at least 2 nodes");
            _spongeBCValue[3]=atof(argv[++i]);
            _spongeBCValue[3]/=100.0; //value must be in %
            assert(ISMID(0.0,_spongeBCValue[3],1.0),
                   "PML Boundry Conditions (%i,%f) value out of range",
                   _spongeBCNodes[3],_spongeBCValue[3]);
            if(_usePML==3 || _usePML==4) {
              float avalue = 0., kvalue = 1.;
              avalue = atof(argv[++i]);
              assert(avalue>=0.,"cPML alpha value %f out of range (>=0)",avalue);
              kvalue = atof(argv[++i]);
              assert(kvalue>=1.,"cPML kappa value %f out of range (>=1)",kvalue);
              _alphaValues[3]=avalue;
              _kValues[3]=kvalue;
              tEprintf(Verbose,
                       "\ty-max cPML absorbing boundry (%i nodes->%.2f%%, alpha=%f, kappa=%f)\n",
                       _spongeBCNodes[3],100.0*_spongeBCValue[3],avalue,kvalue);
            } else {
              tEprintf(Verbose,
                       "\ty-max PML absorbing boundry (%i nodes->%.2f%%)\n",
                       _spongeBCNodes[3],100.0*_spongeBCValue[3]);
              _alphaValues[3] = 0.;
              _kValues[3] = 1.;
            }
            
            //Z min
            _spongeBCNodes[4]=atoi(argv[++i]);
            assert(_spongeBCNodes[4]>1,"PML requires at least 2 nodes");
            _spongeBCValue[4]=atof(argv[++i]);
            _spongeBCValue[4]/=100.0; //value must be in %
            assert(ISMID(0.0,_spongeBCValue[4],1.0),
                   "PML Boundry Conditions (%i,%f) value out of range",
                   _spongeBCNodes[4],_spongeBCValue[4]);
            if(_usePML==3 || _usePML==4) {
              float avalue = 0., kvalue = 1.;
              avalue = atof(argv[++i]);
              assert(avalue>=0.,"cPML alpha value %f out of range (>=0)",avalue);
              kvalue = atof(argv[++i]);
              assert(kvalue>=1.,"cPML kappa value %f out of range (>=1)",kvalue);
              _alphaValues[4]=avalue;
              _kValues[4]=kvalue;
              tEprintf(Verbose,
                       "\tz-min cPML absorbing boundry (%i nodes->%.2f%%, alpha=%f, kappa=%f)\n",
                       _spongeBCNodes[4],100.0*_spongeBCValue[4],avalue,kvalue);
            } else {
              tEprintf(Verbose,
                       "\tz-min PML absorbing boundry (%i nodes->%.2f%%)\n",
                       _spongeBCNodes[4],100.0*_spongeBCValue[4]);
              _alphaValues[4] = 0.;
              _kValues[4] = 1.;
            }
            
            //Z max
            _spongeBCNodes[5]=atoi(argv[++i]);
            assert(_spongeBCNodes[5]>1,"PML requires at least 2 nodes");
            _spongeBCValue[5]=atof(argv[++i]);
            _spongeBCValue[5]/=100.0; //value must be in %
            assert(ISMID(0.0,_spongeBCValue[5],1.0),
                   "PML Boundry Conditions (%i,%f) value out of range",
                   _spongeBCNodes[5],_spongeBCValue[5]);
            if(_usePML==3 || _usePML==4) {
              float avalue = 0., kvalue = 1.;
              avalue = atof(argv[++i]);
              assert(avalue>=0.,"cPML alpha value %f out of range (>=0)",avalue);
              kvalue = atof(argv[++i]);
              assert(kvalue>=1.,"cPML kappa value %f out of range (>=1)",kvalue);
              _alphaValues[5]=avalue;
              _kValues[5]=kvalue;
              tEprintf(Verbose,
                       "\tz-max cPML absorbing boundry (%i nodes->%.2f%%, alpha=%f, kappa=%f)\n",
                       _spongeBCNodes[5],100.0*_spongeBCValue[5],avalue,kvalue);
            } else {
              tEprintf(Verbose,
                       "\tz-max PML absorbing boundry (%i nodes->%.2f%%)\n",
                       _spongeBCNodes[5],100.0*_spongeBCValue[5]);
              _alphaValues[5] = 0.;
              _kValues[5] = 1.;
           }
            if(_usePML==2 || _usePML==4) {
              _xpmlFac = atof(argv[++i]);
              assert(ISMID(0.0,_xpmlFac,1.0),
                     "PML cross-factor Boundry Conditions (%f) value out of range",
                     _xpmlFac);
              tEprintf(Verbose,"\t cross-factor = %f\n",_xpmlFac);
            }
          }else{
            return FALSE;
          }
          if(_usePML==4) {  //kappa must be one for MPML
            bool testFirstKappaWarn = true;
            for(int iii=0;iii<6;++iii) {
              if(_kValues[iii]!=1.) {
                tEprintf(Verbose && testFirstKappaWarn,"MPML must have kappa==1...resetting\n");
                testFirstKappaWarn = false;
                _kValues[iii] = 1.;
              }
            }
          }
          break;
        }
        case 'V':
        case 'W':
          _surfaceBCMode=SURFACE_VEL_FREE;
          tEprintf(Verbose,
                   "Velocity free surface at z=%.2f\n",
                   _modelDef.minZ+
                   ((_mediaType==FIXED_MEDIA)?1.0:2.0)*_modelDef.dz);
          break;
        case 'F':
        case 'P':
          _surfaceBCMode=SURFACE_PRESS_FREE;
          tEprintf(Verbose,
                   "Pressure free surface at z=%.2f\n",
                   _modelDef.minZ+
                   2.0*_modelDef.dz);
          break;
          
        default:
          assert(FALSE,USAGE,argv[0]);
      }
      break;
      
      //The -Q family of flags deals with anacoustic (attuatation using standard linear fluid mechs)
    case 'Q':
      if(strlen(argv[i])>2){
        _qSelectors->parseSelectorArguments(i,2,argc,argv);
      }else{
        assert(argc>i+3,USAGE,argv[0]);
        int numMechs=atoi(argv[++i]);
        float fInfVp=atof(argv[++i]);
        if(fInfVp>10) fInfVp/=100.0; //must be in %
        if(!_qSelectors){
          DEF_MODEL_SIZE(&_modelDef);
          DEF_MODEL_LIMITS(&_modelDef);
          _qSelectors=new selector(NX,NY,NZ,minX,dx,minY,dy,minZ,dz);
          _qSelectors->addField("C",NULL);
          _qSelectors->addField("Rho",NULL);
        }
        
        _qSelectors->addSelector();
        _qSelectors->addIntVar(numMechs);
        _qSelectors->addFloatVar(fInfVp);
        
        float* decayRates,*ampP;
        assert((decayRates=(float*)malloc(numMechs*sizeof(float)))!=NULL,
               "Unable to allocate %i decay rates for Q2Selector %i",
               numMechs,_qSelectors->size()+1);
        assert((ampP=(float*)malloc(numMechs*sizeof(float)))!=NULL,
               "Unable to allocate %i P amp factors for Q2Selector %i",
               numMechs,_qSelectors->size()+1);
        
        tEprintf(Verbose,"Selector #%i; %i mechanisms; C(inf) %.2f%%;\n",
                 _qSelectors->size(),numMechs,100.0*fInfVp);
        for(int iii=0;iii<numMechs;iii++){
          assert(argc>i+3,USAGE,argv[0]);
          decayRates[iii]=atof(argv[++i])*2*PI;
          ampP[iii]=atof(argv[++i]);
          tEprintf(Verbose,"\t%8.3f\t%8.3f\n",
                   decayRates[iii],ampP[iii]);
        }
        _qSelectors->addFloatPtrVar(numMechs,decayRates);
        _qSelectors->addFloatPtrVar(numMechs,ampP);
      }
      break;
      
    default:
      return FALSE;
  }
  return TRUE;
}
//use the acoustic-specific output group
extraOutputGroup* parallelAcoustiCommandParser::allocateExtraOutput(){
  return new acousticOutputGroup(modelDef());
}

//setup defaults for acoustic-specific flags
void parallelAcoustiCommandParser::setDefaultValues(){
  sgfdCommandParser::setDefaultValues();
  
  addFlagDescription(FLAG_OPTIONAL,"-bp","","PML.");
  addFlagDescription(FLAG_OPTIONAL,"-b(FP)","",
                     "Free surface.");
  addFlagDescription(FLAG_OPTIONAL,"-b(WV)","",
                     "Zero Vz surface (implies 0 dP/dz).");
  
  addFlagDescription(FLAG_EXPERT,"-bZ","i fff",
                     "Use Aldridge implementation Z-K impedence BC; thickness, relaxation, omega, Q.");
  addFlagDescription(FLAG_EXPERT,"-bK","i fff",
                     "Use Wilson Z-K impedence BC; thickness, Flow Resistivity, Tortuosity, Porosity.");
  addFlagDescription(FLAG_EXPERT,"-bK3","i s",
                     "Use Wilson 3D Z-K impedence BC; thickness, filename.");
  
  addFlagDescription(FLAG_OPTIONAL,"","","");
  
  
  addFlagDescription(FLAG_EXPERT,"-MF","",
                     "Use FORTRAN interior update subroutines.");
  
  
  addFlagDescription(FLAG_OPTIONAL,"","","");
  addExtraOutputFlags();
  addExtraReceiverFlags();
  addExtraSourceFlags();
  
  //Options for initial conditions: -i[?]
  _initialConditionName[0]='\0'; //-i
  
  _usePML = 0; //-bp
  _xpmlFac = 0.;  //-bpm
  _alphaValues[0]=_alphaValues[1]=_alphaValues[2]=_alphaValues[3]=_alphaValues[4]=_alphaValues[5]=0.;
  _kValues[0]=_kValues[1]=_kValues[2]=_kValues[3]=_kValues[4]=_kValues[5]=1.;
  
  //And properties of internal state variables that are replacing the globals.
  // parallel processes
  _doCowork=FALSE;
  _gridMultiplier=1;
  
  //And of course, flags having to do with media type.
  _mediaType=FIXED_MEDIA;
  
  //True attenuation is off by default
  _qSelectors = NULL;
}

//Initialization routine for the slave processes.  This sets up the
//slave model, dependent, boundary, receiver, source classes.  It
//then moves the pressures to time step zero by inserting any pressure
//source for that time step.
void acousticSlaveProcess::doInit(){
  _dependent=NULL;
  
  //Get the media and pivot type from the next init message.
  int mediaType,pivotType;
  getMessage(_parent,MESSAGE_INITIALIZE,"ii",
             &mediaType,&pivotType);
  
  //Get the intial message and initialize the model.
  getMessage(_parent,MESSAGE_INITIALIZE,NULL);

  //Need some extra pointers for specialized models.
  //slaveAcousticNoInterpModel* noInterpModel=NULL;
  //slaveAcousticPseudoAttenuateModel* pseudoAttenuateModel=NULL;
  slaveAcousticAttenModel* attenModel=NULL;
  switch(mediaType) {
    case FIXED_MEDIA:
      _model=new slaveAcousticModel();
      break;
    case FIXED_MEDIA_ATTEN:
      _model = attenModel = new slaveAcousticAttenModel();
      break;
    default:
      assert(FALSE,"[%d] Unknown media type %d\n",messageRank(),mediaType);
  }

  initSend();
  sendMessage(_parent,MESSAGE_INITIALIZE,"i ff ff ff ff ff",
              FIXED_MEDIA,
              _model->_vMin,_model->_vMax,
              0.0,0.0,
              0.0,0.0,0.0,0.0,0.0,0.0);
  
  //Build the dependent variable class.
  switch(mediaType){
    case FIXED_MEDIA:
      _dependent=new slaveAcousticDependent(_model);
      break;
    case FIXED_MEDIA_ATTEN:
      _dependent=new slaveAcousticAttenDependent(attenModel);
      break;
    default:
      assert(FALSE,
             "[%i] unknown media type %i, must be %i",
             messageRank(),
             mediaType,
             FIXED_MEDIA);
  }
  
  //Check for a change in the time vector.
  char buffer[512];
  getMessage(Parent,MESSAGE_GENERAL,"s",buffer);
  if(!strcmp(buffer,"MESSAGE_RESET_TIME_VECTOR")){
    _model->resetTimeVector();
  }
  
  //Define the major variables to simplify later operations
  DEF_MODEL_SIZE(_model->modelDef());
  //Now reallocate the message buffer to final size
  setMessageBuffer(6*4*sizeof(float)*MAX(NX,NY)*MAX(NY,NZ));
  
  //Forward problem will also need to go inside the switch,
  // before this will work.
  _dependent->setInitialConditions();
  
  //should be receiving additional messages from the parent containing
  // more setup information, these routines receive and process these messages
  setBoundaryConditions(mediaType);
  _dependent->setBoundaryConditions(_boundaries);
  
  //set up sources and receivers
  _dependent->setSources(_model->_rho,FALSE);
  _dependent->setReceivers(TRUE);
  if(_dependent->_checkpointID>0)
    _dependent->_receivers->readCheckpoint(_dependent->_checkpointIteration,
                                           _dependent->_checkpointDir,
                                           _dependent->_checkpointID,0);
  
  //Insert initial source and load initial seismogram.
  _dependent->_sources->
  applyPressure(_dependent->_model->modelDef(),0,
                _dependent->_pressure);
  //load initial seismogram for initial index.
  _dependent->loadSeismogram(0);
  
  //And set the _model and _dependent in the parent class.
  sgfdSlaveProcess::doInit(_model,_dependent);
}

//process some acoustic-specific messages.  These are rarely, if ever,
//used.
int acousticSlaveProcess::processMessage(int msgtag){
  if(!msgtag){
    //Get control message from the parent.
    msgtag=getMessage(_parent,AnyMessageTag,NULL);
  }
  
  //Parse and process the message
  char buffer[512];
  switch(msgtag){
#ifdef SGFD_INVERSION
    case MESSAGE_WRITE_INVERSION_RECEIVERS:
    {
      int iteration;
      char irDirRoot[512];
      unpackMessage("is",&iteration,irDirRoot);
      _dependent->_receivers->
      writeBinaryInversionTraces(_model->modelDef(),_model->parallelDef(),
                                 iteration,
                                 _model->
                                 uniqueFilename("inversionReceiverData",
                                                irDirRoot));
      initSend();
      sendMessage(Parent,MESSAGE_WRITE_INVERSION_RECEIVERS,
                  "i",_dependent->_receivers->inversionSize());
    }
      return TRUE;
#endif
    case MESSAGE_GENERAL:
      unpackMessage("s",buffer);
      if(!strcmp(buffer,"DO_CHECKPOINT") ||
	       !strcmp(buffer,"READ_CHECKPOINT")){
        int iteration,cpID;
        char cpDir[1024];
        unpackMessage("isi",&iteration,cpDir,&cpID);
        if(*buffer=='D'){
          _dependent->doCheckpoint(iteration,cpDir,cpID,0);
        }else{
          _dependent->readCheckpoint(iteration,cpDir,cpID,0);
        }
      }else{
        tEprintf(Verbose,"PID[%i]: unknown MESSAGE_GENERAL option %s\n",
                 _tid,buffer);
        initSend();
        sendMessage(_parent,MESSAGE_FAIL,"P",_model->procLim());
        return FALSE;
      }
      return TRUE;
  }
  
  return sgfdSlaveProcess::processMessage(msgtag);
}

//Set up the boundary conditions (CPML or sponge) and free-surface
//condition (if any).
void acousticSlaveProcess::setBoundaryConditions(int mediaType){
  //read the boundary condition message from the master
  int surfaceBCMode, usePML;
  int spongeNodes[6];
  float spongeValue[6];
  float alphaValues[6];
  float kValues[6];
  float xpmlFac;
  getMessage(Parent,MESSAGE_SET_BC,"iif IFFF",
             &surfaceBCMode,&usePML,&xpmlFac,
             spongeNodes,6,spongeValue,6,alphaValues,6,kValues,6);
  
  //first check the type of model
  if(mediaType==FIXED_MEDIA_ATTEN) {
    //allocate the structure for the boundarys
    if(usePML && !spongeNodes) usePML = 0;
    _boundaries=new acousticSgfdBoundaries(_model,_dependent,
                                           surfaceBCMode,usePML);
    _dependent->setBoundaryConditions(_boundaries);
    if(usePML) {
      if(surfaceBCMode==SURFACE_PRESS_FREE) {  //Zmin should have 2 nodes and no damping if explicit free surface is used
        assert(FALSE,"Not implemented\n");
        spongeNodes[4] = 2;
        spongeValue[4] = 1.0;
      }
      //This call initializes the CPML conditions.
      if(usePML==3)
        assert(FALSE,"Not implemented\n");
        /*setupAcousticCPMLBoundsFull(spongeNodes[0],spongeValue[0],alphaValues[0],kValues[0],spongeNodes[1],spongeValue[1],alphaValues[1],kValues[1],
                                    spongeNodes[2],spongeValue[2],alphaValues[2],kValues[2],spongeNodes[3],spongeValue[3],alphaValues[3],kValues[3],
                                    spongeNodes[4],spongeValue[4],alphaValues[4],kValues[4],spongeNodes[5],spongeValue[5],alphaValues[5],kValues[5],
                                    _dependent->kstart(),_dependent->kstart()-1,_model,surfaceBCMode);*/
      else if(usePML==4)
        setupAcousticAttenMPMLBoundsFull(spongeNodes[0],spongeValue[0],alphaValues[0],kValues[0],spongeNodes[1],spongeValue[1],alphaValues[1],kValues[1],
                                    spongeNodes[2],spongeValue[2],alphaValues[2],kValues[2],spongeNodes[3],spongeValue[3],alphaValues[3],kValues[3],
                                    spongeNodes[4],spongeValue[4],alphaValues[4],kValues[4],spongeNodes[5],spongeValue[5],alphaValues[5],kValues[5],
                                    _dependent->kstart(),_dependent->kstart()-1,_model,surfaceBCMode,xpmlFac);
      spongeValue[0]=spongeValue[1]=spongeValue[2]=spongeValue[3]=spongeValue[4]=spongeValue[5]=1.f;
    } else {
      assert(FALSE,"Not implemented\n");
      //This call initializes the sponge conditions
      setupAcousticSpongeBoundsFull(spongeNodes[0],spongeValue[0],spongeNodes[1],spongeValue[1],
                                    spongeNodes[2],spongeValue[2],spongeNodes[3],spongeValue[3],
                                    spongeNodes[4],spongeValue[4],spongeNodes[5],spongeValue[5],
                                    _dependent->kstart(),_dependent->kstart()-1,_model,surfaceBCMode);
      spongeValue[0]=spongeValue[1]=spongeValue[2]=spongeValue[3]=spongeValue[4]=spongeValue[5]=1.f;
    }
  } else {
    //allocate the structure for the boundarys
    if(usePML && !spongeNodes) usePML = 0;
    _boundaries=new acousticSgfdBoundaries(_model,_dependent,
                                           surfaceBCMode,usePML);
    _dependent->setBoundaryConditions(_boundaries);
    if(usePML) {
      if(surfaceBCMode==SURFACE_PRESS_FREE) {  //Zmin should have 2 nodes and no damping if explicit free surface is used
        spongeNodes[4] = 2;
        spongeValue[4] = 1.0;
      }
      //This call initializes the CPML conditions.
      if(usePML==1 || usePML==3)
        setupAcousticCPMLBoundsFull(spongeNodes[0],spongeValue[0],alphaValues[0],kValues[0],spongeNodes[1],spongeValue[1],alphaValues[1],kValues[1],
                                    spongeNodes[2],spongeValue[2],alphaValues[2],kValues[2],spongeNodes[3],spongeValue[3],alphaValues[3],kValues[3],
                                    spongeNodes[4],spongeValue[4],alphaValues[4],kValues[4],spongeNodes[5],spongeValue[5],alphaValues[5],kValues[5],
                                    _dependent->kstart(),_dependent->kstart()-1,_model,surfaceBCMode);
      else if(usePML==4)
        setupAcousticMPMLBoundsFull(spongeNodes[0],spongeValue[0],alphaValues[0],kValues[0],spongeNodes[1],spongeValue[1],alphaValues[1],kValues[1],
                                    spongeNodes[2],spongeValue[2],alphaValues[2],kValues[2],spongeNodes[3],spongeValue[3],alphaValues[3],kValues[3],
                                    spongeNodes[4],spongeValue[4],alphaValues[4],kValues[4],spongeNodes[5],spongeValue[5],alphaValues[5],kValues[5],
                                    _dependent->kstart(),_dependent->kstart()-1,_model,surfaceBCMode,xpmlFac);
      spongeValue[0]=spongeValue[1]=spongeValue[2]=spongeValue[3]=spongeValue[4]=spongeValue[5]=1.f;
    } else {
      //This call initializes the sponge conditions
      setupAcousticSpongeBoundsFull(spongeNodes[0],spongeValue[0],spongeNodes[1],spongeValue[1],
                                  spongeNodes[2],spongeValue[2],spongeNodes[3],spongeValue[3],
                                  spongeNodes[4],spongeValue[4],spongeNodes[5],spongeValue[5],
                                  _dependent->kstart(),_dependent->kstart()-1,_model,surfaceBCMode);
      spongeValue[0]=spongeValue[1]=spongeValue[2]=spongeValue[3]=spongeValue[4]=spongeValue[5]=1.f;
    }
  }

  DEF_PARALLEL(_model->parallelDef());
  DEF_MODEL_SIZE(_model->modelDef());
  //set the acousticSgfdBoundaries flags that indicate whether this domain
  //has any nodes in one of the boundary (CPML or sponge) zones.  These
  //can be thicker than one single domain.
  if(procLim[0]<spongeNodes[0]) _boundaries->minXb = true;
  if(procLim[2]<spongeNodes[2]) _boundaries->minYb = true;
  if(procLim[4]<spongeNodes[4] && !surfaceBCMode) _boundaries->minZb = true;
  if(procLim[1]>globalNX-spongeNodes[1]+1) _boundaries->maxXb = true;
  if(procLim[3]>globalNY-spongeNodes[3]+1) _boundaries->maxYb = true;
  if(procLim[5]>globalNZ-spongeNodes[5]+1) _boundaries->maxZb = true;

  initSend();
  
  //send a completion message and return
  sendMessage(Parent,MESSAGE_SET_BC,NULL);
}
