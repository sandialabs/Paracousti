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
 *  selector.cc
 *
 *
 *  Defines class functions used to select/set criteria for general data types that
 *  are in an array of values.  This is primarily used for models with attenuation.
 *
 *  Defines the following class functions:
 *  selector::addCriteria
 *  selector::addFloatPtrVar
 *  selector::generateIntVar
 *  selector::generateFloatPtrVar
 *  selector::generateIndexArray
 *  selector::parseSelectorArguments
 *  selector::unpackSelectorMessage
 *  selector::packSelectorMessage
 *  selector::findField
 *
 *  Defines the following macros:
 *  INTERNAL_USAGE_DEFINE
 *  USAGE
 *
 */

#include "selector.hh"

#ifndef USAGE
#define INTERNAL_USAGE_DEFINE 1
#define USAGE "IMPROPER CALL\nUSAGE: %s UNKNOWN"
#endif

float selector::addCriteria(int field,int direction,float value,int index){
  if(index<0) index=_lastMatch;
  criteria* newCriteria=(criteria*)malloc(sizeof(criteria));
  assert(newCriteria!=NULL,
   "selector::addCriteria--unable to allocate %ith criteria for %ith selector",
   1+(*_selections)[index]->criterion->size(),1+_selections->size());
  newCriteria->field=field;
  newCriteria->direction=direction;
  newCriteria->value=value;
  (*_selections)[index]->criterion->Add(newCriteria);
  return value;
}

float* selector::addFloatPtrVar(int n,float* values,int index){
  if(index<0) index=_lastMatch;
  if(values){
    (*_selections)[index]->floatPtrVars->Add(values);
  }else{
    float* newVals=(float*)malloc(n*sizeof(float));
    assert(newVals!=NULL,
     "selector::addFloatPtrVar--unable to allocate %i floats for floatVars in selector %i",
     n,index);
    for(int i=0;i<n;newVals[i++]=0.0);
    (*_selections)[index]->floatPtrVars->Add(newVals);
  }
  (*_selections)[index]->floatPtrSizes->Add(n);
  return (*_selections)[index]->floatPtrVars->last();
}

int* selector::generateIntVar(int vindex){
  int* toReturn=(int*)malloc(size()*sizeof(int));
  assert(toReturn!=NULL,
   "selector::generateIntVar--unable to allocate %i ints",
   size());
  for(int i=0;i<size();i++)
    toReturn[i]=intVar(vindex,i);
  return toReturn;
}
float** selector::generateFloatPtrVar(int vindex){
  float** toReturn=(float**)malloc(size()*sizeof(float*));
  assert(toReturn!=NULL,
   "selector::generateIntVar--unable to allocate %i float*",
   size());
  for(int i=0;i<size();i++)
    toReturn[i]=floatPtrVar(vindex,i);
  return toReturn;
}

//Generate an array of indicies to the correct selector.
int* selector::generateIndexArray(){
  int nxy=_nx*_ny,nxyz=_nx*_ny*_nz;
  int* toReturn=(int*)malloc(nxyz*sizeof(int));
  assert(toReturn!=NULL,
   "selector::generateIndexArray--unable to allocate array of %i ints",
   nxyz);
  for(int kk=0;kk<_nz;kk++){
    for(int jj=0;jj<_ny;jj++){
      for(int ii=0;ii<_nx;ii++){
        int index=ii+jj*_nx+kk*nxy;
        toReturn[index]=-1; //Preset to no selection.

        for(int j=0;j<_selections->size();j++){
          selection* currS=(*_selections)[j];
          int isValid=TRUE;
          for(int k=0;k<currS->criterion->size();k++){
            criteria* currC=(*currS->criterion)[k];
            //Calculate the value for this field.
            float value;
            switch(currC->field){
            case 0:
        value=_minX+ii*_dx;
        break;
            case 1:
        value=_minY+jj*_dy;
        break;
            case 2:
        value=_minZ+kk*_dz;
        break;
            default:
        value=(*_fields)[currC->field-3]->field[index];
            }

            if((currC->direction==-1 && value>currC->value) ||
         (currC->direction== 1 && value<currC->value)){
        //Failure of this selector to fit the criteria.
        isValid=FALSE;
        break;
            }
          }
          if(isValid){
            //This was a successful match.
            toReturn[index]=j;
            break;
          }
        }

        assert(toReturn[index]!=-1,
         "selector::generateIndexArray--%i,%i,%i; no successful matches from %i selectors",
         ii,jj,kk,_selections->size());
      }
    }
  }
  return toReturn;
}

void selector::parseSelectorArguments(int& i,int position,int argc,char *argv[]){
  switch(argv[i][position]){
  case 'x':
  case 'X':
    assert(argc>i+2,USAGE,argv[0]);
    tEprintf(Verbose,"   Selector: X [");
    if(isdigit(argv[++i][0])){
tEprintf(Verbose,"%.4g,",
   addCriteria(0,-1,atof(argv[i])));
    }else{
tEprintf(Verbose,"--,");
    }
    if(isdigit(argv[++i][0])){
tEprintf(Verbose,"%.4g]\n",
   addCriteria(0, 1,atof(argv[i])));
    }else{
tEprintf(Verbose,"--]\n");
    }
    break;

  case 'Y':
  case 'y':
    assert(argc>i+2,USAGE,argv[0]);
    tEprintf(Verbose,"   Selector: Y [");
    if(isdigit(argv[++i][0])){
tEprintf(Verbose,"%.4g,",
   addCriteria(1,-1,atof(argv[i])));
    }else{
tEprintf(Verbose,"--,");
    }
    if(isdigit(argv[++i][0])){
tEprintf(Verbose,"%.4g]\n",
   addCriteria(1, 1,atof(argv[i])));
    }else{
tEprintf(Verbose,"--]\n");
    }
    break;

  case 'Z':
  case 'z':
    assert(argc>i+2,USAGE,argv[0]);
    tEprintf(Verbose,"   Selector: Z [");
    if(isdigit(argv[++i][0])){
tEprintf(Verbose,"%.4g,",
   addCriteria(2,-1,atof(argv[i])));
    }else{
tEprintf(Verbose,"--,");
    }
    if(isdigit(argv[++i][0])){
tEprintf(Verbose,"%.4g]\n",
   addCriteria(2, 1,atof(argv[i])));
    }else{
tEprintf(Verbose,"--]\n");
    }
    break;

  default:
    {
      int fieldIndex=findField(&argv[i][position]);

      assert(argc>i+2,USAGE,argv[0]);
      tEprintf(Verbose,"   Selector: %s [",&argv[i][position]);
      if(isdigit(argv[++i][0])){
        tEprintf(Verbose,"%.4g,",
           addCriteria(fieldIndex,-1,atof(argv[i])));
      }else{
        tEprintf(Verbose,"--,");
      }
      if(isdigit(argv[++i][0])){
        tEprintf(Verbose,"%.4g]\n",
           addCriteria(fieldIndex, 1,atof(argv[i])));
      }else{
        tEprintf(Verbose,"--]\n");
      }
    }
  }
}    

#ifdef PARALLEL_PROGRAM
#include "message_passing.h"
int selector::unpackSelectorMessage(int index){
  int nFields,nSelections;
  
  if(index<1) {
    unpackMessage("ii",&nFields,&nSelections);

    if(nFields){
      //Remember the field names where packed, read them and check that that
      // feild has been loaded.
      char buffer[512];
      for(int i=0;i<nFields;i++){
        unpackMessage("s",buffer);
        findField(buffer);
      }
    }
  }
  
  int startIndex = index<0?0:index;
  int endIndex = index<0?nSelections:index+1;

  //Need to read quite a bit more info for the selections.
  for(int i=startIndex;i<endIndex;i++){
    selection* curr=(selection*)malloc(sizeof(selection));
    _selections->Add(curr);
    assert(curr!=NULL,
     "selector::readSelectorMessage--unable to allocate selection %i",
     i+1);

    int nCriterion,nFloatVars,nIntVars,nFloatPtrVars,nIntPtrVars;
    unpackMessage("i ii ii",&nCriterion,
      &nFloatVars,&nIntVars,&nFloatPtrVars,&nIntPtrVars);

    curr->criterion=new criteriaArray;
    for(int j=0;j<nCriterion;j++){
      criteria* currC=(criteria*)malloc(sizeof(criteria));
      assert(currC!=NULL,
       "selector::readSelectorMessage--unable to allocate criteria %i in selection %i",
       j+1,i+1);

      char buffer[512];
      unpackMessage("si f",buffer,
        &currC->direction,&currC->value);
      currC->field=findField(buffer);
      curr->criterion->Add(currC);
    }

    curr->floatVars=new floatArray;
    for(int j=0;j<nFloatVars;j++){
      float currV;
      unpackMessage("f",&currV);
      curr->floatVars->Add(currV);
    }

    curr->intVars=new intArray;
    for(int j=0;j<nIntVars;j++){
      int currV;
      unpackMessage("i",&currV);
      curr->intVars->Add(currV);
    }
    
    curr->floatPtrVars=new floatPtrArray;
    curr->floatPtrSizes=new intArray;
    for(int j=0;j<nFloatPtrVars;j++){
      int n;
      unpackMessage("i",&n);
      float* currV=(float*)malloc(n*sizeof(float));
      assert(currV!=NULL,
       "selector::readSelectorMessage--unable to allocate floatPtr[%i] %i in selection %i",
       n,j+1,i+1);
       
      unpackMessage("F",currV,n);
      curr->floatPtrVars->Add(currV);
      curr->floatPtrSizes->Add(n);
    }

    curr->intPtrVars=new intPtrArray;
    curr->intPtrSizes=new intArray;
    for(int j=0;j<nIntPtrVars;j++){
      int n;
      unpackMessage("i",&n);
      int* currV=(int*)malloc(n*sizeof(int));
      assert(currV!=NULL,
       "selector::readSelectorMessage--unable to allocate intPtr[%i] %i in selection %i",
       n,j+1,i+1);
       
      unpackMessage("I",currV,n);
      curr->intPtrVars->Add(currV);
      curr->intPtrSizes->Add(n);
    }
  } 
  return _selections->size();
}

//Calculate how big a message we must be able to send for the selectors
int selector::calcMessageSize(int index) {
  int nfields = _fields->size();
  int nsel = _selections->size();
  int maxMessSize = 0;
  for(int i=0;i<_selections->size();++i) {
    selection* curr = (*_selections)[i];
    int nCrit = curr->criterion->size();
    int nIntV=curr->intVars->size();
    int nFloatV=curr->floatVars->size();
    int nFloatPV=curr->floatPtrVars->size();
    int nIntPV=curr->intPtrVars->size();
    int messSize = 100+nCrit*520;  //actually it is onlt 20+..., but leaving a bit of extra
    messSize += nIntV*4;
    messSize += nFloatV*4;
    for(int j=0;j<nFloatPV;++j)
      messSize += 4+(*curr->floatPtrSizes)[j]*4;
    for(int j=0;j<nIntPV;++j)
      messSize += 4+(*curr->intPtrSizes)[j]*4;
    maxMessSize = messSize>maxMessSize?messSize:maxMessSize;
  }
  if(index==-1) maxMessSize *= nsel;
  maxMessSize += 8+nfields*512;
  return maxMessSize;
}

void selector::packSelectorMessage(int index){

  //Just pack the names of the fields, we are assuming that the other end
  // adds the same fields (we will pack field names for the critera so the 
  // order can be different).
  //Only pack the field names if we are doing all selections at once(index=-1)
  //or if index==0
  if(index<1) {
    packMessage("ii",_fields->size(),_selections->size());
    for(int i=0;i<_fields->size();i++)
      packMessage("s",(*_fields)[i]->name);
  }
  
  int startIndex = index<0?0:index;
  int endIndex = index<0?_selections->size():index+1;

  //Need to pack quite a bit more info for the selections.
  for(int i=startIndex;i<endIndex;i++){
    selection* curr=(*_selections)[i];
    int nCriterion=curr->criterion->size(),
    nFloatV=curr->floatVars->size(),nIntV=curr->intVars->size(),
    nFloatPV=curr->floatPtrVars->size(),nIntPV=curr->intPtrVars->size();
    packMessage("i ii ii",nCriterion,
    nFloatV,nIntV,
    nFloatPV,nIntPV);

    for(int j=0;j<curr->criterion->size();j++){
      criteria* currC=(*curr->criterion)[j];
      char buffer[512];
      fieldName(currC->field,buffer);
      packMessage("si f",buffer,currC->direction,currC->value);
    }

    for(int j=0;j<curr->floatVars->size();j++)
      packMessage("f",(*curr->floatVars)[j]);
    for(int j=0;j<curr->intVars->size();j++)
      packMessage("i",(*curr->intVars)[j]);
    
    for(int j=0;j<curr->floatPtrVars->size();j++)
      packMessage("iF",(*curr->floatPtrSizes)[j],
                  (*curr->floatPtrVars)[j],(*curr->floatPtrSizes)[j]);

    for(int j=0;j<curr->intPtrVars->size();j++)
      packMessage("iI",(*curr->intPtrSizes)[j],
                  (*curr->intPtrVars)[j],(*curr->intPtrSizes)[j]);
  }
}
#endif //#if PARALLEL_PROGRAM

int selector::findField(char *name){
  //Check for the coordinate axis fields.
  if(!strcmp(name,"x") || !strcmp(name,"X"))
    return 0;
  if(!strcmp(name,"y") || !strcmp(name,"Y"))
    return 1;
  if(!strcmp(name,"z") || !strcmp(name,"Z"))
    return 2;
  
  for(int i=0;i<_fields->size();i++){
    if(!strcmp(name,(*_fields)[i]->name))
      return i+3;
  }
  assert(FALSE,"selector::findField--unable to match field %s",name);
  return FALSE;
}

#ifdef INTERNAL_USAGE_DEFINE
#undef USAGE
#undef INTERNAL_USAGE_DEFINE
#endif
