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
 *  selector.hh
 *
 *
 *  Declares the class that handles selecting and adding criteria to general data types
 *  in arrays.  This is mostly used for models with attenuation.
 *
 *  Declares the following classes:
 *  selector
 *
 *  Defines the following structures:
 *  fieldStruct
 *  criteria
 *  selection
 *
 */

#ifndef _selector_hh_
#define _selector_hh_

#include <stdio.h>
#include <stdlib.h>

#include "array.hh"

typedef struct _fieldStruct_{
  char name[512];
  float* field;
} fieldStruct,*fieldStructPtr;
typedef arrayI<fieldStructPtr> fieldArray;

typedef struct _criteria_{
  int field,direction;
  float value;
} criteria,*criteriaPtr;
typedef arrayI<criteriaPtr> criteriaArray,*criteriaArrayPtr;

typedef struct _selection_{
  criteriaArray* criterion;

  floatArray* floatVars;
  intArray* intVars;
  floatPtrArray* floatPtrVars;
  intArray* floatPtrSizes;
  intPtrArray* intPtrVars;
  intArray* intPtrSizes;
} selection,*selectionPtr;
typedef arrayI<selectionPtr> selectionArray;

class selector{
protected:
  int _nx,_ny,_nz;
  float _minX,_dx,_minY,_dy,_minZ,_dz;

  //These are additional fields on which to make the selection.
  // This is an array of pointers to arrays which should be of size nxyz.
  fieldArray* _fields;

  //These are the possible selections.
  selectionArray* _selections;
  int _lastMatch;
public:
  //Here are arrays of the variables that hold the parameters that
  // will be used for a valid selection.
  selector(int nx,int ny,int nz,
	   float minX,float dx,
	   float minY,float dy,
	   float minZ,float dz){
    _nx=nx;
    _ny=ny;
    _nz=nz;

    _minX=minX;
    _dx=dx;
    _minY=minY;
    _dy=dy;
    _minZ=minZ;
    _dz=dz;

    _selections=new selectionArray;
    _lastMatch=-1;

    _fields=new fieldArray;
  }

  ~selector(){
    for(int i=0;i<_selections->size();i++){
      (*_selections)[i]->criterion->DeleteElems();
      delete (*_selections)[i]->criterion;

      delete (*_selections)[i]->floatVars;
      delete (*_selections)[i]->intVars;
      
      (*_selections)[i]->floatPtrVars->DeleteElems();
      delete (*_selections)[i]->floatPtrVars;
      delete (*_selections)[i]->floatPtrSizes;
      (*_selections)[i]->intPtrVars->DeleteElems();
      delete (*_selections)[i]->intPtrVars;
      delete (*_selections)[i]->intPtrSizes;
    }
    _selections->DeleteElems();
    delete _selections;

    _fields->DeleteElems();
    delete _fields;
  }

  int size(){return _selections->size();}
  selector* addField(const char* fieldName,float *field){
    fieldStruct* newField=(fieldStruct*)malloc(sizeof(fieldStruct));
    assert(newField!=NULL,
	   "selector::addField--unable to allocate field %i",
	   1+_fields->size());
    strcpy(newField->name,fieldName);
    newField->field=field;
    _fields->Add(newField);
    return this;
  }
  selector* addSelector(){
    selection* newSelector=newSelection();
    _lastMatch=_selections->size();
    _selections->Add(newSelector);
    return this;
  }
  float addCriteria(int field,int direction,float value,int index=-1);
  float& addFloatVar(float value=0,int index=-1){
    if(index<0) index=_lastMatch;
    return (*_selections)[index]->floatVars->Add(value);
  }
  int& addIntVar(int value=0,int index=-1){
    if(index<0) index=_lastMatch;
    (*_selections)[index]->intVars->Add(value);
    return (*_selections)[index]->intVars->last();
  }
  float* addFloatPtrVar(int n,float* values=NULL,int index=-1);

  int& intVar(int vindex,int sindex=-1){
    if(sindex<0) sindex=_lastMatch;
    return (*(*_selections)[sindex]->intVars)[vindex];
  }
  float& floatVar(int vindex,int sindex=-1){
    if(sindex<0) sindex=_lastMatch;
    return (*(*_selections)[sindex]->floatVars)[vindex];
  }
  float* floatPtrVar(int vindex,int sindex=-1){
    if(sindex<0) sindex=_lastMatch;
    return (*(*_selections)[sindex]->floatPtrVars)[vindex];
  }

  int* generateIntVar(int vindex);
  float** generateFloatPtrVar(int vindex);

  //Generate an array of indicies to the correct selector.
  int* generateIndexArray();

  void parseSelectorArguments(int& i,int position,int argc,char *argv[]);

#ifdef PARALLEL_PROGRAM
  int unpackSelectorMessage(int index=-1);
  int calcMessageSize(int index=-1);
  void packSelectorMessage(int index=-1);
#endif //#if PARALLEL_PROGRAM
  
protected:
  char* fieldName(int index,char* buffer){
    switch(index){
    case 0:
      strcpy(buffer,"X");
      break;
    case 1:
      strcpy(buffer,"Y");
      break;
    case 2:
      strcpy(buffer,"Z");
      break;
    default:
      strcpy(buffer,(*_fields)[index-3]->name);
    }
    return buffer;
  }
  int findField(char *name);
  selection* newSelection(){
    selection* newSelection=(selection*)malloc(sizeof(selection));
    assert(newSelection!=NULL,
	   "selector::newSelection--unable to allocate new selection");
    newSelection->criterion=new criteriaArray;
    
    newSelection->floatVars=new floatArray;
    newSelection->intVars=new intArray;

    newSelection->floatPtrVars=new floatPtrArray;
    newSelection->floatPtrSizes=new intArray;
    newSelection->intPtrVars=new intPtrArray;
    newSelection->intPtrSizes=new intArray;
    return newSelection;
  }
};

#endif //#ifndef _selector_hh_
