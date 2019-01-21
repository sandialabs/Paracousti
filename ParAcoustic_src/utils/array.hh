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
 *  array.hh
 *
 *
 *  Declares and defines a series of template classes that handle basic vector-type container
 *  classes.  This should eventually be replaced by the STL vector class.
 *
 *  Declares the following template classes:
 *  array
 *  arrayI
 *  arrayD
 *  sortedArray
 *  sortedSet
 *  sortedArrayI
 *  sortedArrayD
 *  sortedSetI
 *  sortedSetD
 *
 *  Defines the following template class functions (excludes constructors/destructors):
 *  EmptyArray
 *  setNullReturn
 *  setErrorOnEOB
 *  isInBounds
 *  operator []
 *  arrayIndex
 *  ensureRoomToAdd
 *  updateIndexForAdded
 *  Add
 *  Remove
 *  isEmpty
 *  size
 *  maxSize
 *  last
 *  first
 *  current
 *  next
 *  prev
 *  operator ++
 *  operator --
 *  forAll
 *  find
 *  findText
 *  sort
 *  vals
 *  setsize
 *  DeleteElems
 *  FreeElems
 *  SetSortFunc
 *  reset
 *  middle
 *  setBotHalf
 *  setTopHalf
 *  Sort
 *
 */
/*Automatic limit-checking expandable arrays

    These are useful utility structures, they should NOT need to be
    modified or understood by other programmers.
*/

#include <stdio.h>
#include <stdlib.h>

#include "nstdutil.hh"

#ifndef _array
#define _array

///The basis for all array types is Array this is a pure virtual class.
template <class Type> class array{
public:
  array(int low,int size,int delta){
    _low= low;
    _topIndex= low-1; //no elements yet
    _size = size;
    _delta= delta;
    _errorOnEOB=1;

    //_vals = new Type[1+size-low];
    assert((_vals= (Type*)malloc((1+_size-_low)*sizeof(Type)))!=NULL,
     "array::array--unable to allocate array with %i members of size %i (%i bytes)",
     1+_size-_low,sizeof(Type),(1+_size-_low)*sizeof(Type));
    _currIndex = _topIndex;
  }
  array(Type* values,int size,int low,int delta){
    _low=low;
    _size=size;
    _topIndex=_size+_low-1;
    _delta=delta;
    _errorOnEOB=1;

    _vals=values;
    _currIndex=_topIndex;
  }

  virtual ~array(){
    free((void*)_vals);
  }
  virtual int DeleteElems()=0;

  array* EmptyArray(){
    _topIndex=_low-1;
    return this;
  }

  //array indexing
  array* setNullReturn(Type null){
    _null=null;
    return this;
  }
  int setErrorOnEOB(){
    _errorOnEOB=!_errorOnEOB;
    return _errorOnEOB;
  }
  int isInBounds(int index){return ISMID(_low,index,_topIndex);}
  //this sets current
  Type& operator [](int toIndex){
    if(!ISMID(_low,toIndex,_topIndex)){
      deassert(_errorOnEOB,
        "array::[%i] is out of bounds %i to %i",
        toIndex,_low,_topIndex);

      _nullReturn=_null;
      return _nullReturn;
    }else{
      _currIndex=toIndex;
      return _vals[toIndex-_low];
    }
  }
  //this does not set current
  Type& arrayIndex(int toIndex){
    if(!ISMID(_low,toIndex,_topIndex)){
      deassert(_errorOnEOB,
        "array::[%i] is out of bounds %i to %i",
        toIndex,_low,_topIndex);
      return _nullReturn;
    }else{
      return _vals[toIndex-_low];
    }
  }

  //Two routines to allow the user to add elements directly to the _vals space.
  Type* ensureRoomToAdd(int n){
    if(_topIndex+n >= _size){
      _size = _topIndex+n+_delta;
      assert((long)(_vals=(Type*)realloc((void*)_vals,(1+_size-_low)*sizeof(Type))),
       "array::ensureRoomToAdd--unable to re-allocate array to %i members of size %i (%i bytes)",
       1+_size-_low,sizeof(Type),(1+_size-_low)*sizeof(Type));
    }

    return &_vals[_topIndex]+1;
  }
  Type& updateIndexForAdded(int n){
    _topIndex+=n;
    return _vals[_topIndex];
  }

  //Routine for conventional element addition.
  Type& Add(Type elem){
    if(_topIndex == _size){
      _size += _delta;
      //Type* _temp = new Type[1+_size-_low];
      assert((long)(_vals=(Type*)realloc((void*)_vals,(1+_size-_low)*sizeof(Type))),
       "array::Add(Type)--unable to re-allocate array to %i members of size %i (%i bytes)",
       1+_size-_low,sizeof(Type),(1+_size-_low)*sizeof(Type));
    }
    _topIndex++;
    _vals[_topIndex]= elem;
    _currIndex= _topIndex;
    return _vals[_topIndex];
  }
  //might want to add another entire array
  Type& Add(array* source,int toDelete=0){
    source->first();
    for(int i=0;i<source->array<Type>::size();i++){
      Add((*source)++);
    }
    if(toDelete) delete source;
    return array<Type>::last();
  }
  //or add several elements at once
  Type& Add(int n,Type elems[]){
    for(int i=0;i<n;Add(elems[i++]));
    return array<Type>::last();
  }

  void Remove(int index){
    if(index==_low && index == _topIndex)
      _topIndex= -1;
    else if(index==_topIndex)
      _topIndex--;
    else{
      (*this)[index]= (*this)[_topIndex];
      _topIndex--;
    }
  }
  void Remove(Type* elem){
    int index=elem-_vals;
    Remove(index);
  }

  int isEmpty(){
    return !this->array<Type>::size();
  }

  int size(){return _topIndex-_low+1;}
  int maxSize(){return _size;}

  //Iterator functions
  //these set current
  Type& last(){return (*this)[_topIndex];}
  Type& first(){return (*this)[_low];}
  //these don't set current
  Type& current(){return this->arrayIndex(_currIndex);}
  Type& next(){return this->arrayIndex(_currIndex+1);}
  Type& prev(){return this->arrayIndex(_currIndex-1);};
  Type& operator ++(int){return this->arrayIndex(_currIndex++);}
  Type& operator ++(){return this->arrayIndex(++_currIndex);}
  Type& operator --(int){return this->arrayIndex(_currIndex--);}
  Type& operator --(){return this->arrayIndex(--_currIndex);}

  int forAll(void (*func)(Type)){
    first();
    for(int i=0;i<array<Type>::size();i++){
      func((*this)++);
    }
    return array<Type>::size();
  }
  int find(Type key,Type& value,int (*cmpFunc)(Type,Type)){
    first();
    for(int i=0;i<array<Type>::size();i++){
      if(cmpFunc(key,current())){
         value=current();
         return 1;
      }
    }
    return 0;
  }
  //This can be used to search for a text string
  Type& findText(char* key,char *(*extractor)(const Type)){
    for(int i=_low;i<=_topIndex;i++){
      char* tstr=extractor(_vals[i]);
      if(!strcmp(tstr,key))
        return _vals[i];
    }
    return _nullReturn;
  }
  array<Type> *sort(int (*sortFunc)(const void*,const void*)){
    qsort((void*)_vals,_topIndex-_low+1,sizeof(Type),sortFunc);
    return this;
  }

  //dangerous functions for the direct manipulation of array internals
  // sometimes make coding easier/faster
  Type* vals(){return _vals;}
  int setsize(int size){
    _size=size;
    _topIndex=_size+_low-1;
    return _size;
  }
protected:
  int _low,_topIndex,_size,_delta;
  Type* _vals,_nullReturn,_null;
  int _currIndex;
  int _errorOnEOB;
};

/// Real class for indirectly referenced arrays (strings, classes, etc)
template <class Type> class arrayI:public array<Type>{
public:
  arrayI(int low=0,int size=100,int delta=100):
      array<Type>(low,size,delta){}
  arrayI(Type* values,int size,int low=0,int delta=100):
      array<Type>(values,size,low,delta){}

  virtual ~arrayI(){}
  //since this is the indirect array type, it is sometimes useful
  // to delete all the elements
  virtual int DeleteElems(){
    int num=0;
    for(int i=array<Type>::_low;i<=array<Type>::_topIndex;i++){
      num++;
      delete array<Type>::_vals[i];
    }
    array<Type>::_topIndex=array<Type>::_low-1; //now has no elements
    return num;
  }
  int FreeElems(){
    int num=0;
    for(int i=array<Type>::_low;i<=array<Type>::_topIndex;i++){
      num++;
      free(array<Type>::_vals[i]);
    }
    array<Type>::_topIndex=array<Type>::_low-1; //now has no elements
    return num;
  }
};

/// Real class for directly referenced arrays (floats, ints, etc.)
template <class Type> class arrayD:public array<Type>{
public:
  arrayD(int low=0,int size=100,int delta=100):
      array<Type>(low,size,delta){}
  arrayD(Type* values,int size,int low=0,int delta=100):
      array<Type>(values,size,low,delta){}

  virtual ~arrayD(){}
  //just set _topIndex to make the array empty
  virtual int DeleteElems(){
    array<Type>::_topIndex=array<Type>::_low-1; //now has no elements
    return 0;
  }
  Type& Add(){
    if(array<Type>::_topIndex == array<Type>::_size){
      array<Type>::_size += array<Type>::_delta;
      //Type* _temp = new Type[1+_size-_low];
      assert((array<Type>::_vals=(Type*)realloc((void*)array<Type>::_vals,(1+array<Type>::_size-array<Type>::_low)*sizeof(Type)))!=NULL,
	     "arrayD::Add--unable to re-allocate array to %i members of size %i (%i bytes)",
	     1+array<Type>::_size-array<Type>::_low,sizeof(Type),(1+array<Type>::_size-array<Type>::_low)*sizeof(Type));
    }
    array<Type>::_topIndex++;
    array<Type>::_currIndex= array<Type>::_topIndex;
    return array<Type>::_vals[array<Type>::_topIndex];
  }    
  Type& Add(Type elem){
    if(array<Type>::_topIndex == array<Type>::_size){
      array<Type>::_size += array<Type>::_delta;
      //Type* _temp = new Type[1+_size-_low];
      assert((array<Type>::_vals=(Type*)realloc((void*)array<Type>::_vals,(1+array<Type>::_size-array<Type>::_low)*sizeof(Type)))!=NULL,
	     "arrayD::Add--unable to re-allocate array to %i members of size %i (%i bytes)",
	     1+array<Type>::_size-array<Type>::_low,sizeof(Type),(1+array<Type>::_size-array<Type>::_low)*sizeof(Type));
    }
    array<Type>::_topIndex++;
    array<Type>::_vals[array<Type>::_topIndex]= elem;
    array<Type>::_currIndex= array<Type>::_topIndex;
    return array<Type>::_vals[array<Type>::_topIndex];
  }    
  //might want to add another entire array
  Type& Add(array<Type>* source,int toDelete=0){
    source->first();
    for(int i=0;i<source->array<Type>::size();i++){
      Add((*source)++);
    }
    if(toDelete) delete source;
    return array<Type>::last();
  }
  //or add several elements at once
  Type& Add(int n,Type elems[]){
    for(int i=0;i<n;Add(elems[i++]));
    return array<Type>::last();
  }
  Type& Add(int n,Type value){
    for(int i=0;i<n;Add(value),i++);
    return array<Type>::last();
  }

};

/// Array of sorted data; pure virtual class.
template <class Type> class sortedArray:public array<Type>{
 public:
  sortedArray(int(* cmpFunc)(const Type * ,const Type * ),
      int low,int size,int delta):
      array<Type>(low,size,delta){
    _cmpFunc=(int (*)(const void*,const void*))cmpFunc;
    _sorted= 1;
  }
  virtual ~sortedArray(){}
  //normal array maintanance functions
  void SetSortFunc(int (*cmpFunc)(const Type *,const Type *)){
    _cmpFunc=(int (*)(const void*,const void*))cmpFunc;
    _sorted=0;
    Sort();
  }
  Type& Add(Type elem){
    array<Type>::Add(elem);
    if(array<Type>::size()>1){
      _sorted=_sorted&&(_cmpFunc(&array<Type>::prev(),&array<Type>::current())<=0);
    }
    return array<Type>::last();
  }

  void Remove(int index){
    array<Type>::Remove(index);
    _sorted=0;
  }
  Type& operator [](int toIndex){
    Sort();
    array<Type>::_currIndex=toIndex;
    return array<Type>::arrayIndex(toIndex);
  }

  //functions to iterate through an array doing binary type things
  Type& reset(){
    _top=array<Type>::_topIndex;
    _bot=array<Type>::_low;
    return middle();
  }

  int find(Type key,Type& value){
    if(array<Type>::size()<=0) return 0;
    else if(array<Type>::size()==1){
      value=array<Type>::first();
      return !_cmpFunc(&key,&array<Type>::first());
    }
    Type middle=reset();
    for(;;){
      int cmpVal=_cmpFunc(&key,&middle);
      if(cmpVal==0){
        value=array<Type>::current();
        return 1;
      }else if(cmpVal>0){
        if(!this->isInBounds(array<Type>::_currIndex+1)) return 0;
        Type next;
        next=this->next();
        if(next==array<Type>::_nullReturn){
          return 0;
        }
        int cmp2=_cmpFunc(&key,&next);
        if(cmp2==0){
          value=this->next();
          return 1;
        }else if(cmp2<0){
          return 0;
        }
      }else{
        if(!this->isInBounds(array<Type>::_currIndex-1)) return 0;
        Type prev;
        if(!(prev=this->array<Type>::prev())) return 0;
        int cmp2=_cmpFunc(&key,&prev);
        if(cmp2==0){
          value=this->array<Type>::prev();
          return 1;
        }else if(cmp2>0){
          return 0;
        }
      }
      middle=(cmpVal<0)?setBotHalf():setTopHalf();
    }
  }

  Type& middle(){ //returns the middle element of top-bot
    Sort();
    array<Type>::_currIndex = (_top+_bot)/2;
    return (*this)[array<Type>::_currIndex];
  }
  Type& setBotHalf(){ //sets top to curr
    _top = array<Type>::_currIndex-1;
    return this->middle();
  }
  Type& setTopHalf(){    //sets bot to curr
    _bot= array<Type>::_currIndex+1;
    return this->middle();
  }
  int forAll(void (*func)(Type)){
    Sort();
    return array<Type>::forAll(func);
  }
  void Sort(){
    if(!_sorted&&(array<Type>::size()>1)){
      qsort(array<Type>::_vals,array<Type>::_topIndex-array<Type>::_low+1,sizeof(Type),
      (int (*)(const void *,const void *))_cmpFunc);
      _sorted=1;
    }
  }
protected:
  int _sorted,_bot,_top;
  int (*_cmpFunc)(const void*,const void*);
};

///Set of sorted data (no repeats) pure virtual class.
template <class Type> class sortedSet:public sortedArray<Type>{
public:
  sortedSet(int(* cmpFunc)(const Type * ,const Type * ),
      int low,int size,int delta):
      sortedArray<Type>(cmpFunc,low,size,delta){}
  virtual ~sortedSet(){}
  //normal array maintanance functions
  Type& Add(Type elem){
    Type value;
    if(this->find(elem,value)){
      return array<Type>::current();
    }
    return sortedArray<Type>::Add(elem);
  }
};

///Real indirect implementation of sortedArray.
template <class Type> class sortedArrayI:
    public sortedArray<Type>{
public:
  sortedArrayI(int(* cmpFunc)(const Type * ,const Type * ),
      int low=0,int size=50,int delta=100):
      sortedArray<Type>(cmpFunc,low,size,delta){}
  virtual ~sortedArrayI(){}
  //since this is the indirect array type, it is sometimes useful
  // to delete all the elements
  virtual int DeleteElems(){
    int num=0;
    for(int i=array<Type>::_low;i<=array<Type>::_topIndex;i++){
      num++;
      delete array<Type>::_vals[i];
    }
    array<Type>::_topIndex=array<Type>::_low-1; //now has no elements
    return num;
  }
};

///Real direct implementation of sortedArray.
template <class Type> class sortedArrayD:
    public sortedArray<Type>{
public:
  sortedArrayD(int(* cmpFunc)(const Type * ,const Type * ),
      int low=0,int size=50,int delta=100):
      sortedArray<Type>(cmpFunc,low,size,delta){}
  virtual ~sortedArrayD(){}
  //just set _topIndex to make the array empty
  virtual int DeleteElems(){
    array<Type>::_topIndex=array<Type>::_low-1; //now has no elements
    return 0;
  }
};

///Real indirect implementation of sortedSet.
template <class Type> class sortedSetI:
    public sortedSet<Type>{
public:
  sortedSetI(int(* cmpFunc)(const Type * ,const Type * ),
      int low=0,int size=50,int delta=100):
      sortedSet<Type>(cmpFunc,low,size,delta){}
  virtual ~sortedSetI(){
    free((void*)array<Type>::_vals);
  }
  //since this is the indirect array type, it is sometimes useful
  // to delete all the elements
  virtual int DeleteElems(){
    int num=0;
    for(int i=array<Type>::_low;i<=array<Type>::_topIndex;i++){
      num++;
      delete array<Type>::_vals[i];
    }
    array<Type>::_topIndex=array<Type>::_low-1; //now has no elements
    return num;
  }
};

///Real direct implementation of sortedSet.
template <class Type> class sortedSetD:
    public sortedSet<Type>{
public:
  sortedSetD(int(* cmpFunc)(const Type * ,const Type * ),
      int low=0,int size=50,int delta=100):
      sortedSet<Type>(cmpFunc,low,size,delta){}
  virtual ~sortedSetD(){}
  //just set _topIndex to make the array empty
  virtual int DeleteElems(){
    array<Type>::_topIndex=array<Type>::_low-1; //now has no elements
    return 0;
  }
};

//
// And some common types of arrays
//
typedef arrayD<double> doubleArray,*doubleArrayPtr;
typedef sortedArrayD<double> sortedDoubleArray;
typedef sortedSetD<double> sortedDoubleSet;
typedef arrayI<doublePtr> doublePtrArray;

typedef arrayD<float> floatArray,*floatArrayPtr;
typedef sortedArrayD<float> sortedFloatArray;
typedef sortedSetD<float> sortedFloatSet;
typedef arrayI<floatPtr> floatPtrArray;

typedef arrayD<int> intArray,*intArrayPtr;
typedef sortedArrayD<int> sortedIntArray;
typedef sortedSetD<int> sortedIntSet;
typedef arrayI<intPtr> intPtrArray;

typedef arrayI<charPtr> stringArray,*stringArrayPtr;
typedef sortedArrayI<charPtr> sortedStringArray;
typedef sortedSetI<charPtr> sortedStringSet;
#endif
