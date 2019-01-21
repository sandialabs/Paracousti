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
 *  readfile.cc
 *
 *
 *  Defines class functions used for reading a basic file.
 *
 *  Defines the following class functions:
 *  ReadFile::getWord
 *  ReadFile::checkChar
 *  ReadFile::getChar
 *  ReadFile::isEOF
 *  ReadFile::discardWhite
 *
 */
#include "readfile.hh"
//
//routines for class ReadFile
//
char* ReadFile::getWord(){
  if(_bufPos>=strlen(_buffer)){
    _bufPos=0;
    if(!fgets(_buffer,511,_infile)){
      _isEOF=1;
      return NULL;
    }
  }
  discardWhite();
  if(_isEOF) return NULL;

  int wordLen=0;
  for(int i=_bufPos;i<strlen(_buffer);i++){
    if(_buffer[i]=='\0' ||
       isspace(_buffer[i]) ||
       _buffer[i]==',' ||
       _buffer[i]==';' ||
       _buffer[i]=='{' ||
       _buffer[i]=='}') break;
    wordLen++;
  }
  assert(wordLen,"ReadFile::getWord-- no word:buffer=\"%s\"",_buffer+_bufPos);
  char* toReturn=new char[wordLen+2];
  strncpy(toReturn,_buffer+_bufPos,wordLen);
  toReturn[wordLen]='\0';
  _bufPos+=wordLen;
  return toReturn;
}

char ReadFile::checkChar(){
  if(_isEOF) return '\0';
  if(_bufPos>=strlen(_buffer)){
    if(!fgets(_buffer,511,_infile)){
      _isEOF=1;
      return '\0';
    }
    _bufPos=0;
  }
  return _buffer[_bufPos];
}

char ReadFile::getChar(){
  if(_bufPos>=strlen(_buffer)){
    _bufPos=0;
    if(!fgets(_buffer,511,_infile)){
      _isEOF=1;
      return '\0';
    }
  }
  return _buffer[_bufPos++];

}

int ReadFile::isEOF(){
  discardWhite();
  return _isEOF;
}

char ReadFile::discardWhite(){
  if(_isEOF) return '\0';
  char next;
  for(;;){
    next=checkChar();
    if(_isEOF||!isspace(next)) break;
    getChar();
  }
  return next;
}
