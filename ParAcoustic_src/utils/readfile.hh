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
 *  readfile.hh
 *
 *
 *  Declares class used in reading from a basic file.
 *
 *  Declares the following classes:
 *  ReadFile
 *
 */
#ifndef _readfile
#define _readfile

#include <stdio.h>
#include "nstdutil.hh"

class ReadFile{
 public:
  ReadFile(char* name,char* directory=NULL,char* extension=NULL){
    int nameLen=strlen(name)+3;
    if(directory) nameLen+=strlen(directory);
    if(extension) nameLen+=strlen(extension);
    
    _name=new char[nameLen];
    if(directory) sprintf(_name,"%s/",directory);
    sprintf(_name,"%s%s",
	    name,
	    extension?extension:"");

    _buffer[0]='\0';
    _bufPos=0;
    _isEOF=0;

    if(!(_infile=fopen(_name,"r"))){
      _isEOF=1;
    }
  }
  ~ReadFile(){
    if(_name)
      delete _name;
    if(_infile)
      fclose(_infile);
  }

  char* name(){return _name;}

  char* getWord();

  //get or check the next character respectively
  char discardWhite();

  char checkChar();
  char getChar();
  int isEOF();

protected:

  char* _name;
  FILE* _infile;
  char _buffer[512];
  int _bufPos;
  int _isEOF;

};
#endif
