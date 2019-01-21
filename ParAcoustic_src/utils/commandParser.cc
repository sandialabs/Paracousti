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
 *  commandParser.cc
 *
 *
 *  Defines base class for parsing command line arguments.
 *
 *  Defines the following class functions:
 *  flagDescription::printUsage
 *  commandParser::commandParser
 *  commandParser::processArgs
 *  commandParser::recursiveCallProcessArgs
 *  commandParser::addFlagDescription
 *  commandParser::buildUsageString
 *  commandParser::printUsage
 *  commandParser::printHistory
 *  commandParser::repeatCommand
 *
 */

#include "commandParser.hh"
#include "readfile.hh"

//Define a new class for parsing command line arguments. Subclass this to customize
// acceptable args.

void flagDescription::printUsage(char flagTypeDescription[5][256],int printExpert){
  if(!printExpert && _flagType==FLAG_EXPERT)
    return;
  if(!strlen(_flag)){
    fprintf(stderr,"\n");
  }else{
    fprintf(stderr," %s  %s (%s)\n",
            _flag,_fields,
            flagTypeDescription[_flagType]);
    for(int i=0;i<_flagDescriptions->size();i++){
      fprintf(stderr,"    %s\n",
              (*_flagDescriptions)[i]);
    }
  }
}

///Only 1 constructor, set the default values, this includes the usage string.
commandParser::commandParser(int argc,char* argv[]){
  strcpy(_executable,argv[0]);
  
  if(argc<2){
    //If only 1 argument don't need to do any more.
    _flags=NULL;
    
    return;
  }
  //Allocate the command line and add the 0th argument or add a bracket.
  _commandLineSize=1024;
  assert((CommandLine=(char*)malloc(_commandLineSize*sizeof(char)))!=NULL,
         "commandParser--unable to allocate %i chars for command line",
         _commandLineSize);
#ifdef COMPILE_DATE
  sprintf(CommandLine,"%s(Compiled %s) ",argv[0],COMPILE_DATE);
#else
  sprintf(CommandLine,"%s ",argv[0]);
#endif
  
  //Set up stuff for a reasonable help and usage system.
  _usageString=NULL;
  
  strcpy(_flagTypeDescription[1],"REQUIRED");
  strcpy(_flagTypeDescription[2],"Optional");
  strcpy(_flagTypeDescription[3],"Optional may be repeated");
  strcpy(_flagTypeDescription[4],"EXPERT USER ONLY: Optional");
  
  _flags=new flagDescriptionArray();
  addFlagDescription(FLAG_OPTIONAL,"-v","","Toggle verbose model (TRUE).");
  addFlagDescription(FLAG_OPTIONAL,"-help","","Print help.");
#ifdef STORE_RUN_HISTORY
  addFlagDescription(FLAG_OPTIONAL,"-history","[i]",
                     "Print history of program calls.");
  addFlagDescription(FLAG_OPTIONAL,"-repeat","[i]",
                     "Rerun last or index command from history.");
#endif
  addFlagDescription(FLAG_OPTIONAL,"-help_expert","",
                     "Print help with expert options.");
}

///Process the arguments, modifications to behaviour can easily be made
/// through methods processNextArgument or checkProcessFlag.
int commandParser::processArgs(int argc,char* argv[],
                               int recursionDepth,int startArg){
  if(recursionDepth){
    char* endptr=CommandLine+strlen(CommandLine);
    sprintf(endptr,"{");
  }else{
    buildUsageString(argv[0]);
  }
  
  //Step through the arguments and add them to the command line.
  for(int i=startArg;i<argc;i++){
    if(strlen(CommandLine)+strlen(argv[i])+2>_commandLineSize){
      _commandLineSize*=2;
      assert((CommandLine=(char*)realloc(CommandLine,_commandLineSize*sizeof(char)))!=NULL,
             "processArgs--unable to reallocate %i chars for command line",
             _commandLineSize);
    }
    
    char* endptr=CommandLine+strlen(CommandLine);
    sprintf(endptr,"%s ",argv[i]);
  }
  
  //Now step through the arguments.
  for(int i=startArg;i<argc;i++){
    if(argv[i][0]!='-'){
      //check for argument of form ARGUMENT_FILE:*
      if(!strncmp(argv[i],"ARGUMENT_FILE:",14)){
        recursiveCallProcessArgs(argv[i]+14,argv[0],recursionDepth);
      }else if(!strncmp(argv[i],"AF:",3)){ //Allow a short form also
        recursiveCallProcessArgs(argv[i]+3,argv[0],recursionDepth);
      }else{
        processNextArgument(argc,argv,i);
      }
    }else{
      //Check for the -help or --help flags, if there is a match then print some help and
      // exit.
      if(!strcmp(argv[i],"-help") || !strcmp(argv[i],"--help")){
        printUsage(TRUE);
      }else if(!strcmp(argv[i],"-help_expert") || !strcmp(argv[i],"--help_expert")){
        printUsage(TRUE,TRUE);
      }else
#ifdef STORE_RUN_HISTORY
        if(!strcmp(argv[i],"-history") || !strcmp(argv[i],"--history")){
          int nhistory=10;
          if(argc>i+1 && isdigit(argv[i+1][0]))
            nhistory=atoi(argv[i+1]);
          printHistory(nhistory);
        }else if(!strcmp(argv[i],"-repeat") || !strcmp(argv[i],"--repeat")){
          int index=0;
          if(argc>i+1 && isdigit(argv[i+1][0]))
            index=atoi(argv[i+1]);
          repeatCommand(index);
        }
#endif //#ifdef STORE_RUN_HISTORY
      if(!checkProcessFlag(argc,argv,i)){
        //Check for flags that are used in all my codes, verbose is the only one so far.
        switch(argv[i][1]){
          case 'v':
            Verbose=!Verbose;
            if(strlen(argv[i])>2 && isdigit(argv[i][2]))
              Verbose=atoi(argv[i]+2);
            tEprintf(Verbose,"Verbose flag set to %i\n",Verbose);
            break;
            
          default:
            assert(FALSE,_usageString,i,argv[i],argv[0]);
        }
      }
    }
  }
  
  //And finish the command line.
  if(!recursionDepth){
    fprintf(stderr,"\n*******\n%s\n*******\n\n",CommandLine);
    CDFSourceName=(char*)malloc((strlen(CommandLine)+1)*sizeof(char));
    strcpy(CDFSourceName,CommandLine);
  }else{
    char* endptr=CommandLine+strlen(CommandLine);
    sprintf(endptr,"} ");
  }
  
#ifdef STORE_RUN_HISTORY
  char hfn[1024];
  strcpy(hfn,STORE_RUN_HISTORY);
  if(recursionDepth==0 && hfn){
    FILE* historyFile=fopen(hfn,"r");
    if(historyFile){
      //History file exists, copy the contents to a new history file;
      // and then add the new line.
      char nhfn[1024];
      sprintf(nhfn,"%s.new",hfn);
      FILE* newHistoryFile=fopen(nhfn,"w");
      if(!newHistoryFile){
        //Could not open the file for writing.
        fclose(historyFile);
      }else{
        //Count the number of lines in the existing history file.
        int lines=0;
        char buffer[1024]="  ";
        while(!feof(historyFile) && strlen(buffer)>1){
          char* flag=fgets(buffer,1023,historyFile);
          if(!flag || *flag==EOF) break;
          if(strlen(buffer)>1)
            lines++;
        }
        
        //Copy the last 1000 lines of the existing history file to the new
        // history file.
        rewind(historyFile);
        for(int i=0;i<lines;i++){
          fgets(buffer,1023,historyFile);
          if(lines-i<1000)
            fprintf(newHistoryFile,"%s",buffer);
        }
        
        //Write the current command to the end of the new history file.
        fprintf(newHistoryFile,"%s\n",CommandLine);
        //And move the new history file to the old history file.
        fclose(historyFile);
        fclose(newHistoryFile);
        sprintf(buffer,"/bin/mv %s %s",
                nhfn,hfn);
        system(buffer);
      }
    }
  }
#endif //#ifdef STORE_RUN_HISTORY
  
  return argc;
}
//Do a recursive call after reading arguments out of the given file.
int commandParser::recursiveCallProcessArgs(char* filename,char* arg0,
                                            int recursionDepth){
  ReadFile* infile=new ReadFile(filename);
  int fileArgc=1;
  char** fileArgv=(char **)malloc(2*sizeof(char*));
  fileArgv[0]=arg0;
  for(;!infile->isEOF();){
    fileArgv[fileArgc]=infile->getWord();
    if(strncmp(fileArgv[fileArgc],"\\",2) && strncmp(fileArgv[fileArgc],"#",1)){
      fileArgc++;
      fileArgv=(char**)realloc(fileArgv,(fileArgc+1)*sizeof(char*));
    }
  }
  delete infile;
  tEprintf(Verbose,"Processing %i additional arguments in %s\n",
           fileArgc,filename);
  int nargs=processArgs(fileArgc,fileArgv,
                        recursionDepth+1); //recursive call
  
  for(int i=1;i<fileArgc;i++)
    delete fileArgv[i];
  free(fileArgv);
  return nargs;
}

//
//Add a new flag description, assume the strings are all volatile so we need to make
// local copies.
int commandParser::addFlagDescription(int flagType,
                                      const char* flag,const char* flagField,
                                      const char* description,
                                      int extraDescriptionCount,...){
  flagDescription* newFlag;
  newFlag=new flagDescription(flagType,flag,flagField);
  newFlag->addDescription(description);
  
  //Add any additional descriptions.
  if(extraDescriptionCount){
    va_list args;
    va_start(args,extraDescriptionCount);
    for(int i=0;i<extraDescriptionCount;i++){
      description=va_arg(args,char*);
      newFlag->addDescription(description);
    }
    va_end(args);
  }
  _flags->Add(newFlag);
  return _flags->size();
}

//Here is a method to build a usage string from the flag descriptions. This should is called
// at the beginning of process args.
char* commandParser::buildUsageString(char* argv0){
  int usLen=1024;
  assert((_usageString=(char*)malloc((usLen+1)*sizeof(char)))!=NULL,
         "buildUsageString--unable to allocate %i chars for usage string",
         usLen);
  
  sprintf(_usageString,"%s (failure at argument %%i, %%s) [-help]",argv0);
  for(int i=0;i<_flags->size();i++){
    if(!strlen((*_flags)[i]->_flag)) continue;
    
    char newString[1024];
    
    char start[2],end[2];
    switch((*_flags)[i]->_flagType){
      case FLAG_REQUIRED:
        strcpy(start," ");
        strcpy(end," ");
        break;
      case FLAG_OPTIONAL:
      case FLAG_EXPERT:
        strcpy(start,"[");
        strcpy(end,"]");
        break;
      default:
        strcpy(start,"{");
        strcpy(end,"}");
    }
    
    if(!strlen((*_flags)[i]->_fields)){
      sprintf(newString,"%s%s%s",start,(*_flags)[i]->_flag,end);
    }else{
      sprintf(newString,"%s%s %s%s",
              start,(*_flags)[i]->_flag,(*_flags)[i]->_fields,end);
    }
    
    if(strlen(_usageString)+strlen(newString) > usLen){
      usLen*=2;
      assert((_usageString=(char*)realloc(_usageString,(usLen+1)*sizeof(char)))!=NULL,
             "buildUsageString--unable to reallocate %i chars for usage string",
             usLen);
    }
    
    char* endptr=_usageString+strlen(_usageString);
    sprintf(endptr,newString);
  }
  
  return _usageString;
}

///Print the usage instructions on to stderr.
void commandParser::printUsage(int doExit,int printExpert){
  fprintf(stderr,"%s\n",_executable);
  for(int i=0;i<_flags->size();i++){
    (*_flags)[i]->printUsage(_flagTypeDescription,printExpert);
  }
  if(doExit)
    exit(0);
}
#ifdef STORE_RUN_HISTORY
///Write a history of calls.
void commandParser::printHistory(int nhistory){
  char* hfn=getenv("TDAAPS_HISTORY_FILE");
  FILE* historyFile=fopen(hfn,"r");
  assert(hfn && historyFile,
         "printHistory--No \"TDAAPS_HISTORY_FILE\" environment variable or cannot open file.");
  
  //Count the number of lines in the existing history file.
  int lines=0;
  char buffer[1024]="";
  while(!feof(historyFile)){
    fgets(buffer,1023,historyFile);
    if(strlen(buffer)>1)
      lines++;
  }
  
  //Print the last nhistory lines of the history file.
  rewind(historyFile);
  int i=0,curri=1;
  for(i=0;i<lines-nhistory;i++,fgets(buffer,1023,historyFile));
  while(!feof(historyFile)){
    fgets(buffer,1023,historyFile);
    if(strlen(buffer)>1)
      fprintf(stderr,"%04i: %s",i+curri,buffer);
    curri++;
  }
  
  exit(0);
}
///Repeat a command that is read from the history.
void commandParser::repeatCommand(int index){
  //Find the line number that corresponds to index (if index is 0 use the
  // last line).
  char* hfn=getenv("TDAAPS_HISTORY_FILE");
  FILE* historyFile=fopen(hfn,"r");
  assert(hfn && historyFile,
         "repeatCommand--No \"TDAAPS_HISTORY_FILE\" environment variable or cannot open file.");
  
  
  char buffer[1024]="",command[1024];
  if(!index){
    while(!feof(historyFile)){
      fgets(buffer,1023,historyFile);
      if(strlen(buffer)>1)
        strcpy(command,buffer);
    }
  }else{
    //Count the number of lines in the existing history file.
    int lines=0;
    while(!feof(historyFile)){
      fgets(buffer,1023,historyFile);
      if(strlen(buffer)>1)
        lines++;
    }
    
    //Find the correct index.
    rewind(historyFile);
    int i=0;
    for(i=0;i<=index;i++,fgets(buffer,1023,historyFile));
    strcpy(command,buffer);
  }
  //Copy the command to the buffer so I have a copy to safely modify.
  strcpy(buffer,command);
  
  //Now extract the arguments. Note that the executable may have a compiled on
  // parentetical comment to remove.
  char *argv[512];
  int argc=0,lastCharWhite=TRUE,argStartIndex;
  
  for(int i=0;i<strlen(command);i++){
    //Check for an open paren, if present burn charaters until a close paren.
    if(command[i]=='('){
      buffer[i]='\0';
      for(;i<strlen(command);i++){
        if(command[i]==')'){
          argStartIndex=i;
          break;
        }
      }
    }
    
    //Check for white space, if present change to an end of string in the buffer
    // version of the string.
    if(command[i]==' '){
      buffer[i]='\0';
      lastCharWhite=TRUE;
    }else if(lastCharWhite){
      argv[argc++]=&buffer[i];
      lastCharWhite=FALSE;
    }
  }
  argv[argc++]=NULL;
  
  //And run the job.
  //     execv(argv[0],&argv[1]);
  strcpy(command,"");
  for(int i=0;i<argc-1;i++)
    sprintf(&command[strlen(command)],"%s ",argv[i]);
  
  system(command);
  exit(0);
}
#endif //#ifdef STORE_RUN_HISTORY

