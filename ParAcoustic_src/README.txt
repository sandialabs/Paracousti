This directory contains the source files for ParAcousti separated into three subdirectories.  The 'acoustic' directory contains
the Makefile and is where you should be to compile the entire code.  Once in the acoustic/ directory, modify the Makefile_local
file to match your system's specific needs for the compilers and options and locations of include and libraries, if necessary,
for netcdf and mpi.  In the Makefile_local file, you will also specify the location of source and binary files.
The Makefile_local file contains some prompts for SRCHOME and BINDIR that should be modified to match your specific setup.
In general, the hope is that you should not have to modify Makefile itself, but instead all specifics should go into
Makefile_local.  However, it is possible that Makefile itself will have to be modified to get it to compile.  Once 
Makefile_local is set for your specific system, just run 'make' on the command line within the acoustic/ directory.  It will
compile all of the source files and create an executable called 'ParAcousti' if successful.  Netcdf 
(https://www.unidata.ucar.edu/software/netcdf/) at least v. 4.0 includes and libraries and openMPI (https://www.open-mpi.org)
at least version 1.8 are required for compilation.
