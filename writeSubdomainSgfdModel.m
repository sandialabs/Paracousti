function [filenames]=writeSubdomainSgfdModel(filename,x,y,z,t,domains,varargin)
%function writeSubdomainSgfdModel(filename,x,y,z,t,varargin)
%Write a Symons style netcdf file for sgfd modeling. Can be used for either
% elastic or acoustic models. This version will create smaller files that a
% portion of the full model.
% Example:
% > A=@(x,y,z,f) f*ones(length(x),length(y),length(z));
% > writeSubdomainSgfdModel('TestSD.cdf',-100:100,-100:100,-100:100,0:0.0001:1,...
%     [1 201 1 201 1 100;1 201 1 201 101 201],...
%     'vp','func',A,{350},'rho','func',A,{1.2});
% This creates a model consisting of 2 subdomains split along z. Variables
% vp and rho are defined and filled with values 350 and 1.2 respectively.
%
%Neill Symons; Sandia National Laboratories; 1/21/09
%Updated by Leiph Preston, Sandia National Laboratories: 11/13/14, to use
%the internal Matlab netcdf calls.
%Arguments: filename--file to write
%           x, y, z--vectors with spatial position of node centers for
%           entire model
%           t--time vector.
%           domains--breakdown of files by indices.  This is a Nx6 array,
%           with N being the number of domains (this is independent of the
%           number of domains used in mpi).  Each row of this array gives
%           [minXi maxXi minYi maxYi minZi maxZi] for a single domain.  The
%           model limits for this domain will be x(minXi:maxXi),
%           y(minYi:maxYi), z(minZi:maxZi).
%
%Medium Parameter Arguments:
%  These come in groups of 4 arguments:
%  1) medium parameter name: either 'vp' or 'rho'
%  2) how it is defined: only option now is 'func'
%  3) function handle
%  4) arguments to function following required x, y, z arguments in a cell
%  array.
%  In the above example, vp and rho use the same function handle, A, but
%  that is not required.  The function the handle refers to can have any
%  number of arguments, but the first 3 must be x, y, z.  All the remaining
%  arguments are placed in a cell array.  Use {} if there are none.  As
%  another example, suppose you have defined a Matlab function called
%  myfuncVp(x,y,z) for vp, then 1) would be 'vp', 2) is 'func', 3) is
%  @myfuncVp, 4) is {}.
%
%Optional Arguments:
%                    COrder--the function handle is expected to
%                     output in the following ordering:
%                         size(output,1)==length(z)
%                         size(output,2)==length(y)
%                         size(output,3)==length(x)
%                    MatlabOrder--the function handle is
%                     expected to output in the following order:
%                         size(output,1)==length(x)
%                         size(output,2)==length(y)
%                         size(output,3)==length(z)
%                     MatlabOrder is the default and should require less 
%                     Matlab processing/memory which could be important for
%                     large models.  However, the netCDF file stores the
%                     data in COrder.  All the Matlab routines read and
%                     write in MatlabOrder.


%% Check varargin for modifiers to the default arguments.
%   Default values
fileVars={};
varfunc=struct('func',{},'args',{});
nVarfunc=0;
corder=0;

%   Process the variable inputs.
i=1;
while i<=length(varargin)
  currArg=varargin{i};
  i=i+1;
  argType=whos('currArg');
  if ~strcmp(argType.class,'char')
    error('Optional argument %i, type %s must be char',i,argType.class);
  end
  
  switch lower(currArg)
    case {'varfunc' 'vf' 'function' 'func'}
      nVarfunc=nVarfunc+1;
      varfunc(nVarfunc).func=varargin{i};i=i+1;
      varfunc(nVarfunc).args=varargin{i};i=i+1;
    case 'corder'
      corder=1;   
    case 'matlaborder'
      corder=0;
    otherwise
      fileVars=cat(2,fileVars,{currArg});
  end
end

%% Define useful variables.
nDomain=size(domains,1);
filenameRoot=strrep(filename,'.cdf','');

%% Write the base model, then open for writing more stuff.
writeSgfdModel(filename,x,y,z,t);
out=netcdf.open(filename,'write');
netcdf.reDef(out);

%% Define the new variables.
netcdf.defDim(out,'nSubdomains',nDomain);
netcdf.defDim(out,'six',6);

filenames=cell(1,nDomain);
filenames{1}=sprintf('%s/%s_%i_%i_%i_%i_%i_%i.cdf',filenameRoot,filenameRoot,...
  domains(1,1)-1,domains(1,2),domains(1,3)-1,domains(1,4),domains(1,5)-1,domains(1,6));
filenameLen=20+length(filenames{1});
netcdf.defDim(out,'subdomainFileNameLength',filenameLen);
netcdf.close(out);

nccreate(filename,'subdomainLimit','Dimensions',{'six','nSubdomains'},'DataType','int32');
ncwrite(filename,'subdomainLimit',[domains(:,1)-1 domains(:,2) domains(:,3)-1 domains(:,4) domains(:,5)-1 domains(:,6)]');

nccreate(filename,'subdomainFilename','Dimensions',{'subdomainFileNameLength','nSubdomains'},'DataType','char');

for i=1:nDomain
  filenames{i}=sprintf('%s/%s_%i_%i_%i_%i_%i_%i.cdf',filenameRoot,filenameRoot,...
    domains(i,1)-1,domains(i,2),domains(i,3)-1,domains(i,4),domains(i,5)-1,domains(i,6));
  cfnL=length(filenames{i});
  ncwrite(filename,'subdomainFilename',filenames{i}',[1 i]);
end

%% Now create each of the files; define variables if fileVars has been set.
[s,mess]=mkdir(filenameRoot);
if s~=1
  error('Unable to create directory %s; %s',filenameRoot,mess);
end

wb=waitbar(0,sprintf('Filling %i Domains',nDomain));
for i=1:nDomain
  writeSgfdModel(filenames{i},x(domains(i,1):domains(i,2)),...
    y(domains(i,3):domains(i,4)),z(domains(i,5):domains(i,6)),t);
  if ~isempty(fileVars)
    for j=1:length(fileVars)
      nccreate(filenames{i},fileVars{j},'Dimensions',{'NX','NY','NZ'},'DataType','single');
      if nVarfunc>=j
        if corder
          ncwrite(filenames{i},fileVars{j},permute(varfunc(j).func(...
            x(domains(i,1):domains(i,2)),y(domains(i,3):domains(i,4)),z(domains(i,5):domains(i,6)),...
            varfunc(j).args{:}),[3 2 1]));
        else
          ncwrite(filenames{i},fileVars{j},varfunc(j).func(...
            x(domains(i,1):domains(i,2)),y(domains(i,3):domains(i,4)),z(domains(i,5):domains(i,6)),...
            varfunc(j).args{:}));
        end
      end
    end
  end
  waitbar(i/nDomain,wb);
end
try
  close(wb);
catch
end
