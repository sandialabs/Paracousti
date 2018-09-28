function writeSgfdModel(filename,x,y,z,t,varargin)
%function writeSgfdModel(filename,x,y,z,t,varargin)
%Write a Symons style netcdf file for sgfd modeling. Can be used for either
% elastic or acoustic models.
%Neill Symons; Sandia National Laboratories; 4/24/03
%Updated by Leiph Preston, Sandia National Laboratories: 5/22/14, to use
%the internal Matlab netcdf calls.
%Arguments: filename--file to write
%           x, y, z--vectors with spatial position of node centers
%           t--time vector.
%Optional Arguments: vp, vs, rho--these are done optionally for some
%                      flexability in what is actually defined.
%                    comment--add a comment to the file
%                    noclobber--add variables to an exisiting file, useful
%                     for large models.
%                    COrder--input multidimensional arrays are expected to
%                     have the following ordering:
%                         size(vp,1)==length(z)
%                         size(vp,2)==length(y)
%                         size(vp,3)==length(x)
%                    MatlabOrder--input multidimensional arrays are
%                     expected to have the following format:
%                         size(vp,1)==length(x)
%                         size(vp,2)==length(y)
%                         size(vp,3)==length(z)
%                     MatlabOrder is the default and should require less 
%                     Matlab processing/memory which could be important for
%                     large models.  However, the netCDF file stores the
%                     data in COrder.  All the Matlab routines read and
%                     write in MatlabOrder.
%

%Check varargin for modifiers to the default arguments.
corder = 0; %default order of input is MatlabOrder
i=1;
while i<=length(varargin)
  currArg=varargin{i};
  i=i+1;
  argType=whos('currArg');
  if ~strcmp(argType.class,'char')
    error('Optional argument %i, type %s must be char',i,argType.class);
  end
  
  switch lower(currArg)
    case {'velocity' 'vel' 'vp'}
      vp=varargin{i};
      i=i+1;
    case {'alpha' 'c'}
      C=varargin{i};
      i=i+1;
    case {'vs' 'beta'}
      vs=varargin{i};
      i=i+1;
    case {'density' 'rho'}
      rho=varargin{i};
      i=i+1;
    case 'var3d'
      var3dNam=varargin{i+0};
      var3dDat=varargin{i+1};
      i=i+2;
    case {'var1' 'var'}
      var1Nam=varargin{i+0};
      var1Dat=varargin{i+1};
      i=i+2;
      
    case {'xzmodel' 'xzm' 'xz'}
      vpXZ=varargin{i};i=i+1;
      vsXZ=varargin{i};i=i+1;
      rhoXZ=varargin{i};i=i+1;
      
    case {'indexmodel' 'index'}
      vpI=varargin{i+0};
      vsI=varargin{i+1};
      rhoI=varargin{i+2};
      indicies=varargin{i+3};
      i=i+4;
    case {'layermodel' 'layer' 'onedmodel' 'oned' '1d'}
      vpL=varargin{i+0};
      vsL=varargin{i+1};
      rhoL=varargin{i+2};
      i=i+3;
      
    case 'wind'
      vx=varargin{i+0};
      vy=varargin{i+1};
      vz=varargin{i+2};
      i=i+3;
      
    case 'slice'
      sliceName=varargin{i};
      if ~isa(sliceName,'cell')
        sliceName={sliceName};
      end
      
      sliceTimes=varargin{i+1};
      slicePos=varargin{i+2};
      if ~isa(slicePos,'cell')
        slicePos={slicePos};
      end
      
      if length(slicePos)~=length(sliceName)
        error('Must have a position (%i) for each slice type (%i).',...
          length(slicePos),length(sliceName));
      end
      
      i=i+3;
      if exist('slices','var')
        sliceStart=length(slices);
      else
        sliceStart=0;
      end
      
      numSlices=sliceStart+length(sliceTimes)*length(sliceName);
      for jjj=1:length(sliceName)
        slices{jjj+sliceStart}={sliceName{jjj},sliceTimes,slicePos{jjj}};
      end
      clear sliceName sliceTimes slicePos sliceStart;
    case 'receivers'
      receiverType=varargin{i};
      receiverAmp=varargin{i+1};
      receiverX=varargin{i+2};
      receiverY=varargin{i+3};
      receiverZ=varargin{i+4};
      receiverBx=varargin{i+5};
      receiverBy=varargin{i+6};
      receiverBz=varargin{i+7};
      receiverIntegrate=varargin{i+8};
      i=i+9;
    case '3c'
      tReceiverX=varargin{i+0};
      tReceiverY=varargin{i+1};
      tReceiverZ=varargin{i+2};
      i=i+3;
      
      receiverX=zeros(3*length(tReceiverX),1);
      receiverX(1:3:end)=tReceiverX;
      receiverX(2:3:end)=tReceiverX;
      receiverX(3:3:end)=tReceiverX;
      
      receiverY=zeros(3*length(tReceiverY),1);
      receiverY(1:3:end)=tReceiverY;
      receiverY(2:3:end)=tReceiverY;
      receiverY(3:3:end)=tReceiverY;
      
      receiverZ=zeros(3*length(tReceiverZ),1);
      receiverZ(1:3:end)=tReceiverZ;
      receiverZ(2:3:end)=tReceiverZ;
      receiverZ(3:3:end)=tReceiverZ;
            
      receiverType=1+0*receiverX;
      receiverAmp=1+0*receiverX;
      
      receiverBx=0*receiverX;
      receiverBx(1:3:end)=1;
      
      receiverBy=0*receiverX;
      receiverBy(2:3:end)=1;
      
      receiverBz=0*receiverX;
      receiverBz(3:3:end)=1;
      
      receiverIntegrate=0*receiverX;
      
    case {'pressurereceivers' 'pressure' 'press' 'pres'}
      receiverX=varargin{i+0};
      receiverY=varargin{i+1};
      receiverZ=varargin{i+2};
      i=i+3;
      
      receiverType=2+0*receiverX;
      receiverAmp=1+0*receiverX;
      
      receiverBx=0*receiverX;
      receiverBy=0*receiverX;
      receiverBz=0*receiverX;
      receiverIntegrate=0*receiverX;
      
    case {'source' 'pressuresource' 'explosion'}
      pressureSource=varargin{i};
      sourceWaveform=varargin{i+1};
      i=i+2;
      
    case 'comment'
      comment=varargin{i};
      i=i+1;
    case 'history'
      history=varargin{i};
      i=i+1;
      
    case {'noclobber' 'addvar' 'add'}
      %Add a new time plane to an existing file.
      noclobber=1;
      
    case 'corder'
      corder=1;
      
    case 'matlaborder'
      corder=0;
      
    otherwise
      error('Unknown option %s',currArg);
  end
end

%Check that the sizes match up.
NX=length(x);
NY=length(y);
NZ=length(z);
NT=length(t);

if corder
  if exist('vp','var') && (...
      size(vp,1)~=NZ || size(vp,2)~=NY || size(vp,3)~=NX)
    fprintf('COrder requires: (try MatlabOrder option?)\n');
    fprintf(' size(vp,1)==length(z)\n');
    fprintf(' size(vp,2)==length(y)\n');
    fprintf(' size(vp,3)==length(x)\n');
    error('Can not write file');
  end
  if exist('C','var') && (...
      size(C,1)~=NZ || size(C,2)~=NY || size(C,3)~=NX)
    fprintf('COrder requires: (try MatlabOrder option?)\n');
    fprintf(' size(C,1)==length(z)\n');
    fprintf(' size(C,2)==length(y)\n');
    fprintf(' size(C,3)==length(x)\n');
    error('Can not write file');
  end

  if exist('vs','var') && (...
      size(vs,1)~=NZ || size(vs,2)~=NY || size(vs,3)~=NX)
    fprintf('COrder requires: (try MatlabOrder option?)\n');
    fprintf(' size(vs,1)==length(z)\n');
    fprintf(' size(vs,2)==length(y)\n');
    fprintf(' size(vs,3)==length(x)\n');
    error('Can not write file');
  end

  if exist('rho','var') && (...
      size(rho,1)~=NZ || size(rho,2)~=NY || size(rho,3)~=NX)
    fprintf('COrder requires: (try MatlabOrder option?)\n');
    fprintf(' size(rho,1)==length(z)\n');
    fprintf(' size(rho,2)==length(y)\n');
    fprintf(' size(rho,3)==length(x)\n');
    error('Can not write file');
  end
  if exist('var3dDat','var') && (...
      size(var3dDat,1)~=NZ || size(var3dDat,2)~=NY || size(var3dDat,3)~=NX)
    fprintf('COrder requires: (try MatlabOrder option?)\n');
    fprintf(' size(3dVar,1)==length(z)\n');
    fprintf(' size(3dVar,2)==length(y)\n');
    fprintf(' size(3dVar,3)==length(x)\n');
    error('Can not write file');
  end

  if exist('vx','var') && ...
      (size(vx,1)~=NZ || size(vx,2)~=NY || size(vx,3)~=NX) && ...
      (size(vy,1)~=NZ || size(vy,2)~=NY || size(vy,3)~=NX) && ...
      (size(vz,1)~=NZ || size(vz,2)~=NY || size(vz,3)~=NX)
    fprintf('COrder requires: (try MatlabOrder option?)\n');
    fprintf(' size(vx,1)==size(vy,1)==size(vz,1)==length(z)\n');
    fprintf(' size(vx,2)==size(vy,2)==size(vz,2)==length(y)\n');
    fprintf(' size(vx,3)==size(vy,3)==size(vz,3)==length(x)\n');
    error('Can not write file');
  end
else
  if exist('vp','var') && (...
      size(vp,1)~=NX || size(vp,2)~=NY || size(vp,3)~=NZ)
    fprintf('matlabOrder requires: (try Corder option?)\n');
    fprintf(' size(vp,1)==length(x)\n');
    fprintf(' size(vp,2)==length(y)\n');
    fprintf(' size(vp,3)==length(z)\n');
    error('Can not write file');
  end
  if exist('C','var') && (...
      size(C,1)~=NX || size(C,2)~=NY || size(C,3)~=NZ)
    fprintf('matlabOrder requires: (try Corder option?)\n');
    fprintf(' size(C,1)==length(x)\n');
    fprintf(' size(C,2)==length(y)\n');
    fprintf(' size(C,3)==length(z)\n');
    error('Can not write file');
  end

  if exist('vs','var') && (...
      size(vs,1)~=NX || size(vs,2)~=NY || size(vs,3)~=NZ)
    fprintf('matlabOrder requires: (try Corder option?)\n');
    fprintf(' size(vs,1)==length(x)\n');
    fprintf(' size(vs,2)==length(y)\n');
    fprintf(' size(vs,3)==length(z)\n');
    error('Can not write file');
  end

  if exist('rho','var') && (...
      size(rho,1)~=NX || size(rho,2)~=NY || size(rho,3)~=NZ)
    fprintf('matlabOrder requires: (try Corder option?)\n');
    fprintf(' size(rho,1)==length(x)\n');
    fprintf(' size(rho,2)==length(y)\n');
    fprintf(' size(rho,3)==length(z)\n');
    error('Can not write file');
  end
  if exist('var3dDat','var') && (...
      size(var3dDat,1)~=NX || size(var3dDat,2)~=NY || size(var3dDat,3)~=NZ)
    fprintf('matlabOrder requires: (try Corder option?)\n');
    fprintf(' size(3dVar,1)==length(x)\n');
    fprintf(' size(3dVar,2)==length(y)\n');
    fprintf(' size(3dVar,3)==length(z)\n');
    error('Can not write file');
  end

  if exist('vx','var') && ...
      (size(vx,1)~=NX || size(vx,2)~=NY || size(vx,3)~=NZ) && ...
      (size(vy,1)~=NX || size(vy,2)~=NY || size(vy,3)~=NZ) && ...
      (size(vz,1)~=NX || size(vz,2)~=NY || size(vz,3)~=NZ)
    fprintf('matlabOrder requires: (try Corder option?)\n');
    fprintf(' size(vx,1)==size(vy,1)==size(vz,1)==length(x)\n');
    fprintf(' size(vx,2)==size(vy,2)==size(vz,2)==length(y)\n');
    fprintf(' size(vx,3)==size(vy,3)==size(vz,3)==length(z)\n');
    error('Can not write file');
  end
end

%Open the file.  This just opens/creates the files and then closes it.
%This was necessary so that Matlab would not search its entire path for the
%later commands and possibly modify a same-named file elsewhere on the
%search path.
if exist('noclobber','var') && noclobber
  out=netcdf.open(filename,'write');
else
  out=netcdf.create(filename,'clobber');
  
  ncglob = netcdf.getConstant('NC_GLOBAL');
  %Set some global attribute describing how this file was created.
  netcdf.putAtt(out,ncglob,'title','Staggered Grid Finite-Difference Model Input File');
  if exist('comment','var')
    netcdf.putAtt(out,ncglob,'comment',comment);
  end
  if exist('history','var')
    netcdf.putAtt(out,ncglob,'history',history);
  else
    netcdf.putAtt(out,ncglob,'history','Created with matlab writeSgfdModel.m');
  end
  
  %Set the dimensions.
  netcdf.defDim(out,'numCoord',4);
  netcdf.defDim(out,'NX',NX);
  netcdf.defDim(out,'NY',NY);
  netcdf.defDim(out,'NZ',NZ);
  netcdf.defDim(out,'NT',NT);
  netcdf.close(out);
  
  %Define and fill the increment variables.
  nccreate(filename,'minima','Dimensions',{'numCoord'},'DataType','single');
  nccreate(filename,'increments','Dimensions',{'numCoord'},'DataType','single');
  ncwrite(filename,'minima',[x(1) y(1) z(1) t(1)]);
  ncwrite(filename,'increments',[x(2)-x(1) y(2)-y(1) z(2)-z(1) t(2)-t(1)]);
  
  %Define and fill the position variables.
  nccreate(filename,'x','Dimensions',{'NX'},'DataType','single');
  nccreate(filename,'y','Dimensions',{'NY'},'DataType','single');
  nccreate(filename,'z','Dimensions',{'NZ'},'DataType','single');
  nccreate(filename,'time','Dimensions',{'NT'},'DataType','single');
  ncwrite(filename,'x',x);
  ncwrite(filename,'y',y);
  ncwrite(filename,'z',z);
  ncwrite(filename,'time',t);
end

%Define and fill extra scalar variables.
if exist('var1Dat','var')
  nccreate(filename,var1Nam,'DataType','single');
  ncwrite(filename,var1Nam,var1Dat);
end

%Write the defined variables.
if exist('vpL','var') && exist('vsL','var') && exist('rhoL','var')
  nccreate(filename','oneDModelVp','Dimensions',{'NZ'},'DataType','single');
  nccreate(filename','oneDModelVs','Dimensions',{'NZ'},'DataType','single');
  nccreate(filename','oneDModelRho','Dimensions',{'NZ'},'DataType','single');

  ncwrite(filename,'oneDModelVp',vpL);
  ncwrite(filename,'oneDModelVs',vsL);
  ncwrite(filename,'oneDModelRho',rhoL);
end

if exist('indicies','var')
  try
    info = ncinfo(filename,'indexModelValues');
    currIndex = 1+info.Size(2);
  catch
    nccreate(filename,'indexModelValues','Dimensions',{'three',3,'indexModelDim',inf},'DataType','single');
    currIndex = 1;
  end
  ncwrite(filename,'indexModelValues',[vpI vsI rhoI]',[1 currIndex]);
  
  if length(indicies)
    dimName=sprintf('indexModel%iDim',currIndex);
    varName=sprintf('indexModel%iIndicies',currIndex);
    nccreate(filename,varName,'Dimensions',{dimName,size(indicies,1)},'DataType','int');
    ncwrite(filename,varName,indicies);
  end
else
  if exist('vpXZ','var')
    if corder,vpXZ=vpXZ';end
    nccreate(filename,'vpXZ','Dimensions',{'NX','NZ'},'DataType','single');
    ncwrite(filename,'vpXZ',vpXZ);
  end
  if exist('vsXZ','var')
    if corder,vsXZ=vsXZ';end
    nccreate(filename,'vsXZ','Dimensions',{'NX','NZ'},'DataType','single');
    ncwrite(filename,'vsXZ',vsXZ);
  end
  if exist('rhoXZ','var')
    if corder,rhoXZ=rhoXZ';end
    nccreate(filename,'rhoXZ','Dimensions',{'NX','NZ'},'DataType','single');
    ncwrite(filename,'rhoXZ',rhoXZ);
  end
  if exist('vp','var')
    if corder,vp=permute(vp,[3 2 1]);end
    nccreate(filename,'vp','Dimensions',{'NX','NY','NZ'},'DataType','single');
    ncwrite(filename,'vp',vp);
  end
  if exist('C','var')
    if corder,C=permute(C,[3 2 1]);end
    nccreate(filename,'c','Dimensions',{'NX','NY','NZ'},'DataType','single');
    ncwrite(filename,'c',C);
  end
  if exist('vs','var')
    if corder,vs=permute(vs,[3 2 1]);end
    nccreate(filename,'vs','Dimensions',{'NX','NY','NZ'},'DataType','single');
    ncwrite(filename,'vs',vs);
  end
  if exist('rho','var')
    if corder,rho=permute(rho,[3 2 1]);end
    nccreate(filename,'rho','Dimensions',{'NX','NY','NZ'},'DataType','single');
    ncwrite(filename,'rho',rho);
  end

  if exist('vx','var')
    if corder,vx=permute(vx,[3 2 1]);end
    nccreate(filename,'vx','Dimensions',{'NX','NY','NZ'},'DataType','single');
    ncwrite(filename,'vx',vx);
  end
  if exist('vy','var')
    if corder,vy=permute(vy,[3 2 1]);end
    nccreate(filename,'vy','Dimensions',{'NX','NY','NZ'},'DataType','single');
    ncwrite(filename,'vy',vy);
  end
  if exist('vz','var')
    if corder,vz=permute(vz,[3 2 1]);end
    nccreate(filename,'vz','Dimensions',{'NX','NY','NZ'},'DataType','single');
    ncwrite(filename,'vz',vz);
  end
  
  if exist('var3dDat','var')
    if corder,var3dDat=permute(var3dDat,[3 2 1]);end
    nccreate(filename,var3dNam,'Dimensions',{'NX','NY','NZ'},'DataType','single');
    ncwrite(filename,var3dNam,var3dDat);
  end
end


%
%Write extra defined stuff.
%

%Slices.
if exist('slices','var')
  nccreate(filename,'sliceTime','Dimensions',{'numSlices',numSlices},'DataType','single');
  nccreate(filename,'sliceComp','Dimensions',{'numSlices'},'DataType','int');
  nccreate(filename,'slicePlane','Dimensions',{'numSlices'},'DataType','int');
  nccreate(filename,'sliceCoord','Dimensions',{'numSlices'},'DataType','single');
     
  startSlice=0;
  for i=1:length(slices)
    currSlice=slices{i};
    sliceName=currSlice{1};
    sliceTimes=currSlice{2};
    sliceCoord=currSlice{3};
    
    ncwrite(filename,'sliceTime',sliceTimes,startSlice+1);
    for i=1:length(sliceTimes)
      ncwrite(filename,'sliceCoord',sliceCoord,startSlice+i);
      switch lower(sliceName(1:2))
        case 'yz'
          ncwrite(filename,'slicePlane',1,startSlice+i);
          if x(1)>sliceCoord || sliceCoord>x(end)
            error('Slice %s: out of bounds %f<%f<%f',...
              sliceName,x(1),sliceCoord,x(end));
          end
        case 'xz'
          ncwrite(filename,'slicePlane',2,startSlice+i);
          if y(1)>sliceCoord || sliceCoord>y(end)
            error('Slice %s: out of bounds %f<%f<%f',...
              sliceName,y(1),sliceCoord,y(end));
          end
        case 'xy'
          ncwrite(filename,'slicePlane',3,startSlice+i);
          if z(1)>sliceCoord || sliceCoord>z(end)
            error('Slice %s: out of bounds %f<%f<%f',...
              sliceName,z(1),sliceCoord,z(end));
          end
        otherwise
          error('Unknown slice plane %s',sliceName(1:2));
      end
    
      switch lower(sliceName(3:end))
        case 'vx'
          ncwrite(filename,'sliceComp',1,startSlice+i);
        case 'vy'
          ncwrite(filename,'sliceComp',2,startSlice+i);
        case 'vz'
          ncwrite(filename,'sliceComp',3,startSlice+i);
        case 'pressure'
          ncwrite(filename,'sliceComp',4,startSlice+i);
        otherwise
          error('Unknown slice component %s',sliceName(3:end));
      end
      
    end
    startSlice=startSlice+length(sliceTimes);
  end
end

%Receivers.
if exist('receiverType','var')
  nccreate(filename,'receiverType','Dimensions',{numReceivers,length(receiverType)},'DataType','int');
  nccreate(filename,'receiverAmp','Dimensions',{numReceivers},'DataType','float');
  nccreate(filename,'receiverX','Dimensions',{numReceivers},'DataType','float');
  nccreate(filename,'receiverY','Dimensions',{numReceivers},'DataType','float');
  nccreate(filename,'receiverZ','Dimensions',{numReceivers},'DataType','float');
  nccreate(filename,'receiverBx','Dimensions',{numReceivers},'DataType','float');
  nccreate(filename,'receiverBy','Dimensions',{numReceivers},'DataType','float');
  nccreate(filename,'receiverBz','Dimensions',{numReceivers},'DataType','float');
  nccreate(filename,'receiverIntegrate','Dimensions',{numReceivers},'DataType','int');

  ncwrite(filename,'receiverType',receiverType);
  ncwrite(filename,'receiverAmp',receiverAmp);
  ncwrite(filename,'receiverX',receiverX);
  ncwrite(filename,'receiverY',receiverY);
  ncwrite(filename,'receiverZ',receiverZ);
  ncwrite(filename,'receiverBx',receiverBx);
  ncwrite(filename,'receiverBy',receiverBy);
  ncwrite(filename,'receiverBz',receiverBz);
  ncwrite(filename,'receiverIntegrate',receiverIntegrate);
end

%Sources.
if exist('pressureSource','var')
  nccreate(filename,'mSourcesXs','Dimensions',{'numMSources',1},'DataType','single');
  nccreate(filename,'mSourcesYs','Dimensions',{'numMSources'},'DataType','single');
  nccreate(filename,'mSourcesZs','Dimensions',{'numMSources'},'DataType','single');
    
  nccreate(filename,'mSourcesSamp','Dimensions',{'numMSources'},'DataType','single');
  
  nccreate(filename,'mSourcesXxS','Dimensions',{'numMSources'},'DataType','single');
  nccreate(filename,'mSourcesYyS','Dimensions',{'numMSources'},'DataType','single');
  nccreate(filename,'mSourcesZzS','Dimensions',{'numMSources'},'DataType','single');
  
  nccreate(filename,'mSourcesXyS','Dimensions',{'numMSources'},'DataType','single');
  nccreate(filename,'mSourcesXzS','Dimensions',{'numMSources'},'DataType','single');
  nccreate(filename,'mSourcesYzS','Dimensions',{'numMSources'},'DataType','single');
  
  nccreate(filename,'mSourcesXyA','Dimensions',{'numMSources'},'DataType','single');
  nccreate(filename,'mSourcesXzA','Dimensions',{'numMSources'},'DataType','single');
  nccreate(filename,'mSourcesYzA','Dimensions',{'numMSources'},'DataType','single');
  
  ncwrite(filename,'mSourcesXs',pressureSource(1));
  ncwrite(filename,'mSourcesYs',pressureSource(2));
  ncwrite(filename,'mSourcesZs',pressureSource(3));

  ncwrite(filename,'mSourcesSamp',1);
  
  ncwrite(filename,'mSourcesXxS',1);
  ncwrite(filename,'mSourcesYyS',1);
  ncwrite(filename,'mSourcesZzS',1);
  
  ncwrite(filename,'mSourcesXyS',0);
  ncwrite(filename,'mSourcesXzS',0);
  ncwrite(filename,'mSourcesYzS',0);
  
  ncwrite(filename,'mSourcesXyA',0);
  ncwrite(filename,'mSourcesXzA',0);
  ncwrite(filename,'mSourcesYzA',0);
  
  nccreate(filename,'mSourcesData','Dimensions',{'NT','numMSources'},'DataType','single');
  ncwrite(filename,'mSourcesData',sourceWaveform);
end  
