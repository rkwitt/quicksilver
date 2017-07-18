function [I,origin,spacing] = loadMETA2(filename, outDataType, verbose)
%
% loads 2 and 3D META images and vector fields
%
if nargin < 3
  verbose = false;
end

if nargin < 2
  outDataType = 'double';
end

dataFilename='';
dataType = '';
imSize = [0 0 0];
numChannels = 1;
origin = [0 0 0];
spacing = [1 1 1];

% parse header file
%fprintf('Loading image: %s...\n',filename);
fidh = fopen(filename);
while 1
  fline = fgetl(fidh);
  if fline == -1
    break
  end
  key = sscanf(fline,'%s %*256c');
  val = sscanf(fline,'%*s %*1c %256c');
  switch key
   case 'ObjectType'
     if verbose
       fprintf('  Object Type: %s\n', val);
     end
   case 'NDims'
     if verbose     
       fprintf('  Number of Dimensions: %s\n', val);
     end
   case 'DimSize'
     if verbose
       fprintf('  Size: %s\n', val);
     end
    imSize = str2num(val);
   case 'ElementType'
     if verbose
       fprintf('  Element Type: %s\n', val);
     end
     dataType = decideMETADataType(val);
   case 'ElementDataFile'
     if verbose
       fprintf('  DataFile: %s\n', val);
     end
    dataFilename = val;
   case 'ElementNumberOfChannels'
     if verbose
       fprintf('  Number of Channels: %s\n', val);
     end
    numChannels = str2num(val);
   case 'Offset'
     if verbose
       fprintf('  Offset: %s\n', val);
     end
    origin = str2num(val);
   case 'ElementSpacing'
     if verbose
       fprintf('  Spacing: %s\n', val);
     end
    spacing = str2num(val);
    otherwise
     if verbose
       fprintf('  Unknown Key/Val: %s/%s\n',key,val);
     end
  end
end
fclose(fidh);

% load image data
path = fileparts(filename);
fid = fopen([path filesep dataFilename],'r');
if (fid==-1)
  error(sprintf('Can''t read file: %s\n', [path filesep dataFilename]));
end
if isOctave
  [I,count] = fread(fid,prod(imSize)*numChannels,dataType);
else
  [I,count] = fread(fid,prod(imSize)*numChannels,[dataType '=>' outDataType]);
end
fclose(fid);
I=squeeze(reshape(I,[numChannels imSize]));

