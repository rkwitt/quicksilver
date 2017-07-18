function writeMETA(im,filename,dataType,origin,spacing, extraTags)
%
% writes 3D META images and vector fields
%
if nargin < 6
  extraTags = {};
end
if nargin < 5
  spacing = ones(1,ndims(im));
end
if nargin < 4
  origin = zeros(1,ndims(im));
end
if nargin < 3
  dataType = 'MET_FLOAT';
end

[path,file,ext]=fileparts(filename);
%fprintf('path: %s\n', path);
%fprintf('file: %s\n', file);
%fprintf('ext: %s\n', ext);
if strcmp(ext,'') | strcmp(ext,'.mhd')  
  if length(path) > 0
    headerFilename = [path filesep file '.mhd'];
    dataFilename = [path filesep file '.raw'];
    dataFilenameHeaderRelative = [file '.raw'];
  else
    headerFilename = [file '.mhd'];
    dataFilename = [file '.raw'];
    dataFilenameHeaderRelative = [file '.raw'];
  end
else
  headerFilename = [filename '.mhd'];
  dataFilename = [filename '.raw'];  
  dataFilenameHeaderRelative = [path filesep file '.' ext '.raw'];
end

%fprintf('header filename: %s\n', headerFilename);
%fprintf('data filename: %s\n', dataFilename);
%fprintf('data relative filename: %s\n', dataFilenameHeaderRelative);

field = false;
nImDims = ndims(im);
imSize = size(im);
numChannels = 1;
if (ndims(im) == 4)
  field = true;
  nImDims = 3;
  numChannels = 3;
  imSize = imSize(2:end);
end

% write header file
%fprintf('writing header file: %s...\n',headerFilename);
fid = fopen(headerFilename,'w');
fprintf(fid,'ObjectType = %s\n','Image');
if field
  fprintf(fid,'ObjectSubType = %s\n','HField');
end
fprintf(fid,'NDims = %d\n',nImDims);
fprintf(fid,'BinaryData = %s\n','True');
fprintf(fid,'BinaryDataByteOrderMSB = %s\n','False');
fprintf(fid,'Offset = '); fprintf(fid,'%g ',origin); fprintf(fid,'\n');
fprintf(fid,'ElementSpacing = '); fprintf(fid,'%g ',spacing); fprintf(fid,'\n');
fprintf(fid,'DimSize = '); fprintf(fid,'%g ',imSize); fprintf(fid,'\n');
fprintf(fid,'ElementNumberOfChannels = %d\n',numChannels);
fprintf(fid,'ElementType = %s\n',dataType);
fprintf(fid,'ElementDataFile = %s\n',dataFilenameHeaderRelative);
for tagidx = 1:2:length(extraTags)
  fprintf(fid,'%s = %s\n', extraTags{tagidx}, extraTags{tagidx+1});
end

fclose(fid);

fid = fopen(dataFilename,'w');
%fprintf('writing data file: %s...\n',dataFilename);
count = fwrite(fid,im,decideMETADataType(dataType));
fclose(fid);


