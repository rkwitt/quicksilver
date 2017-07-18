function writeMETALandmark(pts,filename,dataType)
%
% writes points to META landmark data file
%
if nargin < 3
  dataType = 'MET_FLOAT';
end

% generate file names
[path,file,ext]=fileparts(filename);
fprintf('path: %s\n', path);
fprintf('file: %s\n', file);
fprintf('ext: %s\n', ext);
if strcmp(ext,'') | strcmp(ext,'.mhd')  
  if length(path) > 0
    headerFilename = [path filesep file '.mhd'];
    dataFilename = [path filesep file '.pts'];
    dataFilenameHeaderRelative = [file '.pts'];
  else
    headerFilename = [file '.mhd'];
    dataFilename = [file '.pts'];
    dataFilenameHeaderRelative = [file '.pts'];
  end
else
  headerFilename = [filename '.mhd'];
  dataFilename = [filename '.pts'];  
  dataFilenameHeaderRelative = [path filesep file '.' ext '.pts'];
end

fprintf('header filename: %s\n', headerFilename);
fprintf('data filename: %s\n', dataFilename);
fprintf('data relative filename: %s\n', dataFilenameHeaderRelative);

[d npts] = size(pts);

% write header file
fprintf('writing header file: %s...\n',headerFilename);
fid = fopen(headerFilename,'w');
fprintf(fid,'ObjectType = %s\n','Image');
fprintf(fid,'PointDim = %d\n',d);
fprintf(fid,'NPoints = %d\n',npts);
fprintf(fid,'Points = %s\n',dataFilenameHeaderRelative);
fclose(fid);

fprintf('writing data file: %s...\n',dataFilename);
dlmwrite(dataFilename,pts,' ');


