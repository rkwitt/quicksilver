function [A,t] = loadAffinePLUNC(filename)
%
% load as matrix and translation
%
% matrix should act on left of vector
%
fid = fopen(filename,'r');
topLine = fgetl(fid);
if ~strcmpi(topLine,'matrix')
  error('Invalid PLUNC matrix file.');
end
c1 = str2num(fgetl(fid));
c2 = str2num(fgetl(fid));
c3 = str2num(fgetl(fid));
t = str2num(fgetl(fid));
A = [c1(1:3)' c2(1:3)' c3(1:3)'];
t = t(1:3)';
fclose(fid);