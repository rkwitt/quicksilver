function IS = downsampleByN(I, N, sigma)

if nargin < 3
  sigma = N;
end

oldDims = size(I);
newDims = ceil(oldDims/N);

% smooth image
ISmooth = filterImage(I,'gaussian',sigma);

% downsample
switch ndims(I)
  case 3,
    IS(1:newDims(1),1:newDims(2),1:newDims(3)) = ...
      ISmooth(1:N:oldDims(1),1:N:oldDims(2),1:N:oldDims(3));
  case 2,
    IS(1:newDims(1),1:newDims(2)) = ISmooth(1:N:oldDims(1),1:N:oldDims(2));    
end