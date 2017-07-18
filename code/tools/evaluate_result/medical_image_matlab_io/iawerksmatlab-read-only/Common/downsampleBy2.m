function IS = downsampleBy2(I, sigma)

oldSize = size(I);
newSize = ceil(oldSize/2);

% smooth image
ISmooth = filterImage(I,'gaussian',sigma);

% downsample
switch ndims(I)
  case 3,
    IS(1:newDims(1),1:newDims(2),1:newDims(3)) = ...
      ISmooth(1:2:oldSize(1),1:2:oldSize(2),1:2:oldSize(3));
  case 2,
    IS(1:newDims(1),1:newDims(2)) = ISmooth(1:2:oldSize(1),1:2:oldSize(2));    
end