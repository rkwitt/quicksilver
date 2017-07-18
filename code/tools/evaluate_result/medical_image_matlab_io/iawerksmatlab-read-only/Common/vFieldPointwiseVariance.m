function v = vFieldPointwiseVariance(vFields)
%
% compute componentwise mean and variance of vectors at each point
% assumes vector components are independent
%

% this took too much memory
v = var(vFields,0,ndims(vFields));
return

s = size(vFields);
dims = ndims(vFields)-2;
fieldSize = s(2:end-1);
v = zeros([3 fieldSize]);

switch dims
 case 3
  for z = 1:fieldSize(3)
    for y = 1:fieldSize(2)
      for x = 1:fieldSize(1)
	v(:,x,y,z) = var(squeeze(vFields(:,x,y,z,:))');
      end
    end
  end
 case 2
  for y = 1:fieldSize(2)
    for x = 1:fieldSize(1)
      v(:,x,y) = var(squeeze(vFields(:,x,y,:))');
    end
  end  
 otherwise
  error('fields must be dim 2 or 3');
end
