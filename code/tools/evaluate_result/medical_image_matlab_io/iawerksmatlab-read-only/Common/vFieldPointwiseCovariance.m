function c = vFieldPointwiseCovariance(vFields)
%
% compute covariance matrix of vectors at each point
%
s = size(vFields);
dims = ndims(vFields)-2;
fieldSize = s(2:end-1);
c = zeros([3 3 fieldSize]);

switch dims
  case 3
    for z = 1:fieldSize(3)
      for y = 1:fieldSize(2)
        for x = 1:fieldSize(1)
          c(:,:,x,y,z) = cov(squeeze(vFields(:,x,y,z,:))');
        end
      end
    end
  case 2
    for y = 1:fieldSize(2)
      for x = 1:fieldSize(1)
        c(:,:,x,y) = cov(squeeze(vFields(:,x,y,:))');
      end
    end
  otherwise
    error('fields must be dim 2 or 3');
end
