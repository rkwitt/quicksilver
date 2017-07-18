function [v,tc] = vFieldPointwiseTraceCovariance(vFields)
%
% compute trace of covariance matrix of vectors at each point
%
s = size(vFields);
dims = ndims(vFields)-2;
fieldSize = s(2:end-1);
tc = zeros(fieldSize);
v  = zeros([3 fieldSize]);

switch dims
 case 3
  for z = 1:fieldSize(3)
    for y = 1:fieldSize(2)
      for x = 1:fieldSize(1)
        d = diag(cov(squeeze(vFields(:,x,y,z,:))'));
	tc(x,y,z) = sum(d);
        v(:,x,y,z) = d;
      end
    end
  end
 case 2
  for y = 1:fieldSize(2)
    for x = 1:fieldSize(1)
      d = cov(squeeze(vFields(:,x,y,:))');
      tc(x,y) = trace(d);
      v(:,x,y) = d;
    end
  end  
 otherwise
  error('fields must be dim 2 or 3');
end
