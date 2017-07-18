function h = composeAffineH(A,t,g,method)
% returns h(x) = Ag(x) + t
% needs to be vectorized
if nargin < 4
  method = 'linear';
end
h = zeros(size(g),class(g));
s = size(g);
dim = ndims(g)-1;

for i=1:prod(s(2:end))
  idx = dim*i-(dim-1);
  h(idx:idx+(dim-1)) = (A*g(idx:idx+(dim-1))' + t(:))';
end

