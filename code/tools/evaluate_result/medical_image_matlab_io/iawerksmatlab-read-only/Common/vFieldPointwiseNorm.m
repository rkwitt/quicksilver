function M = vFieldPointwiseNorm(v,p)
if nargin < 2
  p = 2;
end

switch p
 case 1,
  M = squeeze(sum(abs(v)));
 case 2,
  M = squeeze(sqrt(sum(v.*v)));
 case inf,
  M = squeeze(max(abs(v)));
 case -inf,
  M = squeeze(min(abs(v)));
 otherwise,
  error('p must be 1, 2, inf, or -inf');
end  