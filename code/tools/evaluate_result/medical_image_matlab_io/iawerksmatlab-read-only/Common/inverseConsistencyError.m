function E = inverseConsistencyError(v,vinv)
%  computes inverse consistency error of v \circ vinv
sv = size(v);
hid = eyeHField(sv(2:end));

% first compute error for v \circ vinv
vvinv = composeHFields(v+hid,vinv+hid);
vvinvError = vFieldPointwiseNorm(vvinv-hid);

E = ...
  [min(vvinvError(:)),max(vvinvError(:)),mean(vvinvError(:)),...
  var(vvinvError(:)),median(vvinvError(:))];

