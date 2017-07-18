function h = computeInverseHFieldFromVField(v,inverseMethod)

if nargin < 2
  inverseMethod = 'zero';
  % inverseMethod = 'iter';
end

sv = size(v);
% number of timesteps
N  = sv(end);
hID = eyeHField(sv(2:end-1));
h = hID;

switch ndims(v)-2
 case 3,
  for i=1:N
    fprintf('%d',i);
    % invert the velocity field at t
    vinv = invertVelocityField(v(:,:,:,:,i),inverseMethod);

    % compute h(x) = h(id(x)+vinv(x))
    h = composeHFields(h,hID+vinv);
  end
 case 2,
  for i=1:N
    fprintf('%d',i);
    % invert the velocity field at t
    vinv = invertVelocityField(v(:,:,:,i),inverseMethod);
    
    % compute h(x) = h(id(x)+vinv(x))
    h = composeHFields(h,hID+vinv);
  end
 otherwise,
  error('image must be 2 or 3 dimensional');
end
