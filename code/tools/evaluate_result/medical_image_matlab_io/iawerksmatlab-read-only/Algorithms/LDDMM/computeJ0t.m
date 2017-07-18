function J0t = computeJ0t(I, v,inverseMethod)

if isMatlab
  stderr=2;
end

if nargin < 3
  inverseMethod = 'zero';
  % inverseMethod = 'iter';
end

sv = size(v);
N  = sv(end);
J0t = zeros([sv(2:end-1) N+1], class(v));
hID = eyeHField(sv(2:end-1));
h = hID;

if strcmp(inverseMethod,'iter') 
  fprintf(stderr,'\n');
end

switch ndims(I)
 case 3,
  fprintf(stderr,'1');
  J0t(:,:,:,1) = I;
  for i=2:N+1
    fprintf(stderr,'%d',i);
    % invert the velocity field at t
    vinv = invertVelocityField(v(:,:,:,:,i-1),inverseMethod);

    % compute h(x) = (id+v)(h(x))
    h = composeHFields(h,hID+vinv);

    % pull image back
    J0t(:,:,:,i) = deformImage(I,h);
  end
 case 2,
  fprintf(stderr,'1');
  J0t(:,:,1) = I;
  for i=2:N+1
    fprintf(stderr,'%d',i);
    % invert the velocity field at t
    vinv = invertVelocityField(v(:,:,:,i-1),inverseMethod);
    
    % compute h(x) = h(id+v(x))
    h = composeHFields(h,hID+vinv);

    % pull image back
    J0t(:,:,i) = deformImage(I,h);
  end
 otherwise,
  error('image must be 2 or 3 dimensional');
end
