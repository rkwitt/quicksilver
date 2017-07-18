function J = jacobianHField(h)

J = zeros([size(h,1) size(h)],class(h));
switch ndims(h)-1
 case 3,
  J(1,:,:,:,:) = imageGradient(squeeze(h(1,:,:,:)));
  J(2,:,:,:,:) = imageGradient(squeeze(h(2,:,:,:)));
  J(3,:,:,:,:) = imageGradient(squeeze(h(3,:,:,:)));
 case 2,
  J(1,:,:,:)   = imageGradient(squeeze(h(1,:,:)));
  J(2,:,:,:)   = imageGradient(squeeze(h(2,:,:)));
 otherwise,
  error('fields must be 2 or 3 dimensional');
end
