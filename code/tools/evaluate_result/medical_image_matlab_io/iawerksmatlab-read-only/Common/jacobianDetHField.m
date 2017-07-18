function dJ = jacobianDetHField(h)

switch ndims(h)-1
 case 3,
  d1 = imageGradient(squeeze(h(1,:,:,:)));
  d2 = imageGradient(squeeze(h(2,:,:,:)));
  d3 = imageGradient(squeeze(h(3,:,:,:)));
  dJ = zeros(size(h,2),size(h,3),size(h,4),class(h));

  % write out the determinant of a 3x3 matrix
  dJ = squeeze(d1(1,:,:,:).*d2(2,:,:,:).*d3(3,:,:,:) + ...
    d1(3,:,:,:).*d2(1,:,:,:).*d3(2,:,:,:) + ...
    d1(2,:,:,:).*d2(3,:,:,:).*d3(1,:,:,:) - ...
    d1(3,:,:,:).*d2(2,:,:,:).*d3(1,:,:,:) - ...
    d1(1,:,:,:).*d2(3,:,:,:).*d3(2,:,:,:) - ...
    d1(2,:,:,:).*d2(1,:,:,:).*d3(3,:,:,:));

 case 2,
  d1 = imageGradient(squeeze(h(1,:,:)));
  d2 = imageGradient(squeeze(h(2,:,:)));
  dJ = squeeze(d1(1,:,:).*d2(2,:,:) - d1(2,:,:).*d2(1,:,:));
 otherwise,
  error('fields must be 2 or 3 dimensional');
end
