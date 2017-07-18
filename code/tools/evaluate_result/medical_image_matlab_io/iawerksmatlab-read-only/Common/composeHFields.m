function h = composeHFields(f,g,method)
% returns h = f \circ g
if nargin < 3
  method = 'linear';
end
h = zeros(size(f),class(f));

switch ndims(f)-1
 case 3,
  h(1,:,:,:) = interp3(squeeze(f(1,:,:,:)),...
		       g(2,:,:,:), g(1,:,:,:), g(3,:,:,:),...
		       method);
  h(2,:,:,:) = interp3(squeeze(f(2,:,:,:)),...
		       g(2,:,:,:), g(1,:,:,:), g(3,:,:,:),...
		       method);
  h(3,:,:,:) = interp3(squeeze(f(3,:,:,:)),...
		       g(2,:,:,:), g(1,:,:,:), g(3,:,:,:),...
		       method);
 case 2,
  h(1,:,:) = interp2(squeeze(f(1,:,:)),g(2,:,:),g(1,:,:),method);
  h(2,:,:) = interp2(squeeze(f(2,:,:)),g(2,:,:),g(1,:,:),method);
 otherwise,
  error('fields must be 2 or 3 dimensional');
end
