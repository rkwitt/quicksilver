function J = deformImage(im,h,extrapVal,method)
% compute J = im \circ h
% extrapVal is the value that will be assigned if h(x) is outside of the
% image
if nargin < 4
  method = 'linear';
end
if nargin < 3
  extrapVal = 0;
end

switch ndims(im)
 case 3,
  J = interp3(im,...
	      squeeze(h(2,:,:,:)),...
	      squeeze(h(1,:,:,:)),...
	      squeeze(h(3,:,:,:)),...
	      method, extrapVal);
 case 2,
  J = interp2(im,...
	      squeeze(h(2,:,:)),...
	      squeeze(h(1,:,:)),...
	      method, extrapVal);
  otherwise,
   error('im must be a 2 or 3 dimensional array');
end  
