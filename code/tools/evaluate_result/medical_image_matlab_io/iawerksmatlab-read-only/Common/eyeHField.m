function eyeh = eyeHField(dim,dtype)
if nargin < 2
dtype = 'double';
end
switch length(dim)
 case 2,
  [ey,ex] = meshgrid(cast(1,dtype):cast(dim(2),dtype),...
		     cast(1,dtype):cast(dim(1),dtype));
  eyeh(1,:,:) = ex;
  eyeh(2,:,:) = ey;
 case 3,
  [ey,ex,ez]=meshgrid(cast(1,dtype):cast(dim(2),dtype),...
		      cast(1,dtype):cast(dim(1),dtype),...
		      cast(1,dtype):cast(dim(3),dtype));
  eyeh(1,:,:,:) = ex;
  eyeh(2,:,:,:) = ey;
  eyeh(3,:,:,:) = ez;
 otherwise,
  error('dim must be be length 2 or 3 vector');
end
