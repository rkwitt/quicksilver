function G = imageGradient(im)

switch ndims(im)
  case 3,
   [GY,GX,GZ] = gradient(im);
   G(1,:,:,:) = GX;
   G(2,:,:,:) = GY;
   G(3,:,:,:) = GZ;
 case 2,
  [GY,GX] = gradient(im);
  G(1,:,:) = GX;
  G(2,:,:) = GY;
 otherwise,
  error('im must be a 2 or 3 dimensional array');
end

