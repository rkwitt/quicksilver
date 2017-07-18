function fim = filterImage(im,type,params)

switch type
 case 'gaussian'
   % create 3D gaussian kernel
   p = 5;
   [y,x,z] = meshgrid(1:p,1:p,1:p);
   x = x-(p+1)/2;
   y = y-(p+1)/2;
   z = z-(p+1)/2;
   f = zeros(p,p,p);
   S = eye(3,3)*params(1);
   for i=1:p
     for j=1:p
       for k=1:p
         uvw = [x(i,j,k);y(i,j,k);z(i,j,k)];
         f(i,j,k) = 1/((2*pi)^(3/2) * sqrt(det(S))) * exp(- uvw' * (S\uvw) / 2);
       end
     end
   end
   switch ndims(im)
     case 2,
       % take the central slice of the 3D kernel---which is a 2D kernel
       f = squeeze(f(:,:,3));
       % normalize
       f = f * 1/(sum(f(:)));
       fim = conv2(im,f,'same');
       % compat w/ octave fim = convn(im,f,'same');
     case 3,
       % normalize
       f = f * 1/(sum(f(:)));
       fim = convn(im,f,'same');
     otherwise,
       error('im must be a 2 or 3 dimensional array');
   end
end
