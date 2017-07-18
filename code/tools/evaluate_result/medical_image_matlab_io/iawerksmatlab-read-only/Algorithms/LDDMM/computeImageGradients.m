function G = computeImageGradients(I)
% computes the gradient of a collection of images

if isMatlab
  stderr=2;
end

isize = size(I);

switch ndims(I)-1
 case 3,
  G = zeros([3 isize], class(I));
  for i=1:isize(end)
    fprintf(stderr,'%d',i);
    G(:,:,:,:,i) = imageGradient(I(:,:,:,i));
  end
 case 2,
  G = zeros([2 isize], class(I));
  for i=1:isize(end)
    fprintf(stderr,'%d',i);
    G(:,:,:,i) = imageGradient(I(:,:,i));
  end
 otherwise,
  error('Must be 2 or 3 dimensional image collection');
end
