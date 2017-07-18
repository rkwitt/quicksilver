function im = makeTestImage(dim,type,params,dtype)
if nargin < 3
  params = [];
end
if nargin < 4
  dtype = 'double';
end

switch type
  case 'sphere',
   % params = {radius, noise sigma, smoothing sigma}
   dvec = repmat((dim/2)',[1 dim]) - eyeHField(dim);
   d = vFieldPointwiseNorm(dvec,2);
   im = zeros(dim,dtype);
   r = min(dim)/3;
   if length(params) > 0
     r = params(1);
   end
   im(d<r) = 1;
   im(d>r) = 0;
   
   % add noise
   noiseSigma = 0;
   if length(params > 1)
     noiseSigma = params(2);
   end
   im = im + noiseSigma * randn(size(im));
   
   smoothingSigma = 1;
   if length(params) > 2
     smoothingSigma = params(3);
   end
   im = filterImage(im,'gaussian',smoothingSigma);
 case 'bull',
  im = zeros(dim,dtype);
  
  % compute distance from center
  dvec = repmat((dim/2)',[1 dim]) - eyeHField(dim);
  d = vFieldPointwiseNorm(dvec,2);
  
  % minimum and maximum radius
  minr = 1/7;
  maxr = 1/3;
  
  % compute a random radius from minr to maxr * the image size
  r = (minr + rand()*(maxr-minr))*min(dim(:));
  
  % outer ring
  im(d<r)  = 0.5;
  
  % middle ring
  rm = 3*r/5 + rand()*r/5;
  im(d<rm) = 1.0;
  
  % inner ring
  ri = r/5 + rand*r/5;
  im(d<ri) = 0.5;
  
  % add noise
  noiseSigma = 0;
  if length(params > 0)
    noiseSigma = params(1);
  end
  im = im + noiseSigma * randn(size(im));
  
  % smooth image
  smoothingSigma = 1;
  if length(params) > 1
    smoothingSigma = params(2);
  end
  im = filterImage(im,'gaussian',smoothingSigma);
 case 'ellipse',
end
