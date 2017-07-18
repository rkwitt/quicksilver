function I = MakeBullImage(Dim, R, sigmaNoise, sigmaSmoothing, dtype)

if nargin < 5
  dtype = 'double';
end

if nargin < 4
  sigmaSmoothing = 0.0;
end

if nargin < 3
  sigmaNoise = 0.0;
end

if nargin < 2
  minr = 1/7;
  maxr = 1/3;
  R(3) = minr + rand() * (maxr-minr);
  R(2) = 3*R(3)/5 + rand()*R(3)/5;
  R(1) = R(3)/5 + rand()*R(3)/5;
end

if nargin < 1
  Dim = [256 256]'
end

% create image and compute distance from center
I    = zeros(Dim, dtype);
dvec = repmat((Dim/2)',[1 Dim]) - eyeHField(Dim);
d    = vFieldPointwiseNorm(dvec,2);

% fill in image
mdim = min(Dim(:));
I(d<R(3)*mdim) = 0.5;
I(d<R(2)*mdim) = 1.0;
I(d<R(1)*mdim) = 0.5;

% add noise
if sigmaNoise > 0.0
  I = I + sigmaNoise * randn(size(I));
end

% smooth
if sigmaSmoothing > 0.0
  I = filterImage(I, 'gaussian', sigmaSmoothing);
end

I(I<0)=0;
I(I>1)=1;
