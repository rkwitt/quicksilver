
d = 128;
T = 1000;
numCircles = 200;
errorSigma = 0.12;

clf;

% plot ground truth function
xcts = 1:T;
plot(xcts,0.15*d+0.2*d*sin(pi*xcts/(2*T)));

% choose random time points
t = randperm(length(xcts));
t = sort(t(1:numCircles));
r = 0.15*d + 0.2*d*sin(pi*t/(2*T)) + randn(size(t))*errorSigma*sqrt(d);

% plot examples
hold on;
scatter(t,r);

% generate images with these radii
I = zeros(d,d,numCircles);
dim = [d d];
dtype = 'single';
for i=1:numCircles
    dvec = repmat((dim/2)',[1 dim]) - eyeHField(dim);
    dist = vFieldPointwiseNorm(dvec,2);
    im = zeros(dim,dtype);
    im(dist<r(i)) = 1;

    % add noise
    noiseSigma = 0.05;
    im = im + noiseSigma * randn(size(im));

    smoothingSigma = 1;
    im = filterImage(im,'gaussian',smoothingSigma);
    
    I(:,:,i) = im;
end
