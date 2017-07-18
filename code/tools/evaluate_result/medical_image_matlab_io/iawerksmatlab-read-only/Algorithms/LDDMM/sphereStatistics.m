clear all;
imageDim = 128;
T = 1;
numCircles = 250;
radiusErrorSigma = 0.12;
numCtsSteps = 1000;
numKernelCenters = 20;
kernelCenters = linspace(T/10,9*T/10,numKernelCenters);
kernelSigmaValues = 0.02;
kernelCutOff = 0.02;

% plot ground truth function
clf;
xcts = 0:(T/numCtsSteps):T;
plot(xcts,0.15*imageDim+0.2*imageDim*sin(pi*xcts/(2*1)));

% generate random examples
ti = randperm(length(xcts));
ti = sort(ti(1:numCircles));
t = xcts(ti);
r = 0.15*imageDim + 0.2*imageDim*sin(pi*t/(2*T)) + randn(size(t))*radiusErrorSigma*sqrt(imageDim);

% plot examples
hold on;
scatter(t,r);

%
% compute density at each age given the current kernel
%
kernelDensities = zeros(length(kernelSigmaValues),length(kernelCenters),2);
weights = zeros(length(r),length(kernelCenters));

for s=1:length(kernelSigmaValues)
  currentSigma = kernelSigmaValues(s);

  for j=1:length(kernelCenters)
    kernelCenter = kernelCenters(j);
    currentWeights = zeros(length(r));
    % compute theoretical weights
    for i = 1:length(r)
      currentWeights(i) = regkernel(t(i),kernelCenter,currentSigma);
      %fprintf('t=%g, k=%g, w=%g\n', t(i), kernelCenter, weights(i));
    end

    % normalize weights
    currentWeights = currentWeights / sum(currentWeights);
    % plot original kernel weights
    plot(t, 50*currentWeights,'m');
    origWeights = currentWeights;
    
    % truncate kernel and renormalize
    currentWeights(currentWeights < kernelCutOff) = 0;
    currentWeights = currentWeights / sum(currentWeights);
    % plot truncated kerenl weights
    plot(t(currentWeights>0), 50*currentWeights(currentWeights>0),'c');

    kernelDensities(s,j,1) = sum(origWeights .* r');
    kernelDensities(s,j,2) = sum(currentWeights .* r');

    weights(:,j) = currentWeights;
  end
end

%
% plot the densities
%
plot(kernelCenters, kernelDensities(:,:,1),'r');
plot(kernelCenters, kernelDensities(:,:,2),'g');

%
% generate images
%
I = zeros(imageDim, imageDim, numCircles);
dim = [imageDim imageDim];
dtype = 'single';
for i=1:numCircles
    dvec = repmat((dim/2)',[1 dim]) - eyeHField(dim);
    dist = vFieldPointwiseNorm(dvec,2);
    im = zeros(dim,dtype);
    im(dist<r(i)) = 1;

    % add noise
    noiseSigma = 0.05;
    im = im + noiseSigma * randn(size(im));

    % smooth the image
    smoothingSigma = 1;
    im = filterImage(im,'gaussian',smoothingSigma);
    
    I(:,:,i) = im;
end

save Ir.mat I r t weights;