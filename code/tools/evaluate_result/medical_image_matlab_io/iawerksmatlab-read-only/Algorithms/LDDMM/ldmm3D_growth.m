%
% a template image is given for time 0
% images are given at time 0, 1, ..., T
%

%
% output switches
%
writeDebugImages = false;
writeImages = false;
dispMinMaxL2 = true;
askForInput = false;

% the data type used for storage and computation
% use 'single' for 32-bit floating point
% use 'double' for 64-bit floating point
dtype = 'single';

%
% image names
%
imageFilenames = {}
%ages = 40:2:60;
%sig = 1;
%for i = 1:length(ages)
%  imageFilenames{i} = sprintf('Atlases3D/mean_age%03d_sigma%d_0.mhd',ages(i), sig);
%end

%
% basic parameters
%

% the number of iterations
maxIter = 3;
% the number of timepoints with known images
numTimepoints = 11;
%numTimepoints = length(imageFilenames);
% the number of velocity fields is one less
numVelocityFields = numTimepoints - 1;
% the number of dimensions (of each image)
numDimensions = 3;
% the extent of each image
imageExtent = [64 64 64];
minImageExtent = min(imageExtent);

%
% test image parameters
%

% standard deviation of random noise applied to image
imageNoiseSigma = 0.00; % 0.05
% standard deviation of smoothing kernel applied to image
imageSmoothingSigma = 1.0; % 1

%
% algorithm parameters
%

% increasing sigma decreases the importance of exact image matching
sig = 0.1; % 0.08
% epsilon is the step-size for the optimization
epsilon = 0.001; % 0.001
% how to multiply epsilon if sse goes up
annealingFactor = 0.5;
% alpha, beta, and gamma determine the smoothing of the velocity fields
alpha = 0.5;
beta = 0;
gamma = 1;
% inverse method for inverting velocity fields
inverseMethod = 'zero';

%
% load images 
%
fprintf('Loading images...');
t = cputime;
I = zeros([imageExtent numTimepoints],dtype);
for t=1:numTimepoints
  fprintf('%d',t);
  if length(imageFilenames) > 1
    tmpI = loadMETA(imageFilenames{t});
    size(tmpI)
    I(:,:,:,t) = loadMETA(imageFilenames{t});
  else
    % generated spheres as test images
    I(:,:,:,t) = makeTestImage(imageExtent,'sphere',...
      [minImageExtent*((2*abs(t-numTimepoints/2)+1)/numTimepoints+1)/8,...
      imageNoiseSigma, imageSmoothingSigma],dtype);
  end
end
fprintf(' DONE (%g sec)\n',cputime-t);
if writeDebugImages
  % write original images
  for t = 1:numTimepoints
    writeMETA(squeeze(I(:,:,:,t)),sprintf('debug/input_image_%d.mhd',t));
  end
end

%
% initialize memory for deformation and computation
%
t = cputime;
fprintf('Allocating memory: ');
fprintf('v');
v = zeros([numDimensions imageExtent numTimepoints-1],dtype);
%%%%%%%%%  V,X,Y,Z,t, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% V:  the element of the vector (x, y, or z)
% X, Y, Z: the x, y, and z position of the voxel
% t: the time point (velocity is from t to t+1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% deformed template image from time 0 to each time t
fprintf(',J0t');
J0t    = zeros([imageExtent numTimepoints],dtype);
% (spatial) gradient of the deformed template images
fprintf(',gJ0t');
gJ0t   = zeros([numDimensions imageExtent numTimepoints],dtype); 

% deformed (given) images from t_j to t_i for every i, j s.t. i < j
fprintf(',JTjTi');
for i = 1:numTimepoints-1
  JTjTi{i}    = zeros([imageExtent numTimepoints-i+1]);
end

% determinant of the Jacobian of the mappings \phi_{t_i} \rightarrow
% \phi_{t_j} for every i, j s.t. i < j
fprintf(',dPhiTiTj');
for i = 1:numTimepoints-1
  dPhiTiTj{i} = ones([imageExtent numTimepoints-i+1]);
end

% body force
fprintf(',b');
b      = zeros([numDimensions imageExtent numTimepoints-1],dtype); 
% regularized body force
fprintf(',Linvb');
Linvb  = zeros([numDimensions imageExtent numTimepoints-1],dtype); 

% gradient energy
fprintf(',gE');
gE     = zeros([numDimensions imageExtent numTimepoints-1],dtype); 
fprintf(' DONE (%g sec)\n',cputime-t);

%
% report memory usage
%
memUsageImage = ...
  whos('I','v','J0t','gJ0t','JTjTi','dPhiTiTj','b','Linvb','gE');
memUsageImageBytes = 0;
for i=1:length(memUsageImage)
  memUsageImageBytes = memUsageImageBytes + memUsageImage(i).bytes;
end
fprintf('Memory usage: %d = %0.2f MB\n',...
  memUsageImageBytes, memUsageImageBytes/10^6);

%
% compute initial sum of squared error
%
fprintf('Computing initial error...');
t - cputime;
J0t = computeJ0t(I(:,:,:,1), v, inverseMethod);
diffImage = I - J0t;
sse(1) = sum(diffImage(:).*diffImage(:));
fprintf(' Done (%g sec), sse = %g \n', cputime-t, sse(1));

%
% set up visualization
%
figHandle_Images = figure;
figHandle_Convergence = figure;
figHandle_ErrorHistogram = figure;
figHandles = [figHandle_Images figHandle_Convergence figHandle_ErrorHistogram];
vizSlice = 35;
vizLDMMState_growth(figHandles, v, I, J0t, dPhiTiTj, sse, sig, epsilon, 0, maxIter, vizSlice);

if askForInput 
  x = input('press enter to continue...');
  if x == -1
    askForInput = false;
  end
end

%
% start main loop
%
startTime = cputime;
for iter=1:maxIter
  iterStartTime = cputime;

  fprintf('============== iter: %d, elapsed time: %g ============== \n',...
    iter, cputime-startTime);
  
  if length(sse) > 1 & sse(end) > sse(end-1)
    epsilon(iter+1) = epsilon(iter) * annealingFactor;
  else
    epsilon(iter+1) = epsilon(iter);
  end
  
  %
  % compute template image deformed to each timepoint
  % last dim=numTimePoints
  t = cputime;
  fprintf('Deforming template image to each time point...');
  J0t = computeJ0t(I(:,:,:,1), v, inverseMethod);
  fprintf(' DONE (%g sec)\n', cputime-t);
  
  %
  % compute gradient of these deformed template images
  % last dim=numTimePoints
  t = cputime;
  fprintf('Computing gradients of deformed template images...');
  gJ0t = computeImageGradients(J0t);
  fprintf(' DONE (%g sec)\n',cputime-t);
  
  %
  % for each time point, pull back images from later time points and
  % compute the determinant of the jacobian of this deformation
  % last dim = numTimepoints - n
  t = cputime;  
  fprintf('Deforming images and computing determinant jacobians...');  
  [JTjTi, dPhiTiTj] = computeJTjTi(I, v);
  fprintf('DONE (%g sec)\n', cputime-t);

  %
  % compute "body force" b 
  %
  t = cputime;
  fprintf('Computing body force...');
  for i = 1:numTimepoints-1
    % what about sigma, what about sign
    weights = ...
      - mean((JTjTi{i} - repmat(J0t(:,:,:,i),[1,1,1,numTimepoints-i+1])) .*...
      dPhiTiTj{i}, ...
      numDimensions+1);
    b(1,:,:,:,i) = squeeze(gJ0t(1,:,:,:,i)) .* weights;
    b(2,:,:,:,i) = squeeze(gJ0t(2,:,:,:,i)) .* weights;
    b(3,:,:,:,i) = squeeze(gJ0t(3,:,:,:,i)) .* weights;
  end
  fprintf(' DONE (%g sec)\n',cputime-t);
  if dispMinMaxL2
    fprintf('NUM NAN = %g\n',sum(isnan(b(:))));
    fprintf('MAX = %g\n',max(sqrt(b(:).*b(:))));
    fprintf('MIN = %g\n',min(sqrt(b(:).*b(:))));
  end
    
  %
  % apply Greens function to regularize b
  %
  t = cputime;
  fprintf('Applying Green''s function...');
  for i = 1:numTimepoints-1
    fprintf('%d',i);
    Linvb(:,:,:,:,i) = greensFunction(squeeze(b(:,:,:,:,i)),...
      'laplace',[alpha beta gamma]);
  end
  fprintf(' DONE (%g sec)\n',cputime-t);
  if dispMinMaxL2
    fprintf('NUM NAN = %g\n',sum(isnan(Linvb(:))));
    fprintf('MAX = %g\n',max(sqrt(Linvb(:).*Linvb(:))));
    fprintf('MIN = %g\n',min(sqrt(Linvb(:).*Linvb(:))));
  end

  % zero boundary condition
  Linvb(:,[1 end],:,:,:) = 0;
  Linvb(:,:,[1 end],:,:) = 0;
  Linvb(:,:,:,[1 end],:) = 0;

  %
  % compute gradient energy gE = 2*v - 2/sigma^2 * Linvb
  %
  t = cputime;
  fprintf('Computing gradient energy...');
  gE = 2*v - (2/sig^2)*Linvb;
  fprintf(' DONE (%g sec)\n',cputime-t);
  if dispMinMaxL2
    fprintf('NUM NAN = %g\n',sum(isnan(gE(:))));
    fprintf('MAX = %g\n',max(sqrt(gE(:).*gE(:))));
    fprintf('MIN = %g\n',min(sqrt(gE(:).*gE(:))));
  end

  %
  % update velocity fields v = v - epsilon * gE
  %
  t = cputime;
  fprintf('Updating velocity fields...');
  v = v - epsilon(end)*gE;
  fprintf(' DONE (%g sec)\n',cputime-t);
  if dispMinMaxL2
    fprintf('NUM NAN = %g\n',sum(isnan(v(:))));
    fprintf('MAX = %g\n',max(sqrt(v(:).*v(:))));
    fprintf('MIN = %g\n',min(sqrt(v(:).*v(:))));
  end
    
  %
  % update convergence measure
  %
  t = cputime;
  fprintf('Computing current image match error...');
  J0t = computeJ0t(I(:,:,:,1), v, inverseMethod);
  diffImage = I - J0t;
  sse(iter+1) = sum(diffImage(:).*diffImage(:));
  fprintf(' DONE (%g sec)\n',cputime-t);
  fprintf('################################### sse = %g #### %g%%\n',...
  sse(iter+1), sse(iter+1)/sse(1) * 100);

  fprintf('Iteration Time: %g\n',cputime-iterStartTime);
  vizLDMMState_growth(figHandles, v, I, J0t, dPhiTiTj,...
    sse, sig, epsilon, iter, maxIter, vizSlice);
  
  if askForInput
    x = input('press enter to continue...');
    if x == 0
      askForInput = false;
    end
  end
end

fprintf('Total Time: %g\n',cputime-startTime);

