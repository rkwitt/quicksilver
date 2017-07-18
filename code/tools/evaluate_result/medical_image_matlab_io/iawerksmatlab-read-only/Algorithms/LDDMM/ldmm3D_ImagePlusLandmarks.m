function ldmm3D_ImagesPlusLandmarks(imageFilenames, imageWeights, maxIter)
%
% Ihat is at time T
% input images are at time 0
%

%
% output switches
%
writeDebugImages = false;
writeDebugPoints = false;
writeImages = false;
writePoints = false;
writeFinalAtlasImage = true;
writeDeformationFields = true;
dispMinMaxL2 = true;
extraOutput = false;
displayPlotOverview = true;

writeDebugEverything = false;
if writeDebugEverything
  writeDebugImages = true;
  writeDebugPoints = true;
  writeImages = true;
  writePoints = true;
  dispMinMaxL2 = true;
  extraOutput = true;
  displayPlotOverview = true;
end

skipPointSetMatching  = true;
skipImageMatching     = false;
skipAnnealingSchedule = false;
skipPointMatchSigmaReduction = true;

% the data type used for storage and computation
% use 'single' for 32-bit floating point
% use 'double' for 64-bit floating point
dtype = 'single';

% the number of LDMM (outer) iterations
%maxIter = 50; % 50 % in function arg
% the number of input images
M = length(imageFilenames); % 4 % bull images 9
%M = 4;
% image weights
%imageWeights = [0.25 0.25 0.25 0.25]; % in function arg
% the number of timesteps (i.e. velocity fields)
N = 15; % 10 % bull images 15
% the number of dimensions for the image (2 or 3)
d = 3;
% the extent of the image (number of pixels on a side)
% NB: the image should be a cube, and it's best if s is a power of 2
switch d
  case 2,
    s = 256; % 128
  case 3,
    s = 64;  % 80
end
% this will be [s s] for 2d and [s s s] for 3d
svec = repmat(s,1,d);
% the minimum and maximum number of feature points per image
% the final number of points per image will be in [min, max]
minmaxpts = d*[16 16]; % [16 16]
% standard deviation of random noise applied to the point positions
ptNoiseSigma = 0.0; % 0.1
% standard deviation of random noise applied to image
imageNoiseSigma = 0.025; % 0.05
% standard deviation of smoothing kernel applied to image
imageSmoothingSigma = 1.0; % 1

% sigma weights the smoothed velocity field in the gradient energy computation  
% increasing sigma will decrease the relative weight of the new velocity
% field
% NB: don't use the actual word 'sigma' because it is a built in matlab
% function 
%sigma_image  = 0.06; % 0.06 % for spheres
sigma_image = 1000;
sigma_points = 0.07; % 0.04

% weight the distance measure between points for point matching (this is
% the sigma used in the point set metric)
sigma_ptmatch = 20;

% epsilon weights the gradient energy in the velocity update equation 
% eps=0.5 means replace current velocity with Linvb
% eps=0 means v stays the same
% initial development with 0.001
epsilon(1)  = 0.0005; % for spheres
%epsilon(1) = 0.000001;

% alpha, beta, and gamma determine the smoothing of the velocity fields
operatorType='laplace';
alpha = 0.5;
beta = 0;
gamma = 1;

alpha_pts = 10;
beta_pts = 0;
gamma_pts = 1;

inverseMethod = 'zero'; % 'zero', 'iter'

%
% load images and point sets
%
fprintf('Loading images...');
t = cputime;
I = zeros([svec M],dtype);
for m=1:M
  fprintf('%d',m);
  radius = s*(m/M+1)/8;
  switch d
    case 2,
      I(:,:,m)            = makeTestImage([s s],...
        'sphere',[radius imageNoiseSigma imageSmoothingSigma],dtype);
      %I(:,:,m)            = makeTestImage([s s],...
      %  'bull',[imageNoiseSigma imageSmoothingSigma],dtype);      

      %[P{m}, PWeights{m}] = makeTestPoints([s s],'sphere',...
      %  [radius ptNoiseSigma],minmaxpts,dtype);
      [P{m}, PWeights{m}] = makeTestPoints([s s],'bull',...
        [radius ptNoiseSigma],minmaxpts,dtype);
      %I(:,:,m) = makeTestImage([s s],'sphere',s*(rand+1)/8,dtype);
    case 3,
      currentImage        = loadMETA(imageFilenames{m});
      I(:,:,:,m)          = currentImage;
      mean(currentImage(:))
      std(currentImage(:))
      %I(:,:,:,m)          = makeTestImage([s s s],...
      %  'sphere',[radius imageNoiseSigma imageSmoothingSigma],dtype);
      [P{m}, PWeights{m}] = makeTestPoints([s s s],'sphere',...
        [radius ptNoiseSigma],minmaxpts,dtype);
      %I(:,:,:,m) = makeTestImage([s s s],'sphere',s*(rand+1)/8,dtype);
  end
end
fprintf(' DONE (%g sec)\n',cputime-t);
if writeImages
  % write original images
  for q = 1:M
    switch d
      case 2,
        writeMETA(squeeze(I(:,:,q)),sprintf('debug/input_image_%d.mhd',q));
      case 3,
        writeMETA(squeeze(I(:,:,:,q)),sprintf('debug/input_image_%d.mhd',q));
    end
  end
end
if writePoints
  % write original points
  for q = 1:M
    writeMETALandmark([P{q}; PWeights{q}],...
      sprintf('debug/input_pts_%d.mhd',q));
  end
end

%
% compute initial average image
%
J0T = I;
t = cputime;
fprintf('Computing initial average image...');
% Ihat = mean(J0T,ndims(I));
% weighted sum
Ihat = sum(repmat(reshape(imageWeights,[1,1,1,M]),[s,s,s,1]) .* J0T,ndims(I));

Ihat0 = Ihat;
fprintf(' DONE (%g sec)\n',cputime-t);
% write initial average image
if writeImages
  writeMETA(Ihat,'debug/Ihat_k0.mhd');
end
fprintf('Computing initial variance...');
t = cputime;
Ivar = var(J0T,0,ndims(I));
Ivar0 = Ivar;
avar(1) = sum(Ivar(:));
fprintf(' DONE (%g sec)\n',cputime-t);
fprintf(' === Sum of voxelwise variance %g ::: %g%% === \n',avar(1),100*avar(1)/avar(1));
% write initial variance image
if writeImages
  writeMETA(Ivar,'debug/Ihat_ptwise_variance_k0.mhd');
end

% initial points average (in this case the union of the points). normalize
% weights.
JP0T = P;
t = cputime;
fprintf('Computing initial average points...');
IhatP = cell2mat(JP0T);
IhatP0 = IhatP;
IhatPWeights = cell2mat(PWeights);
IhatPWeights = IhatPWeights / sum(IhatPWeights(:));
fprintf(' DONE (%g sec)\n',cputime-t);
if writePoints
  % write initial points average
  writeMETALandmark([IhatP; IhatPWeights],sprintf('debug/IhatP_k0.mhd'));
end

%
% display initial image set
%
if displayPlotOverview
  figHandle_Overview = figure;
  figHandle_Mean = figure;
  figHandle_DeformingSample = figure;
  figHandles = [figHandle_Overview figHandle_Mean figHandle_DeformingSample];
  switch d
    case 2,
      vizLDMMState(figHandles,...
        J0T, Ihat0, Ihat, Ivar0, Ivar, JP0T, IhatP0, IhatP, avar, epsilon, maxIter);
    case 3,
      vizLDMMState(figHandles,...
        squeeze(J0T(:,:,round(end/2),:)),...
        squeeze(Ihat0(:,:,round(end/2))),...
        squeeze(Ihat(:,:,round(end/2))),...
        squeeze(Ivar0(:,:,round(end/2))),...
        squeeze(Ivar(:,:,round(end/2))),...
        JP0T, IhatP0, IhatP, avar, epsilon, maxIter);
  end
end

%
% initialize memory 
%
t = cputime;
fprintf('Allocating memory for image matching: ');
fprintf('v');
v = zeros([d svec N M],dtype);
%%%%%%%%%  V,X,Y(,Z),t,m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% V:  the element of the vector (x or y)
% X, Y (,Z): the x y z position of the voxel
% t: the time point
% m: the input image index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% image data
% deformed image from time T to each time t
fprintf(',JTt');
JTt    = zeros([svec N],dtype);
% det of Jacobian of h-field mapping \Omega at t to T
fprintf(',dPhitT');
dPhitT = zeros([svec N],dtype);
% deformed image from time 0 to each time t
fprintf(',J0t');
J0t    = zeros([svec N+1],dtype);
% (spatial) gradient of J0t images
fprintf(',gJ0t');
gJ0t   = zeros([d svec N+1],dtype); 
% body force
fprintf(',b');
b      = zeros([d svec N],dtype); 
% regularized body force
fprintf(',Linvb');
Linvb  = zeros([d svec N],dtype); 
% gradient energy term from images
fprintf(',gE');
gE     = zeros([d svec N],dtype); 
fprintf(' DONE (%g sec)\n',cputime-t);
memUsageImage = ...
  whos('I','J0T','Ihat','v','JTt','dPhitT','J0t','gJ0t','b','Linvb','gE');
memUsageImageBytes = 0;
for i=1:length(memUsageImage)
  memUsageImageBytes = memUsageImageBytes + memUsageImage(i).bytes;
end

% landmark data
t = cputime;
fprintf('Allocating memory for point matching: ');
% points deformed from time 0 to each time t
fprintf('JP0t');
JP0t          = zeros([d minmaxpts(2) N+1],dtype);
% velocity field implied by pairwise differences between this pt set and
% the average point set
fprintf(',bIhatPJPDiffs');
bIhatPJPDiffs = zeros([d svec N],dtype);
fprintf(',vIhatPJPDiffs');
vIhatPJPDiffs = zeros([d svec N],dtype);
% velocity field implied by pairwise differences between this pt set and
% itself
fprintf(',bJPJPDiffs');
bJPJPDiffs    = zeros([d svec N],dtype);
fprintf(',vJPJPDiffs');
vJPJPDiffs    = zeros([d svec N],dtype);
% gradient energy term from point sets
fprintf(',gEP');
gEP           = zeros([d svec N],dtype);
fprintf(' DONE (%g sec)\n',cputime-t);
memUsagePts = ...
  whos('P','JP0T','IhatP','JP0t','bIhatPJPDiffs','vIhatPJPDiffs',...
  'bJPJPDiffs','vJPJPDiffs','gEP');
memUsagePtsBytes = 0;
for i=1:length(memUsagePts)
  memUsagePtsBytes = memUsagePtsBytes + memUsagePts(i).bytes;
end
fprintf('Memory usage: (image %d), (points %d), (total %d = %0.2f MB)\n',...
  memUsageImageBytes, memUsagePtsBytes,...
  memUsageImageBytes+memUsagePtsBytes,...
  (memUsageImageBytes+memUsagePtsBytes)/10^6); 

input('Press enter to begin matching...');

%
% start main interation loop
%
startTime = cputime;
for k=1:maxIter
  % decrease step size if energy increased during any of the last set of
  % sub-iterations
  if ~skipAnnealingSchedule && k > 1 && any(diff(avar(end-M:end)) > 0)
    diff(avar(end-M:end))
    epsilon = [epsilon epsilon(end)/2];
  else
    epsilon = [epsilon epsilon(end)];
  end
  
  % reduce point matching radius every iteration
  if ~skipPointMatchSigmaReduction && sigma_ptmatch > 2
    sigma_ptmatch = sigma_ptmatch * 2 / 3;
  end
  
  %
  % start sub-iteration loop
  %
  for m=1:M
    iterStartTime = cputime;

    fprintf('================== iter: %d, image: %d, elapsed time: %g ================== \n',...
	    k,m,cputime-startTime);
    
    %
    % First address image match term
    %
    if ~skipImageMatching
      fprintf('Image matching...\n');
    
      %
      % compute JTt (image (Ihat) at T deformed to each timepoint t)
      % and dPhitT (det jacobian of phi_t, the hfield from t to T)
      %
      t = cputime;
      fprintf('Computing (backward, T-->t) def. images & det. Jacobians...');
      switch d
        case 2,
          [JTt, dPhitT] = computeJTt(Ihat, v(:,:,:,:,m));
        case 3,
          [JTt, dPhitT] = computeJTt(Ihat, v(:,:,:,:,:,m));
      end
      fprintf(' DONE (%g sec)\n',cputime-t);
      if writeDebugImages
        % write debug images
        for q = 1:N
          switch d
            case 2,
              writeMETA(squeeze(JTt(:,:,q)),sprintf('debug/JTt_k%d_m%d_t%d.mhd',k,m,q));
              writeMETA(squeeze(dPhitT(:,:,q)),sprintf('debug/dPhitT_k%d_m%d_t%d.mhd',k,m,q));
            case 3,
              writeMETA(squeeze(JTt(:,:,:,q)),sprintf('debug/JTt_k%d_m%d_t%d.mhd',k,m,q));
              writeMETA(squeeze(dPhitT(:,:,:,q)),sprintf('debug/dPhitT_k%d_m%d_t%d.mhd',k,m,q));
          end
        end
      end
      
      %
      % compute J0t (image at 0 deformed to each timepoint t)
      %
      t = cputime;
      fprintf('Computing (forward, 0-->t) deformed images...');
      switch d
        case 2,
          J0t = computeJ0t(I(:,:,m), v(:,:,:,:,m), inverseMethod);
        case 3,
          J0t = computeJ0t(I(:,:,:,m), v(:,:,:,:,:,m), inverseMethod);
      end
      fprintf(' DONE (%g sec)\n',cputime-t);
      if writeDebugImages
        % write debug images
        for q = 1:N+1
          switch d
            case 2,
              writeMETA(squeeze(J0t(:,:,q)),sprintf('debug/J0t_k%d_m%d_t%d.mhd',k,m,q));
            case 3,
              writeMETA(squeeze(J0t(:,:,:,q)),sprintf('debug/J0t_k%d_m%d_t%d.mhd',k,m,q));
          end
        end
      end
      
      %
      % compute gJ0t (gradient of each deformed (0->t) image)
      %
      t = cputime;
      fprintf('Computing gradients of (forward) deformed images...');
      gJ0t = computeImageGradients(J0t);
      fprintf(' DONE (%g sec)\n',cputime-t);
      if writeDebugImages
        % write debug images
        for q = 1:N+1
          switch d
            case 2,
              writeMETA(squeeze(gJ0t(:,:,:,q)),sprintf('debug/gJ0t_k%d_m%d_t%d.mhd',k,m,q));
            case 3,
              writeMETA(squeeze(gJ0t(:,:,:,:,q)),sprintf('debug/gJ0t_k%d_m%d_t%d.mhd',k,m,q));
          end
        end
      end
      
      %
      % compute "body force" b = dPhitT*(J0t-JTt) * gJ0t at each time step
      %
      t = cputime;
      fprintf('Computing body force...');
      switch d
        case 2,
          b(1,:,:,:) = dPhitT.*(J0t(:,:,1:N)-JTt) .* squeeze(gJ0t(1,:,:,1:N));
          b(2,:,:,:) = dPhitT.*(J0t(:,:,1:N)-JTt) .* squeeze(gJ0t(2,:,:,1:N));
        case 3,
          b(1,:,:,:,:) = dPhitT.*(J0t(:,:,:,1:N)-JTt) .* squeeze(gJ0t(1,:,:,:,1:N));
          b(2,:,:,:,:) = dPhitT.*(J0t(:,:,:,1:N)-JTt) .* squeeze(gJ0t(2,:,:,:,1:N));
          b(3,:,:,:,:) = dPhitT.*(J0t(:,:,:,1:N)-JTt) .* squeeze(gJ0t(3,:,:,:,1:N));
      end
      fprintf(' DONE (%g sec)\n',cputime-t);
      if dispMinMaxL2
        fprintf('NUM NAN = %g\n',sum(isnan(b(:))));
        fprintf('MAX = %g\n',max(b(:)));
        fprintf('MIN = %g\n',min(b(:)));
      end
      if writeDebugImages
        % write debug images
        for q = 1:N
          switch d
            case 2,
              writeMETA(squeeze(b(:,:,:,q)),sprintf('debug/b_k%d_m%d_t%d.mhd',k,m,q));
            case 3,
              writeMETA(squeeze(b(:,:,:,:,q)),sprintf('debug/b_k%d_m%d_t%d.mhd',k,m,q));
          end
        end
      end

      %
      % apply Greens function to regularize b
      %
      t = cputime;
      fprintf('Applying Green''s function...');
      for tm = 1:size(b,ndims(b))
        fprintf('%d',tm);
        switch d
          case 2,
            Linvb(:,:,:,tm) = greensFunction(squeeze(b(:,:,:,tm)),...
              operatorType,[alpha beta gamma]);
          case 3,
            Linvb(:,:,:,:,tm) = greensFunction(squeeze(b(:,:,:,:,tm)),...
              operatorType,[alpha beta gamma]);
        end
      end
      fprintf(' DONE (%g sec)\n',cputime-t);
      if dispMinMaxL2
        fprintf('NUM NAN = %g\n',sum(isnan(Linvb(:))));
        fprintf('MAX = %g\n',max(Linvb(:)));
        fprintf('MIN = %g\n',min(Linvb(:)));
      end
      if writeDebugImages
        % write debug images
        for q = 1:N
          switch d
            case 2,
              writeMETA(squeeze(Linvb(:,:,:,q)),sprintf('debug/Linvb_k%d_m%d_t%d.mhd',k,m,q));
            case 3,
              writeMETA(squeeze(Linvb(:,:,:,:,q)),sprintf('debug/Linvb_k%d_m%d_t%d.mhd',k,m,q));
          end
        end
      end
    
      % zero boundary condition
      switch d
        case 2,
          Linvb(:,[1 end],:,:) = 0;
          Linvb(:,:,[1 end],:) = 0;
        case 3,
          Linvb(:,[1 end],:,:,:) = 0;
          Linvb(:,:,[1 end],:,:) = 0;
          Linvb(:,:,:,[1 end],:) = 0;
      end
      
      %
      % compute gradient energy term gE = - 2/sigma^2 * Linvb
      %
      t = cputime;
      fprintf('Computing gradient energy...');
      gE = -(2/sigma_image^2)*Linvb;
      fprintf(' DONE (%g sec)\n',cputime-t);
      if dispMinMaxL2
        fprintf('NUM NAN = %g\n',sum(isnan(gE(:))));
        fprintf('MAX = %g\n',max(gE(:)));
        fprintf('MIN = %g\n',min(gE(:)));
      end
      if writeDebugImages
        for q = 1:N
          switch d
            case 2,
              writeMETA(squeeze(gE(:,:,:,q)),sprintf('debug/gE_k%d_m%d_t%d.mhd',k,m,q));
            case 3,
              writeMETA(squeeze(gE(:,:,:,:,q)),sprintf('debug/gE_k%d_m%d_t%d.mhd',k,m,q));
          end
        end
      end
    else
      gE(:) = 0;
    end
    
    %
    %
    % now compute velocity field that will minimize energy based on point
    % sets
    %
    %
    if ~skipPointSetMatching
      fprintf('Pointset matching...\n');
      
      %
      % deform image points from time 0 to each time point t according to the
      % current deformation
      %
      t = cputime;
      fprintf('Computing (forward, 0-->t) deformed points...');
      switch d
        case 2,
          JP0t = computeJP0t(P{m}, v(:,:,:,:,m));
        case 3,
          JP0t = computeJP0t(P{m}, v(:,:,:,:,:,m));
      end
      fprintf(' DONE (%g sec)\n',cputime-t);
      if writeDebugPoints
        for q = 1:N+1
          writeMETALandmark(squeeze(J0t(:,:,q)),sprintf('debug/J0t_k%d_m%d_t%d.mhd',k,m,q));
        end
      end
      
      %
      % compute the gradient energy term related to weighted, pairwise
      % differences between Ihat landmarks and these landmarks at time T.
      % The differences are pulled back to each time point t and a velocity
      % field is generated via the green's function.
      %
      t = cputime;
      fprintf('Computing pairwise difference field IhatP-JP...');
      switch d
        case 2,
          bIhatPJPDiffs = computePairwiseWeightedDifferencesField(IhatP,...
            JP0t, v(:,:,:,:,m), sigma_ptmatch);
        case 3,
          bIhatPJPDiffs = computePairwiseWeightedDifferencesField(IhatP,...
            JP0t, v(:,:,:,:,:,m), sigma_ptmatch);
      end
      fprintf(' DONE (%g sec)\n',cputime-t);
      if writeDebugPoints
        for q = 1:N
          switch d
            case 2,
              writeMETA(squeeze(bIhatPJPDiffs(:,:,:,q)),...
                sprintf('debug/bIhatJP_k%d_m%d_t%d.mhd',k,m,q));
            case 3,
              writeMETA(squeeze(bIhatPJPDiffs(:,:,:,:,q)),...
                sprintf('debug/bIhatJP_k%d_m%d_t%d.mhd',k,m,q));
          end
        end
      end
      
      %
      % only need to do this once since everything is linear
      %
      
      t = cputime;
      fprintf('Applying Green''s function...');
      for tm = 1:size(bIhatPJPDiffs,ndims(bIhatPJPDiffs))
        fprintf('%d',tm);
        switch d
          case 2,
            vIhatPJPDiffs(:,:,:,tm) = greensFunction(squeeze(bIhatPJPDiffs(:,:,:,tm)),...
              operatorType,[alpha_pts beta_pts gamma_pts]);
          case 3,
            vIhatPJPDiffs(:,:,:,:,tm) = greensFunction(squeeze(bIhatPJPDiffs(:,:,:,:,tm)),...
              operatorType,[alpha_pts beta_pts gamma_pts]);
        end
      end
      fprintf(' DONE (%g sec)\n',cputime-t);
      if dispMinMaxL2
        fprintf('NUM NAN = %g\n',sum(isnan(vIhatPJPDiffs(:))));
        fprintf('MAX = %g\n',max(vIhatPJPDiffs(:)));
        fprintf('MIN = %g\n',min(vIhatPJPDiffs(:)));
      end
      if writeDebugPoints
        % write debug images
        for q = 1:N
          switch d
            case 2,
              writeMETA(squeeze(vIhatPJPDiffs(:,:,:,q)),sprintf('debug/vIhatJP_k%d_m%d_t%d.mhd',k,m,q));
            case 3,
              writeMETA(squeeze(vIhatPJPDiffs(:,:,:,:,q)),sprintf('debug/vIhatJP_k%d_m%d_t%d.mhd',k,m,q));
          end
        end
      end
      
      % zero boundary condition
      switch d
        case 2,
          vIhatPJPDiffs(:,[1 end],:,:) = 0;
          vIhatPJPDiffs(:,:,[1 end],:) = 0;
        case 3,
          vIhatPJPDiffs(:,[1 end],:,:,:) = 0;
          vIhatPJPDiffs(:,:,[1 end],:,:) = 0;
          vIhatPJPDiffs(:,:,:,[1 end],:) = 0;
      end
      
      %
      % compute the gradient energy term related to weighted, pairwise
      % differences between these landmarks at time T. The differences are
      % pulled back to each time point t and a velocity field is generated
      % via the green's function.
      %
      t = cputime;
      fprintf('Computing pairwise difference field JP-JP...');
      switch d
        case 2,
          bJPJPDiffs = computePairwiseWeightedDifferencesField(JP0t(:,:,end),...
            JP0t, v(:,:,:,:,m), sigma_ptmatch);
        case 3,
          bJPJPDiffs = computePairwiseWeightedDifferencesField(JP0t(:,:,end),...
            JP0t, v(:,:,:,:,:,m), sigma_ptmatch);
      end
      fprintf(' DONE (%g sec)\n',cputime-t);
      if writeDebugPoints
        for q = 1:N
          switch d
            case 2,
              writeMETA(squeeze(bIhatPJPDiffs(:,:,:,q)),...
                sprintf('debug/bJPJP_k%d_m%d_t%d.mhd',k,m,q));
            case 3,
              writeMETA(squeeze(bIhatPJPDiffs(:,:,:,:,q)),...
                sprintf('debug/bJPJP_k%d_m%d_t%d.mhd',k,m,q));
          end
        end
      end
      
      t = cputime;
      fprintf('Applying Green''s function...');
      for tm = 1:size(bJPJPDiffs,ndims(bJPJPDiffs))
        fprintf('%d',tm);
        switch d
          case 2,
            vJPJPDiffs(:,:,:,tm) = ...
              greensFunction(squeeze(bJPJPDiffs(:,:,:,tm)),...
              operatorType,[alpha_pts beta_pts gamma_pts]);
          case 3,
            vJPJPDiffs(:,:,:,:,tm) = ...
              greensFunction(squeeze(bJPJPDiffs(:,:,:,:,tm)),...
              operatorType,[alpha_pts beta_pts gamma_pts]);
        end
      end
      fprintf(' DONE (%g sec)\n',cputime-t);
      if dispMinMaxL2
        fprintf('NUM NAN = %g\n',sum(isnan(vJPJPDiffs(:))));
        fprintf('MAX = %g\n',max(vJPJPDiffs(:)));
        fprintf('MIN = %g\n',min(vJPJPDiffs(:)));
      end
      if writeDebugPoints
        % write debug images
        for q = 1:N
          switch d
            case 2,
              writeMETA(squeeze(vJPJPDiffs(:,:,:,q)),sprintf('debug/vJPJP_k%d_m%d_t%d.mhd',k,m,q));
            case 3,
              writeMETA(squeeze(vJPJPDiffs(:,:,:,:,q)),sprintf('debug/vJPJP_k%d_m%d_t%d.mhd',k,m,q));
          end
        end
      end
      
      % zero boundary condition
      switch d
        case 2,
          vJPJPDiffs(:,[1 end],:,:) = 0;
          vJPJPDiffs(:,:,[1 end],:) = 0;
        case 3,
          vJPJPDiffs(:,[1 end],:,:,:) = 0;
          vJPJPDiffs(:,:,[1 end],:,:) = 0;
          vJPJPDiffs(:,:,:,[1 end],:) = 0;
      end
    
      %
      % compute the gradient energy
      %
      t = cputime;
      fprintf('Computing gradient energy...');
      gEP = (2/sigma_points^2)*vJPJPDiffs-(2/sigma_points^2)*vIhatPJPDiffs;
      fprintf(' DONE (%g sec)\n',cputime-t);
      if dispMinMaxL2
        fprintf('NUM NAN = %g\n',sum(isnan(gEP(:))));
        fprintf('MAX = %g\n',max(gEP(:)));
        fprintf('MIN = %g\n',min(gEP(:)));
      end
      if writeDebugPoints
        for q = 1:N
          switch d
            case 2,
              writeMETA(squeeze(gEP(:,:,:,q)),sprintf('debug/gEP_k%d_m%d_t%d.mhd',k,m,q));
            case 3,
              writeMETA(squeeze(gEP(:,:,:,:,q)),sprintf('debug/gEP_k%d_m%d_t%d.mhd',k,m,q));
          end
        end
      end
    else 
      gEP(:) = 0;
    end
    
    %
    % update velocity fields v = v - epsilon_i * gE_i - epsilon_l * gE_l
    %
    t = cputime;
    fprintf('Updating velocity fields...');
    switch d
      case 2,
        v(:,:,:,:,m) = v(:,:,:,:,m) - epsilon(end)*(2*v(:,:,:,:,m)+gE+gEP);
      case 3,
        v(:,:,:,:,:,m) = v(:,:,:,:,:,m) - epsilon(end)*(2*v(:,:,:,:,:,m)+gE+gEP);        
    end
    fprintf(' DONE (%g sec)\n',cputime-t);
    if dispMinMaxL2
      fprintf('NUM NAN = %g\n',sum(isnan(v(:))));
      fprintf('MAX = %g\n',max(v(:)));
      fprintf('MIN = %g\n',min(v(:)));
    end
    if writeDebugImages    
      for q = 1:N
        switch d
          case 2,
            writeMETA(squeeze(v(:,:,:,q)),sprintf('debug/v_k%d_m%d_t%d.mhd',k,m,q));
          case 3,
            writeMETA(squeeze(v(:,:,:,:,q)),sprintf('debug/v_k%d_m%d_t%d.mhd',k,m,q));
        end
      end
    end

    %
    % update Ihat: need to compute J0T for this image
    %
    t = cputime;
    fprintf('Updating average image...');
    switch d
      case 2,
        J0t = computeJ0t(I(:,:,m), v(:,:,:,:,m), inverseMethod);
        J0T(:,:,m) = J0t(:,:,N+1);
      case 3,
        J0t = computeJ0t(I(:,:,:,m), v(:,:,:,:,:,m), inverseMethod);
        J0T(:,:,:,m) = J0t(:,:,:,N+1);        
    end
    % Ihat = mean(J0T,ndims(I));
    % weighted mean
    Ihat = sum(repmat(reshape(imageWeights,[1,1,1,M]),[s,s,s,1]) .* J0T,ndims(I));
    fprintf(' DONE (%g sec)\n',cputime-t);
    if writeDebugImages
      writeMETA(Ihat,sprintf('debug/Ihat_k%d_m%d.mhd',k,m));
    end
    
    %
    % update Ihatp: need to redeform points for this image
    %
    t = cputime;
    fprintf('Updating average points...');
    switch d
      case 2,
        JP0t = computeJP0t(P{m}, v(:,:,:,:,m));
      case 3,
        JP0t = computeJP0t(P{m}, v(:,:,:,:,:,m));        
    end
    JP0T{m} = JP0t(:,:,N+1);
    IhatP = cell2mat(JP0T);
    fprintf(' DONE (%g sec)\n',cputime-t);
    if writePoints
      % write initial points average
      writeMETALandmark([IhatP; IhatPWeights],sprintf('debug/IhatP_k%d_m%d.mhd',k,m));
    end
    
    %
    % update convergence measure (TODO: add landmarks convergence measure)
    %
    t = cputime;
    fprintf('Computing variance...');
    Ivar = var(J0T,0,ndims(I));
    avar(M*(k-1) + m + 1) = sum(Ivar(:));
    fprintf(' DONE (%g sec)\n',cputime-t);
    if writeDebugImages
      writeMETA(Ivar,sprintf('debug/Ihat_ptwise_variance_k%d_m%d.mhd',k,m));
    end
    
    fprintf('::::::::::::::: Sum of voxelwise variance %g ::: %g%% :::::::::::::::\n',...
	    avar(M*(k-1) + m + 1),...
	    100*avar(M*(k-1) + m + 1)/avar(1));
    fprintf('Iteration Time: %g\n',cputime-iterStartTime);
    
    if displayPlotOverview
      switch d
        case 2,
          vizLDMMState(figHandles,...
            J0T, Ihat0, Ihat, Ivar0, Ivar, JP0T, IhatP0, IhatP,...
            avar, epsilon,maxIter);
        case 3,
          vizLDMMState(figHandles,...
            squeeze(J0T(:,:,round(end/2),:)),...
            squeeze(Ihat0(:,:,round(end/2))),...
            squeeze(Ihat(:,:,round(end/2))),...
            squeeze(Ivar0(:,:,round(end/2))),...
            squeeze(Ivar(:,:,round(end/2))),...
            JP0T, IhatP0, IhatP, avar, epsilon, maxIter);
      end
    end
    %input('Press enter to continue...');
  end
  
  % write images at end of iteration
  if writeImages
    writeMETA(Ihat,sprintf('debug/Ihat_k%d.mhd',k));
  end
  if writeImages
    writeMETA(Ivar,sprintf('debug/Ihat_ptwise_variance_k%d.mhd',k));
  end
end

if writeFinalAtlasImage
  % write final atlas image
  writeMETA(Ihat,'debug/Ihat_final.mhd');
end

if writeDeformationFields
  % write final deformation field
  for m = 1:M
    defField = computeInverseHFieldFromVField(v(:,:,:,:,:,m), inverseMethod);
    writeMETA(defField,sprintf('debug/h_final_%d.mhd',m));
  end
end

fprintf('Total Time: %g (sec)\n',cputime-startTime);
fprintf('%0.3g, \n',avar/avar(1));
