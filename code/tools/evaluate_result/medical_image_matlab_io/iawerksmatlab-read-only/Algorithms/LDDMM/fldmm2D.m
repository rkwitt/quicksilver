function [Ihat, Ivar, VarLog] = fldmm2D(I, K, optParams, fluidParams, writeParams)

% octave vs. matlab
if isMatlab
  stderr = 2;
end

% optParams is 
%   maxIterations
%   numTimesteps
%   sigma
%   epsilon
%   inverseMethod

% fluidParams is
%   alpha
%   beta
%   gamma

% Ihat is at time T
% input images are at time 0

% "global" variables

% the number of images
M = size(I,3);

% image dimensions
numDim = 2;
Dim = [size(I,1) size(I,2)];

% the data type used for storage and computation
% use 'single' for 32-bit floating point
% use 'double' for 64-bit floating point
dtype = 'single'

%
% output
%
if (nargin < 5)
  writeParams.writeDebugImages = false;
  writeParams.writeImages      = false;
  writeParams.dispMinMaxL2     = false;
  writeParams.extraOutput      = false;
  writeParams.filePrefix       = '';
end

if (nargin < 2)
  K = ones(M,1)/M;
end

KWeights = zeros(size(I));
for ii = 1:M
  KWeights(:,:,ii) = K(ii);
end
squeeze(KWeights(1,1,:))
sum(squeeze(KWeights(1,1,:)))

if (nargin < 3)
  optParams.maxIterations = 50;
  optParams.numTimesteps  = 10;
  optParams.sigma         = 0.08;
  optParams.epsilon       = 0.01;
  optParams.inverseMethod = 'zero';
end

if (nargin < 4)
  fluidParams.alpha = 0.5;
  fluidParams.beta  = 0.0;
  fluidParams.gamma = 1.0;
end

%
% Debug: write original images
if writeParams.writeImages
  fprintf(stderr,'DEBUG: writing original images...');
  t = cputime;
  for q = 1:M
    writeMETA(squeeze(I(:,:,q)),sprintf('debug/%sinput_image_%d.mhd',writeParams.filePrefix,q));
  end
  fprintf(stderr,' DONE (%g sec)\n',cputime-t);
end

%
% compute initial average image
J0T = I;
fprintf(stderr,'Computing initial average image...');
t = cputime;
Ihat = sum(J0T .* KWeights, ndims(I));
fprintf(stderr,' DONE (%g sec)\n',cputime-t);

%
% write initial average image
if writeParams.writeImages
  fprintf(stderr,'DEBUG: writing initial average image...');
  t = cputime;
  writeMETA(Ihat,sprintf('debug/%sIhat_k0.mhd',writeParams.filePrefix));
  fprintf(stderr,' DONE (%g sec)\n',cputime-t);
end

%
% compute initial variance
fprintf(stderr,'Computing initial variance...');
t = cputime;
Ivar = sum((J0T - repmat(Ihat,[1 1 M])).^2 .* KWeights, ndims(I));
avar(1) = sum(Ivar(:));
fprintf(stderr,' DONE (%g sec)\n',cputime-t);
fprintf(stderr,' === Sum of voxelwise variance %g ::: %g%% === \n',avar(1),100*avar(1)/avar(1));

%
% write initial variance image
if writeParams.writeImages
  fprintf(stderr,'DEBUG: writing initial variance image...');
  t = cputime;
  writeMETA(Ivar,sprintf('debug/%sIhat_ptwise_variance_k0.mhd',writeParams.filePrefix));
  fprintf(stderr,' DONE (%g sec)\n',cputime-t);
end

%
% initialize memory 
%
t = cputime;
fprintf(stderr,'Allocating memory: ');

% the velocity fields
%%%%%%%%%  V,X,Y,t,m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% V:  the element of the vector (x or y)
% X, Y: the x and y position of the voxel
% t: the time point
% m: the input image index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(stderr,'v');
v = zeros([numDim Dim optParams.numTimesteps M],dtype);

% deformed image from time T to each time t
fprintf(stderr,',JTt');
JTt    = zeros([Dim optParams.numTimesteps],dtype);

% det of Jacobian of h-field mapping \Omega at t to T
fprintf(stderr,',dPhitT');
dPhitT = zeros([Dim optParams.numTimesteps],dtype);

% deformed image from time 0 to each time t
fprintf(stderr,',J0t');
J0t    = zeros([Dim optParams.numTimesteps+1],dtype);

% (spatial) gradient of J0t images
fprintf(stderr,',gJ0t');
gJ0t   = zeros([numDim Dim optParams.numTimesteps+1],dtype); 

% body force
fprintf(stderr,',b');
%b      = zeros([d s s N],dtype); 

% regularized body force
fprintf(stderr,',Linvb');
bLinvb  = zeros([numDim Dim optParams.numTimesteps],dtype); 

% gradient energy
fprintf(stderr,',gE');
gE      = zeros([numDim Dim optParams.numTimesteps],dtype); 
fprintf(stderr,' DONE (%g sec)\n',cputime-t);

startTime = cputime;
for k=1:optParams.maxIterations
  for m=1:M
    iterStartTime = cputime;

    fprintf(stderr,'================== iter: %d/%d, image: %d/%d, elapsed time: %g ================== \n',...
	    k,optParams.maxIterations,m,M,cputime-startTime);
    
    %
    % compute JTt (image (Ihat) at T deformed to each timepoint t)
    % and dPhitT (det jacobian of phi_t, the hfield from t to T)
    %
    t = cputime;
    fprintf(stderr,'Computing (backward, T-->t) def. images & det. Jacobians...');
    [JTt, dPhitT] = computeJTt(Ihat, v(:,:,:,:,m)); 
    fprintf(stderr,' DONE (%g sec)\n',cputime-t);
    if writeParams.dispMinMaxL2
      fprintf(stderr,'NUM NAN = %g\n',sum(isnan(JTt(:))));
      fprintf(stderr,'MAX = %g\n',max(JTt(:)));
      fprintf(stderr,'MIN = %g\n',min(JTt(:)));
    end
    if writeParams.dispMinMaxL2
      fprintf(stderr,'NUM NAN = %g\n',sum(isnan(dPhitT(:))));
      fprintf(stderr,'MAX = %g\n',max(dPhitT(:)));
      fprintf(stderr,'MIN = %g\n',min(dPhitT(:)));
    end
    if writeParams.writeDebugImages
      % write debug images
      for q = 1:optParams.numTimesteps
        writeMETA(squeeze(JTt(:,:,q)),sprintf('debug/%sJTt_k%d_m%d_t%d.mhd',writeParams.filePrefix,k,m,q));
        writeMETA(squeeze(dPhitT(:,:,q)),sprintf('debug/%sdPhitT_k%d_m%d_t%d.mhd',writeParams.filePrefix,k,m,q));
      end
    end

    %
    % compute J0t (image at 0 deformed to each timepoint t)
    % 
    t = cputime;
    fprintf(stderr,'Computing (forward, 0-->t) deformed images...');
    J0t = computeJ0t(I(:,:,m), v(:,:,:,:,m), optParams.inverseMethod); 
    fprintf(stderr,' DONE (%g sec)\n',cputime-t);
    if writeParams.dispMinMaxL2
      fprintf(stderr,'NUM NAN = %g\n',sum(isnan(J0t(:))));
      fprintf(stderr,'MAX = %g\n',max(J0t(:)));
      fprintf(stderr,'MIN = %g\n',min(J0t(:)));
    end
    if writeParams.writeDebugImages
      % write debug images
      for q = 1:optParams.numTimesteps+1
        writeMETA(squeeze(J0t(:,:,q)),sprintf('debug/%sJ0t_k%d_m%d_t%d.mhd',writeParams.filePrefix,k,m,q));
      end
    end
    
    %
    % compute gJ0t (gradient of each deformed (0->t) image)
    %
    t = cputime;
    fprintf(stderr,'Computing gradients of (forward) deformed images...');
    gJ0t = computeImageGradients(J0t);
    fprintf(stderr,' DONE (%g sec)\n',cputime-t);
    if writeParams.dispMinMaxL2
      fprintf(stderr,'NUM NAN = %g\n',sum(isnan(gJ0t(:))));
      fprintf(stderr,'MAX = %g\n',max(gJ0t(:)));
      fprintf(stderr,'MIN = %g\n',min(gJ0t(:)));
    end
    if writeParams.writeDebugImages
      % write debug images
      for q = 1:optParams.numTimesteps+1
        writeMETA(squeeze(gJ0t(:,:,:,q)),sprintf('debug/%sgJ0t_k%d_m%d_t%d.mhd',writeParams.filePrefix,k,m,q));
      end
    end

    %
    % compute "body force" b = dPhitT*(J0t-JTt) * gJ0t at each time step
    %
    t = cputime;
    fprintf(stderr,'Computing body force...');    
    bLinvb(1,:,:,:) = dPhitT.*(J0t(:,:,1:optParams.numTimesteps)-JTt) .* squeeze(gJ0t(1,:,:,1:optParams.numTimesteps));
    bLinvb(2,:,:,:) = dPhitT.*(J0t(:,:,1:optParams.numTimesteps)-JTt) .* squeeze(gJ0t(2,:,:,1:optParams.numTimesteps));
    fprintf(stderr,' DONE (%g sec)\n',cputime-t);
    if writeParams.dispMinMaxL2
      fprintf(stderr,'NUM NAN = %g\n',sum(isnan(bLinvb(:))));
      fprintf(stderr,'MAX = %g\n',max(bLinvb(:)));
      fprintf(stderr,'MIN = %g\n',min(bLinvb(:)));
    end
    if writeParams.writeDebugImages
      % write debug images
      for q = 1:optParams.numTimesteps
        writeMETA(squeeze(bLinvb(:,:,:,q)),sprintf('debug/%sb_k%d_m%d_t%d.mhd',writeParams.filePrefix,k,m,q));
      end
    end

    %
    % apply Greens function to regularize b
    %
    t = cputime;
    fprintf(stderr,'Applying Green''s function...');
    for tm = 1:size(bLinvb,ndims(bLinvb))
      fprintf(stderr,'%d',tm);
      bLinvb(:,:,:,tm) = greensFunction(squeeze(bLinvb(:,:,:,tm)),...
					 'laplace',[fluidParams.alpha fluidParams.beta fluidParams.gamma]); 
    end
    fprintf(stderr,' DONE (%g sec)\n',cputime-t);
    if writeParams.dispMinMaxL2
      fprintf(stderr,'NUM NAN = %g\n',sum(isnan(bLinvb(:))));
      fprintf(stderr,'MAX = %g\n',max(bLinvb(:)));
      fprintf(stderr,'MIN = %g\n',min(bLinvb(:)));
    end
    if writeParams.writeDebugImages
      % write debug images
      for q = 1:optParams.numTimesteps
        writeMETA(squeeze(bLinvb(:,:,:,q)),sprintf('debug/%sLinvb_k%d_m%d_t%d.mhd',writeParams.filePrefix,k,m,q));
      end
    end
    
    % zero boundary condition
    bLinvb(:,[1 end],:,:) = 0;
    bLinvb(:,:,[1 end],:) = 0;

    %
    % compute gradient energy gE = 2*v - 2/sigma^2 * Linvb
    %
    t = cputime;
    fprintf(stderr,'Computing gradient energy...');
    gE = 2*v(:,:,:,:,m) - (2/optParams.sigma^2)*bLinvb;
    fprintf(stderr,' DONE (%g sec)\n',cputime-t);
    if writeParams.dispMinMaxL2
      fprintf(stderr,'NUM NAN = %g\n',sum(isnan(gE(:))));
      fprintf(stderr,'MAX = %g\n',max(gE(:)));
      fprintf(stderr,'MIN = %g\n',min(gE(:)));
    end
    if writeParams.writeDebugImages
      for q = 1:optParams.numTimesteps
        writeMETA(squeeze(gE(:,:,:,q)),sprintf('debug/%sgE_k%d_m%d_t%d.mhd',writeParams.filePrefix,k,m,q));
      end
    end

    %
    % update velocity fields v = v - epsilon * gE
    %
    t = cputime;
    fprintf(stderr,'Updating velocity fields...');
    v(:,:,:,:,m) = v(:,:,:,:,m) - optParams.epsilon*gE;
    fprintf(stderr,' DONE (%g sec)\n',cputime-t);
    if writeParams.dispMinMaxL2
      fprintf(stderr,'NUM NAN = %g\n',sum(isnan(v(:))));
      fprintf(stderr,'MAX = %g\n',max(v(:)));
      fprintf(stderr,'MIN = %g\n',min(v(:)));
    end
    if writeParams.writeDebugImages    
      for q = 1:optParams.numTimesteps
        writeMETA(squeeze(v(:,:,:,q)),sprintf('debug/%sv_k%d_m%d_t%d.mhd',writeParams.filePrefix,k,m,q));
      end
    end

    %
    % update Ihat: need to compute J0T for this image
    %
    t = cputime;
    fprintf(stderr,'Updating average image...');
    J0t = computeJ0t(I(:,:,m), v(:,:,:,:,m), optParams.inverseMethod); 
    J0T(:,:,m) = J0t(:,:,optParams.numTimesteps+1);
    Ihat = sum(J0T .* KWeights, ndims(I));    
    fprintf(stderr,' DONE (%g sec)\n',cputime-t);
    if writeParams.writeDebugImages
      writeMETA(Ihat,sprintf('debug/%sIhat_k%d_m%d.mhd',writeParams.filePrefix,k,m));
    end
    
    %
    % update convergence measure
    %
    t = cputime;
    fprintf(stderr,'Computing variance...');
    Ivar = sum((J0T - repmat(Ihat,[1 1 M])).^2 .* KWeights, ndims(I));
    avar(M*(k-1) + m + 1) = sum(Ivar(:));
    fprintf(stderr,' DONE (%g sec)\n',cputime-t);
    if writeParams.writeDebugImages
      writeMETA(Ivar,sprintf('debug/%sIhat_ptwise_variance_k%d_m%d.mhd',writeParams.filePrefix,k,m));
    end
    
    fprintf(stderr,'::::::::::::::: Sum of voxelwise variance %g ::: %g%% :::::::::::::::\n',...
	    avar(M*(k-1) + m + 1),...
	    100*avar(M*(k-1) + m + 1)/avar(1));
    fprintf(stderr,'Iteration Time: %g\n',cputime-iterStartTime);
  end
  
  % write images at end of iteration
  if writeParams.writeImages
    writeMETA(Ihat,sprintf('debug/%sIhat_k%d.mhd',writeParams.filePrefix,k));
  end
  if writeParams.writeImages
    writeMETA(Ivar,sprintf('debug/%sIhat_ptwise_variance_k%d.mhd',writeParams.filePrefix,k));
  end
  
end

% write final deformed images
if writeParams.writeImages
  fprintf(stderr,'DEBUG: writing final deformed images...');
  t = cputime;
  for q = 1:M
    writeMETA(squeeze(J0T(:,:,q)),sprintf('debug/%sdeformed_image_%d.mhd',writeParams.filePrefix,q));
  end
  fprintf(stderr,' DONE (%g sec)\n',cputime-t);
end

fprintf(stderr,'Total Time: %g\n',cputime-startTime);
fprintf(stderr,'%0.3g, \n',avar/avar(1));
VarLog = avar;
