%
% Ihat is at time T
% input images are at time 0
%

%
% todo: 
%

%
% output
%
writeDebugImages = false;
writeImages = false;
dispMinMaxL2 = true;
extraOutput = true;

% the data type used for storage and computation
% use 'single' for 32-bit floating point
% use 'double' for 64-bit floating point
dtype = 'single';

% the number of LDMM iterations
maxIter = 50;
% the number of input images
M = 4;
% the number of timesteps (i.e. velocity fields)
N = 10;
% the number of dimensions 
d = 2;
% the extent of the test image
s = 128;

% sigma weights the smoothed velocity field in the gradient energy computation 
% don's use the actual word 'sigma' because it is a built in matlab function
sig = 0.08;

% standard deviation of random noise applied to image
imageNoiseSigma = 0.00;
% standard deviation of smoothing kernel applied to image
imageSmoothingSigma = 1;

% epsilon weights the gradient energy in the velocity update equation 
% eps=0.5 means replace current velocity with Linvb
% eps=0 means v stays the same
% initial development with 0.001
epsilon = 0.01;

% alpha, beta, and gamma determine the smoothing of the velocity fields
alpha = 0.5;
beta = 0;
gamma = 1;

%
% load images
%
fprintf('Loading images...');
t = cputime;
I = zeros([s s M],dtype);
for m=1:M
  fprintf('%d',m);
  %ds!!!
  I(:,:,m) = makeTestImage([s s],'sphere',[s*(m/M+1)/8 imageNoiseSigma imageSmoothingSigma],dtype);
  %ds!!!
  %I(:,:,m) = makeTestImage([s s],'sphere',s*(rand+1)/8,dtype);
end
fprintf(' DONE (%g sec)\n',cputime-t);
if writeDebugImages
  % write original images
  for q = 1:M
    %ds!!!
    writeMETA(squeeze(I(:,:,q)),sprintf('debug/input_image_%d.mhd',q));
  end
end

%
% compute initial average image
%
J0T = I;
t = cputime;
fprintf('Computing initial average image...');
Ihat = mean(J0T,ndims(I));
fprintf(' DONE (%g sec)\n',cputime-t);
% write initial average image
if writeImages
  writeMETA(Ihat,'debug/Ihat_k0.mhd');
end
fprintf('Computing initial variance...');
t = cputime;
Ivar = var(J0T,0,ndims(I));
% write initial variance image
if writeImages
  writeMETA(Ivar,'debug/Ihat_ptwise_variance_k0.mhd');
end
avar(1) = sum(Ivar(:));
fprintf(' DONE (%g sec)\n',cputime-t);
fprintf(' === Sum of voxelwise variance %g ::: %g%% === \n',avar(1),100*avar(1)/avar(1));

%
% initialize memory 
%
t = cputime;
fprintf('Allocating memory: ');
fprintf('v');
v = zeros([d s s N M],dtype);
%%%%%%%%%  V,X,Y,t,m  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% V:  the element of the vector (x or y)
% X, Y: the x and y position of the voxel
% t: the time point
% m: the input image index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% deformed image from time T to each time t
fprintf(',JTt');
JTt    = zeros([s s N],dtype);
% det of Jacobian of h-field mapping \Omega at t to T
fprintf(',dPhitT');
dPhitT = zeros([s s N],dtype);
% deformed image from time 0 to each time t
fprintf(',J0t');
J0t    = zeros([s s N+1],dtype);
% (spatial) gradient of J0t images
fprintf(',gJ0t');
gJ0t   = zeros([d s s N+1],dtype); 
% body force
fprintf(',b');
%b      = zeros([d s s N],dtype); 
% regularized body force
fprintf(',Linvb');
bLinvb  = zeros([d s s N],dtype); 
% gradient energy
fprintf(',gE');
gE     = zeros([d s s N],dtype); 
fprintf(' DONE (%g sec)\n',cputime-t);

startTime = cputime;
for k=1:maxIter
  for m=1:M
    iterStartTime = cputime;

    fprintf('================== iter: %d, image: %d, elapsed time: %g ================== \n',...
	    k,m,cputime-startTime);
    
    %
    % compute JTt (image (Ihat) at T deformed to each timepoint t)
    % and dPhitT (det jacobian of phi_t, the hfield from t to T)
    %
    t = cputime;
    fprintf('Computing (backward, T-->t) def. images & det. Jacobians...');
    [JTt, dPhitT] = computeJTt(Ihat, v(:,:,:,:,m)); 
    fprintf(' DONE (%g sec)\n',cputime-t);
    if writeDebugImages
      % write debug images
      for q = 1:N
        writeMETA(squeeze(JTt(:,:,q)),sprintf('debug/JTt_k%d_m%d_t%d.mhd',k,m,q));
        writeMETA(squeeze(dPhitT(:,:,q)),sprintf('debug/dPhitT_k%d_m%d_t%d.mhd',k,m,q));
      end
    end

    %
    % compute J0t (image at 0 deformed to each timepoint t)
    % 
    t = cputime;
    fprintf('Computing (forward, 0-->t) deformed images...');
    J0t = computeJ0t(I(:,:,m), v(:,:,:,:,m)); 
    fprintf(' DONE (%g sec)\n',cputime-t);
    if writeDebugImages
      % write debug images
      for q = 1:N+1
        writeMETA(squeeze(J0t(:,:,q)),sprintf('debug/J0t_k%d_m%d_t%d.mhd',k,m,q));
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
        writeMETA(squeeze(gJ0t(:,:,:,q)),sprintf('debug/gJ0t_k%d_m%d_t%d.mhd',k,m,q));
      end
    end

    %
    % compute "body force" b = dPhitT*(J0t-JTt) * gJ0t at each time step
    %
    t = cputime;
    fprintf('Computing body force...');    
    bLinvb(1,:,:,:) = dPhitT.*(J0t(:,:,1:N)-JTt) .* squeeze(gJ0t(1,:,:,1:N));
    bLinvb(2,:,:,:) = dPhitT.*(J0t(:,:,1:N)-JTt) .* squeeze(gJ0t(2,:,:,1:N));
    fprintf(' DONE (%g sec)\n',cputime-t);
    if dispMinMaxL2
      fprintf('NUM NAN = %g\n',sum(isnan(bLinvb(:))));
      fprintf('MAX = %g\n',max(bLinvb(:)));
      fprintf('MIN = %g\n',min(bLinvb(:)));
    end
    if writeDebugImages
      % write debug images
      for q = 1:N
        writeMETA(squeeze(bLinvb(:,:,:,q)),sprintf('debug/b_k%d_m%d_t%d.mhd',k,m,q));
      end
    end

    %
    % apply Greens function to regularize b
    %
    t = cputime;
    fprintf('Applying Green''s function...');
    for tm = 1:size(bLinvb,ndims(bLinvb))
      fprintf('%d',tm);
      bLinvb(:,:,:,tm) = greensFunction(squeeze(bLinvb(:,:,:,tm)),...
					 'laplace',[alpha beta gamma]); 
    end
    fprintf(' DONE (%g sec)\n',cputime-t);
    if dispMinMaxL2
      fprintf('NUM NAN = %g\n',sum(isnan(bLinvb(:))));
      fprintf('MAX = %g\n',max(bLinvb(:)));
      fprintf('MIN = %g\n',min(bLinvb(:)));
    end
    if writeDebugImages
      % write debug images
      for q = 1:N
        writeMETA(squeeze(bLinvb(:,:,:,q)),sprintf('debug/Linvb_k%d_m%d_t%d.mhd',k,m,q));
      end
    end
    
    % zero boundary condition
    bLinvb(:,[1 end],:,:) = 0;
    bLinvb(:,:,[1 end],:) = 0;

    %
    % compute gradient energy gE = 2*v - 2/sigma^2 * Linvb
    %
    t = cputime;
    fprintf('Computing gradient energy...');
    gE = 2*v(:,:,:,:,m) - (2/sig^2)*bLinvb;
    fprintf(' DONE (%g sec)\n',cputime-t);
    if dispMinMaxL2
      fprintf('NUM NAN = %g\n',sum(isnan(gE(:))));
      fprintf('MAX = %g\n',max(gE(:)));
      fprintf('MIN = %g\n',min(gE(:)));
    end
    if writeDebugImages
      for q = 1:N
        writeMETA(squeeze(gE(:,:,:,q)),sprintf('debug/gE_k%d_m%d_t%d.mhd',k,m,q));
      end
    end

    %
    % update velocity fields v = v - epsilon * gE
    %
    t = cputime;
    fprintf('Updating velocity fields...');
    v(:,:,:,:,m) = v(:,:,:,:,m) - epsilon*gE;
    fprintf(' DONE (%g sec)\n',cputime-t);
    if dispMinMaxL2
      fprintf('NUM NAN = %g\n',sum(isnan(v(:))));
      fprintf('MAX = %g\n',max(v(:)));
      fprintf('MIN = %g\n',min(v(:)));
    end
    if writeDebugImages    
      for q = 1:N
        writeMETA(squeeze(v(:,:,:,q)),sprintf('debug/v_k%d_m%d_t%d.mhd',k,m,q));
      end
    end

    %
    % update Ihat: need to compute J0T for this image
    %
    t = cputime;
    fprintf('Updating average image...');
    J0t = computeJ0t(I(:,:,m), v(:,:,:,:,m)); 
    J0T(:,:,m) = J0t(:,:,N+1);
    Ihat = mean(J0T,ndims(I));    
    fprintf(' DONE (%g sec)\n',cputime-t);
    if writeDebugImages
      writeMETA(Ihat,sprintf('debug/Ihat_k%d_m%d.mhd',k,m));
    end
    
    %
    % update convergence measure
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
  end
  
  % write images at end of iteration
  if writeImages
    writeMETA(Ihat,sprintf('debug/Ihat_k%d.mhd',k));
  end
  if writeImages
    writeMETA(Ivar,sprintf('debug/Ihat_ptwise_variance_k%d.mhd',k));
  end
  
end
fprintf('Total Time: %g\n',cputime-startTime);
fprintf('%0.3g, \n',avar/avar(1));
