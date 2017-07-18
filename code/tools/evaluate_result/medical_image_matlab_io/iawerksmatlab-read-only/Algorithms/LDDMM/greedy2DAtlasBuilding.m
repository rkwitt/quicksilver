%
% 2-D Greedy Deformable reg
%
% bcd 2005
%

%
% todo: search for todo
%

clear all; close all;

%
% algorithm global variables
%
M = 16;            % number of images
maxIter = 500;    % maximum number of iterations
maxPert = 0.5;    % the most someting can move during one iteration

%
% debug and visualization flags
%
debug   = true;
display = true;
skip    = 1;        % how many iterations to process before pausing
ps = ceil(sqrt(M)); % size of side for subplotting
sprange = reshape(1:ps*ps*2,ps*2,ps)';
avrange = sprange(end/2+1:end,1:end/2);
plrange = sprange(end/2+1:end,end/2+1:end);

%
% create images
%
I = createImageSample(M,'disk2',64);
imHeight = size(I,1);
imWidth = size(I,2);

%
% setup displays
%
if (display)
  InitialImages_handel  = figure;  % shows original images
  InitialAverage_handel = figure;  % original average image
  DeformedImages_handel = figure;  % shows deformed images
  
  figure(InitialImages_handel);
  for i=1:M
    subplot(ps,ps,i,'align'); imshow(I(:,:,i)); colormap('gray'); axis off; axis image; %title(sprintf('I%d',i));
  end

  figure(DeformedImages_handel);
  for i=1:M
    subplot(ps,ps*2,i,'align'); imshow(I(:,:,i)); colormap('gray'); title(sprintf('Def I%d',i)); axis off; axis image;
  end
  subplot(ps,ps*2,avrange(:),'align'); imshow(mean(I,3)); colormap('gray'); title('Average'); axis off; axis image;
  subplot(ps,ps*2,plrange(:),'align'); plot(0,0); axis([0 maxIter 0 1]); title('Var/(Initial Var) at Ihat');

  figure(InitialAverage_handel);
  imagesc(mean(I,3)); colormap gray; axis image; axis off;

  input('press enter to begin...');
end

%
% initialize algorithm variables
%
k = 0;                                   % iteration number
hx = zeros(imHeight, imWidth, M);        % transformation
hy = zeros(imHeight, imWidth, M);
vx = zeros(imHeight, imWidth);           % velocity field
vy = zeros(imHeight, imWidth);     
dhx = zeros(imHeight, imWidth);          % v(h(x))
dhy = zeros(imHeight, imWidth);     
bx = zeros(imHeight, imWidth);           % body force
by = zeros(imHeight, imWidth);      
gIx = zeros(imHeight, imWidth,M);        % gradients of moving image
gIy = zeros(imHeight, imWidth,M);         
J = zeros(imHeight, imWidth, M);         % deformed images
[idx,idy] = meshgrid(1:imWidth,1:imHeight);

%
% initialize deformed images and transformations
%
J = I;
Ihat = mean(J,3);
hx = repmat(idx, [1,1,M]);
hy = repmat(idy, [1,1,M]);
%figure; quiver(hx(:,:,1), hy(:,:,1));

%
% compute gradients
%
[gIx,gIy] = imageGradient(J);

%
% compute original variation
%
origAtlasVariance = sum(var(J,0,3));
atlasVariance = [];

for k=1:maxIter       
  % iterate over images 
  for m=1:M 
    fprintf('iter: %d/%d, image: %d/%d\n',k,maxIter,m,M);
    
    %
    % compute body force
    %
    bx = (Ihat - J(:,:,m)) .* gIx(:,:,m);
    by = (Ihat - J(:,:,m)) .* gIy(:,:,m);
    %figure; quiver(bx,by); title('body force'); axis equal;

    %
    % compute velocity Lv=b
    %
    [vx, vy] = greensFunction(bx,by,'fluid');    
    %figure; quiver(vx,vy); title('velocity'); axis equal;

    if k == 1 && m == 1
      %
      % enforce maximum perurbation
      %
      mpert = sqrt(max(max(vx .* vx + vy .* vy)));
      delta = maxPert / mpert;
      fprintf('delta = %g\n',delta);       
    end

    % 
    % update transformation for this image
    % hx = vx(hx)
    [dhx,dhy] = composeH(hx(:,:,m),hy(:,:,m),delta*vx,delta*vy);
    hx(:,:,m) = hx(:,:,m) + dhx;
    hy(:,:,m) = hy(:,:,m) + dhy;
    %figure; quiver(hx(:,:,m),hy(:,:,m)); title('xform'); axis equal;

    %
    % update deformed image
    %
    J(:,:,m) = deformImage(I(:,:,m), hx(:,:,m), hy(:,:,m));

    %
    % update Ihat
    %
    Ihat = mean(J,3);

    if m==M
      %
      % update error
      %
      newAtlasVariance = sum(var(J,0,3));
      atlasVariance(k) = newAtlasVariance/origAtlasVariance;

      %
      % update figures
      %
      figure(DeformedImages_handel);
      for i=1:M
        subplot(ps,ps*2,i,'align'); imshow(J(:,:,i)); 
        colormap('gray'); title(sprintf('Def I%d',i)); axis off; axis image;
      end
      subplot(ps,ps*2,avrange(:),'align'); imshow(Ihat); 
      colormap('gray'); title('Average'); axis off; axis image;
      subplot(ps,ps*2,plrange(:),'align'); plot(0:k,[1 atlasVariance]); 
      axis([0 maxIter 0 1]); title('Var/(Initial Var) at Ihat');
    end

    %
    % pause and wait for user
    %
    if (m==M & mod(k,skip) == 0)
      newSkip = input('Press enter to continue...');
      if (size(newSkip) == 1)
        skip = newSkip;
      end
    end
  end
end



