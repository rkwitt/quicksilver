clear all; close all;

% image size
imageDims     = [256 256];
minDim        = min(imageDims(:));

% granularity of time scale
numTimePoints = 10000;

% how many images should be produced
sampleSize    = 100

% noise in pixels added to radii
rNoiseSigmaCombined     = 6.0;
rNoiseSigmaIndividual   = 3.0;

% image intensity noise and smoothing
imageNoiseSigma = 0.3;
imageSmoothingSigma = 1.0;

%
% create time values
t             = linspace(0,1,numTimePoints);

%
% compute the inner, middle, and outer radii multipliers (which
% interpolate between min and max radius)
r = tToRadiiPixels_ICCV2008(t',imageDims);

%
% collect the sample and add noise
iSample = randperm(length(t));
iSample = sort(iSample(1:sampleSize));
tSample = t(iSample);
rSample = r(iSample,:);
rSample = rSample ...
          + repmat(rNoiseSigmaCombined * randn(size(rSample(:,1))), [1 size(rSample,2)]);
rSample = rSample ...
          + rNoiseSigmaIndividual * randn(size(rSample));

%
% plot the radii 
figure; hold on;
plot(t, r(:,1),'b','LineWidth',4);
plot(t, r(:,2),'g','LineWidth',4);
plot(t, r(:,3),'r','LineWidth',4);
plot(tSample, rSample(:,1),'xb','MarkerSize',2,'LineWidth',2);
plot(tSample, rSample(:,2),'og','MarkerSize',2,'LineWidth',2);
plot(tSample, rSample(:,3),'<r','MarkerSize',2,'LineWidth',2);
title('Bullseye Disk Radii');
xlabel('t');
ylabel('Radius in Pixels');
print(['RadiiCurvesNoise' num2str(rNoiseSigmaCombined) '_' num2str(rNoiseSigmaIndividual) '.png'],'-dpng');

%
% show a histogram of the sample t values to make sure they look
% uniform across [0,1]
figure; hist(tSample,10);

%
% generate the images
I = zeros([imageDims sampleSize]);
for imageIndex = 1:sampleSize
  I(:,:,imageIndex) = MakeBullImage(imageDims, rSample(imageIndex,:)/minDim, imageNoiseSigma, imageSmoothingSigma, 'single');
  writeMETA(squeeze(I(:,:,imageIndex)),sprintf('Bull%04d.mhd',imageIndex), 'MET_FLOAT',[0 0],[1 1], {'RNoiseSigmaCombined',sprintf('%g',rNoiseSigmaCombined),'RNoiseSigmaIndividual',sprintf('%g',rNoiseSigmaIndividual),'R1',sprintf('%g',rSample(imageIndex,1)),'R2',sprintf('%g',rSample(imageIndex,2)),'R3',sprintf('%g',rSample(imageIndex,3))});  
end

%
% save radii data
save(['RadiiDataNoiseIJCV2008.mat'], 't','r','tSample','rSample','rNoiseSigmaCombined','rNoiseSigmaIndividual');

