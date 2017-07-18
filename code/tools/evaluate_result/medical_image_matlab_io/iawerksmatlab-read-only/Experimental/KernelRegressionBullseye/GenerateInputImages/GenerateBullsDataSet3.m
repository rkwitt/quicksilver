clear all; close all;

% image size
imageDims     = [256 256];
minDim        = min(imageDims(:));

% granularity of time scale
numTimePoints = 10000;

% how many images should be produced
sampleSize    = 200;

% noise in pixels added to radii
rNoiseSigma   = 3.0;

imageNoiseSigma = 0.3;
imageSmoothingSigma = 1.0;

%
% create time values
t             = linspace(0,1,numTimePoints);

%
% compute the inner, middle, and outer radii multipliers (which
% interpolate between min and max radius)
r = tToRadiiPixels(t',imageDims);

%
% collect the sample and add noise
iSample = randperm(length(t));
iSample = sort(iSample(1:sampleSize));
tSample = t(iSample);
rSample = r(iSample,:);
rSample = rSample + rNoiseSigma * randn(size(rSample));

%
% plot the radii 
figure; hold on;
plot(t, r(:,1),'b','LineWidth',3);
plot(t, r(:,2),'g','LineWidth',3);
plot(t, r(:,3),'r','LineWidth',3);
plot(tSample, rSample(:,1),'xb','MarkerSize',10,'LineWidth',2);
plot(tSample, rSample(:,2),'xg','MarkerSize',10,'LineWidth',2);
plot(tSample, rSample(:,3),'xr','MarkerSize',10,'LineWidth',2);
title('Bullseye Disk Radii');
xlabel('t');
ylabel('Radius in Pixels');
print(['RadiiCurvesNoise' num2str(rNoiseSigma) '.png'],'-dpng');

%
% show a histogram of the sample t values to make sure they look
% uniform across [0,1]
figure; hist(tSample);

%
% generate the images
I = zeros([imageDims sampleSize]);
for imageIndex = 1:sampleSize
  I(:,:,imageIndex) = MakeBullImage(imageDims, rSample(imageIndex,:)/minDim, imageNoiseSigma, imageSmoothingSigma, 'single');
  writeMETA(squeeze(I(:,:,imageIndex)),sprintf('Bull%04d.mhd',imageIndex), 'MET_FLOAT',[0 0],[1 1], {'RNoiseSigma',sprintf('%g',rNoiseSigma),'R1',sprintf('%g',rSample(imageIndex,1)),'R2',sprintf('%g',rSample(imageIndex,2)),'R3',sprintf('%g',rSample(imageIndex,3))});  
end

%
% save radii data
save(['RadiiDataNoise' num2str(rNoiseSigma) '.mat'], 'tSample','rSample');

