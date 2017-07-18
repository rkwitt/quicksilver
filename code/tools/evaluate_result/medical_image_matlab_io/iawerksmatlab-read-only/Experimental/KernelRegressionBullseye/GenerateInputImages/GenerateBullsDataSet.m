clear all; close all;

% image size
imageDims     = [256 256];
minDim        = min(imageDims(:));

% granularity of time scale
numTimePoints = 10000;

% how many images should be produced
sampleSize    = 200;

% what range of radii (subset of [0,1]) should be used for each band
rMultRange    = [0.25 0.75];
% smallest band min radius (ratio of image width)
rInnerMin     = 1/8;
% largest band max radius (ratio of image width)
rOuterMax     = 4/9;

% noise in pixels added to radii
rNoiseSigma   = 2.5;

imageNoiseSigma = 0.05;
imageSmoothingSigma = 1.0;

%
% create time values
t             = linspace(0,1,numTimePoints);

%
% compute the inner, middle, and outer radii multipliers (which
% interpolate between min and max radius)
rMult(:,1) = f1(t');
rMult(:,2) = f2(t');
rMult(:,3) = f3(t');
rMult      = rMultRange(1) + rMult*(rMultRange(2)-rMultRange(1));

%
% plot the radius multipliers
%figure; plot(t,rMult);

%
% compute the actual radii in voxels
rMiddleMin    = rInnerMin+(rOuterMax-rInnerMin)/3;
rOuterMin     = rInnerMin+(rOuterMax-rInnerMin)*2/3;

innerMinR = rInnerMin*minDim;
innerMaxR = rMiddleMin*minDim;
r(:,1) = innerMinR + rMult(:,1)*(innerMaxR-innerMinR);

middleMinR = rMiddleMin*minDim;
middleMaxR = rOuterMin*minDim;
r(:,2) = middleMinR + rMult(:,2)*(middleMaxR-middleMinR);

outerMinR = rOuterMin*minDim;
outerMaxR = rOuterMax*minDim;
r(:,3) = outerMinR + rMult(:,3)*(outerMaxR-outerMinR);

%
% collect the sample and add noise
iSample = randperm(length(t));
iSample = sort(iSample(1:sampleSize));
tSample = t(iSample);
rSample = r(iSample,:);
rSample = rSample + rNoiseSigma * randn(size(rSample));

%
% plot the radii 
figure; plot(t, r);
hold on;
plot(tSample, rSample(:,1),'xb');
plot(tSample, rSample(:,2),'xg');
plot(tSample, rSample(:,3),'xr');
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
save('-binary',['RadiiDataNoise' num2str(rNoiseSigma) '.odf'], 'tSample','rSample');

