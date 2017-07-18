clear all; close all;

%
% load images
sampleStride  = 1; % !!!
sampleSize    = 200 / sampleStride;
imageDims     = [256 256];
I = zeros([imageDims sampleSize],'single');
for i = 1:sampleSize
  I(:,:,i) = loadMETA2(sprintf('../InputImagesNoise25/Bull%04d.mhd',1+(i-1)*sampleStride), 'single');
end

%
% load time information
load('../InputImagesNoise25/RadiiDataNoise25.mat');
tSample = tSample(1:sampleStride:end);
rSample = rSample(1:sampleStride:end,:);


for ti = 1:length(tSample)
  %
  % display image with circles drawn on top
  clf; hold on;
  imagesc(I(:,:,ti)); axis square; axis xy; axis off; colormap gray;
  %r = tToRadiiPixels(tSample(ti), imageDims);
  %circleT = linspace(0,2*pi,500);
  %plot(imageDims(1)/2+r(1)*cos(circleT),imageDims(2)/2+r(1)*sin(circleT),'b','LineWidth',3);
  %plot(imageDims(1)/2+r(2)*cos(circleT),imageDims(2)/2+r(2)*sin(circleT),'g','LineWidth',3);
  %plot(imageDims(1)/2+r(3)*cos(circleT),imageDims(2)/2+r(3)*sin(circleT),'r','LineWidth',3);
  plotTitle = sprintf('BullInput_t%03d.png',ti);
  print(gcf,'-dpng',plotTitle);
end
