clear all; close all;

jobParcel1 = false;
jobParcel2 = true;

%
% this was estimated from the data
bandwidth = 0.0682;
kernelCutoff = 2;

%
% optimization parameters
ldmmOptParams.maxIterations = 20;    % !!!
ldmmOptParams.numTimesteps  = 10;
ldmmOptParams.sigma         = 0.08;
ldmmOptParams.epsilon       = 0.01;

%
% load images
sampleStride  = 2; % !!!
sampleSize    = 200 / sampleStride;
imageDims     = [256 256];
I = zeros([imageDims sampleSize],'single');
for i = 1:sampleSize
  I(:,:,i) = loadMETA2(sprintf('../InputImagesNoise25/Bull%04d.mhd',1+(i-1)*sampleStride), 'single');
end

%
% load time information
%load('-binary','../InputImagesNoise25/RadiiDataNoise25.odf');
load('../InputImagesNoise25/RadiiDataNoise25.mat');
tSample = tSample(1:sampleStride:end);
rSample = rSample(1:sampleStride:end,:);

%
% estimation timepoints
numTimePoints = 20; % !!!
t             = linspace(0,1,numTimePoints);

%
% if running 2 jobs, choose piece of the work
tIndices  = 1:length(t);
tIndices1 = tIndices(1:round(end/2));
tIndices2 = tIndices(1+round(end/2):end);
if jobParcel1 & jobParcel2
  tIndices = tIndices;
elseif jobParcel1
  tIndices = tIndices1;
elseif jobParcel2
  tIndices = tIndices2;
end

%
% do the regression
figure; hold on;
for ti = tIndices
  fprintf('Generating mean %d/%d\n', ti, length(t));
  
  %
  % select data that fits within n standard deviations
  tMin = t(ti) - kernelCutoff*bandwidth;
  tMax = t(ti) + kernelCutoff*bandwidth;
  
  cohortMask   = (tSample > tMin) & (tSample < tMax);
  cohortTimes  = tSample(cohortMask);
  cohortRadii  = rSample(cohortMask,:);
  cohortImages = I(:,:,cohortMask);
  cohortSize   = sum(cohortMask(:));

  %
  % compute normalized weights for cohort data
  cohortWeights = zeros(1,cohortSize);
  for ci = 1:cohortSize
    cohortWeights(ci) = regkernel(t(ti), cohortTimes(ci), bandwidth);
  end
  cohortWeights = cohortWeights / sum(cohortWeights(:));
  plot(cohortTimes,cohortWeights);

  %
  % run lddmm average
  [Ihat, Ivar, VarLog] = fldmm2D(cohortImages, cohortWeights, ldmmOptParams);

  %
  % save average image
  r = tToRadiiPixels(t(ti), imageDims);
  writeMETA(Ihat,sprintf('Mean_ti%03d.mhd',ti),'MET_FLOAT',[0 0],[1 1],{'T',sprintf('%g',t(ti)),'CohortSize',sprintf('%d',cohortSize),'R1',sprintf('%g',r(1)),'R2',sprintf('%g',r(2)),'R3',sprintf('%g',r(3))});
  writeMETA(Ivar,sprintf('Var_ti%03d.mhd',ti),'MET_FLOAT',[0 0],[1 1],{'T',sprintf('%g',t(ti)),'CohortSize',sprintf('%d',cohortSize),'R1',sprintf('%g',r(1)),'R2',sprintf('%g',r(2)),'R3',sprintf('%g',r(3))});

  %
  % display image with circles drawn on top
  clf; hold on;
  imagesc(Ihat); axis square; axis xy; axis off; colormap gray;
  circleT = linspace(0,2*pi,500);
  plot(imageDims(1)/2+r(1)*cos(circleT),imageDims(2)/2+r(1)*sin(circleT),'b','LineWidth',3);
  plot(imageDims(1)/2+r(2)*cos(circleT),imageDims(2)/2+r(2)*sin(circleT),'g','LineWidth',3);
  plot(imageDims(1)/2+r(3)*cos(circleT),imageDims(2)/2+r(3)*sin(circleT),'r','LineWidth',3);
  plotTitle = sprintf('Mean_ti%03d.png',ti);
  print(gcf,'-dpng',plotTitle);
  
  %
  % display image with circles drawn on top
  clf; hold on;
  imagesc(Ivar); axis square; axis xy; axis off; colormap gray;
  circleT = linspace(0,2*pi,500);
  plot(imageDims(1)/2+r(1)*cos(circleT),imageDims(2)/2+r(1)*sin(circleT),'b','LineWidth',3);
  plot(imageDims(1)/2+r(2)*cos(circleT),imageDims(2)/2+r(2)*sin(circleT),'g','LineWidth',3);
  plot(imageDims(1)/2+r(3)*cos(circleT),imageDims(2)/2+r(3)*sin(circleT),'r','LineWidth',3);
  plotTitle = sprintf('Var_ti%03d.png',ti);
  print(gcf,'-dpng',plotTitle);  
  
  %
  % plot and save optimization record
  clf;
  plot(VarLog);
  title(sprintf('ISE Log t=%g',t(ti)));
  xlabel('Iteration');
  ylabel('ISE');
  plotTitle = sprintf('VarLog_ti%03d.png',ti);
  print(gcf,'-dpng',plotTitle);
end
