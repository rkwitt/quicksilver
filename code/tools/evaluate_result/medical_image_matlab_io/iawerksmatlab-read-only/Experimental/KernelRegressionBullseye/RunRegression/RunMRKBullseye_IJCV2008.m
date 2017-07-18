function RunMRKBullseye_IJCV2008(partID, maxJobs, outFilePrefix)

if nargin < 2
  partID = 1;
  maxJobs = 1;
end

if nargin < 3
  outFilePrefix = '';
end

if isMatlab
  stderr = 2;
end

jobParcel  = zeros(1,maxJobs);
jobParcel(partID) = 1;

%
% this was estimated from the data
bandwidth = 0.045;                        % 0.045
kernelCutoff = 2;                       % 2

%
% optimization parameters
ldmmOptParams.maxIterations = 25;       % 20/30
ldmmOptParams.numTimesteps  = 12;       % 10
ldmmOptParams.sigma         = 0.08;    % 0.08
ldmmOptParams.epsilon       = 0.01;    % 0.01
ldmmOptParams.inverseMethod = 'zero';  % 'zero'/'iter'

%
% fluid params
ldmmFluidParams.alpha = 0.5;
ldmmFluidParams.beta  = 0.0;
ldmmFluidParams.gamma = 1.0;

%
% algorithm write params
ldmmWriteParams.writeDebugImages = false;
ldmmWriteParams.writeImages      = false;
ldmmWriteParams.dispMinMaxL2     = false;
ldmmWriteParams.extraOutput      = false;
ldmmWriteParams.filePrefix       = outFilePrefix;

%
% load images
sampleStride  = 1; 
sampleSize    = 100 / sampleStride;
imageDims     = [256 256];
I = zeros([imageDims sampleSize],'single');
fprintf(stderr, 'Loading images...');
for i = 1:sampleSize
  I(:,:,i) = loadMETA2(sprintf('./Bull%04d.mhd',1+(i-1)*sampleStride), 'single');
end
fprintf(stderr, 'DONE\n');

%
% load time information
%load('-binary','../InputImagesNoise25/RadiiDataNoise25.odf');
fprintf(stderr, 'Loading time and radii data...');
load('RadiiDataNoiseIJCV2008_mat.mat');
tSample = tSample(1:sampleStride:end);
rSample = rSample(1:sampleStride:end,:);
fprintf(stderr, 'DONE\n');

%
% estimation timepoints
numTimePoints = 8;                    % 30
t             = linspace(0,1,numTimePoints);

%
% if running multiple jobs, choose piece of the work
tIndices  = 1:length(t);
tIndicesJob = [];
for ji = 1:maxJobs
  if jobParcel(ji)
    tIndicesJob = [tIndicesJob tIndices(1+round((ji-1)*end/maxJobs):round(ji*end/maxJobs))]
  end
end
tIndices = tIndicesJob;

fprintf(stderr, 'Job indices: ');
fprintf(stderr, '%d ', tIndices);
fprintf(stderr, '\n');

%
% do the regression
figure; hold on;
for ti = tIndices
  outFilePrefixCurrent = '';
  if length(outFilePrefix)
    outFilePrefixCurrent = sprintf('%s%03d',outFilePrefix, ti);
    ldmmWriteParams.filePrefix = outFilePrefixCurrent;
  end

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

  fprintf(stderr, 'Generating mean %d/%d; Cohort size: %d\n', ti, length(t), cohortSize);

  %
  % run lddmm average
  try 
    [Ihat, Ivar, VarLog] = fldmm2D(cohortImages, cohortWeights, ldmmOptParams, ldmmFluidParams,ldmmWriteParams);
    fprintf(stderr, 'Mean generated successfully.');
  catch
    fprintf(stderr, 'Error generating mean!!!');
    me = lasterror
    continue;
  end

  %
  % save average image
  r = tToRadiiPixels_IJCV2008(t(ti), imageDims);
  writeMETA(Ihat,sprintf('%sMean_ti%03d.mhd',outFilePrefixCurrent,ti),'MET_FLOAT',[0 0],[1 1],{'T',sprintf('%g',t(ti)),'CohortSize',sprintf('%d',cohortSize),'R1',sprintf('%g',r(1)),'R2',sprintf('%g',r(2)),'R3',sprintf('%g',r(3))});
  writeMETA(Ivar,sprintf('%sVar_ti%03d.mhd',outFilePrefixCurrent,ti),'MET_FLOAT',[0 0],[1 1],{'T',sprintf('%g',t(ti)),'CohortSize',sprintf('%d',cohortSize),'R1',sprintf('%g',r(1)),'R2',sprintf('%g',r(2)),'R3',sprintf('%g',r(3))});

  %
  % display image with circles drawn on top
  clf; hold on;
  imagesc(Ihat); axis square; axis xy; axis off; colormap gray;
  circleT = linspace(0,2*pi,500);
  plot(imageDims(1)/2+r(1)*cos(circleT),imageDims(2)/2+r(1)*sin(circleT),'b','LineWidth',3);
  plot(imageDims(1)/2+r(2)*cos(circleT),imageDims(2)/2+r(2)*sin(circleT),'g','LineWidth',3);
  plot(imageDims(1)/2+r(3)*cos(circleT),imageDims(2)/2+r(3)*sin(circleT),'r','LineWidth',3);
  plotTitle = sprintf('%sMean_ti%03d.png',outFilePrefixCurrent,ti);
  print(gcf,'-dpng',plotTitle);
  
  %
  % display image with circles drawn on top
  clf; hold on;
  imagesc(Ivar); axis square; axis xy; axis off; colormap gray;
  circleT = linspace(0,2*pi,500);
  plot(imageDims(1)/2+r(1)*cos(circleT),imageDims(2)/2+r(1)*sin(circleT),'b','LineWidth',3);
  plot(imageDims(1)/2+r(2)*cos(circleT),imageDims(2)/2+r(2)*sin(circleT),'g','LineWidth',3);
  plot(imageDims(1)/2+r(3)*cos(circleT),imageDims(2)/2+r(3)*sin(circleT),'r','LineWidth',3);
  plotTitle = sprintf('%sVar_ti%03d.png',outFilePrefixCurrent,ti);
  print(gcf,'-dpng',plotTitle);  

  %
  % plot and save optimization record
  save(sprintf('ISELog%03d.mat',ti),'VarLog');
  clf;
  plot(VarLog);
  title(sprintf('ISE Log t=%g',t(ti)));
  xlabel('Iteration');
  ylabel('ISE');
  plotTitle = sprintf('%sVarLog_ti%03d.png',outFilePrefixCurrent,ti);
  print(gcf,'-dpng',plotTitle);
end
