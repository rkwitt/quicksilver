function vizLDMMState(handles, v, I, J0t, dPhiTiTj, sse, sig, epsilon, currentIter, maxIter, vizSlice)

imagesHandle = handles(1);
convergenceHandle = handles(2);
errorHistogramHandle = handles(3);

numImages = size(I,length(size(I)));

%
% plot convergence of integrated pixelwise variance at Ihat
%
figure(convergenceHandle);
subplot(2,2,1:2,'align');
hold off;
plot((0:length(sse)-1),sse/sse(1),'b.:');
hold on;
plot(0:length(epsilon)-1,epsilon/epsilon(1),'gs:');
axis([-1 maxIter+1 -0.1 1.1]);
xlabel('Interation Number');
ylabel('Normalized SSE');
title({'Sum of Squared Errors';sprintf('Current Iter=%03d/%03d Sigma=%g Epsilon=%g  SV=%g (%g%%)',...
  currentIter, maxIter, sig, epsilon(end),sse(end),100*sse(end)/sse(1))});

subplot(2,2,3,'align');
origDiff = I - repmat(J0t(:,:,:,1),[1,1,1,numImages]);
origDiff = origDiff(origDiff > 0.0001 | origDiff < -0.0001);
currentDiff = I(:)-J0t(:);
currentDiff = currentDiff(currentDiff > 0.0001 | currentDiff < -0.0001);

maxDiffOrig  = max(abs(origDiff));
meanDiffOrig = mean(origDiff);
%meanDiffCurrent = mean(currentDiff);

stdDiffOrig = std(origDiff);
%stdDiffCurrent = std(currentDiff);

%hist(origDiff, 300, 'r');
%hold on;
%hist(currentDiff, 300, 'b');
randspots = randperm(min([length(origDiff) length(currentDiff)]'));
numRandSpots = min(2000, length(randspots))
scatter(origDiff(randspots(1:numRandSpots)), currentDiff(randspots(1:numRandSpots)), '.');
axis equal;
hold on;
plot([-maxDiffOrig maxDiffOrig], [-maxDiffOrig maxDiffOrig], 'r:');
plot([maxDiffOrig -maxDiffOrig], [maxDiffOrig -maxDiffOrig], 'r:');
xlim([meanDiffOrig-5*stdDiffOrig meanDiffOrig+5*stdDiffOrig]);
ylim([meanDiffOrig-5*stdDiffOrig meanDiffOrig+5*stdDiffOrig]);
title('Intensity Error');
xlabel('Initial Error');
ylabel('Current Error');

subplot(2,2,4,'align')
imagesc(squeeze(dPhiTiTj{1}(end:-1:1,:,vizSlice,end))');
axis image; axis off; colormap jet; colorbar;
title('Det. Jacobian');

printFilename = sprintf('images2D/convergence_%03d_iter%03d.png', size(I,1), currentIter);
print('-dpng', printFilename);

%
% plot images along with deformed template
%
figure(imagesHandle);
for i = 1:numImages
  % initial image
  subplot(numImages, 5, 5*i-4, 'align');
  imagesc(squeeze(I(end:-1:1,:,vizSlice,i))', [0.4 0.9]);
  axis image; axis off;
  
  % deformed template
  subplot(numImages, 5, 5*i-3, 'align');
  imagesc(squeeze(J0t(end:-1:1,:,vizSlice,i))', [0.4 0.9]);
  axis image; axis off;
  
  % difference between image and static template
  subplot(numImages, 5, 5*i-2, 'align');
  imagesc(squeeze(J0t(end:-1:1,:,vizSlice,1) - I(end:-1:1,:,vizSlice,i))', [meanDiffOrig-3*stdDiffOrig meanDiffOrig+3*stdDiffOrig]);
  axis image; axis off; colorbar;

  % difference between image and deformed template
  subplot(numImages, 5, 5*i-1, 'align');
  imagesc(squeeze(J0t(end:-1:1,:,vizSlice,i) - I(end:-1:1,:,vizSlice,i))', [meanDiffOrig-3*stdDiffOrig meanDiffOrig+3*stdDiffOrig]);
  axis image; axis off; colorbar;
  
  % determinant jacobian
  if i < numImages
    subplot(numImages, 5, 5*i-0, 'align');
    imagesc(squeeze(dPhiTiTj{i}(end:-1:1,:,vizSlice,2))');
    axis image; axis off; colorbar;
    colormap jet;
  end
  
  if i < numImages
    vdist = v(:,:,:,:,i) .* v(:,:,:,:,i);
    maxVDist = max(vdist(:));
  else
    maxVDist = 0;
  end
  
  imageDiffs = (J0t(:,:,:,i) - I(:,:,:,i)).^2;
  currentSSE = sum(imageDiffs(:));
  origImageDiffs = (J0t(:,:,:,1) - I(:,:,:,i)).^2;
  origSSE = sum(origImageDiffs(:));
  
  if i > 1
    ssePercent = 100*currentSSE/origSSE;
  else
    ssePercent = 0;
  end
    
  title(sprintf('maxv=%g, sse=%g / %g%%', ...
    maxVDist, currentSSE, ssePercent));
end

printFilename = sprintf('images2D/images_%03d_iter%03d.png', size(I,1), currentIter);
print('-dpng', printFilename);

figure(errorHistogramHandle);
clf;
