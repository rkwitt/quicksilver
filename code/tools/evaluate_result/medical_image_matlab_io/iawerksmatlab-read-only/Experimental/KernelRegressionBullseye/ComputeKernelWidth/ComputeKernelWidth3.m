clear all; close all;

%
% load radii data
load('RadiiDataNoise3.mat');

%
% possible bandwidths
bandwidths     = 0.04:0.0001:0.05;
%bandwidths     = [0.0539];
%bandwidths    = [0.0682];
%bandwidths    = [0.1];

%
% where estimator will be applied
numTimePoints = 250;
t             = linspace(0,1,numTimePoints);

%
% needed for error checking
imageDims     = [256 256];
minDim        = min(imageDims(:));
rMultRange    = [0.25 0.75];
rInnerMin     = 1/8;
rOuterMax     = 4/9;
rMiddleMin    = rInnerMin+(rOuterMax-rInnerMin)/3;
rOuterMin     = rInnerMin+(rOuterMax-rInnerMin)*2/3;
innerMinR     = rInnerMin*minDim;
innerMaxR     = rMiddleMin*minDim;
middleMinR    = rMiddleMin*minDim;
middleMaxR    = rOuterMin*minDim;
outerMinR     = rOuterMin*minDim;
outerMaxR     = rOuterMax*minDim;

%
% do the regression
Estimates     = zeros(numTimePoints,3,length(bandwidths));
sse           = zeros(size(bandwidths));
for bi = 1:length(bandwidths)
  for ti = 1:length(t)
    numerator = zeros(1,3);
    denominator = 0;
    for oi = 1:length(tSample)
      k = regkernel(tSample(oi), t(ti), bandwidths(bi));
      numerator = numerator + k * rSample(oi,:);
      denominator = denominator + k;
    end
    Estimates(ti,:,bi) = numerator / denominator;

    gs1 = rMultRange(1) + f1(t(ti))*(rMultRange(2)-rMultRange(1));
    gs1 = innerMinR + gs1*(innerMaxR-innerMinR);

    gs2 = rMultRange(1) + f2(t(ti))*(rMultRange(2)-rMultRange(1));
    gs2 = middleMinR + gs2*(middleMaxR-middleMinR);

    gs3 = rMultRange(1) + f3(t(ti))*(rMultRange(2)-rMultRange(1));
    gs3 = outerMinR + gs3*(outerMaxR-outerMinR);

    sse1 = (Estimates(ti,1,bi) - gs1)^2;
    sse2 = (Estimates(ti,2,bi) - gs2)^2;
    sse3 = (Estimates(ti,3,bi) - gs3)^2;
    sse(bi) = sse(bi) + sse1 + sse2 + sse3;
  end
end

%
% plot results
figure; hold on;
[dummy, minsseindex] = min(sse);
plot(tSample, rSample(:,1),'xb');
plot(tSample, rSample(:,2),'xg');
plot(tSample, rSample(:,3),'xr');
plot(t,squeeze(Estimates(:,:,minsseindex)), 'k');
title(sprintf('Best Bandwidth: %g',bandwidths(minsseindex)));
xlabel('Bandwidth');
ylabel('ISE');
print('-dpng','BestRegression.png');

%
% plot errors
figure;
plot(bandwidths,sse);
title(sprintf('Best Bandwidth: %g',bandwidths(minsseindex)));
xlabel('t');
ylabel('Radii in Pixels');
print('-dpng','BestBandwidth.png');
