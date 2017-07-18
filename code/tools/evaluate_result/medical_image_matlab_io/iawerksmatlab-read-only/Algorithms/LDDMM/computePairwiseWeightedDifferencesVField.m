function pvf = computePairwiseWeightedDifferencesField(targetPts,basePts,...
  v, sigma, method)

vsize = size(v);
pvf = zeros(vsize,class(v));

% need to find a better way (vectorize)
for tPt = 1:size(targetPts,2)
  for bPt = 1:size(basePts,2)
    d = targetPts(:,tPt) - basePts(:,bPt,end);
    gamma_diff = 2/sigma * exp(-);
    for t = 1:vsize(end)
      appPoint = basePts(:,bPt,t);
      % nearest neighbor needs to change
      % need to pull d back by jacobian
      % need to multiply by distance weighting factor
      pvf(:,round(appPoint(1)),round(appPoint(2)),round(appPoint(3)),t) = d;
    end
  end
end

% be careful because there will be one more set of points than velocity
% fields
