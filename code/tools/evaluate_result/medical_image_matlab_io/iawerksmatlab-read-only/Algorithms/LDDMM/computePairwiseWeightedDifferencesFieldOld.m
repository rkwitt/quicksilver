function pvf = computePairwiseWeightedDifferencesField(targetPts,basePts,...
  v, sigma)

pvf = zeros(size(v),class(v));
eps_cutoff = 0.001;
addressedTotal = 0;
addressedSkipped = 0;

% uncomment these when we are ready for optomization
%gamma_diff_precomp_a = 2/sigma;
%gamma_diff_precomp_b = exp(-1/sigma);
%gamma_diff = 0;
% need to find a better way (vectorize)
for tPt = 1:size(targetPts,2)
  for bPt = 1:size(basePts,2)
    tPoint = targetPts(:,tPt);
    bPoint = basePts(:,bPt,end);
    d = tPoint - bPoint;
    % should do some precomputation outside of the inner loop
    %gamma_diff = gamma_diff_precomp_a * gamma_diff_precomp_b^(d'*d);
    gamma_diff = 2/sigma * exp(-d'*d/sigma);
    %gamma_diff = exp(-d'*d/sigma^2);
    addressedTotal = addressedTotal + 1;
    % this should speed it up but get it working first
    %if gamma_diff < eps_cutoff || d'*d < eps_cutoff
    %  addressedSkipped = addressedSkipped + 1;
    %  continue;
    %end
    for t = 1:size(v,ndims(v))
      appPoint = basePts(:,bPt,t);
      % nearest neighbor needs to change
      % need to pull d back by jacobian!!!!!!!
      %weights!!!
      pvf(:,round(appPoint(1)),round(appPoint(2)),t) = ...
        pvf(:,round(appPoint(1)),round(appPoint(2)),t) + ...
        gamma_diff * d;
      %appPoint
      %d
      %gamma_diff
      %tPoint
      %bPoint
      %gamma_diff * d
    end
  end
end

fprintf('Ratio of skipped to addressed = %d/%d\n',addressedSkipped,addressedTotal);
% be careful because there will be one more set of points than velocity
% fields
