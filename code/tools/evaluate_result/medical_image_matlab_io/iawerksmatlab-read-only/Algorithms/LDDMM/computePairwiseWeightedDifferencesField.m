function pvf = computePairwiseWeightedDifferencesField(targetPts,basePts,...
  v, sigma)

dimsv = ndims(v)-2;
sv = size(v);
pvf = zeros(sv,class(v));
hID = eyeHField(sv(2:end-1));
phi_t = hID;

% ptwise jacobian matrix of phi_t 
Jphi_t = zeros([sv(1) sv],class(v));

for t = size(v,ndims(v)):-1:1
  fprintf('%d',t);
  % compose v_t to phi_{t+1,T} to get phi_{t,T}  
  switch dimsv
    case 2,
      phi_t = composeHFields(phi_t,hID+v(:,:,:,t));
    case 3,
      phi_t = composeHFields(phi_t,hID+v(:,:,:,:,t));      
  end

  % compute jacobian of phi_{t,T}
  % it may be better to use the chain rule here
  Jphi_t = jacobianHField(phi_t);
  
  % loop over all point pairs
  for tPt = 1:size(targetPts,2)
    targetPoint = targetPts(:,tPt);
    
    for bPt = 1:size(basePts,2)
      basePoint   = basePts(:,bPt,end);
      appPoint    = basePts(:,bPt,t);
      rAppPoint   = round(appPoint);
      
      % difference of pts
      d = targetPoint - basePoint;
      
      % compute nearness measure \gamma in my writeup
      % maybe we pull some of this out of the inner loop
      gamma_diff = 2/sigma * exp(-d'*d/sigma);

      % get jacobian of phi_{t,T} at the application point
      % use nearest neighbor for now---better idea?
      Jphi_t_appPoint = ...
        Jphi_t(:,:,rAppPoint(1),rAppPoint(2));
      
      % transform d by the jacobian and nearness measure
      pvf(:,rAppPoint(1),rAppPoint(2),t) =...
        pvf(:,rAppPoint(1),rAppPoint(2),t) + ...
        gamma_diff * Jphi_t_appPoint' *d;
    end
  end
end
