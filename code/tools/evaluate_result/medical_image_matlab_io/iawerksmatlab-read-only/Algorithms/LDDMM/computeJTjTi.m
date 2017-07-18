function [JTjTi, dPhiTiTj] = computeJTjTi(I,v)

sizeI = size(I);
numTimepoints = sizeI(end);
hID = eyeHField(sizeI(1:end-1));
phi = hID;

switch ndims(I)
  case 4,
    for i = 1:numTimepoints-1
      JTjTi{i}             = zeros([sizeI(1:end-1) numTimepoints-i+1]);
      dPhiTiTj{i}          = zeros([sizeI(1:end-1) numTimepoints-i+1]);
      JTjTi{i}(:,:,:,1)    = I(:,:,:,i);
      dPhiTiTj{i}(:,:,:,1) = 1;
    end
    
    for j = numTimepoints:-1:2
      fprintf('%d',j);
      phi = hID;
      
      for i = j-1:-1:1
        % compose incremental deformation and pull back image
        phi = composeHFields(phi, hID + v(:,:,:,:,i));
        JTjTi{i}(:,:,:,j-i+1)    = deformImage(I(:,:,:,j), phi);
        
        % compute jacobian of this incremental deformation
        dJv = jacobianDetHField(hID + v(:,:,:,:,i));
        % apply composition rule to get det jacobian of phi
        if i == j-1
          dPhiTiTj{i}(:,:,:,j-i+1) = dJv;
        else
          dPhiTiTj{i}(:,:,:,j-i+1) = ...
            dJv .* deformImage(dPhiTiTj{i+1}(:,:,:,j-i), hID + v(:,:,:,:,i));
        end
      end
    end
  case 3,
    for i = 1:numTimepoints-1
      JTjTi{i}             = zeros([sizeI(1:end-1) numTimepoints-i+1]);
      dPhiTiTj{i}          = zeros([sizeI(1:end-1) numTimepoints-i+1]);
      JTjTi{i}(:,:,1)    = I(:,:,i);
      dPhiTiTj{i}(:,:,1) = 1;
    end
    
    for j = numTimepoints:-1:2
      fprintf('%d',j);
      phi = hID;
      
      for i = j-1:-1:1
        % compose incremental deformation and pull back image
        phi = composeHFields(phi, hID + v(:,:,:,i));
        JTjTi{i}(:,:,j-i+1)    = deformImage(I(:,:,j), phi);
        
        % compute jacobian of this incremental deformation
        dJv = jacobianDetHField(hID + v(:,:,:,i));
        % apply composition rule to get det jacobian of phi
        if i == j-1
          dPhiTiTj{i}(:,:,j-i+1) = dJv;
        else
          dPhiTiTj{i}(:,:,j-i+1) = ...
            dJv .* deformImage(dPhiTiTj{i+1}(:,:,j-i), hID + v(:,:,:,i));
        end
      end
    end
  otherwise,
    error('image must be 2 or 3 dimensional');
end