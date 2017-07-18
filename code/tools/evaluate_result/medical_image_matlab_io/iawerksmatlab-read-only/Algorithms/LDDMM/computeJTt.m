function [JTt,dPhitT] = computeJTt(Ihat, v)

if isMatlab
  stderr=2;
end

sv = size(v);
N  = sv(end);
JTt = zeros(sv(2:end), class(v));
dPhitT = zeros(sv(2:end), class(v));
hID = eyeHField(sv(2:end-1));
h = hID;

switch ndims(Ihat)
 case 3,
  for i=N:-1:1
    fprintf(stderr,'%d',i);
    % compose this v field
    h = composeHFields(h,hID+v(:,:,:,:,i));
    % pull image back
    JTt(:,:,:,i) = deformImage(Ihat,h);

    % compute jacobian of this incremental xform
    dJv = jacobianDetHField(hID+v(:,:,:,:,i));
    % apply composition rule to get jacobian of h (entire xform)
    if (i==N)
      dPhitT(:,:,:,i) = dJv;
    else
      dPhitT(:,:,:,i) = dJv .* deformImage(dPhitT(:,:,:,i+1),hID+v(:,:,:,:,i));
    end
  end
 case 2,
  for i=N:-1:1
    fprintf(stderr,'%d',i);
    % compose this v field
    h = composeHFields(h,hID+v(:,:,:,i));
    % pull image back
    JTt(:,:,i) = deformImage(Ihat,h);

    % compute jacobian of this incremental xform
    dJv = jacobianDetHField(hID+v(:,:,:,i));
    % apply composition rule to get jacobian of h (entire xform)
    if (i==N)
      dPhitT(:,:,i) = dJv;
    else
      dPhitT(:,:,i) = dJv .* deformImage(dPhitT(:,:,i+1),hID+v(:,:,:,i));
    end
  end
 otherwise,
  error('image must be 2 or 3 dimensional');
end
