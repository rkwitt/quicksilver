function JP0t = computeJP0t(pts, v)

sv = size(v);
N  = sv(end);
JP0t = zeros([size(pts) N+1], class(pts));
hID = eyeHField(sv(2:end-1));

switch size(pts,1)
  case 3,
    JP0t(:,:,1) = pts;
    for i=2:N+1
      fprintf('%d',i);
      % push points forward
      JP0t(:,:,i) = deformPointsForward(JP0t(:,:,i-1), hID+v(:,:,:,:,i-1));
    end
  case 2,
    JP0t(:,:,1) = pts;
    for i=2:N+1
      fprintf('%d',i);
      % push points forward
      JP0t(:,:,i) = deformPointsForward(JP0t(:,:,i-1), hID+v(:,:,:,i-1));
    end
  otherwise,
    error('image must be 2 or 3 dimensional');
end
