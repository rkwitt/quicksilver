function [Isample, Iinfo] = createImageSample(sampleSize, type, imSize)

if nargin < 3
  imSize = 50;
end
if nargin < 2
  type = 'disk';
end
if nargin < 1
  sampleSize = 4;
end

Iinfo = [];

if (strcmp(type,'disk2'))
  for i = 1:sampleSize
    I = zeros(imSize);

    % compute distance from center
    center = imSize/2;
    [X,Y] = meshgrid(1:imSize,1:imSize);
    d = sqrt((X-center).*(X-center) + (Y-center).*(Y-center));

    % uniformly sample radius from 1/3 to 1/8 the image size
    % set image to foreground in disk
    r(i) = 1 / (3 + floor(i/sampleSize * 6)) * imSize;
    I(d<r(i)) = 1;

    % smooth image
    h = fspecial('gaussian',5,1);
    Isample(:,:,i) = imfilter(I,h,'circular');
  end
  fprintf('Radius Mean Arithmetic: %g Geometric: %g\n',mean(r),geomean(r));
  Iinfo = [center, center, mean(r), geomean(r)];
elseif (strcmp(type,'disk'))
  for i = 1:sampleSize
    I = zeros(imSize);

    % compute distance from center
    center = imSize/2;
    [X,Y] = meshgrid(1:imSize,1:imSize);
    d = sqrt((X-center).*(X-center) + (Y-center).*(Y-center));

    % get a random radius from 1/3 to 1/8 the image size
    % set image to foreground in disk
    r(i) = 1 / (3 + floor(rand() * 6)) * imSize;
    I(d<r(i)) = 1;

    % smooth image
    h = fspecial('gaussian',5,1);
    Isample(:,:,i) = imfilter(I,h,'circular');
  end
  fprintf('Radius Mean Arithmetic: %g Geometric: %g\n',mean(r),geomean(r));
  Iinfo = [center, center, mean(r), geomean(r)];  
elseif (strcmp(type,'ellipse_rotated'))
  for i = 1:sampleSize
    % smooth image
    I = createEllipseImage(20,10,2*pi*i/sampleSize,imSize,imSize,1,0);
    h = fspecial('gaussian',5,1);
    Isample(:,:,i) = imfilter(I,h,'circular');
  end
elseif (strcmp(type,'ellipse_bent'))
  for i = 1:sampleSize
    % smooth image
    beta = -0.05 + (i-1)/(sampleSize-1) * 0.1;
    I = createBentEllipse(20,10,0,beta,0,imSize,imSize,1,0);
    h = fspecial('gaussian',5,1);
    Isample(:,:,i) = imfilter(I,h,'circular');
  end
elseif (strcmp(type,'bull'))
  for i = 1:sampleSize
    I = zeros(imSize);

    % compute distance from center
    center = imSize/2;
    [X,Y] = meshgrid(1:imSize,1:imSize);
    d = sqrt((X-center).*(X-center) + (Y-center).*(Y-center));

    % get a random radius from 1.2*(1/3 to 1/8) the image size
    r(i) = 1 / (3 + floor(rand() * 6)) * imSize * 1.2;
    I(d<r(i)) = 0.5;
    rm(i) = 3*r(i)/5 + rand()*r(i)/5;
    I(d<rm(i)) = 1;
    rs(i) = r(i)/5 + rand()*r(i)/5;
    I(d<rs(i)) = 0.5;

    % smooth image
    h = fspecial('gaussian',5,0.5);
    Isample(:,:,i) = imfilter(I,h,'circular');
  end 
  fprintf('Radius Mean Arithmetic: %g Geometric: %g\n',mean(r),geomean(r));
  fprintf('Radius Mean Arithmetic: %g Geometric: %g\n',mean(rm),geomean(rm));
  fprintf('Radius Mean Arithmetic: %g Geometric: %g\n',mean(rs),geomean(rs)); 
  Iinfo = [center, center, mean(r), geomean(r), mean(rm), geomean(rm), mean(rs), geomean(rs)];  
elseif (strcmp(type,'stripedisk2'))
  % use stripedisk with some crazy cases
  for i = 1:sampleSize
    I = zeros(imSize);

    % compute distance from center
    center = imSize/2;
    [X,Y] = meshgrid(1:imSize,1:imSize);
    d = sqrt((X-center).*(X-center) + (Y-center).*(Y-center));

    % get a random radius from 1/3 to 1/8 the image size
    % set image to foreground in disk
    r = 1 / (3 + floor(rand() * 6)) * imSize;

    I(d<r*1.2) = 0.5;
    if i == 1
      % do nothing
    elseif i == 2
      I(d<r) = 0;
    else
      I(d<r) = 1;
    end
    I(d<r*0.3) = 0.5;

    % smooth image
    h = fspecial('gaussian',5,0.5);
    Isample(:,:,i) = imfilter(I,h,'circular');
  end
elseif (strcmp(type,'brain'))
  for i = 1:4
    filename = sprintf('images/slice%d.png',i);
    I0 = imread(filename,'png');
    I0 = double(I0(:,:,1));
    I0 = I0 / max(max(I0));
    backgroundSize = 224;
    background = zeros(backgroundSize,backgroundSize);
    insetStart = floor((backgroundSize-size(I0))/2);
    insetEnd   = insetStart + size(I0) - 1;
    background(insetStart(1):insetEnd(1),insetStart(2):insetEnd(2)) = I0;
    Isample(:,:,i) = background;
  end
end
