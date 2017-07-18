function [pts, ptweights] = makeTestPoints(dim,type,params,minmaxpts,dtype)
if nargin < 3
  params = [1 0];
end
if nargin < 4
  minmaxpts = [10 10];
end
if nargin < 5
  dtype = 'double';
end

% decide how many points to generate
% generate random integer in the set minpts:maxpts
n = minmaxpts(2) - minmaxpts(1) + 1;
npts = ceil(n*rand) + minmaxpts(1) - 1;

switch type
  case 'sphere'
    switch length(dim)
      case 2,
        r = params(1);
        tEps = 1/npts;
        t = linspace(0,1-tEps,npts);
        pts = [dim(1)/2+r*cos(2*pi*t); dim(2)/2+r*sin(2*pi*t)];
      case 3,
        r = params(1);
        [X,Y,Z] = ellipsoid(dim(1)/2,dim(2)/2,dim(3)/2,...
          r,r,r,round(sqrt(npts)));
        pts = [X(:)';Y(:)';Z(:)'];        
    end
  case 'bull'
    switch length(dim)
      case 2,
        r = params(1);
        tEps = 1/npts;
        t = linspace(0,1-tEps,npts);
        
        % middle band
        rm = 3*r/5 + rand()*r/5;
        % inner band
        ri = r/5 + rand*r/5;

        pts = [dim(1)/2+[rm*cos(2*pi*t) ri*cos(2*pi*t)]; dim(2)/2+[rm*sin(2*pi*t) ri*sin(2*pi*t)]];

      case 3,
        % not done!!!
        r = params(1);
        [X,Y,Z] = ellipsoid(dim(1)/2,dim(2)/2,dim(3)/2,...
          r,r,r,round(sqrt(npts)));
        pts = [X(:)';Y(:)';Z(:)'];        
    end    
end

%weights
ptweights = ones([1 npts])/npts;

% add random noise
pts = pts + params(2)*randn(size(pts));
pts = cast(pts,dtype);
