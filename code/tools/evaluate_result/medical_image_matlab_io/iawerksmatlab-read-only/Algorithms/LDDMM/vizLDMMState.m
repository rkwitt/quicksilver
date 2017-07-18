function vizLDMMState(handles, J0T, Ihat0, Ihat,...
  Ivar0,Ivar,JP0T, IhatP0, IhatP, avar, epsilon,maxIter)

overallFigureHandle = handles(1);
averageFigureHandle = handles(2);
initialSampleHandle = handles(3);

sJ0T = size(J0T);
M = sJ0T(end);

%
% initial sample
%
figure(initialSampleHandle);
isSide = ceil(sqrt(M));
for i = 1:M
  %isRow = ceil(i/isSide);
  %isCol = i - (isRow-1)*isSide;
  subplot(isSide,isSide,i,'align');
  imagesc(J0T(:,:,i));
  axis image; axis off;
  hold on;
  scatter(JP0T{i}(1,:),JP0T{i}(2,:),'b');
  colormap gray;
end
  
%
% mean view
%
% mean image
figure(averageFigureHandle);
imagesc(Ihat);
colormap('hot');
colorbar;
axis image; axis off;
title('Current Mean');
% current mean points
hold on;
scatter(IhatP(1,:),IhatP(2,:),'b');

%
% overall view
%


figure(overallFigureHandle);

spc = 2;
spr = 3;

%
% initial mean and ptwise variance
%

% initial mean image
subplot(spr,spc,1,'align'); 
imagesc(Ihat0);
colormap('gray');
colorbar;
axis image; axis off;
title('Initial Mean');

% initial mean points
hold on;
scatter(IhatP0(1,:),IhatP0(2,:));

subplot(spr,spc,2,'align'); 
imagesc(Ivar0);
colormap('hot');
colorbar;
axis image; axis off;
title('Initial Ptwise Variance');

%
% current mean and ptwise variance
%

% current mean image
subplot(spr,spc,3,'align');
imagesc(Ihat);
colormap('gray');
colorbar;
axis image; axis off;
title('Current Mean');

% current mean points
hold on;
scatter(IhatP(1,:),IhatP(2,:));

subplot(spr,spc,4,'align'); 
imagesc(Ivar);
colormap('hot');
colorbar;
axis image; axis off;
title('Current Ptwise Variance');

%
% plot convergence of integrated pixelwise variance at Ihat
%
subplot(spr,spc,5:6,'align');
plot((0:length(avar)-1)/M,avar/avar(1),'b.:'); 
hold on;
plot(0:length(avar(1:M:end))-1,avar(1:M:end)/avar(1),'r*-');
plot(0:length(epsilon)-1,epsilon/epsilon(1),'gs:');
axis([-1 maxIter+1 -0.1 1.1]);
xlabel('Interation Number');
ylabel('Normalized Variance');
title({'Integrated Pointwise Variance';sprintf('Epsilon=%g  SV=%g (%g%%)',...
  epsilon(end),avar(end),100*avar(end)/avar(1))});

